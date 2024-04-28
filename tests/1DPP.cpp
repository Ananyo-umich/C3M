// @sec3{Include files}
// Author: Ananyo Bhattacharya
// Affiliation: University of Michigan
// Email: ananyo@umich.edu
// This code represents a 1-d photochemistry model for planetary atmosphere.
// Note: The equations are solved for mass transport
#include <cantera/base/Solution.h>

// ThermoPhase object stores the thermodynamic state
#include <cantera/thermo.h>

// Kinetics object stores the chemical kinetics information
#include <cantera/kinetics.h>
#include <cantera/kinetics/Reaction.h>

// output stream
#include <cantera/numerics/eigen_dense.h>
#include <cantera/numerics/eigen_sparse.h>

#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>
// Athena++ header
#include <parameter_input.hpp>

// C3M header
#include <CustomRate.hpp>
#include <CustomTransport.hpp>
#include <PhotoChemistry.hpp>
#include <RadTran.hpp>
#include <interpolation.hpp>

// NetCDF Output
#if NETCDFOUTPUT
#include <netcdf.h>
#endif

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace Cantera;
using namespace std;

int main(int argc, char** argv) {
  // Reading input file
  IOWrapper infile;
  infile.Open("test_athena.inp", IOWrapper::FileMode::read);
  ParameterInput* pinput = new ParameterInput();

  pinput->LoadFromFile(infile);
  infile.Close();

  // Loading the input parameters for atmospheric profile and reaction network
  std::string atm_file = pinput->GetString("problem", "planet");
  std::string PlanetName = pinput->GetString("problem", "planetName");
  std::string network_file = pinput->GetString("problem", "network");
  std::string profile_file =
      pinput->GetOrAddString("problem", "chemInp", "nan");

  // Reading the chemical kinetics network
  auto sol = newSolution(network_file);
  auto gas = sol->thermo();
  auto gas_kin = sol->kinetics();
  int nsp = gas->nSpecies();
  int nrxn = gas_kin->nReactions();
  std::cout << "Number of reactions: " << nrxn << std::endl;
  std::cout << "Network imported" << std::endl;

  // Initial condition for mole fraction
  VectorXd mole_fractions = VectorXd::Zero(nsp);

  // Initial condition for boundary fluxes
  VectorXd flux_lower = VectorXd::Zero(nsp);
  VectorXd flux_upper = VectorXd::Zero(nsp);

  // Loading the input parameters for initial condition
  // Homogeneous input condition
  std::string init_species_list = pinput->GetString("init", "species");
  std::regex pattern("[a-zA-Z0-9_+-]+");
  std::smatch m;
  int species_inx;
  while (std::regex_search(init_species_list, m, pattern)) {
    for (auto x : m) {
      std::string species_init_condition = pinput->GetString("init", x);
      species_inx = gas->speciesIndex(x);
      mole_fractions(species_inx) = atof(species_init_condition.c_str());
    }
    init_species_list = m.suffix().str();
  }

  // Loading the input for upper boundary condition
  std::string ub_species_list =
      pinput->GetString("upperboundaryflux", "species");
  while (std::regex_search(ub_species_list, m, pattern)) {
    for (auto x : m) {
      std::string species_upperboundary_condition =
          pinput->GetString("upperboundaryflux", x);
      species_inx = gas->speciesIndex(x);
      flux_upper(species_inx) = atof(species_upperboundary_condition.c_str());
    }
    ub_species_list = m.suffix().str();
  }

  // Loading the input for lower boundary condition
  std::string lb_species_list =
      pinput->GetString("lowerboundaryflux", "species");
  while (std::regex_search(lb_species_list, m, pattern)) {
    for (auto x : m) {
      std::string species_lowerboundary_condition =
          pinput->GetString("lowerboundaryflux", x);
      species_inx = gas->speciesIndex(x);
      flux_lower(species_inx) = atof(species_lowerboundary_condition.c_str());
    }
    lb_species_list = m.suffix().str();
  }
  std::cout << "Boundary conditions loaded" << std::endl;
  // Integrator
  std::string time_step = pinput->GetString("integrator", "dt");
  std::string max_time = pinput->GetString("integrator", "Tmax");
  double dt = atof(time_step.c_str());
  double Tmax = atof(max_time.c_str());
  double Ttot = 0.0;
  std::cout << "Time step: " << dt << std::endl;

  // Atmosphere Properties, indices for I/O storage
  fstream InFile;
  int nSize = 0;
  string data1, data2, data3, data4, data5;
  InFile.open(atm_file);
  getline(InFile, data1);
  getline(InFile, data1);
  while (getline(InFile, data1)) nSize++;
  InFile.close();

  // Atmospheric Profile Data

  MatrixXd AtmData(5, nSize);
  int iTemp = 0;
  int iPress = 1;
  int iKzz = 2;
  int iAlt = 3;
  int iNd = 4;

  // Input from txt file
  int inx = 0;
  InFile.open(atm_file);
  getline(InFile, data1);
  getline(InFile, data1);
  double kmax;
  while (InFile >> data1 >> data2 >> data3 >> data4 >> data5) {
    AtmData(iPress, inx) = atof(data1.c_str()) * 1E2;  // Pressure (mbar) -> Pa
    AtmData(iTemp, inx) = atof(data2.c_str());         // Temperature (K)
    AtmData(iKzz, inx) = atof(data3.c_str()) * 1E-4;   // Kzz (cm^2/s) -> m^2/s
    AtmData(iAlt, inx) = atof(data4.c_str()) * 1E3;    // Altitude (m)
    AtmData(iNd, inx) = atof(data5.c_str()) * 1E6 /
                        (6.022E23 * 1E3);  // Concentration (1/cm^3) -> kmol/m^3
    if (inx == 0) {
      kmax = AtmData(iKzz, inx);
    }
    if (inx != 0) {
      if (kmax < AtmData(iKzz, inx)) {
        kmax = AtmData(iKzz, inx);
      }
    }

    inx++;
  }
  InFile.close();

  // Setting initial condition for chemical species
  std::vector<Eigen::MatrixXd> ChemMole;
  MatrixXd ChemMoleFrac(nsp, nSize);
  MatrixXd ChemConc(nsp, nSize);
  MatrixXd TChem(nsp, nSize);
  MatrixXd TDyn(nsp, nSize);
  MatrixXd a(nsp, nSize);
  MatrixXd n_conc(nsp, nSize);
  MatrixXd conv(nsp, nSize);
  MatrixXd ProdRates(nsp, nSize);
  MatrixXd DiffRates(nsp, nSize);
  VectorXd Time(1);

  for (int i = 0; i < nSize; i++) {
    gas->setState_TP(AtmData(iTemp, i),
                     (AtmData(iPress, i) / 1.0132E5) * OneAtm);
    ChemMoleFrac.col(i) = mole_fractions;
  }

  // Chemical species profile from input file
  if (profile_file != "nan") {
    std::string input_species_list = pinput->GetString("profile", "species");
    std::regex pattern("[a-zA-Z0-9_+-]+");
    int chemP = 0;
    while (std::regex_search(input_species_list, m, pattern)) {
      for (auto x : m) {
        chemP++;
      }
      input_species_list = m.suffix().str();
    }
    InFile.open(profile_file);
    std::string inp1, inp2;
    int i = 0;
    getline(InFile, inp1);
    while (i < nSize) {
      int chem_inx = 0;
      input_species_list = pinput->GetString("profile", "species");
      while (std::regex_search(input_species_list, m, pattern)) {
        for (auto x : m) {
          species_inx = gas->speciesIndex(x);
          if (chem_inx < chemP - 1) {
            getline(InFile, inp2, ',');
            ChemMoleFrac(species_inx, i) = atof(inp2.c_str());
          }

          if (chem_inx == chemP - 1) {
            getline(InFile, inp2, '\n');
            ChemMoleFrac(species_inx, i) = atof(inp2.c_str());
          }

          chem_inx++;
        }

        input_species_list = m.suffix().str();
      }
      i++;
    }
  }
  InFile.close();

  for (int i = 0; i < nSize; i++) {
    gas->setState_TP(AtmData(iTemp, i),
                     (AtmData(iPress, i) / 1.0132E5) * OneAtm);
    gas->setMoleFractions(&ChemMoleFrac.col(i)[0]);
    gas->getConcentrations(&ChemConc.col(i)[0]);
  }

  std::cout << "All inputs loaded into C3M " << std::endl;

  // Condition for turning on photochemistry
  std::string photochem = pinput->GetString("problem", "photochem");
  if (photochem == "true") {
    std::cout << "Starting photochemistry calculations" << std::endl;
    // Reading the Stellar Irradiance Input
    std::string stellar_input_file = pinput->GetString("radtran", "solar");
    std::string radius = pinput->GetString("radtran", "radius");
    std::string reference = pinput->GetString("radtran", "reference");
    std::string SZA = pinput->GetString("radtran", "SZA");
    double rad = atof(radius.c_str());
    double ref = atof(reference.c_str());
    double sz_angle = atof(SZA.c_str());
    MatrixXd stellar_input =
        ReadStellarRadiationInput(stellar_input_file, rad, ref);
    MatrixXd d_wavelength(stellar_input.row(0).size(), 1);
    MatrixXd wavelength = stellar_input.row(0);

    // Modification of solar flux based on SZA input, along with conversion from
    // deg. to radians
    stellar_input.row(1) = stellar_input.row(1) * cos(sz_angle * 3.14 / 180);

    // Discretization of the wavelength
    for (int dw = 0; dw < (stellar_input.row(0).size()) - 1; dw++) {
      d_wavelength(dw) = wavelength(dw + 1) - wavelength(dw);
    }

    std::cout << "!! Radiation Input loaded in C3M !!" << std::endl;
    // Input Wavelength and Cross Section Data for all the absorbers
    int AbsorberSize = 0;
    std::string absorber_species_list =
        pinput->GetString("abscross", "absorbers");
    while (std::regex_search(absorber_species_list, m, pattern)) {
      for (auto x : m) {
        AbsorberSize++;
      }
      absorber_species_list = m.suffix().str();
    }
    Eigen::MatrixXd absorber_cross_data(stellar_input.row(0).size(),
                                        AbsorberSize);
    std::cout << "!! Initiating absorbers  !!" << std::endl;

    // Reading absorber cross section for all absorbers
    int ab_inx = 0;
    absorber_species_list = pinput->GetString("abscross", "absorbers");
    while (std::regex_search(absorber_species_list, m, pattern)) {
      for (auto x : m) {
        std::string absorber_cross_info = pinput->GetString("abscross", x);
        std::string absorber_cross_database =
            absorber_cross_info.substr(0, absorber_cross_info.find(","));
        std::string absorber_cross_file =
            absorber_cross_info.substr(absorber_cross_info.find(",") + 1,
                                       absorber_cross_info.length() - 1);

        // Reading the absorption cross section for absorber as per the database
        // structure Interpolating the cross sections to reference grid
        if (absorber_cross_database == "VULCAN") {
          MatrixXd cross_info =
              ReadVULCANPhotoAbsCrossSection(absorber_cross_file);
          MatrixXd cross_section_data = InterpolateCrossSection(
              stellar_input.row(0), cross_info.row(0), cross_info.row(1));
          absorber_cross_data.col(ab_inx) = cross_section_data;
        }

        /*
           if(absorber_cross_database == "KINETICS"){
             std::cout << absorber_cross_database << std::endl;

            }
        */

        if (absorber_cross_database == "VPL") {
          MatrixXd cross_info = ReadAtmosCrossSection(absorber_cross_file);
          MatrixXd cross_section_data = InterpolateCrossSection(
              stellar_input.row(0), cross_info.row(0), cross_info.row(1));
          absorber_cross_data.col(ab_inx) = cross_section_data;
        }

        // species_inx = gas->speciesIndex(x);
        ab_inx++;
      }
      absorber_species_list = m.suffix().str();
    }
    std::cout << "!! Obtaining Photochemical Cross Sections !!" << std::endl;

    // Input Wavelength and Cross Section Data for all the absorbers
    int ScatSize = 0;
    std::string scatter_species_list =
        pinput->GetString("scatcross", "scatterers");
    while (std::regex_search(scatter_species_list, m, pattern)) {
      for (auto x : m) {
        ScatSize++;
      }
      scatter_species_list = m.suffix().str();
    }
    Eigen::MatrixXd scat_cross_data(stellar_input.row(0).size(), ScatSize);
    std::cout << "!! Initiating scattering  !!" << std::endl;

    // Reading the scattering cross section for atmospheric constituents as per
    // database structure
    int sc_inx = 0;
    scatter_species_list = pinput->GetString("scatcross", "scatterers");
    while (std::regex_search(scatter_species_list, m, pattern)) {
      for (auto x : m) {
        std::string scat_cross_info = pinput->GetString("scatcross", x);
        std::string scat_cross_database =
            scat_cross_info.substr(0, scat_cross_info.find(","));
        std::string scat_cross_file = scat_cross_info.substr(
            scat_cross_info.find(",") + 1, scat_cross_info.length() - 1);

        // Reading the absorption cross section for absorber as per the database
        // structure Interpolating the cross sections to reference grid
        if (scat_cross_database == "VULCAN") {
          MatrixXd sccross_info = ReadVULCANScatCrossSection(scat_cross_file);
          MatrixXd sccross_section_data = InterpolateCrossSection(
              stellar_input.row(0), sccross_info.row(0), sccross_info.row(1));
          scat_cross_data.col(sc_inx) = sccross_section_data;
        }
        sc_inx++;
      }
      scatter_species_list = m.suffix().str();
    }

    // Storing the photochemical cross section data
    int PhotoRxn = 0;

    for (int irxn = 0; irxn < nrxn; irxn++) {
      auto& rxnObj = *(gas_kin->reaction(irxn));
      std::string rxnEquation = rxnObj.equation();
      std::cout << rxnEquation << std::endl;
      int pos = rxnEquation.find("=");
      rxnEquation.replace(pos, 2, "->");
      std::string photo_cross_info =
          pinput->GetOrAddString("photocross", rxnEquation, "nan");
      // std::cout << photo_cross_info << std::endl;
      if (photo_cross_info != "nan") {
        PhotoRxn++;
      }
    }

    Eigen::MatrixXd photo_cross_data(stellar_input.row(0).size(), PhotoRxn);
    Eigen::MatrixXd qyield_data(stellar_input.row(0).size(), PhotoRxn);
    Eigen::VectorXd RxnIndex(PhotoRxn);
    Eigen::VectorXd Jrate(PhotoRxn);

    // Reading cross section for all photochemical reactions
    int ph_inx = 0;
    for (int irxn = 0; irxn < nrxn; irxn++) {
      auto& rxnObj = *(gas_kin->reaction(irxn));
      std::string rxnEquation = rxnObj.equation();
      int pos = rxnEquation.find("=");
      rxnEquation.replace(pos, 2, "->");
      std::string photo_cross_info =
          pinput->GetOrAddString("photocross", rxnEquation, "nan");
      if (photo_cross_info != "nan") {
        std::cout << rxnEquation << std::endl;
        RxnIndex(ph_inx) = irxn;
        std::string photo_cross_database =
            photo_cross_info.substr(0, photo_cross_info.find(","));
        std::string photo_cross_ = photo_cross_info.substr(
            photo_cross_info.find(",") + 1, photo_cross_info.length() - 1);

        if (photo_cross_database == "VULCAN_Ion") {
          MatrixXd photoXsection_info =
              ReadVULCANPhotoIonCrossSection(photo_cross_);
          photo_cross_data.col(ph_inx) = InterpolateCrossSection(
              stellar_input.row(0), photoXsection_info.row(0),
              photoXsection_info.row(1));
        }

        if (photo_cross_database == "VULCAN_Diss") {
          MatrixXd photoXsection_info =
              ReadVULCANPhotoDissCrossSection(photo_cross_);
          photo_cross_data.col(ph_inx) = InterpolateCrossSection(
              stellar_input.row(0), photoXsection_info.row(0),
              photoXsection_info.row(1));
          // std::cout <<    photo_cross_data.col(ph_inx).transpose() <<
          // std::endl;
        }

        if (photo_cross_database == "VPL") {
          MatrixXd photoXsection_info = ReadAtmosCrossSection(photo_cross_);
          photo_cross_data.col(ph_inx) = InterpolateCrossSection(
              stellar_input.row(0), photoXsection_info.row(0),
              photoXsection_info.row(1));
        }

        if (photo_cross_database == "MPD") {
          MatrixXd photoXsection_info = ReadMPCrossSection(photo_cross_);
          photo_cross_data.col(ph_inx) = InterpolateCrossSection(
              stellar_input.row(0), photoXsection_info.row(0),
              photoXsection_info.row(1));
        }

        if (photo_cross_database == "KINETICS") {
          MatrixXd photoXsection_info = ReadKINETICSCrossSection(
              atoi(photo_cross_.c_str()));  // In this case its not file name
                                            // but the reaction number
          photo_cross_data.col(ph_inx) = InterpolateCrossSection(
              stellar_input.row(0), photoXsection_info.row(0),
              photoXsection_info.row(1));
        }

        //**Space for custom Photochemical cross section database**

        ph_inx++;
      }
    }

    std::cout << "Photochemical cross sections stored in C3M" << std::endl;
    ph_inx = 0;

    // Reading Quantum Yield for photochemical reactions in the network
    for (int irxn = 0; irxn < nrxn; irxn++) {
      auto& rxnObj = *(gas_kin->reaction(irxn));
      std::string rxnEquation = rxnObj.equation();
      int pos = rxnEquation.find("=");
      rxnEquation.replace(pos, 2, "->");
      std::string qyield_info =
          pinput->GetOrAddString("qyield", rxnEquation, "nan");
      if (qyield_info != "nan") {
        std::string QYType = qyield_info.substr(0, qyield_info.find(";"));

        // Quantum Yield for VULCAN database
        if (QYType == "VULCAN") {
          std::string col_str = qyield_info.substr(
              qyield_info.find(";") + 1,
              qyield_info.find(",") - qyield_info.find(";") - 1);
          int col_num = atoi(col_str.c_str());
          std::string qyield_file = qyield_info.substr(
              qyield_info.find(",") + 1,
              qyield_info.length() - qyield_info.find(",") - 1);
          MatrixXd QYield_info = ReadQYield(qyield_file);
          qyield_data.col(ph_inx) =
              InterpolateQYield(stellar_input.row(0), QYield_info.row(0),
                                QYield_info.row(col_num));
        }

        // Quantum Yield for KINETICS7 database (QY == 1)
        if (QYType == "KINETICS") {
          qyield_data.col(ph_inx) =
              MatrixXd::Ones(stellar_input.row(0).size(), 1);
        }

        //**Space for custom Quantum Yield database**

        ph_inx++;
      }
    }
    std::cout << "Quantum Yields stored in C3M" << std::endl;

    // Chemical evolution
    double dh;  // m
    int iPrev, iNext;

    // Generate a opacity profile for the atmosphere. The opacity is calculated
    // at every height
    MatrixXd Stellar_activity =
        MatrixXd::Zero(stellar_input.row(0).size(), nSize);
    VectorXd one = VectorXd::Ones(stellar_input.row(0).size());

    // Initiating Matrices for Time Evolution
    Eigen::SparseMatrix<double> m_wjac;
    m_wjac.resize(nsp, nsp);
    MatrixXd mat1(nsp, nsp);
    mat1 = MatrixXd::Identity(nsp, nsp);
    MatrixXd mat2(nsp, nsp);
    VectorXd Un(nsp);
    VectorXd alpha = VectorXd::Zero(nsp);  // Thermal diffusion coefficient
    VectorXd Unext(nsp);
    VectorXd Uprev(nsp);
    VectorXd flux1(nsp);
    VectorXd dQ(nsp);
    double Keddy_j;
    double Keddy_prev;
    double Keddy_next;
    double Temp, Press, Kzz;
    VectorXd m_wdot(nsp);
    VectorXd mole_frac(nsp);

    double Kb = 1.38E-23;
    int i = 0;

    std::string a_species_list = pinput->GetString("thermaldiff", "species");
    while (std::regex_search(a_species_list, m, pattern)) {
      for (auto x : m) {
        std::string species_thermaldiff_condition =
            pinput->GetString("thermaldiff", x);
        species_inx = gas->speciesIndex(x);
        alpha(species_inx) = atof(species_thermaldiff_condition.c_str());
      }
      a_species_list = m.suffix().str();
    }

    std::string equbm = pinput->GetString("problem", "equilibrate");
    if (equbm == "true") {
      std::cout << "!!  Initiating Equilibrium  !!" << std::endl;

      // Setting Thermochemical Equilibrium
      for (int j = 0; j < nSize; j++) {
        Temp = AtmData(iTemp, j);
        Press = AtmData(iPress, j);
        gas->setState_TP(Temp, (Press / 1.0132E5) * OneAtm);
        mole_frac = ChemMoleFrac.col(j);
        gas->setMoleFractions(&mole_frac[0]);
        gas->equilibrate("TP");
      }
      std::cout << "!! Atmosphere in thermochemical equilibrium !! \n"
                << std::endl;
    }

    // Extracting boundary condition flags from input file
    std::string top = pinput->GetString("problem", "top");
    std::string bot = pinput->GetString("problem", "bot");

    double err = 1e-12;
    double tol = 0.01;
    std::cout << "!! Starting the Simulation !!" << std::endl;
    conv = ChemMoleFrac;
    double old_conv;
    double mm_prev, mm_next;
    // Acceleration due to gravity (m/s)
    std::string grav_str = pinput->GetString("grav", "g");
    double g = atof(grav_str.c_str());

    // Correct for the atmospheric scattering and absorption above the model
    // boundary Using column density of species to determine the absorption and
    // scattering not accounted for
    int clden = 0;
    std::string clden_species_list =
        pinput->GetOrAddString("coldensity", "species", "nan");
    std::cout << clden_species_list << std::endl;
    Eigen::VectorXd cldMag =
        VectorXd::Zero(gas->nSpecies());  // Stores the magnitude of column
                                          // density difference for TOA
    Eigen::VectorXd cldMag_diff =
        VectorXd::Zero(gas->nSpecies());  // Stores the magnitude of column
                                          // density difference for each species
    while (std::regex_search(clden_species_list, m, pattern)) {
      for (auto x : m) {
        // Column density at top of atmosphere [From input file]
        std::string clden_toa = pinput->GetOrAddString("coldensity", x, "nan");
        if (clden_toa != "nan") {
          species_inx = gas->speciesIndex(x);
          cldMag(species_inx) =
              atof(clden_toa.c_str()) * 1e4;  // 1/cm^2 -> 1/m^2
        }
        clden++;
      }
      clden_species_list = m.suffix().str();
    }

    // Initialize the block matrix
    Eigen::MatrixXd BlockMatrix = MatrixXd::Zero(nsp * nSize, nsp * nSize);
    Eigen::MatrixXd Upresent = VectorXd::Zero(nsp * nSize);
    Eigen::MatrixXd Ufuture = VectorXd::Zero(nsp * nSize);
    VectorXd mWt(nsp);  // Molecular weight
    VectorXd B = VectorXd::Zero(nsp);
    VectorXd A = VectorXd::Zero(nsp);
    VectorXd C = VectorXd::Zero(nsp);
    MatrixXd conv(nsp, nSize);
    int counter = 0;

    std::cout << "Time starts" << std::endl;
    while (Ttot < Tmax) {
      MatrixXd Opacity = MatrixXd::Zero(stellar_input.row(0).size(), nSize);
      for (int j = 0; j < nSize; j++) {
        // Set the output of previous step as present condition of the system

        // Setting T, P, X for each grid point
        iPrev = j - 1;
        iNext = j + 1;
        auto sol2 = newSolution(network_file);
        auto gas2 = sol2->thermo();
        auto gas_kin2 = sol2->kinetics();
        Temp = AtmData(iTemp, j);
        Press = AtmData(iPress, j);
        mole_frac = ChemMoleFrac.col(j);
        gas2->setMoleFractions(&mole_frac[0]);
        gas2->setState_TP(Temp, (Press / 1.0132E5) * OneAtm);
        Cantera::Kinetics* gasRawPtr = sol2->kinetics().get();
        Cantera::ThermoPhase* gasThermo = gas2.get();

        // Radiative transfer (Beer-Lambert law) for actinic flux, and rate
        // calculation
        if (j != 0) {
          // For each absorber, find the total absorption
          dh = (AtmData(iAlt, j - 1) - AtmData(iAlt, j));

          // Atomic and Molecular absorption
          int Absorber = 0;
          std::string absorber_species_list =
              pinput->GetString("abscross", "absorbers");
          while (std::regex_search(absorber_species_list, m, pattern)) {
            for (auto x : m) {
              species_inx = gas2->speciesIndex(x);
              double number_density = Press / (Kb * Temp);
              Opacity.col(j) = Opacity.col(j) -
                               (mole_fractions(species_inx) * number_density *
                                absorber_cross_data.col(Absorber) * dh /
                                cos(sz_angle * 3.14 / 180));
              Absorber++;
            }
            absorber_species_list = m.suffix().str();
          }

          // Rayleigh scattering
          int Scatter = 0;
          std::string scat_species_list =
              pinput->GetString("scatcross", "scatterers");
          while (std::regex_search(absorber_species_list, m, pattern)) {
            for (auto x : m) {
              species_inx = gas2->speciesIndex(x);
              double number_density = Press / (Kb * Temp);
              Opacity.col(j) = Opacity.col(j) -
                               (mole_fractions(species_inx) * number_density *
                                scat_cross_data.col(Scatter) * dh /
                                cos(sz_angle * 3.14 / 180));
              Scatter++;
            }
            scat_species_list = m.suffix().str();
          }

          // The transmission coefficient from opacity
          Opacity.col(j) = Opacity.col(j) -
                           (dh * handleCustomOpacity(
                                     PlanetName, gasThermo, Press, Temp,
                                     AtmData(iAlt, j), stellar_input.row(0)));
          // OpacityStorage.col(j) = Opacity.col(j);
          Opacity.col(j) = Opacity.col(j).array().exp().matrix();
          Opacity.col(j) =
              (Opacity.col(j).array() * Opacity.col(j - 1).array()).matrix();
          // Stellar spectrum at each altitude
          Stellar_activity.col(j) =
              (Stellar_activity.col(j - 1).array() * Opacity.col(j).array())
                  .transpose()
                  .matrix();
          // std::cout << stellar_activity.col(j) << std::endl;
        }
        if (j == 0) {
          clden_species_list =
              pinput->GetOrAddString("coldensity", "species", "nan");
          if (clden_species_list != "nan") {
            // Updating the difference in column density
            cldMag_diff = cldMag;  // 1/m^2
            // Atomic and molecular absorption
            int Absorber = 0;
            std::string absorber_species_list =
                pinput->GetString("abscross", "absorbers");
            while (std::regex_search(absorber_species_list, m, pattern)) {
              for (auto x : m) {
                species_inx = gas2->speciesIndex(x);
                double number_density = Press / (Kb * Temp);
                Opacity.col(j) =
                    Opacity.col(j) -
                    (absorber_cross_data.col(Absorber) *
                     cldMag_diff(species_inx) / cos(sz_angle * 3.14 / 180));
                // std::cout <<
                // absorber_cross_data.col(Absorber).transpose()*cldMag_diff(species_inx)/cos(sz_angle*3.14/180)
                // << std::endl;
                // std::cout << Opacity.col(j).transpose() << std::endl;
                Absorber++;
              }
              absorber_species_list = m.suffix().str();
            }

            // std::cout << "Absorption done! " << std::endl;
            // Rayleigh scattering
            int Scatter = 0;
            std::string scat_species_list =
                pinput->GetString("scatcross", "scatterers");
            while (std::regex_search(absorber_species_list, m, pattern)) {
              for (auto x : m) {
                species_inx = gas2->speciesIndex(x);
                double number_density = Press / (Kb * Temp);
                Opacity.col(j) =
                    Opacity.col(j) -
                    (scat_cross_data.col(Scatter) * cldMag_diff(species_inx) /
                     cos(sz_angle * 3.14 / 180));
                Scatter++;
              }
              scat_species_list = m.suffix().str();
            }
            // std::cout << "Scattering done at TOA!" << std::endl;
            // Updating the stellar activity at upper boundary
            // The transmission coefficient from opacity
            //  std::cout << Opacity.col(j) << std::endl;
            //  OpacityStorage.col(j) = Opacity.col(j);
            Opacity.col(j) = Opacity.col(j).array().exp().matrix();
            // Stellar spectrum at each altitude
            Stellar_activity.col(j) =
                (stellar_input.row(1).transpose().array() *
                 Opacity.col(j).array())
                    .transpose()
                    .matrix();
            // std::cout << "TOA" << std::endl;
          }

          if (clden_species_list == "nan") {
            Opacity.col(0) = VectorXd::Ones(stellar_input.row(0).size());
            Stellar_activity.col(j) = (stellar_input.row(1).transpose());
          }
          // std::cout << "Setting TOA opacity" << std::endl;
        }

        for (int rx = 0; rx < PhotoRxn; rx++) {
          double j_rate = QPhotoChemRate(
              stellar_input.row(0), d_wavelength, photo_cross_data.col(rx),
              qyield_data.col(rx), Stellar_activity.col(j));
          gas_kin2->setMultiplier(RxnIndex(rx), j_rate);
          auto& rxnObj = *(gas_kin->reaction(RxnIndex(rx)));
          std::string rxnEquation = rxnObj.equation();
          // std::cout<< j << " " << AtmData(iAlt, j)/1e3 << " " << rxnEquation
          // << " " << j_rate << std::endl;
        }

        // Solving the net production for each species
        gas_kin2->getNetProductionRates(
            &m_wdot[0]);  // Extracting net production rates from Cantera
        m_wjac = gas_kin2->netProductionRates_ddCi();  // Extracting Jacobian
                                                       // from Cantera
        for (int insp = 0; insp < nsp; insp++) {
          mWt(insp) = gas2->molecularWeight(insp);  // Kg/kml
        }

        // Chemical time scales
        TChem.col(j) =
            abs(mole_frac.array() * AtmData(iNd, j) / m_wdot.array()).matrix();
        TChem.col(j) =
            (isnan(abs(TChem.col(j).array()))).select(1e20, TChem.col(j));
        TChem.col(j) =
            (isinf(abs(TChem.col(j).array()))).select(1e20, TChem.col(j));
        TChem.col(j) =
            (TChem.col(j).array().abs() == 0).select(1e20, TChem.col(j));
        //    std::cout << TChem.col(j).transpose() << std::endl;
        // Solve for diffusion terms
        if ((j > 0) && (j < nSize - 1)) {
          mole_frac = ChemMoleFrac.col(j - 1);
          gas2->setMoleFractions(&mole_frac[0]);
          mm_prev = gas2->meanMolecularWeight();  // Mean molcular weight

          mole_frac = ChemMoleFrac.col(j + 1);
          gas2->setMoleFractions(&mole_frac[0]);
          mm_next = gas2->meanMolecularWeight();  // Mean molcular weight

          mole_frac = ChemMoleFrac.col(j);
          gas2->setMoleFractions(&mole_frac[0]);
          gas2->setState_TP(Temp, (Press / 1.0132E5) * OneAtm);

          double K_next = (AtmData(iKzz, j + 1) + AtmData(iKzz, j)) / 2;
          double K_prev = (AtmData(iKzz, j - 1) + AtmData(iKzz, j)) / 2;
          VectorXd D_this = handleCustomMolecularDiffusion(
              PlanetName, gasThermo, AtmData(iPress, j), AtmData(iTemp, j),
              mWt);
          VectorXd D_next = handleCustomMolecularDiffusion(
              PlanetName, gasThermo, AtmData(iPress, j + 1),
              AtmData(iTemp, j + 1), mWt);
          VectorXd D_prev = handleCustomMolecularDiffusion(
              PlanetName, gasThermo, AtmData(iPress, j - 1),
              AtmData(iTemp, j - 1), mWt);
          D_next = (D_next + D_this) / 2;
          D_prev = (D_prev + D_this) / 2;
          // std::cout << AtmData(iAlt,j) << std::endl;
          // std::cout << D_prev.transpose() << std::endl;
          // std::cout << K_prev << std::endl;
          // std::cout << "Averaging diffusion coefficient" << std::endl;
          double T_next = (AtmData(iTemp, j) + AtmData(iTemp, j + 1)) / 2;
          double T_prev = (AtmData(iTemp, j) + AtmData(iTemp, j - 1)) / 2;
          double dz = (AtmData(iAlt, j - 1) - AtmData(iAlt, j + 1)) / 2;
          VectorXd I = VectorXd::Ones(nsp);
          // Dynamic Time Scales
          VectorXd Diff = AtmData(iKzz, j) * I + D_this;
          TDyn.col(j) = ((dz * dz) / Diff.array()).matrix();
          // std::cout << "Identity matrix" << std::endl;
          VectorXd d_next =
              ((mm_next * I - mWt).array() * D_next.array()).matrix() *
              (1 / (2 * dz * dz)) *
              (g * 1e-3 * dz / (6.022E23 * 1.38E-23 * T_next));
          VectorXd d_prev =
              ((mm_prev * I - mWt).array() * D_prev.array()).matrix() *
              (1 / (2 * dz * dz)) *
              (g * 1e-3 * dz / (6.022E23 * 1.38E-23 * T_prev));

          double T_this = AtmData(iTemp, j);
          VectorXd dTdz_prev =
              alpha * (T_this - T_prev) * 2 / ((T_this + T_prev) * dz);
          VectorXd dTdz_next =
              alpha * (T_next - T_this) * 2 / ((T_this + T_next) * dz);
          d_next = d_next + ((dTdz_next.array() * D_next.array()).matrix() *
                             (1 / (2 * dz * dz)));
          d_prev = d_prev + ((dTdz_next.array() * D_next.array()).matrix() *
                             (1 / (2 * dz * dz)));

          // std::cout << "d complete" << std::endl;
          VectorXd k_next = (((K_next * I) + D_next) / (dz * dz));
          VectorXd k_prev = (((K_prev * I) + D_prev) / (dz * dz));
          // std::cout << "k complete" << std::endl;
          double N_next = AtmData(iNd, j + 1);
          double N_this = AtmData(iNd, j);
          double N_prev = AtmData(iNd, j - 1);
          double N_n = (N_this + N_next) / 2;
          double N_p = (N_this + N_prev) / 2;
          // std::cout << "Averaging densities complete" << std::endl;

          // Contribution of U_i
          B = ((k_next * N_n / N_this) + (d_next * N_n / N_this) +
               (k_prev * N_p / N_this) - (d_prev * N_p / N_this));
          //  std::cout << "U_i complete" << std::endl;
          // Contribution of U_i-1
          A = -1 * ((k_prev * N_p / N_prev) + (d_prev * N_p / N_prev));
          //  std::cout << "U_i-1 complete" << std::endl;
          // Contribution of U_i+1
          C = -1 * ((k_next * N_n / N_next) - (d_next * N_n / N_next));
          //  std::cout << "U_i+1 complete" << std::endl;
          // Check for dn/dt
          conv.col(j) =
              m_wdot -
              (B.array() * N_this * ChemMoleFrac.col(j).array()).matrix() -
              (A.array() * N_prev * ChemMoleFrac.col(j - 1).array()).matrix() -
              (C.array() * N_prev * ChemMoleFrac.col(j + 1).array()).matrix();

          // Inserting the terms into block matrix

          // Main diagonal terms - Jacobian from chemistry
          BlockMatrix.block(j * nsp, j * nsp, nsp, nsp) =
              ((mat1) - (m_wjac * dt));
          Un = ((mat1) - (m_wjac * dt)) * N_this * mole_frac;
          //   std::cout << "main diagonal set" << std::endl;
          // Main diagonal - diffusion term
          for (int sp = 0; sp < nsp; sp++) {
            Upresent(j * nsp + sp) = (m_wdot(sp) * dt) + Un(sp);
            BlockMatrix(j * nsp + sp, j * nsp + sp) =
                BlockMatrix(j * nsp + sp, j * nsp + sp) + (B(sp) * dt);

            // Off diagonal terms - diffusion
            BlockMatrix((j)*nsp + sp, (j - 1) * nsp + sp) = A(sp) * dt;
            BlockMatrix((j)*nsp + sp, (j + 1) * nsp + sp) = C(sp) * dt;
          }
          //  std::cout << "Off diagonal terms set" << std::endl;
        }

        // Upper boundary condition
        if (j == 0) {
          mole_frac = ChemMoleFrac.col(j + 1);
          gas2->setMoleFractions(&mole_frac[0]);
          mm_next = gas2->meanMolecularWeight();  // Mean molcular weight

          if (top == "Dirichlet") {
            std::string l_species_list =
                pinput->GetOrAddString("upperboundaryMixRat", "species", "nan");
            if (l_species_list != "nan") {
              while (std::regex_search(l_species_list, m, pattern)) {
                for (auto x : m) {
                  std::string species_lowerboundary_MixRat =
                      pinput->GetString("upperboundaryMixRat", x);
                  species_inx = gas2->speciesIndex(x);
                  mole_frac(species_inx) =
                      atof(species_lowerboundary_MixRat.c_str());
                }
                l_species_list = m.suffix().str();
              }
            }

            gas2->setMoleFractions(&mole_frac[0]);
            mm_prev = gas2->meanMolecularWeight();  // Mean molcular weight
          }

          if (bot == "Neumann") {
            mole_frac = ChemMoleFrac.col(j);
            gas2->setMoleFractions(&mole_frac[0]);
            mm_prev = gas2->meanMolecularWeight();  // Mean molcular weight
          }

          double K_next = (AtmData(iKzz, j + 1) + AtmData(iKzz, j)) / 2;
          double K_prev = (AtmData(iKzz, j) + AtmData(iKzz, j)) / 2;
          VectorXd D_this = handleCustomMolecularDiffusion(
              PlanetName, gasThermo, AtmData(iPress, j), AtmData(iTemp, j),
              mWt);
          VectorXd D_next = handleCustomMolecularDiffusion(
              PlanetName, gasThermo, AtmData(iPress, j + 1),
              AtmData(iTemp, j + 1), mWt);
          VectorXd D_prev = handleCustomMolecularDiffusion(
              PlanetName, gasThermo, AtmData(iPress, j), AtmData(iTemp, j),
              mWt);
          //   std::cout << "Diffusion coefficient complete!" << std::endl;
          D_next = (D_next + D_this) / 2;
          D_prev = (D_prev + D_this) / 2;
          // std::cout << "Averaging diffusion coefficient" << std::endl;
          double T_next = (AtmData(iTemp, j) + AtmData(iTemp, j + 1)) / 2;
          double T_prev = (AtmData(iTemp, j) + AtmData(iTemp, j)) / 2;
          double dz = (AtmData(iAlt, j) - AtmData(iAlt, j + 1));
          VectorXd I = VectorXd::Ones(nsp);
          // Dynamic Time Scales
          VectorXd Diff = AtmData(iKzz, j) * I + D_this;
          TDyn.col(j) = ((dz * dz) / Diff.array()).matrix();
          // std::cout << "Identity matrix" << std::endl;
          VectorXd d_next =
              ((mm_next * I - mWt).array() * D_next.array()).matrix() *
              (1 / (2 * dz * dz)) *
              (g * 1e-3 * dz / (6.022E23 * 1.38E-23 * T_next));
          VectorXd d_prev =
              ((mm_prev * I - mWt).array() * D_prev.array()).matrix() *
              (1 / (2 * dz * dz)) *
              (g * 1e-3 * dz / (6.022E23 * 1.38E-23 * T_prev));

          double T_this = AtmData(iTemp, j);
          VectorXd dTdz_prev =
              alpha * (T_this - T_prev) * 2 / ((T_this + T_prev) * dz);
          VectorXd dTdz_next =
              alpha * (T_next - T_this) * 2 / ((T_this + T_next) * dz);
          d_next = d_next + ((dTdz_next.array() * D_next.array()).matrix() *
                             (1 / (2 * dz * dz)));
          d_prev = d_prev + ((dTdz_next.array() * D_next.array()).matrix() *
                             (1 / (2 * dz * dz)));

          // std::cout << "d complete" << std::endl;
          VectorXd k_next = (((K_next * I) + D_next) / (dz * dz));
          VectorXd k_prev = (((K_prev * I) + D_prev) / (dz * dz));
          // std::cout << "k complete" << std::endl;
          double N_next = AtmData(iNd, j + 1);
          double N_this = AtmData(iNd, j);
          double N_prev = AtmData(iNd, j);
          double N_n = (N_this + N_next) / 2;
          double N_p = (N_this + N_prev) / 2;

          // Contribution of U_i
          B = ((k_next * N_n / N_this) + (d_next * N_n / N_this) +
               (k_prev * N_p / N_this) - (d_prev * N_p / N_this));
          //  std::cout << "U_i complete" << std::endl;
          // Contribution of U_i-1
          A = -1 * ((k_prev * N_p / N_prev) + (d_prev * N_p / N_prev));
          //  std::cout << "U_i-1 complete" << std::endl;
          // Contribution of U_i+1
          C = -1 * ((k_next * N_n / N_next) - (d_next * N_n / N_next));
          //  std::cout << "U_i+1 complete" << std::endl;
          // Check for dn/dt
          conv.col(j) =
              m_wdot - (B.array() * N_this * mole_frac.array()).matrix() -
              (C.array() * N_prev * ChemMoleFrac.col(j + 1).array()).matrix() -
              (A.array() * N_prev * mole_frac.array()).matrix();

          // Inserting the terms into block matrix

          // Dirichlet boundary condition
          if (top == "Dirichlet") {
            // Main diagonal terms - Jacobian from chemistry
            BlockMatrix.block(j * nsp, j * nsp, nsp, nsp) =
                ((mat1) - (m_wjac * dt));
            Un = ((mat1) - (m_wjac * dt)) * N_this * mole_frac;
            // Main diagonal - diffusion term
            for (int sp = 0; sp < nsp; sp++) {
              // Dirichlet conditions only for species for whom boundary
              // conditions have been specified
              std::string speciesName = gas->speciesName(sp);
              std::string speciesD = pinput->GetOrAddString(
                  "upperboundaryMixRat", speciesName, "nan");
              if (speciesD != "nan") {
                Upresent(j * nsp + sp) = (m_wdot(sp) * dt) + Un(sp) -
                                         (A(sp) * dt * N_next * mole_frac(sp));
                BlockMatrix(j * nsp + sp, j * nsp + sp) =
                    BlockMatrix(j * nsp + sp, j * nsp + sp) + (B(sp) * dt);
              }
              // If boundary conditions are not specified, by default its a zero
              // flux boundary condition
              if (speciesD == "nan") {
                Upresent(j * nsp + sp) = (m_wdot(sp) * dt) + Un(sp);
                BlockMatrix(j * nsp + sp, j * nsp + sp) =
                    BlockMatrix(j * nsp + sp, j * nsp + sp) +
                    ((A(sp) + B(sp)) * dt);
              }

              // Off diagonal terms - diffusion
              BlockMatrix(j * nsp + sp, (j + 1) * nsp + sp) = C(sp) * dt;
            }
          }

          // Neumann boundary condition
          if (top == "Neumann") {
            // Main diagonal terms - Jacobian from chemistry
            BlockMatrix.block(j * nsp, j * nsp, nsp, nsp) =
                ((mat1) - (m_wjac * dt));
            Un = ((mat1) - (m_wjac * dt)) * N_this * mole_frac;

            // Main diagonal - diffusion term
            for (int sp = 0; sp < nsp; sp++) {
              Upresent(j * nsp + sp) = (m_wdot(sp) * dt) + Un(sp);
              BlockMatrix(j * nsp + sp, j * nsp + sp) =
                  BlockMatrix(j * nsp + sp, j * nsp + sp) +
                  ((A(sp) + B(sp)) * dt);

              // Off diagonal terms - diffusion
              BlockMatrix(j * nsp + sp, (j + 1) * nsp + sp) = C(sp) * dt;
            }
          }
        }

        // Lower boundary condition
        if (j == nSize - 1) {
          mole_frac = ChemMoleFrac.col(j - 1);
          gas2->setMoleFractions(&mole_frac[0]);
          mm_prev = gas2->meanMolecularWeight();  // Mean molcular weight

          if (bot == "Dirichlet") {
            std::string l_species_list =
                pinput->GetOrAddString("lowerboundaryMixRat", "species", "nan");
            if (l_species_list != "nan") {
              while (std::regex_search(l_species_list, m, pattern)) {
                for (auto x : m) {
                  std::string species_lowerboundary_MixRat =
                      pinput->GetString("lowerboundaryMixRat", x);
                  species_inx = gas2->speciesIndex(x);
                  mole_frac(species_inx) =
                      atof(species_lowerboundary_MixRat.c_str());
                }
                l_species_list = m.suffix().str();
              }
            }

            gas2->setMoleFractions(&mole_frac[0]);
            mm_next = gas2->meanMolecularWeight();  // Mean molcular weight
          }

          if (bot == "Neumann") {
            mole_frac = ChemMoleFrac.col(j);
            gas2->setMoleFractions(&mole_frac[0]);
            mm_next = gas2->meanMolecularWeight();  // Mean molcular weight
          }

          // Setting the values into mass matrix
          double K_next = (AtmData(iKzz, j) + AtmData(iKzz, j)) / 2;
          double K_prev = (AtmData(iKzz, j - 1) + AtmData(iKzz, j)) / 2;
          VectorXd D_this = handleCustomMolecularDiffusion(
              PlanetName, gasThermo, AtmData(iPress, j), AtmData(iTemp, j),
              mWt);
          VectorXd D_next = handleCustomMolecularDiffusion(
              PlanetName, gasThermo, AtmData(iPress, j), AtmData(iTemp, j),
              mWt);
          VectorXd D_prev = handleCustomMolecularDiffusion(
              PlanetName, gasThermo, AtmData(iPress, j - 1),
              AtmData(iTemp, j - 1), mWt);
          // std::cout << "Diffusion coefficient complete!" << std::endl;
          D_next = (D_next + D_this) / 2;
          D_prev = (D_prev + D_this) / 2;
          // std::cout << "Averaging diffusion coefficient" << std::endl;
          double T_next = (AtmData(iTemp, j) + AtmData(iTemp, j)) / 2;
          double T_prev = (AtmData(iTemp, j) + AtmData(iTemp, j - 1)) / 2;
          double dz = (AtmData(iAlt, j - 1) - AtmData(iAlt, j));
          VectorXd I = VectorXd::Ones(nsp);
          // Dynamic Time Scales
          VectorXd Diff = AtmData(iKzz, j) * I + D_this;
          TDyn.col(j) = ((dz * dz) / Diff.array()).matrix();
          // std::cout << "Identity matrix" << std::endl;
          VectorXd d_next =
              ((mm_next * I - mWt).array() * D_next.array()).matrix() *
              (1 / (2 * dz * dz)) *
              (g * 1e-3 * dz / (6.022E23 * 1.38E-23 * T_next));
          VectorXd d_prev =
              ((mm_prev * I - mWt).array() * D_prev.array()).matrix() *
              (1 / (2 * dz * dz)) *
              (g * 1e-3 * dz / (6.022E23 * 1.38E-23 * T_prev));
          // std::cout << (d_next - d_prev).transpose() << std::endl;

          double T_this = AtmData(iTemp, j);
          VectorXd dTdz_prev =
              alpha * (T_this - T_prev) * 2 / ((T_this + T_prev) * dz);
          VectorXd dTdz_next =
              alpha * (T_next - T_this) * 2 / ((T_this + T_next) * dz);
          d_next = d_next + ((dTdz_next.array() * D_next.array()).matrix() *
                             (1 / (2 * dz * dz)));
          d_prev = d_prev + ((dTdz_next.array() * D_next.array()).matrix() *
                             (1 / (2 * dz * dz)));

          // std::cout << "d complete" << std::endl;
          VectorXd k_next = (((K_next * I) + D_next) / (dz * dz));
          VectorXd k_prev = (((K_prev * I) + D_prev) / (dz * dz));
          // std::cout << "k complete" << std::endl;
          double N_next = AtmData(iNd, j);
          double N_this = AtmData(iNd, j);
          double N_prev = AtmData(iNd, j - 1);
          double N_n = (N_this + N_next) / 2;
          double N_p = (N_this + N_prev) / 2;
          // std::cout << "Averaging densities complete" << std::endl;
          // Contribution of U_i
          B = ((k_next * N_n / N_this) + (d_next * N_n / N_this) +
               (k_prev * N_p / N_this) - (d_prev * N_p / N_this));
          //  std::cout << "U_i complete" << std::endl;
          // Contribution of U_i-1
          A = -1 * ((k_prev * N_p / N_prev) + (d_prev * N_p / N_prev));
          //  std::cout << "U_i-1 complete" << std::endl;
          // Contribution of U_i+1
          C = -1 * ((k_next * N_n / N_next) - (d_next * N_n / N_next));
          //  std::cout << "U_i+1 complete" << std::endl;

          // Check for dn/dt
          conv.col(j) =
              m_wdot - (B.array() * N_this * mole_frac.array()).matrix() -
              (A.array() * N_prev * ChemMoleFrac.col(j - 1).array()).matrix() -
              (C.array() * N_prev * mole_frac.array()).matrix();

          // Inserting the terms into block matrix
          if (bot == "Dirichlet") {
            // Main diagonal terms - Jacobian from chemistry
            BlockMatrix.block(j * nsp, j * nsp, nsp, nsp) =
                ((mat1) - (m_wjac * dt));
            Un = ((mat1) - (m_wjac * dt)) * N_this * mole_frac;
            // Main diagonal - diffusion term
            for (int sp = 0; sp < nsp; sp++) {
              // Dirichlet conditions only for species for whom boundary
              // conditions have been specified
              std::string speciesName = gas->speciesName(sp);
              std::string speciesD = pinput->GetOrAddString(
                  "lowerboundaryMixRat", speciesName, "nan");

              if (speciesD != "nan") {
                Upresent(j * nsp + sp) = (m_wdot(sp) * dt) + Un(sp) -
                                         (C(sp) * dt * N_next * mole_frac(sp));
                BlockMatrix(j * nsp + sp, j * nsp + sp) =
                    BlockMatrix(j * nsp + sp, j * nsp + sp) + (B(sp) * dt);
              }

              if (speciesD == "nan") {
                Upresent(j * nsp + sp) = (m_wdot(sp) * dt) + Un(sp);
                BlockMatrix(j * nsp + sp, j * nsp + sp) =
                    BlockMatrix(j * nsp + sp, j * nsp + sp) +
                    ((C(sp) + B(sp)) * dt);
              }

              // Off diagonal terms - diffusion
              BlockMatrix(j * nsp + sp, (j - 1) * nsp + sp) = A(sp) * dt;
            }
          }

          if (bot == "Neumann") {
            // Main diagonal terms - Jacobian from chemistry
            BlockMatrix.block(j * nsp, j * nsp, nsp, nsp) =
                ((mat1) - (m_wjac * dt));
            Un = ((mat1) - (m_wjac * dt)) * N_this * mole_frac;
            // Main diagonal - diffusion term
            for (int sp = 0; sp < nsp; sp++) {
              Upresent(j * nsp + sp) = (m_wdot(sp) * dt) + Un(sp);
              BlockMatrix(j * nsp + sp, j * nsp + sp) =
                  BlockMatrix(j * nsp + sp, j * nsp + sp) +
                  ((C(sp) + B(sp)) * dt);

              // Off diagonal terms - diffusion
              BlockMatrix(j * nsp + sp, (j - 1) * nsp + sp) = A(sp) * dt;
            }
          }
        }
      }

      // Inverting the matrix

      Ufuture = BlockMatrix.inverse() * Upresent;

      // Checking for convergence

      double conv_factor = conv.maxCoeff();
      double tau = 1.1;
      std::cout << "convergence: " << conv_factor << std::endl;
      std::cout << "Chemical time scale: " << TChem.minCoeff() << std::endl;
      std::cout << "Dynamic time scale: " << TDyn.minCoeff() << std::endl;
      // if((conv_factor >= 1E-9) && (dt > TDyn.minCoeff())){
      //  dt = 0.1*TDyn.minCoeff();
      //  counter++;
      // }
      // if((conv_factor < 1E-9)){
      // if(counter < 1){
      //  dt = dt*2;}
      //  }
      if (dt < 0.01 * TDyn.minCoeff()) {
        dt = dt * 2;
      }
      // Updating the solution
      for (int j = 0; j < nSize; j++) {
        for (int sp = 0; sp < nsp; sp++) {
          ChemMoleFrac(sp, j) = Ufuture(j * nsp + sp) / AtmData(iNd, j);
          ;
          ChemConc(sp, j) = Ufuture(j * nsp + sp);
        }
      }
      Ttot = Ttot + dt;
      std::cout << "Simulation completed at time step t = " << Ttot
                << std::endl;

      /*
        if(conv_factor >= tau){
      //Reject the solution, and reduce the step size
        dt = 0.9*dt*pow(tau/conv_factor, 0.5);
        std::cout << "Time step: " << dt << std::endl;
        }

        if(conv_factor < tau){
      //Updating the solution
        for (int j = 0; j < nSize; j++) {
           for(int sp = 0; sp < nsp; sp++){
               ChemMoleFrac(sp, j) = Ufuture(j*nsp + sp)/AtmData(iNd, j);;
               ChemConc(sp, j) = Ufuture(j*nsp + sp);
           } }
        std::cout << "Simulation completed at time step t = " << Ttot <<
      std::endl; dt = dt*1.25; Ttot = Ttot + dt;
      }
      */
    }

    // This is where the photochemistry definition ends (Anything beyond is not
    // defined in the scope of photochemistry)
  }

  // Writing output into NetCDF file
  int num3D = ChemMole.size();

  std::vector<Eigen::MatrixXd> matrix3d(Time.size(),
                                        Eigen::MatrixXd(nsp, nSize));

  VectorXd cPress = AtmData.row(iPress);
  VectorXd cTemp = AtmData.row(iTemp);
  VectorXd cHeight = AtmData.row(iAlt);
  VectorXd cKzz = AtmData.row(iKzz);
  MatrixXd dPRates = ProdRates;
  MatrixXd dDRates = DiffRates;
#if NETCDFOUTPUT
  int ifile;
  string fname = pinput->GetString("output", "out_file");
  int status = nc_create(fname.c_str(), NC_NETCDF4, &ifile);
  int dimids[3];
  int alt, iPres, iTem, iKeddy, iAltz, iChem, isp, iProd, iDiff, iTs;

  // Atmospheric Properties (All in SI units!)
  nc_def_dim(ifile, "Time", Time.size(), &dimids[0]);

  nc_def_dim(ifile, "Pressure", nSize, &dimids[2]);

  nc_def_dim(ifile, "Species", nsp, &dimids[1]);

  nc_def_var(ifile, "Pressure", NC_DOUBLE, 1, &dimids[2], &iPres);
  nc_put_var_double(ifile, iPres, &cPress[0]);

  nc_def_var(ifile, "Keddy", NC_DOUBLE, 1, &dimids[2], &iKeddy);
  nc_put_var_double(ifile, iKeddy, &cKzz[0]);

  nc_def_var(ifile, "Temp", NC_DOUBLE, 1, &dimids[2], &iTem);
  nc_put_var_double(ifile, iTem, &cTemp[0]);

  nc_def_var(ifile, "Altitude", NC_DOUBLE, 1, &dimids[2], &iAltz);
  nc_put_var_double(ifile, iAltz, &cHeight[0]);

  // nc_def_var(ifile, "Series", NC_DOUBLE, 3, dimids, &iTs);
  // nc_put_var_double(ifile, iTs, matrix3d.data());

  init_species_list = pinput->GetString("output", "species");
  while (std::regex_search(init_species_list, m, pattern)) {
    for (auto x : m) {
      std::string species = x;
      std::string prod = "Prod";
      std::string diff = "Diff";
      std::string prod_species = prod + species;
      std::string diff_species = diff + species;
      species_inx = gas->speciesIndex(x);
      VectorXd cChem =
          ChemConc.row(species_inx) * 6.022E23 * 1E3 * 1E-6;  // #/cm^-3
      VectorXd pChem = ProdRates.row(species_inx);  // kmol/m^3s (Cantera units)
      VectorXd dChem = DiffRates.row(species_inx);  // kmol/m^3s (Cantera units)
      const char* ccx = &species[0];
      const char* prodx = &prod_species[0];
      const char* diffx = &diff_species[0];
      nc_def_var(ifile, ccx, NC_DOUBLE, 1, &dimids[2], &iChem);
      nc_put_var_double(ifile, iChem, &cChem[0]);

      nc_def_var(ifile, prodx, NC_DOUBLE, 1, &dimids[2], &iProd);
      nc_put_var_double(ifile, iProd, &pChem[0]);

      nc_def_var(ifile, diffx, NC_DOUBLE, 1, &dimids[2], &iDiff);
      nc_put_var_double(ifile, iDiff, &dChem[0]);
    }
    init_species_list = m.suffix().str();
  }

  nc_close(ifile);
#endif

  delete pinput;
}
