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
#include <cantera/transport/Transport.h>

// output stream
#include <iostream>
#include <fstream>
#include <cantera/numerics/eigen_dense.h>
#include <cantera/numerics/eigen_sparse.h>
#include <vector>
#include <string>
#include <sstream>
#include <regex>

//OpenMP
#include <omp.h>

// Athena++ header
#include <parameter_input.hpp>

// C3M header
#include <PhotoChemistry.hpp>
#include <RadTran.hpp>
#include <interpolation.hpp>
#include <CustomRate.hpp>
#include <CustomTransport.hpp>
// NetCDF Output
#if NETCDFOUTPUT
  #include <netcdf.h>
#endif

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace Cantera;
using namespace std;

//Function for Eddy diffusion [Hu et al., 2012]
MatrixXd Kflux(int inxn, int inxNext, int inxPrev, MatrixXd Sp, MatrixXd Np, int ns, double Kzz, double Kzz_prev, double Kzz_next, double dx1, double dx2){
  VectorXd diff_flux(ns), Xn(ns), Xnext(ns), Xprev(ns);
  double Nn, Nnext, Nprev; 
  double beta1, beta2;
  Xn = Sp.col(inxn); //mole fraction
  Xnext = Sp.col(inxNext); //mole fraction
  Xprev = Sp.col(inxPrev); //mole fraction
  Nn = Np.colwise().sum()(inxn); //concentration (kmol/m^3)
  Nnext = Np.colwise().sum()(inxNext); //concentration (kmol/m^3)
  Nprev = Np.colwise().sum()(inxPrev); //concentration (kmol/m^3)
  beta1 = -1*((Kzz + Kzz_prev)/2)*((Nn + Nprev)/2); //j-1/2
  beta2 = -1*((Kzz + Kzz_next)/2)*((Nn + Nnext)/2); //j+1/2
  double dx = (dx1 + dx2)/2;
  diff_flux =   -1*(((((Xnext - Xn)*beta2)/dx2) - (((Xn - Xprev))*beta1/dx1))/dx);
return diff_flux;
}



//Function for Molecular, and thermodiffusion [Hu et al., 2012]
MatrixXd Dflux(int inxn, int inxNext, int inxPrev, MatrixXd Sp, MatrixXd Np, int ns, VectorXd D_this, VectorXd D_prev, VectorXd D_next, double dx1, double dx2, MatrixXd mW, double T_prev, double T_this, double T_next, double mWav, double ga){
  VectorXd diff_flux(ns), Xn(ns), Xnext(ns), Xprev(ns);
  double Nn, Nnext, Nprev; 
  VectorXd beta1(ns), beta2(ns), Hs_prev(ns), Hs_next(ns);
  VectorXd phi_prev(ns), phi_next(ns);
  double Hs0_prev = 2*mWav*ga*1e-3/(6.022e23*1.38e-23*(T_prev + T_this)); //scale height of atmosphere -> 1/m
  double Hs0_next = 2*mWav*ga*1e-3/(6.022e23*1.38e-23*(T_next + T_this)); //scale height of atmosphere -> 1/m
  double dTdz1 = 2*(T_this - T_prev)/(dx1*(T_prev + T_this)); //Thermodiffusion
  double dTdz2 = 2*(T_next - T_this)/(dx2*(T_next + T_this)); //Thermodiffusion
  double dx = (dx1 + dx2)/2;
  VectorXd I = VectorXd::Ones(ns);
  Hs_prev = 2*mW*ga*1e-3/(6.022e23*1.38e-23*(T_prev + T_this)); //scale height of species -> 1/m
  Hs_next = 2*mW*ga*1e-3/(6.022e23*1.38e-23*(T_next + T_this)); //scale height of species -> 1/m
  Xn = Sp.col(inxn); //mole fraction
  Xnext = Sp.col(inxNext); //mole fraction
  Xprev = Sp.col(inxPrev); //mole fraction
  Nn = Np.colwise().sum()(inxn); //concentration (kmol/m^3)
  Nnext = Np.colwise().sum()(inxNext); //concentration (kmol/m^3)
  Nprev = Np.colwise().sum()(inxPrev); //concentration (kmol/m^3)
  beta1 = -1*((D_this + D_prev)/2)*((Nn + Nprev)/2); //j-1/2
  beta2 = -1*((D_this + D_next)/2)*((Nn + Nnext)/2); //j+1/2
  phi_prev = -0.5*(((I*Hs0_prev) - Hs_prev - (I*dTdz1)).array()*beta1.array()).matrix();
  phi_next = -0.5*(((I*Hs0_next) - Hs_next - (I*dTdz2)).array()*beta2.array()).matrix();
  diff_flux =   -1*(((((Xnext - Xn).array()*beta2.array()).matrix()/dx2) - (((Xn - Xprev).array()*beta1.array()).matrix()/dx1))/dx); 
  diff_flux = diff_flux - ((phi_next - phi_prev)/dx);
return diff_flux;
}




//Function for ambipolar diffusion
MatrixXd Aflux(int inxn, int inxNext, int inxPrev, MatrixXd Sp, MatrixXd Np, int ns, VectorXd D_this, VectorXd D_prev, VectorXd D_next, double dx1, double dx2, MatrixXd mW, double T_prev, double T_this, double T_next, double mWav, double ga, int iE, VectorXd alph){
  VectorXd diff_flux(ns), Xn(ns), Xnext(ns), Xprev(ns);
  double Nn, Nnext, Nprev; 
  VectorXd beta1(ns), beta2(ns), Hs_prev(ns), Hs_next(ns);
  VectorXd phi_prev(ns), phi_next(ns);
  double Ne_prev, Ne_this, Ne_next; //Electron concentrtaion (kmol/m^3)
//Mean scale height
  double Hs0_prev = 2*mWav*ga*1e-3/(6.022e23*1.38e-23*(T_prev + T_this)); //scale height of atmosphere -> 1/m
  double Hs0_next = 2*mWav*ga*1e-3/(6.022e23*1.38e-23*(T_next + T_this)); //scale height of atmosphere -> 1/m
//Lapse rate for thermodiffusion  
  double dTdz1 = 2*(T_this - T_prev)/(dx1*(T_prev + T_this)); //Thermodiffusion
  double dTdz2 = 2*(T_next - T_this)/(dx2*(T_next + T_this)); //Thermodiffusion
  double dx = (dx1 + dx2)/2;
  VectorXd I = VectorXd::Ones(ns);
//Scale unit for all the species
  Hs_prev = 2*mW*ga*1e-3/(6.022e23*1.38e-23*(T_prev + T_this)); //scale height of species -> 1/m
  Hs_next = 2*mW*ga*1e-3/(6.022e23*1.38e-23*(T_next + T_this)); //scale height of species -> 1/m
//Mole fraction array
  Xn = Sp.col(inxn); //mole fraction
  Xnext = Sp.col(inxNext); //mole fraction
  Xprev = Sp.col(inxPrev); //mole fraction
//Electron density
  Ne_prev = Np.col(inxPrev)(iE);
  Ne_this = Np.col(inxn)(iE);
  Ne_next = Np.col(inxNext)(iE);
//Ambipolar diffusion term
  double dNedz1 = 2*((Ne_this*T_this) - (Ne_prev*T_prev))/(dx1*((Ne_prev*T_prev) + (Ne_this*T_this)));
  double dNedz2 = 2*((Ne_next*Ne_next) - (Ne_this*T_this))/(dx2*((Ne_next*T_next) + (Ne_this*T_this)));
  if (std::isnan(std::abs(dNedz1))) {
    dNedz1 = 0.0;
  }
  if (std::isnan(std::abs(dNedz2))) {
    dNedz2 = 0.0;
  }
//Total number density
  Nn = Np.colwise().sum()(inxn); //concentration (kmol/m^3)
  Nnext = Np.colwise().sum()(inxNext); //concentration (kmol/m^3)
  Nprev = Np.colwise().sum()(inxPrev); //concentration (kmol/m^3)
  beta1 = -1*((D_this + D_prev)/2)*((Nn + Nprev)/2); //j-1/2
  beta2 = -1*((D_this + D_next)/2)*((Nn + Nnext)/2); //j+1/2
//Total flux due to diffusion of ions
  phi_prev = -0.5*(((I*Hs0_prev) - Hs_prev - (I*dNedz1) - ((I + alph)*dTdz1)).array()*beta1.array()).matrix();
  phi_next = -0.5*(((I*Hs0_next) - Hs_next - (I*dNedz2) - ((I + alph)*dTdz2)).array()*beta2.array()).matrix();
  diff_flux =   -1*(((((Xnext - Xn).array()*beta2.array()).matrix()/dx2) - (((Xn - Xprev).array()*beta1.array()).matrix()/dx1))/dx); 
  diff_flux = diff_flux - ((phi_next - phi_prev)/dx);
  //std::cout << diff_flux.transpose() << std::endl;
return diff_flux;
}


int main(int argc, char **argv) {
  
//Reading input file
  IOWrapper infile;
  infile.Open("IonRatesJEDIonly.inp", IOWrapper::FileMode::read);
  ParameterInput *pinput = new ParameterInput();

  pinput->LoadFromFile(infile);
  infile.Close();

//Loading the input parameters for atmospheric profile and reaction network
  std::string atm_file = pinput->GetString("problem", "planet");
  std::string PlanetName = pinput->GetString("problem", "planetName");
  std::string network_file = pinput->GetString("problem", "network");
  std::string profile_file = pinput->GetOrAddString("problem", "chemInp", "nan");
  std::string ionkinetics_file = pinput->GetString("problem", "ionkinetics");
  
//Reading the chemical kinetics network
  auto sol = newSolution(network_file);
  auto gas = sol->thermo();
  auto gas_kin = sol->kinetics(); 
  auto gas_tr  = sol->transport(); 
  int nsp = gas->nSpecies();
  int nrxn = gas_kin->nReactions();
  std::cout << "Number of reactions: " << nrxn << std::endl;
  std::cout << "Network imported" << std::endl;
//Initial condition for mole fraction
  VectorXd mole_fractions = VectorXd::Zero(nsp);
  
//Initial condition for boundary fluxes
  VectorXd flux_lower = VectorXd::Zero(nsp);
  VectorXd flux_upper = VectorXd::Zero(nsp);

//Loading the input parameters for initial condition
//Homogeneous input condition
  std::string init_species_list = pinput->GetString("init", "species");
  std::regex pattern ("[a-zA-Z0-9_+-]+");
  std::smatch m;
  int species_inx;
  while (std::regex_search (init_species_list,m,pattern)) {
    for (auto x:m){
    std::string species_init_condition = pinput->GetString("init", x);
    species_inx = gas->speciesIndex(x);
    mole_fractions(species_inx) = atof(species_init_condition.c_str());
    }
    init_species_list = m.suffix().str();
  }


//Loading the input for upper boundary condition
  std::string ub_species_list = pinput->GetString("upperboundaryflux", "species");
  while (std::regex_search (ub_species_list,m,pattern)) {
    for (auto x:m){
    std::string species_upperboundary_condition = pinput->GetString("upperboundaryflux", x);
    species_inx = gas->speciesIndex(x);
    flux_upper(species_inx) = atof(species_upperboundary_condition.c_str());
    }
    ub_species_list = m.suffix().str();
  }

//Loading the input for lower boundary condition 
  std::string lb_species_list = pinput->GetString("lowerboundaryflux", "species");
  while (std::regex_search (lb_species_list,m,pattern)) {
    for (auto x:m){
    std::string species_lowerboundary_condition = pinput->GetString("lowerboundaryflux", x);
    species_inx = gas->speciesIndex(x);
    flux_lower(species_inx) = atof(species_lowerboundary_condition.c_str());
    }
    lb_species_list = m.suffix().str();
  }
std::cout << "Boundary conditions loaded" << std::endl;
//Integrator
  std::string time_step = pinput->GetString("integrator", "dt");
  std::string total_time_steps = pinput->GetString("integrator", "nSteps");
  std::string max_time = pinput->GetString("integrator", "Tmax");
  std::string start_time = pinput->GetString("integrator", "Tstart");
  std::string precp_time = pinput->GetString("integrator", "Tprecp");
  double dt = atof(time_step.c_str());
  double Tmax = atof(max_time.c_str());
  double Tmin = dt;
  double Ttot = 0.0;
  double Tstart = atof(start_time.c_str());
  double Tprecp = atof(precp_time.c_str());
  std::cout << "Time step: " << dt << std::endl;
  int nTime = stoi(total_time_steps.c_str());


//Atmosphere Properties, indices for I/O storage
  fstream InFile;
  int nSize = 0;
  string data1, data2, data3, data4;
  InFile.open(atm_file); 
  getline(InFile, data1);
  getline(InFile, data1);
  while (getline(InFile, data1))
    nSize++;
  InFile.close();
  

//Atmospheric Profile Data
  MatrixXd AtmData(5, nSize);
  int iTemp = 0;
  int iPress =  1;
  int iKzz =  2;
  int iAlt = 3;
  int iMlDf= 4;

//Input from txt file
  int inx = 0;
  InFile.open(atm_file); 
  getline(InFile, data1);
  getline(InFile, data1);
  double kmax;
  while (InFile >> data1 >> data2 >> data3 >> data4){
      AtmData(iPress, inx) = atof(data1.c_str())*1E2; //Pressure (mbar) -> Pa
      AtmData(iTemp,inx) = atof(data2.c_str()); //Temperature (K)
      AtmData(iKzz, inx) = atof(data3.c_str())*1E-4; //Kzz (cm^2/s) -> m^2/s
      AtmData(iAlt, inx) = atof(data4.c_str())*1E3; //Altitude (m) 
      AtmData(iMlDf, inx) = 0.0;
      if(inx == 0)
      {kmax = AtmData(iKzz, inx); }
      if(inx != 0)
      {if(kmax < AtmData(iKzz, inx))
      {
         kmax = AtmData(iKzz, inx);
      }
      
      }

      inx++;
      }
  InFile.close();

//Setting initial condition for chemical species
  MatrixXd ChemMoleFrac(nsp, nSize);
  MatrixXd ChemConc(nsp, nSize);
  VectorXd MolDCoeff(nsp);
  MatrixXd a(nsp, nSize);
  MatrixXd n_conc(nsp, nSize);
  MatrixXd conv(nsp, nSize);
  MatrixXd ProdRates(nsp, nSize);
  MatrixXd DiffRates(nsp, nSize);
  
  for (int i = 0; i < nSize; i++) {
      gas->setState_TP(AtmData(iTemp,i), (AtmData(iPress,i)/1.0132E5)*OneAtm);
      ChemMoleFrac.col(i) = mole_fractions;
    }

//Chemical species profiles from input file
  if(profile_file != "nan"){
  std::cout << "Setting initial condition from input profile" << std::endl;
  std::string input_species_list = pinput->GetString("profile", "species");
  std::regex pattern ("[a-zA-Z0-9_+-]+");
  int chemP = 0;
  while (std::regex_search (input_species_list,m,pattern)) {
    for (auto x:m){
    chemP++;
    }
    input_species_list = m.suffix().str();
    }
  InFile.open(profile_file);
  std::string inp1, inp2;
  int i = 0;
  getline(InFile, inp1);
  while(i < nSize){
    int chem_inx = 0;
    gas->setState_TP(AtmData(iTemp,i), (AtmData(iPress,i)/1.0132E5)*OneAtm);
    input_species_list = pinput->GetString("profile", "species");
    while (std::regex_search (input_species_list,m,pattern)) {
    for (auto x:m){
     species_inx = gas->speciesIndex(x);
     if(chem_inx < chemP-1){
       getline(InFile,inp2,',');
       ChemMoleFrac(species_inx,i) = atof(inp2.c_str());
       }

     if(chem_inx == chemP-1){
       getline(InFile,inp2, '\n');
       ChemMoleFrac(species_inx,i) = atof(inp2.c_str());
       }
   
     chem_inx++;

    }
    
    input_species_list = m.suffix().str();
    }
    gas->setMoleFractions(&ChemMoleFrac.col(i)[0]);
    //std::cout << ChemMoleFrac.col(i).transpose() << std::endl;
    // std::cout << i << " " << nSize << std::endl;
     gas->getConcentrations(&ChemConc.col(i)[0]);
   //std::cout << "Got the concentrations" << std::endl;

    i++;
  }
  
  }
 InFile.close();

std::cout << "Extracting Ion kinetics reaction rates" << std::endl;

//Extracting the reaction number corresponding to electron impact process
  std::string ionkin_list = pinput->GetString("ionkinetics", "rxn");
  int ion_rxn = 0;
  while (std::regex_search (ionkin_list,m,pattern)) {
    for (auto x:m){
    ion_rxn++;
    }
    ionkin_list = m.suffix().str();
    }
  
std::cout << "Number of reactions" << ion_rxn << std::endl;  
  MatrixXd ion_rate(ion_rxn, nSize);
  VectorXd Irate(ion_rxn);
  VectorXd ion_rxnid(ion_rxn);
  ionkin_list = pinput->GetString("ionkinetics", "rxn");


  ion_rxn = 0;

  while (std::regex_search (ionkin_list,m,pattern)) {
    for (auto x:m){
    std::string xtemp = x;
    std::cout << "Reaction number: " <<x << std::endl;
    ion_rxnid(ion_rxn) = atof(xtemp.c_str());
    ion_rxn++;
    }
    ionkin_list = m.suffix().str();
    }

//Reading Input file
  InFile.open(ionkinetics_file);
  std::string inp3, inp4;
  int ih = 0;
  getline(InFile, inp3);
    while(ih < nSize){
    for(int ixn = 0; ixn < ion_rxn; ixn++){
      if(ixn < ion_rxn-1){
       getline(InFile,inp4,',');
      // std::cout << inp4 << std::endl;
       ion_rate(ixn,ih) = atof(inp4.c_str());
       }

     if(ixn == ion_rxn-1){
       getline(InFile,inp4, '\n');
      // std::cout << inp4 << std::endl;
       ion_rate(ixn,ih) = atof(inp4.c_str());} 
    }
    //std::cout << AtmData(iAlt, ih)/1E3  << " " << ion_rate.col(ih).transpose() << std::endl;
    ih++;
   }
  InFile.close();


  std::cout << "All inputs loaded into C3M " << std::endl;
  
  //Reading the Stellar Irradiance Input
  std::string stellar_input_file = pinput->GetString("radtran", "solar");
  std::string radius = pinput->GetString("radtran", "radius");
  std::string reference =   pinput->GetString("radtran", "reference");
  std::string SZA =   pinput->GetString("radtran", "SZA");
  double rad = atof(radius.c_str());
  double ref = atof(reference.c_str());
  double sz_angle = atof(SZA.c_str());
  MatrixXd stellar_input = ReadStellarRadiationInput(stellar_input_file, rad, ref);
  MatrixXd d_wavelength(stellar_input.row(0).size(), 1);
  MatrixXd wavelength = stellar_input.row(0);

  MatrixXd OpacityStorage = MatrixXd::Zero(stellar_input.row(0).size(),nSize);
//Condition to turn on photochemistry
  std::string photochem = pinput->GetString("problem", "photochem");
  if(photochem == "true")
  {
  
  std::cout << "Starting photochemistry calculations" << std::endl;

//Modification of solar flux based on SZA input, along with conversion from deg. to radians
  stellar_input.row(1) = stellar_input.row(1)*cos(sz_angle*3.14/180);
  
  for(int dw = 0;  dw < (stellar_input.row(0).size()) - 1; dw++){
  d_wavelength(dw) = wavelength(dw+1) - wavelength(dw);
  
  }

  std::cout << "!! Radiation Input loaded in C3M !!" << std::endl;
//Input Wavelength and Cross Section Data for all the absorbers 
  int AbsorberSize = 0; 
  std::string absorber_species_list = pinput->GetString("abscross", "absorbers");
  while (std::regex_search (absorber_species_list,m,pattern)) {
    for (auto x:m){
     AbsorberSize++;
    }
    absorber_species_list = m.suffix().str();
  }
  Eigen::MatrixXd absorber_cross_data(stellar_input.row(0).size(), AbsorberSize); 
  std::cout << "!! Initiating absorbers  !!" << std::endl;
//Reading absorber cross section for all absorbers
  int ab_inx = 0;
  absorber_species_list = pinput->GetString("abscross", "absorbers");
  while (std::regex_search (absorber_species_list,m,pattern)) {
    for (auto x:m){
    std::string absorber_cross_info = pinput->GetString("abscross", x);
    std::string absorber_cross_database = absorber_cross_info.substr(0,absorber_cross_info.find(","));  
    std::string absorber_cross_file = absorber_cross_info.substr(absorber_cross_info.find(",")+1,absorber_cross_info.length()-1);
    
//Reading the absorption cross section for absorber as per the database structure
//Interpolating the cross sections to reference grid
   if(absorber_cross_database == "VULCAN"){
     MatrixXd cross_info = ReadVULCANPhotoAbsCrossSection(absorber_cross_file);
     MatrixXd cross_section_data = InterpolateCrossSection(stellar_input.row(0), cross_info.row(0), cross_info.row(1));
     absorber_cross_data.col(ab_inx) = cross_section_data;
   }
    
   if(absorber_cross_database == "VPL"){
     MatrixXd cross_info = ReadAtmosCrossSection(absorber_cross_file);
     MatrixXd cross_section_data = InterpolateCrossSection(stellar_input.row(0), cross_info.row(0), cross_info.row(1));
     absorber_cross_data.col(ab_inx) = cross_section_data;
   } 

    //species_inx = gas->speciesIndex(x);
    ab_inx++;
    }
    absorber_species_list = m.suffix().str();
  }
  std::cout << "!! Obtaining Photochemical Cross Sections !!" << std::endl;

//Input Wavelength and Cross Section Data for all the absorbers 
  int ScatSize = 0;
  std::string scatter_species_list = pinput->GetString("scatcross", "scatterers");
  while (std::regex_search (scatter_species_list,m,pattern)) {
    for (auto x:m){
     ScatSize++;
    }
    scatter_species_list = m.suffix().str();
  }
  Eigen::MatrixXd scat_cross_data(stellar_input.row(0).size(), ScatSize);
  std::cout << "!! Initiating absorbers  !!" << std::endl;

//Reading the scattering cross section for atmospheric constituents as per database
//structure
  int sc_inx = 0;
  scatter_species_list = pinput->GetString("scatcross", "scatterers");
  while (std::regex_search (scatter_species_list,m,pattern)) {
    for (auto x:m){
    std::string scat_cross_info = pinput->GetString("scatcross", x);
    std::string scat_cross_database = scat_cross_info.substr(0,scat_cross_info.find(","));
    std::string scat_cross_file = scat_cross_info.substr(scat_cross_info.find(",")+1,scat_cross_info.length()-1);
    
//Reading the absorption cross section for absorber as per the database structure
//Interpolating the cross sections to reference grid
   if(scat_cross_database == "VULCAN"){
     MatrixXd sccross_info = ReadVULCANScatCrossSection(scat_cross_file);
     MatrixXd sccross_section_data = InterpolateCrossSection(stellar_input.row(0), sccross_info.row(0), sccross_info.row(1));
     scat_cross_data.col(sc_inx) = sccross_section_data;
   }
   sc_inx++;}
   scatter_species_list = m.suffix().str();
}

 
//Storing the photochemical cross section data
  int PhotoRxn = 0;   
  
  for(int irxn = 0; irxn < nrxn; irxn++){
    auto& rxnObj = *(gas_kin->reaction(irxn));
    std::string rxnEquation = rxnObj.equation();
    std::cout << rxnEquation << std::endl;
    int pos = rxnEquation.find("=");
    rxnEquation.replace(pos, 2, "->");
    std::string photo_cross_info = pinput->GetOrAddString("photocross", rxnEquation, "nan");
   // std::cout << photo_cross_info << std::endl;
    if(photo_cross_info != "nan"){
    PhotoRxn++;
    }
    
    }
  Eigen::MatrixXd photo_cross_data(stellar_input.row(0).size(), PhotoRxn);
  Eigen::MatrixXd qyield_data(stellar_input.row(0).size(), PhotoRxn);
  Eigen::VectorXd RxnIndex(PhotoRxn);
  double Jrate;
  double Irate;
  
//Reading cross section for all photochemical reactions
  int ph_inx = 0;
  for(int irxn = 0; irxn < nrxn; irxn++){
    auto& rxnObj = *(gas_kin->reaction(irxn));
    std::string rxnEquation = rxnObj.equation();
    int pos = rxnEquation.find("=");
    rxnEquation.replace(pos, 2, "->");
    std::string photo_cross_info = pinput->GetOrAddString("photocross", rxnEquation, "nan"); 
    if(photo_cross_info != "nan"){
    std::cout << rxnEquation << std::endl;
    RxnIndex(ph_inx) = irxn;
    std::string photo_cross_database = photo_cross_info.substr(0,photo_cross_info.find(","));  
    std::string photo_cross_ = photo_cross_info.substr(photo_cross_info.find(",")+1,photo_cross_info.length()-1);
    
//VULCAN Photoionization cross section
    if(photo_cross_database == "VULCAN_Ion"){
     MatrixXd photoXsection_info = ReadVULCANPhotoIonCrossSection(photo_cross_);
     photo_cross_data.col(ph_inx) = InterpolateCrossSection(stellar_input.row(0), photoXsection_info.row(0), photoXsection_info.row(1)); 
   }
   
   
//VULCAN Photodissociation cross section 
    if(photo_cross_database == "VULCAN_Diss"){
     MatrixXd photoXsection_info = ReadVULCANPhotoDissCrossSection(photo_cross_);
     photo_cross_data.col(ph_inx) = InterpolateCrossSection(stellar_input.row(0), photoXsection_info.row(0), photoXsection_info.row(1));
   }
   
//ATMOS photochemical cross section
    if(photo_cross_database == "VPL"){
     MatrixXd photoXsection_info = ReadAtmosCrossSection(photo_cross_);
     photo_cross_data.col(ph_inx) = InterpolateCrossSection(stellar_input.row(0), photoXsection_info.row(0), photoXsection_info.row(1));
   
   }

//Max Planck Database photochemical cross section
   if(photo_cross_database == "MPD"){
     MatrixXd photoXsection_info = ReadMPCrossSection(photo_cross_);
     photo_cross_data.col(ph_inx) = InterpolateCrossSection(stellar_input.row(0), photoXsection_info.row(0), photoXsection_info.row(1));

   }

//KINETICS7 Photochemical cross section
   if(photo_cross_database == "KINETICS"){
     MatrixXd photoXsection_info = ReadKINETICSCrossSection(atoi(photo_cross_.c_str())); 
     photo_cross_data.col(ph_inx) = InterpolateCrossSection(stellar_input.row(0), photoXsection_info.row(0), photoXsection_info.row(1));
   }
   
   ph_inx++;
    }
    
   }


ph_inx = 0;
//Reading all the quantum yield
  for(int irxn = 0; irxn < nrxn; irxn++){
    auto& rxnObj = *(gas_kin->reaction(irxn));
    std::string rxnEquation = rxnObj.equation();
    int pos = rxnEquation.find("=");
    rxnEquation.replace(pos, 2, "->");
    std::string qyield_info = pinput->GetOrAddString("qyield", rxnEquation, "nan");
    if(qyield_info != "nan"){
      std::string QYType = qyield_info.substr(0,qyield_info.find(";")); 
      
//Quantum Yield for VULCAN database
      if(QYType == "VULCAN"){
      std::string col_str = qyield_info.substr(qyield_info.find(";")+1,qyield_info.find(",") - qyield_info.find(";")-1 );
      int col_num = atoi(col_str.c_str());
      std::string qyield_file = qyield_info.substr(qyield_info.find(",")+1,qyield_info.length()-qyield_info.find(",") -1);
      MatrixXd QYield_info = ReadQYield(qyield_file);
      qyield_data.col(ph_inx) = InterpolateQYield(stellar_input.row(0), QYield_info.row(0), QYield_info.row(col_num));
      }
//Quantum Yield for KINETICS7 database (QY == 1)
      if(QYType == "KINETICS"){
      qyield_data.col(ph_inx) = MatrixXd::Ones(stellar_input.row(0).size(), 1);

      }


      ph_inx++;
      }

    }

//Chemical evolution
double dh; //m
int iPrev, iNext;
 

// Generate a opacity profile for the atmosphere. The opacity is calculated at every height
 MatrixXd Stellar_activity = MatrixXd::Zero(stellar_input.row(0).size(),nSize);
 VectorXd one = VectorXd::Ones(stellar_input.row(0).size());

   
//Initiating Matrices for Time Evolution
Eigen::SparseMatrix<double>  m_wjac;
m_wjac.resize(nsp, nsp);
MatrixXd mat1(nsp, nsp);
mat1 = MatrixXd::Identity(nsp, nsp);
MatrixXd mat2(nsp, nsp);
VectorXd Un(nsp);
VectorXd Unext(nsp);
VectorXd Uprev(nsp);
VectorXd fluxT(nsp); //Total flux due to diffusion
VectorXd fluxM(nsp); //Flux due to molecular diffusion
VectorXd fluxK(nsp); //Flux due to eddy diffusion
VectorXd fluxA(nsp); //Flux due to ambipolar diffusion

VectorXd Source(nsp);
VectorXd Loss(nsp);
VectorXd D_j(nsp);
VectorXd D_prev(nsp);
VectorXd D_next(nsp);
VectorXd mWt(nsp); //Molecular weight
VectorXd Qn(nsp); //Charge number for filtering ions and neutrals
VectorXd Qe(nsp); //Charge number for filtering ions from neutrals and electron
VectorXd alpha = VectorXd::Zero(nsp); //Thermal diffusion coefficient
VectorXd oem = Eigen::VectorXd::Ones(nsp);
VectorXd dQ(nsp);
double Keddy_j, Keddy_prev, Keddy_next;
double Temp, Press, Kzz;
double charge;
double mm;
VectorXd m_wdot(nsp);
VectorXd mole_frac(nsp);
VectorXd m_subcycle(nsp);

std::string a_species_list = pinput->GetString("thermaldiff", "species");
  while (std::regex_search (a_species_list,m,pattern)) {
    for (auto x:m){
    std::string species_thermaldiff_condition = pinput->GetString("thermaldiff", x);
    species_inx = gas->speciesIndex(x);
    alpha(species_inx) = atof(species_thermaldiff_condition.c_str());
    }
    a_species_list = m.suffix().str();
  }

int iE = gas->speciesIndex("E");
//std::cout << "alpha: " << alpha.transpose() << std::endl;

//Filters for charged particles, and neutrals
for (int insp = 0; insp < nsp; insp++) {
  mWt(insp) = gas->molecularWeight(insp); // Kg/kml
  int elementIndex = gas->elementIndex("E");
  charge = gas->nAtoms(insp, elementIndex);
  if(charge == 0){
  Qn(insp) = 0;
  Qe(insp) = 0;
  }
  else if(charge != 0){
  Qn(insp) = 1;
  Qe(insp) = 1;
  if(charge == 1){
  Qe(insp) = 0;
  }

  }

}

double Kb = 1.38E-23;
int i = 0;
int iSource = 0;
std::string equbm = pinput->GetString("problem", "equilibrate");
if(equbm == "true"){
std::cout << "!!  Initiating Equilibrium  !!" << std::endl;
//Equilibrium as initial condition
  for (int j = 0; j < nSize; j++) {
   Temp = AtmData(iTemp,j);
   Press = AtmData(iPress,j);
   gas->setState_TP(Temp, (Press/1.0132E5)*OneAtm);
   mole_frac = ChemMoleFrac.col(j);
   gas->setMoleFractions(&mole_frac[0]);
   gas->equilibrate("TP");        
   }
std::cout << "!! Atmosphere in thermochemical equilibrium !! \n" << std::endl;
   }

int niter = 0;
int tol = 0;
int ftime = 0;
double t_chem;
int num_chem = 1;
std::cout << "!! Starting the Simulation !!" << std::endl;


//Acceleration due to gravity (m/s)
std::string grav_str = pinput->GetString("grav", "g");
double g = atof(grav_str.c_str()); 

//Scale height (km)
std::string scale_height = pinput->GetString("grav", "H0");
double H0 = atof(scale_height.c_str())*1E3; //km -> m
double D_Max;
//Correct for the atmospheric scattering and absorption above the model boundary
//Using column density of species to determine the absorption and scattering not
//accounted for
int clden = 0;
std::string clden_species_list = pinput->GetOrAddString("coldensity", "species", "nan");
std::cout << clden_species_list << std::endl;
Eigen::VectorXd cldMag = VectorXd::Zero(gas->nSpecies()); //Stores the magnitude of column density difference for TOA
Eigen::VectorXd cldMag_diff = VectorXd::Zero(gas->nSpecies()); //Stores the magnitude of column density difference for each species
while (std::regex_search (clden_species_list,m,pattern)) {
     for (auto x:m){
//Column density at top of atmosphere [From input file]
       std::string clden_toa = pinput->GetOrAddString("coldensity", x, "nan");
       if(clden_toa != "nan"){
       species_inx = gas->speciesIndex(x);
       cldMag(species_inx) = atof(clden_toa.c_str())*1e4; //1/cm^2 -> 1/m^2
       }
       clden++;
     }
      clden_species_list = m.suffix().str();
      }

while(Ttot  < Tmax) {
  MatrixXd Opacity = MatrixXd::Zero(stellar_input.row(0).size(),nSize);
  //Opacity.col(0) = VectorXd::Ones(stellar_input.row(0).size());
  a = ChemMoleFrac;
  n_conc = ChemConc;
  for (int j = 0; j < nSize; j++) {
//Setting T, P, X for each grid point
  //  std::cout << "Setting T, P, X" << std::endl;
    Temp = AtmData(iTemp,j);
    Press = AtmData(iPress,j);
    auto sol2 = newSolution(network_file);
    auto gas2 = sol2->thermo();
    //sol2->setTransportModel("mixture-averaged");
    auto gas_kin2 = sol2->kinetics();
    auto gas_tr2 = sol2->transport();
    Cantera::Kinetics *gasRawPtr = sol2->kinetics().get();
    Cantera::ThermoPhase *gasThermo = gas2.get(); 
    mole_frac = ChemMoleFrac.col(j);
    //std::cout << mole_frac.transpose() << std::endl;
    gas2->setMoleFractions(&mole_frac[0]);
    gas2->setState_TP(Temp, (Press/1.0132E5)*OneAtm);
    
//Setting the multiplier for each reaction involving electron impact process
    std::vector<double> fwd_rates(gas_kin2->nReactions());
  //  std::cout << "Setting rates for electron impact process" << std::endl;
      for(int ionx = 0; ionx < ion_rxn; ionx++){
        gas_kin2->setMultiplier(ion_rxnid(ionx), ion_rate(ionx, j));
        gas_kin2->getFwdRatesOfProgress(fwd_rates.data());
       //std::cout << AtmData(iAlt, j)/1e3 << " " <<  ion_rate(ionx, j) << std::endl;
    }
   // std::cout << "Setting rates for electron impact process done" << std::endl;
//Calculating the photochemical reaction rate
    
//Opacity at the given altitude 
   if(j != 0){
//For each absorber, find the total absorption 
   dh = (AtmData(iAlt,j-1) - AtmData(iAlt,j));

//Atomic and Molecular absorption
   int Absorber = 0; 
   std::string absorber_species_list = pinput->GetString("abscross", "absorbers");
   while (std::regex_search (absorber_species_list,m,pattern)) {
     for (auto x:m){
      species_inx = gas2->speciesIndex(x);
      double number_density = Press/(Kb*Temp);
      Opacity.col(j) = Opacity.col(j) - (mole_fractions(species_inx)*number_density*absorber_cross_data.col(Absorber)*dh/cos(sz_angle*3.14/180));
      Absorber++;
     }
      absorber_species_list = m.suffix().str();
    }

//Rayleigh scattering
    int Scatter = 0;
    std::string scat_species_list = pinput->GetString("scatcross", "scatterers");
    while (std::regex_search (absorber_species_list,m,pattern)) {
     for (auto x:m){
      species_inx = gas2->speciesIndex(x);
      double number_density = Press/(Kb*Temp);
      Opacity.col(j) = Opacity.col(j) - (mole_fractions(species_inx)*number_density*scat_cross_data.col(Scatter)*dh/cos(sz_angle*3.14/180));
      Scatter++;
     }
      scat_species_list = m.suffix().str();
      }
    
//The transmission coefficient from opacity
    Opacity.col(j) = Opacity.col(j) - (dh*handleCustomOpacity(PlanetName, gasThermo, Press, Temp, AtmData(iAlt,j), stellar_input.row(0)));
    //OpacityStorage.col(j) = Opacity.col(j);
    Opacity.col(j) = Opacity.col(j).array().exp().matrix();    
    Opacity.col(j) =  (Opacity.col(j).array()*Opacity.col(j-1).array()).matrix();
//Stellar spectrum at each altitude
    Stellar_activity.col(j) = (Stellar_activity.col(j-1).array()*Opacity.col(j).array()).transpose().matrix();
     }
   if(j == 0){
     clden_species_list = pinput->GetOrAddString("coldensity", "species", "nan");
     if(clden_species_list != "nan"){
//Updating the difference in column density
       cldMag_diff = cldMag; //1/m^2
//Atomic and molecular absorption
       int Absorber = 0;
       std::string absorber_species_list = pinput->GetString("abscross", "absorbers");
      while (std::regex_search (absorber_species_list,m,pattern)) {
        for (auto x:m){
         species_inx = gas2->speciesIndex(x);
         double number_density = Press/(Kb*Temp);
         Opacity.col(j) = Opacity.col(j) - (absorber_cross_data.col(Absorber)*cldMag_diff(species_inx)/cos(sz_angle*3.14/180));
         //std::cout << absorber_cross_data.col(Absorber).transpose()*cldMag_diff(species_inx)/cos(sz_angle*3.14/180) << std::endl;
        // std::cout << Opacity.col(j).transpose() << std::endl;
         Absorber++;
     }
         absorber_species_list = m.suffix().str();
    }
    
    //std::cout << "Absorption done! " << std::endl;
//Rayleigh scattering
      int Scatter = 0;
      std::string scat_species_list = pinput->GetString("scatcross", "scatterers");
      while (std::regex_search (absorber_species_list,m,pattern)) {
      for (auto x:m){
        species_inx = gas2->speciesIndex(x);
        double number_density = Press/(Kb*Temp);
        Opacity.col(j) = Opacity.col(j) - (scat_cross_data.col(Scatter)*cldMag_diff(species_inx)/cos(sz_angle*3.14/180));
        Scatter++;
     }
        scat_species_list = m.suffix().str();
      }
   // std::cout << "Scattering done at TOA!" << std::endl;
//Updating the stellar activity at upper boundary
//The transmission coefficient from opacity
   // std::cout << Opacity.col(j) << std::endl;
   // OpacityStorage.col(j) = Opacity.col(j);
    Opacity.col(j) = Opacity.col(j).array().exp().matrix();
//Stellar spectrum at each altitude
    Stellar_activity.col(j) = (stellar_input.row(1).transpose().array()*Opacity.col(j).array()).transpose().matrix(); 
    //std::cout << "TOA" << std::endl;
    }

    if(clden_species_list == "nan"){
    //Opacity.col(0) = VectorXd::Ones(stellar_input.row(0).size());
    Stellar_activity.col(j) = (stellar_input.row(1).transpose());}
   // std::cout << "Setting TOA opacity" << std::endl;
     }


   for(int rx = 0; rx < PhotoRxn; rx++){
     double j_rate = QPhotoChemRate(stellar_input.row(0),d_wavelength, photo_cross_data.col(rx), qyield_data.col(rx), Stellar_activity.col(j));
     gas_kin2->setMultiplier(RxnIndex(rx), j_rate);
    //std::cout<< RxnIndex(rx) << " " << j_rate << std::endl;
      
   }

//Adding custom reaction rates
   handleCustomChemistry(PlanetName, gasRawPtr, Press, Temp);
   
   ftime++;
//Solving the net production for each species
    gas_kin2->getNetProductionRates(&m_wdot[0]); //Extracting net production rates from Cantera 
    m_wjac = gas_kin2->netProductionRates_ddX(); //Extracting Jacobian from Cantera
    mm  = gas->meanMolecularWeight(); //Mean molcular weight
//Upper boundary flux
    if (j == 0){
    fluxT = flux_upper;
    D_Max = handleCustomMolecularDiffusion(PlanetName, gasThermo, AtmData(iPress,j), AtmData(iTemp,j), mWt).maxCoeff();
    }
    
    if ((j > 0) && (j < nSize-1)) {
    iPrev = j-1;
    iNext = j+1;
    dh = 0.5*(AtmData(iAlt,j-1) - AtmData(iAlt,j+1));

    Keddy_j = AtmData(iKzz, j);

//Update the total diffusion coefficient in terms of scale height
    D_j = handleCustomMolecularDiffusion(PlanetName, gasThermo, AtmData(iPress,j), AtmData(iTemp,j), mWt);
    Keddy_prev = AtmData(iKzz, j-1);
    D_prev = handleCustomMolecularDiffusion(PlanetName, gasThermo, AtmData(iPress,j-1), AtmData(iTemp,j-1), mWt);
    Keddy_next = AtmData(iKzz, j+1);
    D_next = handleCustomMolecularDiffusion(PlanetName, gasThermo, AtmData(iPress,j+1), AtmData(iTemp,j+1), mWt);
//Diffusion terms (central difference scheme)

//Molecular diffusion 
    fluxM = Dflux(j,iNext,iPrev,a,n_conc,nsp,D_j,D_prev,D_next,AtmData(iAlt,j-1) - AtmData(iAlt,j), AtmData(iAlt,j) - AtmData(iAlt,j+1), mWt,AtmData(iTemp,j-1),AtmData(iTemp,j), AtmData(iTemp,j+1), mm, g);

//Ambipolar diffusion
    fluxA = Aflux(j, iNext, iPrev, a, n_conc, nsp, D_j, D_prev, D_next, AtmData(iAlt,j-1) - AtmData(iAlt,j), AtmData(iAlt,j) - AtmData(iAlt,j+1), mWt, AtmData(iTemp,j-1), AtmData(iTemp,j), AtmData(iTemp,j+1), mm, g, iE, alpha);

//Eddy diffusion
    fluxK = Kflux(j,iNext,iPrev,a,n_conc,nsp,Keddy_j,Keddy_prev,Keddy_next,AtmData(iAlt,j-1) - AtmData(iAlt,j), AtmData(iAlt,j) - AtmData(iAlt,j+1));

/*
std::cout << "Flux K" << std::endl;
std::cout << fluxK.transpose() << std::endl;
std::cout << "Flux M" << std::endl;
std::cout << fluxM.transpose() << std::endl;
std::cout << "Flux A" << std::endl;
std::cout << fluxA.transpose() << std::endl;

//Total diffusion
    fluxT = ((oem - Qn).array()*(fluxM + fluxK).array()).matrix();
    fluxT = fluxT + (Qe.array()*fluxA.array()).matrix();

std::cout << "Flux T" << std::endl;
std::cout << fluxT.transpose() << std::endl;
*/ 
 }

//Lower boundary flux
    if (j == nSize - 1){
    fluxT = flux_lower;
    }

    t_chem = dt;
    num_chem = 1;
    m_subcycle = mole_frac;

/*
//Semi-implicit scheme (private communication with Aaron Ridley)
 
//Source term 
  gas_kin2->getCreationRates(&Source[0]); //Extracting source term from Cantera

//Adding diffusion to the source term
  Source = Source - fluxT;

//Loss term and normalized loss
  gas_kin2->getDestructionRates(&Loss[0]); //Extracting loss term from Cantera
  ProdRates.col(j) = Source;
  DiffRates.col(j) = Loss;
  Loss = (Loss.array()/n_conc.col(j).array()).matrix(); //Normalizing the loss term
  Source =  (isnan(abs(Loss.array())) ).select(0, Source);
  Loss =  (isnan(abs(Loss.array())) ).select(0, Loss);
//New concentration
  VectorXd Ix = VectorXd::Ones(nsp);
  n_conc.col(j) = (((Source*dt) + n_conc.col(j)).array()/(Ix - (Loss*dt)).array()).matrix();
  double N_this = n_conc.colwise().sum()(j);
  mole_frac = n_conc.col(j)/N_this;
  std::cout << mole_frac.transpose() << std::endl;

*/


//Backward Euler Scheme (Li and Chen, 2019)
    mat2 = ((mat1) - (m_wjac*dt/gas2->molarDensity() ));
    mat2 = mat2.inverse();
    mat2 = dt*mat2*((m_wdot ) + fluxT);
    //std::cout << fluxT.transpose() << std::endl;
    dQ = mat2/gas2->molarDensity();
    ProdRates.col(j) = m_wdot;
    DiffRates.col(j) = fluxT;
    conv.col(j) = (abs(dQ.array()/mole_frac.array())).matrix();
    conv.col(j) =  (isnan(abs(conv.col(j).array())) ).select(0, conv.col(j));   
    mole_frac = mole_frac + dQ;
//    std::cout << mole_frac.transpose() << std::endl;


/*
//Subcycling Mechanism (Ridley et al., 2006)
   do{   
    for(int num_check = 1; num_check <= num_chem; num_check++){
    mat2 = ((mat1) - (m_wjac*t_chem/gas2->molarDensity() ));
    mat2 = mat2.inverse();
    mat2 = t_chem*mat2*((m_wdot ) - fluxT);
    dQ = mat2/gas2->molarDensity();
//    std::cout << "dt: " << t_chem << std::endl;
    ProdRates.col(j) = m_wdot;
    DiffRates.col(j) = fluxT;
    conv.col(j) = (abs(dQ.array()/mole_frac.array())).matrix();
    conv.col(j) =  (isnan(abs(conv.col(j).array())) ).select(0, conv.col(j));   
    m_subcycle = m_subcycle + dQ;
    if(dt == Tmin){
      conv.col(j) = (abs(conv.col(j).array()).isInf()).select(0, conv.col(j));
    }
    
//   std::cout << "outside the update loop: " << conv.col(j).maxCoeff() << std::endl;

    if(conv.col(j).maxCoeff() >= 0.25){
      t_chem = t_chem/10;
      num_chem = num_chem*10;
      num_check =  1;
      m_subcycle = mole_frac;
      std::cout << "Inside the update loop: " << t_chem << std::endl;
    } }
//     std::cout << conv.col(j).maxCoeff() << std::endl;
    } 
    while(conv.col(j).maxCoeff() >= 0.25 );
*/

//Testing the change in mole fraction
// std::cout << "came out of loop" << std::endl;
//Update the mole fractions
//    mole_frac = m_subcycle;
    VectorXd krate(nrxn);
    gas_kin2->getFwdRateConstants(&krate[0]);
    
/*   
//Upper boundary mixing ratio
    if (j == 0){
    std::string u_species_list = pinput->GetOrAddString("upperboundaryMixRat", "species", "nan");
    if(u_species_list != "nan"){
    while (std::regex_search (u_species_list,m,pattern)) {
      for (auto x:m){
        std::string species_upperboundary_MixRat = pinput->GetString("upperboundaryMixRat", x);
        species_inx = gas->speciesIndex(x);
        mole_frac(species_inx) = atof(species_upperboundary_MixRat.c_str());
    }
      u_species_list = m.suffix().str();
  } } }
*/


//Lower boundary mixing ratio
    if (j == nSize-1){
    std::string l_species_list = pinput->GetOrAddString("lowerboundaryMixRat", "species", "nan");
    mole_frac = VectorXd::Ones(nsp);
    mole_frac = mole_frac*1e-30;
    if(l_species_list != "nan"){
    while (std::regex_search (l_species_list,m,pattern)) {
      for (auto x:m){
        std::string species_lowerboundary_MixRat = pinput->GetString("lowerboundaryMixRat", x);
        species_inx = gas2->speciesIndex(x);
        mole_frac(species_inx) = atof(species_lowerboundary_MixRat.c_str());
    }
      l_species_list = m.suffix().str();

  } } 

     conv.col(nSize -1) = VectorXd::Zero(gas2->nSpecies());
  }
 
  
//Setting the mole fractions and concentrations
    gas2->setMoleFractions(&mole_frac[0]);
    gas2->setState_TP(Temp, (Press/1.0132E5)*OneAtm);
    gas2->getMoleFractions(&ChemMoleFrac.col(j)[0]);  
    gas2->getConcentrations(&ChemConc.col(j)[0]);  
    
//Update the pressure in storage
    //AtmData(iPress,j) = gas->pressure();
} 
   
    std::cout << "Simulation completed at time step t = " << Ttot << std::endl;
  //  std::cout << D_Max  << std::endl;
//Limit set by Explicit diffusion

    if(6E8*dt/(dh*dh) < 0.1){
      dt = dt*(1.25);
     }
    Ttot = Ttot + dt;
    OpacityStorage = Opacity;
} 
   

//This is where the photochemistry definition ends (Anything beyond is not defined in the scope of photochemistry)  
  }

  if(photochem != "true"){
   Eigen::SparseMatrix<double>  m_wjac;
   m_wjac.resize(nsp, nsp);
   MatrixXd mat1(nsp, nsp);
   mat1 = MatrixXd::Identity(nsp, nsp);
   MatrixXd mat2(nsp, nsp);
   VectorXd Un(nsp);
   VectorXd Unext(nsp);
   VectorXd Uprev(nsp);
   VectorXd flux1(nsp);
   VectorXd dQ(nsp);
   double Keddy_j, Keddy_next, Keddy_prev;
   double Temp, Press, Kzz;
   VectorXd m_wdot(nsp);
   VectorXd mole_frac(nsp);
   double dh;
   int iPrev, iNext;
   double Kb = 1.38E-23;
   int i = 0;
   int iSource = 0;
   std::cout << "Starting the simulation now (without photochemistry ofc!)" << std::endl;
   while(Ttot  < Tmax){
    for (int j = 0; j < nSize; j++) {
//Setting T, P, X for each grid point
    iPrev = j-1;
    iNext = j+1;
    Temp = AtmData(iTemp,j);
    Press = AtmData(iPress,j);
    auto sol2 = newSolution(network_file);
    auto gas2 = sol2->thermo();
    auto gas_kin2 = sol2->kinetics();  
    mole_frac = ChemMoleFrac.col(j);
    gas2->setMoleFractions(&mole_frac[0]);   
    gas2->setState_TP(Temp, (Press/1.0132E5)*OneAtm);
    Cantera::Kinetics *gasRawPtr = sol2->kinetics().get();
    Cantera::ThermoPhase *gasThermo = gas2.get(); 
//Setting the multiplier for each reaction based Ion kinetics
    for(int ionx = 0; ionx < ion_rxn; ionx++){
    gas_kin2->setMultiplier(ion_rxnid(ionx), ion_rate(ionx, j));

    }

//Adding custom reaction rates
   handleCustomChemistry(PlanetName, gasRawPtr, Press, Temp);

//Solving the net production for each species
    gas_kin->getNetProductionRates(&m_wdot[0]); //Extracting net production rates from Cantera 
    m_wjac = gas_kin->netProductionRates_ddX(); //Extracting Jacobian from Cantera
    
//Upper boundary flux
    if (j == 0){
    flux1 = flux_upper;
    }
    
    if ((j > 0) && (j < nSize-1)) {
    dh = 0.5*(AtmData(iAlt,j-1) - AtmData(iAlt,j+1));
    Keddy_j = AtmData(iKzz, j);
    Keddy_prev = AtmData(iKzz, j-1);
    Keddy_next = AtmData(iKzz, j+1);

//Diffusion terms (central difference scheme)
    flux1 = Kflux(j,iNext,iPrev,a,n_conc,nsp,Keddy_j,Keddy_prev,Keddy_next,AtmData(iAlt,j-1) - AtmData(iAlt,j), AtmData(iAlt,j) - AtmData(iAlt,j+1));
    }

//Lower boundary flux
    if (j == nSize - 1){
    flux1 = flux_lower;
    }

//Backward Euler Scheme (Li and Chen, 2019)
    mat2 = ((mat1) - (m_wjac*dt/gas->molarDensity() ));
    mat2 = mat2.inverse();
    mat2 = dt*mat2*((m_wdot ) - flux1);
    dQ = mat2/gas->molarDensity();
    mole_frac = mole_frac + dQ;
    ProdRates.col(j) = m_wdot;
    DiffRates.col(j) = flux1;
    conv.col(j) = m_wdot - flux1 ;//(dQ.array()/mole_frac.array()).matrix();
    conv.col(j) =  (isnan(abs(conv.col(j).array())) ).select(0, conv.col(j));
    
/*    
//Upper boundary mixing ratio
    if (j == 0){
    std::string u_species_list = pinput->GetOrAddString("upperboundaryMixRat", "species", "nan");
    if(u_species_list != "nan"){
    while (std::regex_search (u_species_list,m,pattern)) {
      for (auto x:m){
        std::string species_upperboundary_MixRat = pinput->GetString("upperboundaryMixRat", x);
        species_inx = gas->speciesIndex(x);
        mole_frac(species_inx) = atof(species_upperboundary_MixRat.c_str());
    }
      u_species_list = m.suffix().str();
  } } }
*/


//Lower boundary mixing ratio
    if (j == nSize-1){
    std::string l_species_list = pinput->GetOrAddString("lowerboundaryMixRat", "species", "nan");
    if(l_species_list != "nan"){
    while (std::regex_search (l_species_list,m,pattern)) {
      for (auto x:m){
        std::string species_lowerboundary_MixRat = pinput->GetString("lowerboundaryMixRat", x);
        species_inx = gas->speciesIndex(x);
        mole_frac(species_inx) = atof(species_lowerboundary_MixRat.c_str());
    }
      l_species_list = m.suffix().str();

  } } 

     conv.col(nSize -1) = VectorXd::Zero(gas->nSpecies());
  }
 
  
//Setting the mole fractions and concentrations
    gas->setMoleFractions(&mole_frac[0]);
    gas2->setState_TP(Temp, (Press/1.0132E5)*OneAtm);
    gas->getMoleFractions(&ChemMoleFrac.col(j)[0]);  
    gas->getConcentrations(&ChemConc.col(j)[0]);
} 

    Ttot = Ttot + dt;

    if(kmax*dt/(dh*dh) < 0.01){
    dt = dt*1.25;
    }
    if((Ttot > Tstart)  && (iSource == 0)){
    dt = 1e-8;
    iSource = 1;
    }
    std::cout << "Simulation completed at time step t = " << Ttot << std::endl;
    i++;
}


}
//Writing output into NetCDF file

//std::vector<Eigen::MatrixXd> matrix3d(Time.size(), Eigen::MatrixXd(nsp, nSize));

VectorXd cPress = AtmData.row(iPress);
VectorXd cTemp = AtmData.row(iTemp);
VectorXd cHeight = AtmData.row(iAlt);
VectorXd cKzz = AtmData.row(iKzz);
MatrixXd dPRates = ProdRates;
MatrixXd dDRates = DiffRates;
MatrixXd OP = OpacityStorage;
VectorXd Wavelength = stellar_input.row(0);
#if NETCDFOUTPUT
int ifile;
string fname = pinput->GetString("output", "out_file");
int status = nc_create(fname.c_str(), NC_NETCDF4, &ifile);
int dimids[3];
int dWav[2], iWav;
int alt, iPres, iTem, iKeddy, iMolDiff, iAltz, iChem, isp, iOpacity, iProd, iDiff, iTs;


// Atmospheric Properties (All in SI units!)
nc_def_dim(ifile, "Pressure", nSize, &dimids[2]);

nc_def_dim(ifile, "Species", nsp, &dimids[1]);

nc_def_dim(ifile, "Wavelength", Wavelength.size() ,  &dWav[0]);
nc_def_dim(ifile, "Altitude", nSize ,  &dWav[1]);

nc_def_var(ifile, "Wavelength", NC_DOUBLE, 1, &dWav[0], &iWav);
nc_put_var_double(ifile, iWav, &Wavelength[0]);

nc_def_var(ifile, "Pressure", NC_DOUBLE, 1, &dimids[2], &iPres);
nc_put_var_double(ifile, iPres, &cPress[0]);

nc_def_var(ifile, "Keddy", NC_DOUBLE, 1, &dimids[2], &iKeddy);
nc_put_var_double(ifile, iKeddy, &cKzz[0]);

nc_def_var(ifile, "Temp", NC_DOUBLE, 1, &dimids[2], &iTem);
nc_put_var_double(ifile, iTem, &cTemp[0]);

nc_def_var(ifile, "Altitude", NC_DOUBLE, 1, &dimids[2], &iAltz);
nc_put_var_double(ifile, iAltz, &cHeight[0]);

nc_def_var(ifile, "Opacity", NC_DOUBLE, 2, dWav, &iOpacity);
nc_put_var_double(ifile, iOpacity, OP.data());

init_species_list = pinput->GetString("output", "species");
  while (std::regex_search (init_species_list,m,pattern)) {
    for (auto x:m){
    std::string species = x;
    std::string prod = "Prod";
    std::string diff = "Diff";
    std::string prod_species = prod + species;
    std::string diff_species = diff + species;
    species_inx = gas->speciesIndex(x);
    VectorXd cChem = ChemConc.row(species_inx)*1E3*6.022E23*1E-6; // #/cm^-3
    VectorXd pChem = ProdRates.row(species_inx); // kmol/m^3s (Cantera units)
    VectorXd dChem = DiffRates.row(species_inx); // kmol/m^3s (Cantera units)
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
