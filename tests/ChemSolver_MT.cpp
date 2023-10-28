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

// output stream
#include <iostream>
#include <fstream>
#include <cantera/numerics/eigen_dense.h>
#include <cantera/numerics/eigen_sparse.h>
#include <vector>
#include <string>
#include <sstream>
#include <regex>

// Athena++ header
#include <parameter_input.hpp>

// C3M header
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

//Function for Central Difference Scheme
MatrixXd diff_flux(int inxn, int inxNext, int inxPrev, MatrixXd Sp, MatrixXd Np, int ns, double Kzz, double Kzz_prev, double Kzz_next, double dx){
  //std::cout << "Lets start the diffusion process" << std::endl;
  VectorXd diff_flux(ns), Xn(ns), Xnext(ns), Xprev(ns), Nn(ns), Nnext(ns), Nprev(ns); 
  VectorXd beta1(ns), beta2(ns);
  Xn = Sp.col(inxn);
  Xnext = Sp.col(inxNext);
  Xprev = Sp.col(inxPrev);
  Nn = Np.col(inxn).colwise().sum();
  Nnext = Np.col(inxPrev).colwise().sum();
  Nprev = Np.col(inxNext).colwise().sum();
  //std::cout << Nn << Nprev << std::endl;
  beta1 = ((Kzz + Kzz_prev)/2)*((Nn + Nprev)/2); //j-1
  beta2 = ((Kzz + Kzz_next)/2)*((Nn + Nnext)/2); //j+1
  //std::cout << Xn*beta2  << std::endl; 
  diff_flux =   (((-1*((Xnext - Xn)*beta2)/dx) + (((Xn - Xprev))*beta1/dx))/dx);
return diff_flux;
}

int main(int argc, char **argv) {
  
//Reading input file
  IOWrapper infile;
  infile.Open("test_athena.inp", IOWrapper::FileMode::read);
  ParameterInput *pinput = new ParameterInput();

  pinput->LoadFromFile(infile);
  infile.Close();

//Loading the input parameters for atmospheric profile and reaction network
  std::string atm_file = pinput->GetString("problem", "planet");
  std::string network_file = pinput->GetString("problem", "network");
  std::string profile_file = pinput->GetOrAddString("problem", "chemInp", "nan");

//Reading the chemical kinetics network
  auto sol = newSolution(network_file);
  //auto sol = newSolution(network_file);
  auto gas = sol->thermo();
  auto gas_kin = sol->kinetics();  
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
  double dt = atof(time_step.c_str());
  double Tmax = atof(max_time.c_str());
  double Ttot = 0.0;
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
  MatrixXd AtmData(4, nSize);
  int iTemp = 0;
  int iPress =  1;
  int iKzz =  2;
  int iAlt = 3;
  
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
  MatrixXd a(nsp, nSize);
  MatrixXd n_conc(nsp, nSize);
  for (int i = 0; i < nSize; i++) {
      gas->setState_TP(AtmData(iTemp,i), (AtmData(iPress,i)/1.0132E5)*OneAtm);
      ChemMoleFrac.col(i) = mole_fractions;
    }


//Chemical profiles from input file
  if(profile_file != "nan"){
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
    input_species_list = pinput->GetString("profile", "species");
    while (std::regex_search (input_species_list,m,pattern)) {
    for (auto x:m){
     species_inx = gas->speciesIndex(x);
     if(chem_inx != chemP){
       getline(InFile,inp2,',');
       std::cout << inp2 << std::endl;
       ChemMoleFrac(species_inx,i) = atof(inp2.c_str());
       }

     if(chem_inx == chemP-1){
       getline(InFile,inp2, '\n');
       std::cout << inp2 << std::endl;
       ChemMoleFrac(species_inx,i) = atof(inp2.c_str());}
       
     chem_inx++;
    }
    
    input_species_list = m.suffix().str();
    }
    i++;
  }}
 InFile.close();
for (int i = 0; i < nSize; i++) {
      gas->setState_TP(AtmData(iTemp,i), (AtmData(iPress,i)/1.0132E5)*OneAtm);
      gas->setMoleFractions(&ChemMoleFrac.col(i)[0]);
      gas->getConcentrations(&ChemConc.col(i)[0]);
    }

std::cout << "All inputs loaded into C3M " << std::endl;

//Photochemistry
  std::string photochem = pinput->GetString("problem", "photochem");
  if(photochem == "true")
  {
  
  std::cout << "Starting photochemistry calculations" << std::endl;
//Radiative transfer
//Reading the Stellar Irradiance Input
  std::string stellar_input_file = pinput->GetString("radtran", "solar");
  std::string radius = pinput->GetString("radtran", "radius");
  std::string reference =   pinput->GetString("radtran", "reference");
  double rad = atof(radius.c_str());
  double ref = atof(reference.c_str());
  MatrixXd stellar_input = ReadStellarRadiationInput(stellar_input_file, rad, ref);
  MatrixXd d_wavelength(stellar_input.row(0).size(), 1);
  MatrixXd wavelength = stellar_input.row(0);
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

/*    
   if(absorber_cross_database == "KINETICS"){
     std::cout << absorber_cross_database << std::endl;

    }
*/
    
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
//Storing the photochemical cross section data
  int PhotoRxn = 0;
  
  
  
  for(int irxn = 0; irxn < nrxn; irxn++){
    std::cout << gas_kin->reactionString(irxn) << std::endl;
    std::string rxnEquation = gas_kin->reactionString(irxn);
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
  Eigen::VectorXd Jrate(PhotoRxn);
  
//Reading cross section for all photochemical reactions
  int ph_inx = 0;
  for(int irxn = 0; irxn < nrxn; irxn++){
    std::string rxnEquation = gas_kin->reactionString(irxn);
    int pos = rxnEquation.find("=");
    rxnEquation.replace(pos, 2, "->");
    std::string photo_cross_info = pinput->GetOrAddString("photocross", rxnEquation, "nan"); 
    if(photo_cross_info != "nan"){
    RxnIndex(ph_inx) = irxn;
    std::string photo_cross_database = photo_cross_info.substr(0,photo_cross_info.find(","));  
    std::string photo_cross_ = photo_cross_info.substr(photo_cross_info.find(",")+1,photo_cross_info.length()-1);
    
    
    if(photo_cross_database == "VULCAN_Ion"){
     MatrixXd photoXsection_info = ReadVULCANPhotoIonCrossSection(photo_cross_);
     photo_cross_data.col(ph_inx) = InterpolateCrossSection(stellar_input.row(0), photoXsection_info.row(0), photoXsection_info.row(1)); 
   }
    
    if(photo_cross_database == "VULCAN_Diss"){
     MatrixXd photoXsection_info = ReadVULCANPhotoDissCrossSection(photo_cross_);
     photo_cross_data.col(ph_inx) = InterpolateCrossSection(stellar_input.row(0), photoXsection_info.row(0), photoXsection_info.row(1));
     //std::cout <<    photo_cross_data.col(ph_inx).transpose() << std::endl;
   }
   
    if(photo_cross_database == "VPL"){
     MatrixXd photoXsection_info = ReadAtmosCrossSection(photo_cross_);
     photo_cross_data.col(ph_inx) = InterpolateCrossSection(stellar_input.row(0), photoXsection_info.row(0), photoXsection_info.row(1));
   
   }

   if(photo_cross_database == "MPD"){
     MatrixXd photoXsection_info = ReadMPCrossSection(photo_cross_);
     photo_cross_data.col(ph_inx) = InterpolateCrossSection(stellar_input.row(0), photoXsection_info.row(0), photoXsection_info.row(1));

   }
   
   if(photo_cross_database == "KINETICS"){
     MatrixXd photoXsection_info = ReadKINETICSCrossSection(atoi(photo_cross_.c_str())); //In this case its not file name but the reaction number
     photo_cross_data.col(ph_inx) = InterpolateCrossSection(stellar_input.row(0), photoXsection_info.row(0), photoXsection_info.row(1));
   }
   
   ph_inx++;
    }
    
//If the string is not equal to nan and then proceed and load the photochemical cross sectikons
   }


ph_inx = 0;
//Reading all the quantum yield
  for(int irxn = 0; irxn < nrxn; irxn++){
    std::string rxnEquation = gas_kin->reactionString(irxn);
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
      //std::cout << qyield_data.col(ph_inx).transpose() << std::endl;
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


   
// Beer-Lambert Solver

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

std::cout << "!! Starting the Simulation !!" << std::endl;
while(Ttot  < Tmax) {
  MatrixXd Opacity = MatrixXd::Zero(stellar_input.row(0).size(),nSize);
  Opacity.col(0) = VectorXd::Ones(stellar_input.row(0).size());
  a = ChemMoleFrac;
  n_conc = ChemConc;
  for (int j = 0; j < nSize; j++) {
//Setting T, P, X for each grid point
    iPrev = j-1;
    iNext = j+1;
    Temp = AtmData(iTemp,j);
    Press = AtmData(iPress,j);
    gas->setState_TP(Temp, (Press/1.0132E5)*OneAtm);
    mole_frac = ChemMoleFrac.col(j);
    gas->setMoleFractions(&mole_frac[0]);
    //std::cout << "Temp: " << Temp << " Press: " << Press << std::endl;
    //std::cout << "mole_frac: " << mole_frac.transpose() << std::endl;
//Calculating the photochemical reaction rate
    
//Opacity at the given altitude 
   if(j != 0){
//For each absorber, find the total absorption 
   dh = (AtmData(iAlt,j-1) - AtmData(iAlt,j));
   int Absorber = 0; 
   std::string absorber_species_list = pinput->GetString("abscross", "absorbers");
   while (std::regex_search (absorber_species_list,m,pattern)) {
     for (auto x:m){
      species_inx = gas->speciesIndex(x);
      double number_density = Press/(Kb*Temp);
      //std::cout << mole_fractions(species_inx)*number_density*absorber_cross_data.col(Absorber)*dh << std::endl;
      //std::cout << number_density << std::endl;
      Opacity.col(j) = Opacity.col(j) - (mole_fractions(species_inx)*number_density*absorber_cross_data.col(Absorber)*dh);
      Absorber++;
     }
      absorber_species_list = m.suffix().str();
    }
    
//The transmission coefficient from opacity
    Opacity.col(j) = Opacity.col(j).array().exp().matrix();    
    Opacity.col(j) =  (Opacity.col(j).array()*Opacity.col(j-1).array()).matrix();
//Stellar spectrum at each altitude
    Stellar_activity.col(j) = (Stellar_activity.col(j-1).array()*Opacity.col(j).array()).transpose().matrix();
     }
   if(j == 0){
     Stellar_activity.col(j) = (stellar_input.row(1).transpose());
     }
     //std::cout << "level " << AtmData(iAlt,j)/1E3 << " km" << std::endl;
//Computing the photochemical reaction rate for each reaction
   for(int rx = 0; rx < PhotoRxn; rx++){
     double j_rate = QPhotoChemRate(stellar_input.row(0),d_wavelength, photo_cross_data.col(rx), qyield_data.col(rx), Stellar_activity.col(j));
     //double j_rate = PhotoChemRate(stellar_input.row(0),photo_cross_data.col(rx), Stellar_activity.col(j));
     //std::cout << AtmData(iAlt,j)/1e3   << " " << j_rate <<  " " << gas_kin->reactionString(RxnIndex(rx)) << std::endl;
//Setting the multiplier for each reaction
     
   if(i == 0){
     //std::cout << gas_kin->multiplier(RxnIndex(rx)) << std::endl;
     gas_kin->setMultiplier(RxnIndex(rx), j_rate);
     //std::cout << gas_kin->reactionString(RxnIndex(rx)) <<" "  <<  j_rate << std::endl;
     //std::cout << "new photochemical rates for reaction " << gas_kin->reactionString(RxnIndex(rx)) << std::endl;
     //std::cout << gas_kin->multiplier(RxnIndex(rx)) << std::endl;
     Jrate(rx) = j_rate;

   }
   if(i != 0){
     if(Jrate(rx) != 0){
       double multiplier = j_rate/Jrate(rx);
       //std::cout << gas_kin->multiplier(RxnIndex(rx)) << std::endl;
       gas_kin->setMultiplier(RxnIndex(rx), multiplier);
       //std::cout << gas_kin->reactionString(RxnIndex(rx)) <<" "  <<  j_rate << std::endl;
       //std::cout << "new photochemical rates for reaction" << gas_kin->reactionString(RxnIndex(rx)) << std::endl;
       //std::cout << gas_kin->multiplier(RxnIndex(rx)) << std::endl;
       
       Jrate(rx) = j_rate;
   }
   i++;
   }
   
   }

//Solving the net production for each species
    gas_kin->getNetProductionRates(&m_wdot[0]); //Extracting net production rates from Cantera
    m_wjac = gas_kin->netProductionRates_ddX(); //Extracting Jacobian from Cantera
    //std::cout << m_wdot.transpose() << std::endl;
//Backward Euler Scheme (Li and Chen, 2019)
    if (j == 0){
    flux1 = flux_upper;
    }
    
    if ((j > 0) && (j < nSize-1)) {
    dh = 0.5*(AtmData(iAlt,j-1) - AtmData(iAlt,j+1));
    Keddy_j = AtmData(iKzz, j);
    Keddy_prev = AtmData(iKzz, j-1);
    Keddy_next = AtmData(iKzz, j+1);
//Diffusion terms (central difference scheme)
    
    flux1 = diff_flux(j,iNext,iPrev,a,n_conc,nsp,Keddy_j,Keddy_prev,Keddy_next,dh);
    }
//Backward Euler Scheme (Li and Chen, 2019)
    mat2 = ((mat1) - (m_wjac*dt/gas->molarDensity() ));
    mat2 = mat2.inverse();
    mat2 = dt*mat2*((m_wdot ) + flux1);
    dQ = mat2/gas->molarDensity();
    mole_frac = mole_frac + dQ;

//Integration for each species



/*    
//Forced boundary condition
//Upper boundary
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
//Lower boundary
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
  } } }
 
  
//Setting the mole fractions and concentrations
    gas->setMoleFractions(&mole_frac[0]);
    gas->getMoleFractions(&ChemMoleFrac.col(j)[0]);  
    gas->getConcentrations(&ChemConc.col(j)[0]);


} 
    std::cout << "Simulation completed at time step t = " << Ttot << std::endl;
    Ttot = Ttot + dt;
    if(kmax*dt/(dh*dh) < 0.1){
    dt = dt*1.5;
    }
    
} 
      
//This is where the photochemistry definition ends (Anything beyond is not defined in the scope of photochemistry)  
  }

  if(photochem != "true"){
//If photochemistry is set equal to false
//Initiating Matrices for Time Evolution
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
std::cout << "Starting the simulation now" << std::endl;
while(Ttot  < Tmax){
  for (int j = 0; j < nSize; j++) {
//Setting T, P, X for each grid point
    iPrev = j-1;
    iNext = j+1;
    Temp = AtmData(iTemp,j);
    Press = AtmData(iPress,j);
    gas->setState_TP(Temp, (Press/1.0132E5)*OneAtm);
    mole_frac = ChemMoleFrac.col(j);
    gas->setMoleFractions(&mole_frac[0]);   
    

//Solving the net production for each species
    gas_kin->getNetProductionRates(&m_wdot[0]); //Extracting net production rates from Cantera
    m_wjac = gas_kin->netProductionRates_ddX()/gas->molarDensity(); //Extracting Jacobian from Cantera
    a = ChemMoleFrac;
    if (j == 0){
    flux1 = flux_upper;
    }
    if (j == nSize-1){
    flux1 = flux_lower;
    }
    if ((j > 0) && (j < nSize-1)) {
    dh = 0.5*(AtmData(iAlt,j+1) - AtmData(iAlt,j-1));
    Keddy_j = AtmData(iKzz, j);
    Keddy_prev = AtmData(iKzz, j-1);
    Keddy_next = AtmData(iKzz, j+1);
//Diffusion terms (central difference scheme)
    flux1 = diff_flux(j,iNext,iPrev,a,n_conc,nsp,Keddy_j,Keddy_prev,Keddy_next,dh);
    }
//Integration for each species
//Backward Euler Scheme (Li and Chen, 2019)
    mat2 = ((mat1) - (dt*m_wjac/gas->molarDensity()));
    mat2 = mat2.inverse();
    mat2 = mat2*((m_wdot/gas->molarDensity() ) - (flux1/gas->molarDensity() ));
    dQ = mat2*dt;
    mole_frac = mole_frac + dQ;

//Forcing the lower boundary condition   
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
  } } }
    
    
    
    gas->setMoleFractions(&mole_frac[0]);
    gas->getMoleFractions(&ChemMoleFrac.col(j)[0]);
    gas->getConcentrations(&ChemConc.col(j)[0]);
} 

    Ttot = Ttot + dt;

    if(kmax*dt/(dh*dh) < 1){
    dt = dt*1.25;
    }
    std::cout << "Simulation completed at time step t = " << Ttot << std::endl;
    i++;
}


}
//Writing output into NetCDF file

VectorXd cPress = AtmData.row(iPress);
VectorXd cTemp = AtmData.row(iTemp);
VectorXd cHeight = AtmData.row(iAlt);
VectorXd cKzz = AtmData.row(iKzz);
#if NETCDFOUTPUT
int ifile;
string fname = pinput->GetString("output", "out_file");
int status = nc_create(fname.c_str(), NC_NETCDF4, &ifile);
int alt, iPres, iTem, iKeddy, iAltz, iChem;
// Atmospheric Properties (All in SI units!)
nc_def_dim(ifile, "Pressure", nSize, &alt);
nc_def_var(ifile, "Pressure", NC_DOUBLE, 1, &alt, &iPres);
nc_put_var_double(ifile, iPres, &cPress[0]);

nc_def_var(ifile, "Keddy", NC_DOUBLE, 1, &alt, &iKeddy);
nc_put_var_double(ifile, iKeddy, &cKzz[0]);

nc_def_var(ifile, "Temp", NC_DOUBLE, 1, &alt, &iTem);
nc_put_var_double(ifile, iTem, &cTemp[0]);

nc_def_var(ifile, "Altitude", NC_DOUBLE, 1, &alt, &iAltz);
nc_put_var_double(ifile, iAltz, &cHeight[0]);



init_species_list = pinput->GetString("output", "species");
  while (std::regex_search (init_species_list,m,pattern)) {
    for (auto x:m){
    std::string species = x;
    species_inx = gas->speciesIndex(x);
    VectorXd cChem = ChemConc.row(species_inx)*6.022E23*1E3*1E-6;
    const char* ccx = &species[0];
    nc_def_var(ifile, ccx, NC_DOUBLE, 1, &alt, &iChem);
    nc_put_var_double(ifile, iChem, &cChem[0]);
    
    }
    init_species_list = m.suffix().str();
  }

nc_close(ifile);
#endif


delete pinput;
}
