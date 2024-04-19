// @sec3{Include files}
// Author: Ananyo Bhattacharya
// Affiliation: University of Michigan
// Email: ananyo@umich.edu
// It is a test code for constructing one dimensional plane parallel atmosphere using cantera functions
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
#include <cantera/numerics/Integrator.h>
#include <cantera/numerics/CVodesIntegrator.h>
#include <vector>
#include <string>
#include <sstream>
#include <regex>
#include <cantera/oneD/Sim1D.h>
#include <cantera/oneD/StFlow.h>
#include <cantera/oneD/Boundary1D.h>
#include <cantera/oneD/DomainFactory.h>
#include <cantera/base/stringUtils.h>
#include <cantera/transport/Transport.h>
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

int main(int argc, char **argv) {
  
//Reading input file
  IOWrapper infile;
  infile.Open("test_athena.inp", IOWrapper::FileMode::read);
  ParameterInput *pinput = new ParameterInput();

  pinput->LoadFromFile(infile);
  infile.Close();

//Loading the input parameters for atmospheric profile and reaction network
  std::string atm_file = pinput->GetString("problem", "planet");
  std::string PlanetName = pinput->GetString("problem", "planetName");
  std::string network_file = pinput->GetString("problem", "network");
  std::string profile_file = pinput->GetOrAddString("problem", "chemInp", "nan");

//Reading the chemical kinetics network
  auto sol = newSolution(network_file);
  auto gas = sol->thermo();
  auto gas_kin = sol->kinetics();  
  auto gas_tp = sol->transport();
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
  std::string max_time = pinput->GetString("integrator", "Tmax");
  double dt = atof(time_step.c_str());
  double Tmax = atof(max_time.c_str());
  double Ttot = 0.0;
  std::cout << "Time step: " << dt << std::endl;


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
  MatrixXd AtmData(6, nSize);
  int iTemp = 0;
  int iPress =  1;
  int iKzz =  2;
  int iAlt = 3;
  int locn = 4;
  int iU = 5;
  
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
      AtmData(locn, inx) = double(inx)/double(nSize-1);
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
  std::vector<Eigen::MatrixXd> ChemMole;
  MatrixXd ChemMoleFrac(nsp, nSize);
  MatrixXd ChemConc(nsp, nSize);
  MatrixXd ChemMassFrac(nsp, nSize);
  MatrixXd a(nsp, nSize);
  MatrixXd n_conc(nsp, nSize);
  MatrixXd conv(nsp, nSize);
  MatrixXd ProdRates(nsp, nSize);
  MatrixXd DiffRates(nsp, nSize);
  VectorXd Time(1);

    
//Chemical species profile from input file
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
    i++;
  }}
 InFile.close();

//Feeding the inputs to Cantera One-D objects

  auto flow = newDomain<StFlow>("gas-flow", sol, "flow");
  flow->setupGrid(AtmData.row(iAlt).size(),  &AtmData.row(iAlt).reverse()[0]);
  double massFlowRate = 0.0; //Inlet flow rate


//Setting up the T, P, Y for the One-D object

for (int i = 0; i < nSize; i++) {
      std::vector<double> data((1 + nsp), 1);
      flow->setTemperature(i, AtmData(iTemp,nSize -1 - i));
    //  flow->setPressure(i, AtmData(iPress,nSize -1 - i));
    //  std::cout << AtmData(iTemp,nSize -1 - i) << " " << AtmData(iPress,nSize -1 - i) << std::endl;
      gas->setMoleFractions(&ChemMoleFrac.col(nSize -1 - i)[0]);
      gas->setState_TP(AtmData(iTemp,nSize -1 - i), (AtmData(iPress,nSize -1 - i)/1.0132E5)*OneAtm);
      gas->getMassFractions(ChemMassFrac.col(nSize -1 - i).data());
      AtmData(iU, nSize -1 - i) = massFlowRate;
    }

///  flow->setMassFrac(ChemMoleFrac);

flow->setKinetics(gas_kin);
flow->setTransport(gas_tp);
flow->doEnergy(0);

/*Setting up the boundary conditions
As the domain starts from TOA, the outer refers to the upper boundary, and inlet will refer to the lower boundary or surface
depending upon the application
*/

//Inlet condition with T, mole fraction and total mass flow rate
//auto inlet = std::make_shared<Inlet1D>();
auto inlet = newDomain<Inlet1D>("inlet", sol);
std::cout << "Inlet set up done! " << std::endl;
inlet->setMoleFractions(&ChemMoleFrac.col(nSize-1)[0]);
inlet->setTemperature(AtmData(iTemp,nSize-1));
inlet->setMdot(massFlowRate);


//Outlet condition
auto outlet = newDomain<Outlet1D>("outlet", sol);
//outlet->setMoleFractions(&ChemMoleFrac.col(0)[0]);
outlet->setTemperature(AtmData(iTemp,0));


//Setting up solution domain
vector<shared_ptr<Domain1D>> domains { inlet, flow, outlet };
//std::vector<std::shared_ptr<Cantera::Domain1D>> domains = {inlet, flow, outlet};
Sim1D sim(domains);

std::cout << "Simulation object has been set up" << std::endl;


//Solving the equations for momentum, continuity, and chemical diffusion
//auto transportModel = newTransportMgr("Multi", &gas);
//gas.setTransport(transportModel);
flow->setTransportModel("mixture-averaged"); // Ignore species diffusion
std::cout << AtmData.row(locn) << std::endl;
std::vector<double> locs;
std::vector<double> Tvalues;
std::vector<double> Uvalues;
std::vector<double> Xvalues;

locs.resize(AtmData.cols());
Tvalues.resize(AtmData.cols());
Uvalues.resize(AtmData.cols());
Xvalues.resize(AtmData.cols());

Eigen::Map<Eigen::RowVectorXd>(&locs[0], AtmData.cols()) = AtmData.row(locn);
Eigen::Map<Eigen::RowVectorXd>(&Tvalues[0], AtmData.cols()) = AtmData.row(iTemp);
Eigen::Map<Eigen::RowVectorXd>(&Uvalues[0], AtmData.cols()) = AtmData.row(iU);


sim.setInitialGuess("T", locs, Tvalues);
sim.setInitialGuess("velocity", locs, Uvalues);
flow->setBounds(flow->componentIndex("T"), 100.0, 2500.0);

for (int ij = 0; ij < gas->nSpecies(); ++ij) {
    int index = flow->componentIndex(gas->speciesName(ij)); // Get index for each species
    flow->setBounds(index, 0.0, 1.0); // Set bounds for each species
    Eigen::Map<Eigen::RowVectorXd>(&Xvalues[0], AtmData.cols()) = ChemMassFrac.row(ij);
    sim.setInitialGuess(gas->speciesName(ij), locs, Xvalues);
}

//gas_kin->disableAllReactions();
flow->solveEnergyEqn();
sim.fixedTemperature();
std::cout << "Transport conditions have been set up" << std::endl;

//sim.setInitialGuess("auto");
//sim.solve(1, true);

//Setting up the integrator
//CVodesIntegrator integrator;
//integrator.setMethod(BDF_Method);
//integrator.setTolerances(1e-9, 1e-12);
//integrator.setMinStepSize(1e-6);
//integrator.setMaxStepSize(1.0e-5);
double t = 10.0;
//integrator.initialize(t, &Cantera::Func);
//integrator.initialize(flow);

//Time-dependent solution
std::cout << "Solver running...." << std::endl;
dt = 3;
//integrator.integrate(t);


//while (t < Tmax) {
//      t +=dt;
 //     integrator.integrate(t);
 //     std::cout << "Time: " << t << "s\n";
       
//       }


//double final_time = integrator.currentTime();
std::cout << "Solution complete" << std::endl;




}
