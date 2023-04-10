// C/C++ headers
#include <iostream>
#include <fstream>

// C3M headers
#include "PhotoChemistry.hpp"
using Eigen::MatrixXd;
using Eigen::VectorXd;

//NetCDF
#if NETCDFOUTPUT
  #include <netcdf.h>
#endif



//Print wavelengths
void printWavelength(Eigen::VectorXd wavelengths_)
{
  for (int i = 0; i < wavelengths_.size(); ++i) {
    std::cout << wavelengths_[i] << std::endl;
  }
}

//Print Cross Sections
void printCrossSection(Eigen::VectorXd crossSection_)
{

for (int i = 0; i < crossSection_.size(); ++i) {
    std::cout << crossSection_[i] << std::endl;
  }
}

// Calculate the photochemical reaction rate based on wavelength, cross section and actinic flux
double PhotoChemRate(Eigen::MatrixXd wavelengths_, Eigen::MatrixXd crossSection_, Eigen::MatrixXd Spectral_radiance)
{
int size = wavelengths_.size();
double sum = 0;
double sum1 = 0;
double pi = 3.141592;
double h = 6.626e-34;
double c = 3e8;
Spectral_radiance = (Spectral_radiance.array()*wavelengths_.array().transpose()/(h*c)).matrix();
MatrixXd integrd = (crossSection_.array()*Spectral_radiance.array()).matrix();
for(int  i = 0; i < size-1; i++){

  sum = sum + (integrd(i)*(wavelengths_(i+1) - wavelengths_(i)));
  sum1 = sum1 + (integrd(i+1)*(wavelengths_(i+1) - wavelengths_(i)));
}
sum = (sum + sum1)/2;
sum = sum*2*pi;
return sum;
}


// Function to read the photoionization cross sections from VULCAN database
// The output file will contain wavelength (nm) and photoionization cross section (cm^2)
Eigen::MatrixXd  ReadVULCANPhotoIonCrossSection(string VULCAN_ID){
  fstream InFile;
  InFile.open(VULCAN_ID); 
  string wavlength;
  string photoabs;
  string photodiss;
  string photoion; 
  int num = 0;
  int rows = 0;
  while (getline(InFile, wavlength))
  rows++;
  InFile.close();
  
  Eigen::MatrixXd Output(2, rows-1);
  InFile.open(VULCAN_ID); 
  getline(InFile,wavlength);
  while(getline(InFile,wavlength, ',')){
  getline(InFile,photoabs,',');
  getline(InFile,photodiss,',');
  getline(InFile,photoion,'\n');
  
  
  Output(0, num) = atof(wavlength.c_str())*1E-9; //nm to m
  Output(1, num) = atof(photoion.c_str())*1E-4; //cm^2 to m^2
  
  num++;
  }
  InFile.close();
  
  return Output;
}


// Function to read the photodissociation cross sections from VULCAN database
// The output file will contain wavelength (nm) and photodissociation cross section (cm^2)
Eigen::MatrixXd  ReadVULCANPhotoDissCrossSection(string VULCAN_ID){
  fstream InFile;
  InFile.open(VULCAN_ID); 
  string wavlength;
  string photoabs;
  string photodiss;
  string photoion; 
  int num = 0;
  int rows = 0;
  while (getline(InFile, wavlength))
  rows++;
  InFile.close();
  
  Eigen::MatrixXd Output(2, rows-1);
  InFile.open(VULCAN_ID); 
  getline(InFile,wavlength);
  while(getline(InFile,wavlength, ',')){
  getline(InFile,photoabs,',');
  getline(InFile,photodiss,',');
  getline(InFile,photoion,'\n');
  
  
  Output(0, num) = atof(wavlength.c_str())*1E-9; //nm to m
  Output(1, num) = atof(photodiss.c_str())*1E-4; //cm^2 to m^2
  
  num++;
  }
  InFile.close();
  return Output;
}

// Function to read cross section from VPL photochemical database
// The output file contains the wavelength (Angstroms) and cross section (cm^2)
Eigen::MatrixXd ReadAtmosCrossSection(string AtmosFileName){
  fstream InFile;
  InFile.open(AtmosFileName); 
  string wavlength;
  string photoabs;
  int num = 0;
  int rows = 0;
  while (getline(InFile, wavlength))
  rows++;
  InFile.close();
  
  
  Eigen::MatrixXd Output(2, rows-4);
  InFile.open(AtmosFileName); 
  getline(InFile,wavlength);
  getline(InFile,wavlength);
  getline(InFile,wavlength);
  getline(InFile,wavlength);
  while(InFile >> wavlength >> photoabs){
  
  Output(0, num) = atof(wavlength.c_str())*1E-10; //Angstrom to m
  Output(1, num) = atof(photoabs.c_str())*1E-4; //cm^2 to m^2

  
  num++;
  }
  InFile.close();
  return Output;
  
  
}


Eigen::MatrixXd ReadMPCrossSection(string MPFileName){
  fstream InFile;
  InFile.open(MPFileName);
  string wavlength;
  string photoabs;
  int num = 0;
  int rows = 0;
  while (getline(InFile, wavlength))
  rows++;
  InFile.close();


  Eigen::MatrixXd Output(2, rows);
  InFile.open(MPFileName);
  while(InFile >> wavlength >> photoabs){

  Output(0, num) = atof(wavlength.c_str())*1E-10; 
  Output(1, num) = atof(photoabs.c_str())*1E-4; //cm^2 to m^2

  num++;
  }
  InFile.close();
  return Output;



}

Eigen::MatrixXd ReadKINETICSCrossSection(string SpeciesName){
#if NETCDFOUTPUT
  Eigen::MatrixXd Output(2, 10);
  //int fileid, dimid, varid, err;
  //char tname[80];
  //fname = "KINETICS7_Bhattacharya.nc";
//Reading the netCDF file
  //nc_open(fname.c_str(), NC_NETCDF4, &fileid);


//Defining the dimensions and lengths of dimensions



//Use the reaction number to get the wavelength and cross section
  /*nc_inq_dimid(fileid, "Wavenumber", &dimid);
  nc_inq_dimlen(fileid, dimid, (size_t*)len_);
  nc_inq_dimid(fileid, "Pressure", &dimid);
  nc_inq_dimlen(fileid, dimid, (size_t*)(len_ + 1));
  strcpy(tname, "T_");
  strcat(tname, name_.c_str());
  nc_inq_dimid(fileid, tname, &dimid);
  nc_inq_dimlen(fileid, dimid, (size_t*)(len_ + 2));

  axis_.resize(len_[0] + len_[1] + len_[2]);

  nc_inq_varid(fileid, "Wavenumber", &varid);
  nc_get_var_double(fileid, varid, axis_.data());
  err = nc_inq_varid(fileid, "Pressure", &varid);
  if (err != NC_NOERR)
    throw std::runtime_error(nc_strerror(err));
  err = nc_get_var_double(fileid, varid, axis_.data() + len_[0]);
  if (err != NC_NOERR)
    throw std::runtime_error(nc_strerror(err));
  nc_inq_varid(fileid, tname, &varid);
  nc_get_var_double(fileid, varid, axis_.data() + len_[0] + len_[1]);

  Real *temp = new Real[len_[1]];
  nc_inq_varid(fileid, "Temperature", &varid);
  nc_get_var_double(fileid, varid, temp);

  refatm_.NewAthenaArray(NHYDRO, len_[1]);
  for (int i = 0; i < len_[1]; i++) {
    refatm_(ipr,i) = axis_[len_[0] + i];
    refatm_(idn,i)  = temp[i];
  }

  kcoeff_.resize(len_[0]*len_[1]*len_[2]);
  nc_inq_varid(fileid, name_.c_str(), &varid);
  nc_get_var_double(fileid, varid, kcoeff_.data());
  nc_close(fileid);
  delete[] temp;

  */
#endif
  return Output;

}

