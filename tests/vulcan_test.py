from pylab import *
from glob import glob
import re
from equation import rate_equation
from config import *
from glob import glob
import os
import netCDF4 as nc



#Reading NETCDF output file
data = nc.Dataset('vulcan.nc', 'r' ,format = 'NETCDF4')
CH4c = data['CH4'][:]
H2Oc = data['H2O'][:]
CO2c = data['CO2'][:]
inx = len(CH4c)
#Standard solution of VULCAN test
std_C_O = 1
std_O_H = 3
model_C_O = 1
model_O_H = 3

#Error analysis
C_O_error = (model_C_O - std_C_O)/std_C_O
O_H_error = (model_O_H - std_O_H)/std_O_H
print("Error for VULCAN C-H-O scheme")
print("C/O = " + str(C_O_error))
print("O/H = " + str(O_H_error))


