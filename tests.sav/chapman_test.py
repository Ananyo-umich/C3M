#Chapman test for C3M box model
#The reference pressure and temperature correspond to 30 km altitude from Seinfield and Pandis
from pylab import *
from glob import glob
import re
import os
import netCDF4 as nc



#Reading NETCDF output file
data = nc.Dataset('chapman.nc', 'r' ,format = 'NETCDF4')
Oc = data['O'][:]
O2c = data['O2'][:] 
O3c = data['O3'][:]
O_1Dc = data['O_1D'][:]
inx = len(Oc)
#Standard solution of Chapman test
#At 30 km
standard = 3e-5 #O/O3 ratio at 30 km
model = Oc[inx-1]/O3c[inx-1]



#Error analysis
err = 100*(model-standard)/standard

print("Error for chapman process")
print("err ="+str(err) +  " %")