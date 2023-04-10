from pylab import *
from glob import glob
from scipy.interpolate import interp1d
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import re, matplotlib
from scipy.interpolate import UnivariateSpline
import netCDF4 as nc


data = nc.Dataset('build/bin/Venus_test.nc', 'r' ,format = 'NETCDF4')

z = (data['Altitude'][:])/1000
x1f = data['SO2'][:]
x2f = data['H2O'][:]


plot(x1f*1E6, z, c = 'r')
plot(x2f*1E6, z, c = 'b')
xlabel('X (ppm)')
ylabel('Altitude (km)')
#xscale('log')
ylim([58, 112])
legend(['SO$_{2}$', 'H$_{2}$O'], loc = 'best')
savefig('Venus_test.png')




