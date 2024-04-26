from pylab import *
from glob import glob
from scipy.interpolate import interp1d
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import re, matplotlib
from scipy.interpolate import UnivariateSpline
import netCDF4 as nc
import pandas as pd

# ta(time, plev, latitude)
# plev(plev)
# zg(time, plev, latitude)
data = nc.Dataset("cira.nc", "r", format="NETCDF4")
temp = data["ta"][0, :, 16]  # K
pres = data["plev"][:]  # mbar
z = data["zg"][0, :, 16] / 1000  # m to km

"""
print("Temp")
print(len(temp))
print("Press")
print(len(pres))
print("z")
print(len(z))
"""
# Eddy diffusion data from (Hu et al., 2012), Earth model
fname = "Earth_Eddy_Diffusion_Hu2012.csv"
kdata = pd.read_csv(fname, usecols=[0, 1], names=["kzz", "z"])
z_kdata = kdata["z"]
kzz_kdata = kdata["kzz"]
# print(z_kdata, kzz_kdata)


# make a standard atmosphere from 0 to 80 km
zgrid = linspace(0.0, 100, 50)
pfunc = interp1d(z, pres)
tfunc = interp1d(z, temp)
kzzfunc = interp1d(z_kdata, kzz_kdata)
pnew = pfunc(zgrid)
tnew = tfunc(zgrid)
kzznew = kzzfunc(zgrid)
inx = len(kzznew)
Ndnew = 1e-6 * pnew * 1e2 / (1.38e-23 * tnew)
for i in range(inx):
    print(
        pnew[inx - i - 1],
        tnew[inx - i - 1],
        kzznew[inx - i - 1],
        zgrid[inx - i - 1],
        Ndnew[inx - i - 1],
    )
