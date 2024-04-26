from pylab import *
import pandas as pd
import netCDF4 as nc
from scipy.interpolate import interp1d

data_file = "/data4/ananyo/models/C3M/data/VULCAN/O2/O2_cross.csv"

data = pd.read_csv(
    data_file, skiprows=2, usecols=[0, 1, 2, 3], names=["wav", "diss", "abs", "ion"]
)

w = data["wav"][:] / 1e-9
ab = data["abs"][:] * 1e-4

planet_file = "cira.nc"

data = nc.Dataset("cira.nc", "r", format="NETCDF4")
temp = data["ta"][0, :, 16]  # K
pres = data["plev"][:] * 100  # mbar
z = data["zg"][0, :, 16] / 1000  # m to km

# make a standard atmosphere from 0 to 80 km
zgrid = linspace(0.0, 100, 1000)
dz = 100 / 1000
pfunc = interp1d(z, pres)
tfunc = interp1d(z, temp)
pnew = pfunc(zgrid)
tnew = tfunc(zgrid)
kb = 1.38e-23
N = pnew / (kb * tnew)
print(len(ab))
inx = len(N)
nn, abb = meshgrid(zgrid, ab)
OD = zeros(len(ab))
Odepth = zeros([len(ab), inx])
for i in range(inx):
    Odepth[:, i] = OD + (N[i] * ab * dz)
    OD = OD + (N[i] * ab * dz)

contourf(flip(zgrid), w, log10(Odepth))
colorbar()
xlabel("Height")
ylabel("wav")
savefig("OD.png")
