from pylab import *
import pandas as pd
import netCDF4 as nc
from scipy.interpolate import interp1d

"""Photochemical cross section from VULCAN database"""

data_file = "/data4/ananyo/models/C3M/data/VULCAN/N2/N2_cross.csv"

data = pd.read_csv(
    data_file, skiprows=2, usecols=[0, 1, 2, 3], names=["wav", "diss", "abs", "ion"]
)

w = data["wav"][:] * 1e-9
ab = data["diss"][:] * 1e-4


"""Photochemical cross section from Caltech/JPL KINETICS database"""


"""Input text file for planet from C3M"""
planet_file = "Venus_Bhattacharya.txt"

"""
SO2 -> 150 ppmv
H2O -> 40 ppmv
CO2 -> 0.965
N2 -> 0.035
H2S -> 1 ppmv
"""

data = genfromtxt(planet_file)
x = 0.035
P = data[2:-1, 0] * 100
T = data[2:-1, 1]
Kzz = data[2:-1, 2]
H = data[2:-1, 3] * 1000
print(P)

"""Computing the optical depth"""
kb = 1.38e-23
N = x * P / (kb * T)
inx = len(N)
nn, abb = meshgrid(H, ab)
OD = zeros(len(ab))
Odepth = zeros([inx, len(ab)])
for i in range(1, inx):
    Odepth[i, :] = Odepth[i, :] + (N[i] * ab * (H[i - 1] - H[i]))

contourf(w / 1e-9, (H - H[-1]) / 1000, log10(Odepth))
colorbar()
ylabel("Height (km)")
xlabel("wavlength (nm)")
savefig("photo_peak_Venus_N2.png")
