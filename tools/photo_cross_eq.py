from pylab import *
import pandas as pd
import netCDF4 as nc
from scipy.interpolate import interp1d

"""Photochemical cross section from VULCAN database"""

data_file1 = "/data4/ananyo/models/C3M/data/VULCAN/N2/N2_cross.csv"
data_file2 = "/data4/ananyo/models/C3M/data/VULCAN/CO2/CO2_cross.csv"
data_file3 = "/data4/ananyo/models/C3M/data/VULCAN/SO2/SO2_cross.csv"
data_file4 = "/data4/ananyo/models/C3M/data/VULCAN/H2O/H2O_cross.csv"
data_file5 = "/data4/ananyo/models/C3M/data/VULCAN/COS/COS_cross.csv"
data_file6 = "/data4/ananyo/models/C3M/data/VULCAN/S2/S2_cross.csv"

data1 = pd.read_csv(
    data_file1, skiprows=2, usecols=[0, 1, 2, 3], names=["wav", "diss", "abs", "ion"]
)
data2 = pd.read_csv(
    data_file2, skiprows=2, usecols=[0, 1, 2, 3], names=["wav", "diss", "abs", "ion"]
)
data3 = pd.read_csv(
    data_file3, skiprows=2, usecols=[0, 1, 2, 3], names=["wav", "diss", "abs", "ion"]
)
data4 = pd.read_csv(
    data_file4, skiprows=2, usecols=[0, 1, 2, 3], names=["wav", "diss", "abs", "ion"]
)
data5 = pd.read_csv(
    data_file5, skiprows=2, usecols=[0, 1, 2, 3], names=["wav", "diss", "abs", "ion"]
)
data6 = pd.read_csv(
    data_file6, skiprows=2, usecols=[0, 1, 2, 3], names=["wav", "diss", "abs", "ion"]
)

w1 = data1["wav"][:]
ab1 = data1["abs"][:]

w2 = data2["wav"][:]
ab2 = data2["abs"][:]

w3 = data3["wav"][:]
ab3 = data3["abs"][:]

w4 = data4["wav"][:]
ab4 = data4["abs"][:]

w5 = data5["wav"][:]
ab5 = data5["abs"][:]

w6 = data6["wav"][:]
ab6 = data6["abs"][:]

plot(w1, ab1, "k-")
plot(w2, ab2, "r-")
plot(w3, ab3, "m-")
plot(w4, ab4, "g-")
plot(w5, ab5, "b-")
plot(w6, ab6, "C7-")
yscale("log")
ylim([1e-30, 1e-15])

xlim([0.0, 400])
legend(["N$_{2}$", "CO$_{2}$", "SO$_{2}$", "H$_{2}$O", "COS", "S$_{2}$"], loc="best")
xlabel("Wavelength (nm)")
ylabel("Cross section (cm$^{2}$)")
savefig("photo_cross_compile_abs.png")
