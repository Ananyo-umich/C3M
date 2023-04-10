from pylab import *
import netCDF4 as nc
from glob import glob
import cantera as ct


T = 200.0
P = 100000.0
#Creatng chemistry network
gas = ct.Solution("vulcanCHOnetwork.yaml")
