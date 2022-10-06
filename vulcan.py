from pylab import *
import netCDF4 as nc
from glob import glob
import cantera as ct


T = 200.0
P = 100000.0
#Creatng chemistry network
gas = ct.Solution("vulcan.yaml")
gas.TPX = 1000, 32*ct.one_atm, 'CH4:0.1,H2:1,CO2:10,H2O:10'
ct.use_legacy_rate_constants(False)

dt = 1e-5
Nsteps = 1000
for nt in range(Nsteps):
  Xi = gas.net_production_rates
  phi = gas.net_production_rates_ddP
  idt = identity(len(Xi))
  A = idt - phi
  B = linalg.inv(A)*Xi
  gas.X = gas.X + B

