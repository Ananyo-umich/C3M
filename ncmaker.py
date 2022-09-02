#The program writes the initial conditions for C3M 

from pylab import *
from glob import glob
import multiprocessing as mp
import concurrent.futures
import re
import netCDF4 as nc


dim = 0

#Making the box model initial condition file

#Making the one-dimensional model initial and boundary conditions
if dim == 1:
  AtmosProfile = 'Venus.txt'
  ChemProfile = 'chemical_data.txt'
  T 
  P
  Kzz
  
