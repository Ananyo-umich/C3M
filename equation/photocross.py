#The program computes the cross sections for photochemical reactions of neutral species
from pylab import *
import re
from glob import glob
from scipy import integrate

def CalcCross_Section(FileName):
  with open(FileName) as f:
    lines = f.readlines()[4:]
    
  wavelength, cross = [],[]
  for line in lines:
    wav = line.split(" ")[0]  
    csc = line.split(" ")[1]  
    wavelength.append(wav)
    cross.append(csc)
    
#Integration scheme for determining cross section
  cross_sec = integrate.trapezoid(cross, wavelength)
  return cross_sec