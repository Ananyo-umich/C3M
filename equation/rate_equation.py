#The code includes functions for calculating rate coefficients from photochemistry databases from JPL, VULCAN and SWRI
from pylab import *
from glob import glob
import multiprocessing as mp
import concurrent.futures
import re
from scipy import integrate


def JPL_cross(FileName):
  File = sort(glob('data/cross/'+FileName))
  with open(File[0]) as f:
    lines = f.readlines()[4:]
  wavelength, cross = [],[]
  for line in lines:
    line_int = re.findall('[+0-9.Ee-]+',line)
    wav = float(line_int[0])
    csc = float(line_int[1]) 
    if(csc <= 1): 
     wavelength.append(wav)
     cross.append(csc)  
#Integration scheme for determining cross section
  print(cross)
  print(File[0])
  cross_sec = integrate.simps(cross, wavelength)
  return cross_sec*1e-40


def VULCAN_cross(CrossFileName,CrossCol,BranchFileName,BranchCol):
  File = sort(glob('data/photo_cross/'+CrossFileName))
  with open(File[0]) as f:
    lines = f.readlines()[1:]
  wavelength, cross = [],[]
  for line in lines:
    line_int = re.findall('[+0-9.Ee-]+',line)
    wav = float(line_int[0]) #nm
    csc = float(line_int[int(CrossCol)]) #cm^2
    wavelength.append(wav)
    cross.append(csc)
  cross_sec = integrate.simps(cross, wavelength)
  return cross_sec
  


def ChemicalReactionRateCoefficient(coeff): 
  Flag = coeff[2]
  if(Flag == 'p_JPL'): #JPL photochemical database
    rate = str(JPL_cross(str(coeff[3])))
  if(Flag == 'p_SWRI'): #SWRI photochemical database
    rate = str(CalcCross_Section(str(coeff[3])))
  if(Flag == 'p_VULCAN'): #SWRI photochemical database
    rate = str(VULCAN_cross(str(coeff[3]), str(coeff[4]),str(coeff[5]),str(coeff[6])))
  if(Flag == 't'): #Thermochemical reactions
    rate = coeff[3]+"*pow("+coeff[4]+"/T,"+coeff[5]+")*exp("+coeff[6]+"/T)"
  if(Flag == 't0'): #Thermochemical reaction rate
    rate = coeff[3]
  if(Flag == 'm0'): #Pressure dependent reactions type 0
    rate = coeff[3]+"*N"
  if(Flag == 'm1'): #Pressure dependent reactions type 1
    rate = coeff[3]+"*pow("+coeff[4]+"/T,"+coeff[5]+")*exp("+coeff[6]+"/T)*N"
  if(Flag == 'm2'): #Pressure dependent reactions type 2
    k0 = "("+coeff[3]+"*pow("+coeff[4]+"/T,"+coeff[5]+")*exp("+coeff[6]+"/T))"
    kinf = "("+coeff[7]+"*pow("+coeff[8]+"/T,"+coeff[9]+")*exp("+coeff[10]+"/T))"
    rate = "("+k0+"*N*"+kinf+")/(("+k0+"*N) +" + kinf + ")"
  if(Flag == 'm3'): #Pressure dependent reactions type 3
    k0 = "("+coeff[3]+"*pow("+coeff[4]+"/T,"+coeff[5]+")*exp("+coeff[6]+"/T)) + " + "("+coeff[7]+"*pow("+coeff[8]+"/T,"+coeff[9]+")*exp("+coeff[10]+"/T))"
    kinf = "("+coeff[11]+"*pow("+coeff[12]+"/T,"+coeff[13]+")*exp("+coeff[14]+"/T))" 
    rate = "("+k0+"*N*"+kinf+")/(("+k0+"*N) +" + kinf + ")"
  return rate