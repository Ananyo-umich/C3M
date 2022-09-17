from pylab import * 
import netCDF4 as nc 
from glob import glob


T = 200.0
P = 100000.0
Kb = 1.38e-23
N = P/(Kb*T) 

def Jacobian(x):
	xjac = zeros([57,57])
	xjac[0,0] = 	xjac[0,0] - (3.2e-11*pow(300/T,0)*exp(70/T)*1.0*pow(x[0],0.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,0] = 	xjac[0,0] - (1.8e-11*pow(300/T,0)*exp(110/T)*1.0*pow(x[0],0.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,0] = 	xjac[0,0] - (7.4e-11*pow(300/T,0)*exp(120/T)*1.0*pow(x[0],0.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,0] = 	xjac[0,0] - (1.2e-10*1.0*pow(x[0],0.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,0] = 	xjac[0,0] - (1.2e-10*1.0*pow(x[0],0.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,0] = 	xjac[0,0] - (1.1e-10*1.0*pow(x[0],0.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,0] = 	xjac[0,0] - (2.2e-10*1.0*pow(x[0],0.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,0] = 	xjac[0,0] - (1.55e-10*1.0*pow(x[0],0.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,0] = 	xjac[0,0] - (5.25e-11*1.0*pow(x[0],0.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,0] = 	xjac[0,0] - (1e-10*1.0*pow(x[0],0.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,0] = 	xjac[0,0] - (3.6e-11*1.0*pow(x[0],0.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,0] = 	xjac[0,0] - (1.35e-11*1.0*pow(x[0],0.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,0] = 	xjac[0,0] - (3.6e-10*1.0*pow(x[0],0.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,0] = 	xjac[0,0] - (3.6e-10*1.0*pow(x[0],0.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,1] = 	xjac[0,1] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*1.0*pow(x[1],0.0))#O_1D + O2 = O + O2
	xjac[0,1] = 	xjac[0,1] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,1] = 	xjac[0,1] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,1] = 	xjac[0,1] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,1] = 	xjac[0,1] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,1] = 	xjac[0,1] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,1] = 	xjac[0,1] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,1] = 	xjac[0,1] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,1] = 	xjac[0,1] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,1] = 	xjac[0,1] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,1] = 	xjac[0,1] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,1] = 	xjac[0,1] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,1] = 	xjac[0,1] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,1] = 	xjac[0,1] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,2] = 	xjac[0,2] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,2] = 	xjac[0,2] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,2] = 	xjac[0,2] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,2] = 	xjac[0,2] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,2] = 	xjac[0,2] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,2] = 	xjac[0,2] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,2] = 	xjac[0,2] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,2] = 	xjac[0,2] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,2] = 	xjac[0,2] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,2] = 	xjac[0,2] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,2] = 	xjac[0,2] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,2] = 	xjac[0,2] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,2] = 	xjac[0,2] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,2] = 	xjac[0,2] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,3] = 	xjac[0,3] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,3] = 	xjac[0,3] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*1.0*pow(x[3],0.0))#O_1D + N2 = O + N2
	xjac[0,3] = 	xjac[0,3] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,3] = 	xjac[0,3] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,3] = 	xjac[0,3] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,3] = 	xjac[0,3] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,3] = 	xjac[0,3] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,3] = 	xjac[0,3] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,3] = 	xjac[0,3] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,3] = 	xjac[0,3] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,3] = 	xjac[0,3] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,3] = 	xjac[0,3] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,3] = 	xjac[0,3] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,3] = 	xjac[0,3] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,4] = 	xjac[0,4] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,4] = 	xjac[0,4] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,4] = 	xjac[0,4] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*1.0*pow(x[4],0.0))#O_1D + CO2 = O + CO2
	xjac[0,4] = 	xjac[0,4] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,4] = 	xjac[0,4] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,4] = 	xjac[0,4] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,4] = 	xjac[0,4] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,4] = 	xjac[0,4] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,4] = 	xjac[0,4] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,4] = 	xjac[0,4] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,4] = 	xjac[0,4] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,4] = 	xjac[0,4] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,4] = 	xjac[0,4] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,4] = 	xjac[0,4] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,5] = 	xjac[0,5] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,5] = 	xjac[0,5] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,5] = 	xjac[0,5] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,5] = 	xjac[0,5] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,5] = 	xjac[0,5] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,5] = 	xjac[0,5] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,5] = 	xjac[0,5] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,5] = 	xjac[0,5] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,5] = 	xjac[0,5] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,5] = 	xjac[0,5] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,5] = 	xjac[0,5] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,5] = 	xjac[0,5] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,5] = 	xjac[0,5] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,5] = 	xjac[0,5] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,6] = 	xjac[0,6] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,6] = 	xjac[0,6] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,6] = 	xjac[0,6] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,6] = 	xjac[0,6] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,6] = 	xjac[0,6] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,6] = 	xjac[0,6] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,6] = 	xjac[0,6] - (2.2e-10*pow(x[0],1.0)*1.0*pow(x[6],0.0))#O_1D + H2O = 2OH
	xjac[0,6] = 	xjac[0,6] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,6] = 	xjac[0,6] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,6] = 	xjac[0,6] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,6] = 	xjac[0,6] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,6] = 	xjac[0,6] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,6] = 	xjac[0,6] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,6] = 	xjac[0,6] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,7] = 	xjac[0,7] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,7] = 	xjac[0,7] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,7] = 	xjac[0,7] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,7] = 	xjac[0,7] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,7] = 	xjac[0,7] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,7] = 	xjac[0,7] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,7] = 	xjac[0,7] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,7] = 	xjac[0,7] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,7] = 	xjac[0,7] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,7] = 	xjac[0,7] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,7] = 	xjac[0,7] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,7] = 	xjac[0,7] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,7] = 	xjac[0,7] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,7] = 	xjac[0,7] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,8] = 	xjac[0,8] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,8] = 	xjac[0,8] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,8] = 	xjac[0,8] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,8] = 	xjac[0,8] - (1.2e-10*pow(x[0],1.0)*1.0*pow(x[8],0.0))#O_1D + O3 = 2O2
	xjac[0,8] = 	xjac[0,8] - (1.2e-10*pow(x[0],1.0)*1.0*pow(x[8],0.0))#O_1D + O3 = 2O + O2
	xjac[0,8] = 	xjac[0,8] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,8] = 	xjac[0,8] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,8] = 	xjac[0,8] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,8] = 	xjac[0,8] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,8] = 	xjac[0,8] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,8] = 	xjac[0,8] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,8] = 	xjac[0,8] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,8] = 	xjac[0,8] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,8] = 	xjac[0,8] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,9] = 	xjac[0,9] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,9] = 	xjac[0,9] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,9] = 	xjac[0,9] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,9] = 	xjac[0,9] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,9] = 	xjac[0,9] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,9] = 	xjac[0,9] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,9] = 	xjac[0,9] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,9] = 	xjac[0,9] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,9] = 	xjac[0,9] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,9] = 	xjac[0,9] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,9] = 	xjac[0,9] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,9] = 	xjac[0,9] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,9] = 	xjac[0,9] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,9] = 	xjac[0,9] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,10] = 	xjac[0,10] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,10] = 	xjac[0,10] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,10] = 	xjac[0,10] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,10] = 	xjac[0,10] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,10] = 	xjac[0,10] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,10] = 	xjac[0,10] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,10] = 	xjac[0,10] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,10] = 	xjac[0,10] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,10] = 	xjac[0,10] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,10] = 	xjac[0,10] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,10] = 	xjac[0,10] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,10] = 	xjac[0,10] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,10] = 	xjac[0,10] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,10] = 	xjac[0,10] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,11] = 	xjac[0,11] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,11] = 	xjac[0,11] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,11] = 	xjac[0,11] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,11] = 	xjac[0,11] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,11] = 	xjac[0,11] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,11] = 	xjac[0,11] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,11] = 	xjac[0,11] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,11] = 	xjac[0,11] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,11] = 	xjac[0,11] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,11] = 	xjac[0,11] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,11] = 	xjac[0,11] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,11] = 	xjac[0,11] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,11] = 	xjac[0,11] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,11] = 	xjac[0,11] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,12] = 	xjac[0,12] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,12] = 	xjac[0,12] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,12] = 	xjac[0,12] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,12] = 	xjac[0,12] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,12] = 	xjac[0,12] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,12] = 	xjac[0,12] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,12] = 	xjac[0,12] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,12] = 	xjac[0,12] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,12] = 	xjac[0,12] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,12] = 	xjac[0,12] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,12] = 	xjac[0,12] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,12] = 	xjac[0,12] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,12] = 	xjac[0,12] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,12] = 	xjac[0,12] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,13] = 	xjac[0,13] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,13] = 	xjac[0,13] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,13] = 	xjac[0,13] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,13] = 	xjac[0,13] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,13] = 	xjac[0,13] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,13] = 	xjac[0,13] - (1.1e-10*pow(x[0],1.0)*1.0*pow(x[13],0.0))#O_1D + H2 = H + OH
	xjac[0,13] = 	xjac[0,13] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,13] = 	xjac[0,13] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,13] = 	xjac[0,13] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,13] = 	xjac[0,13] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,13] = 	xjac[0,13] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,13] = 	xjac[0,13] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,13] = 	xjac[0,13] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,13] = 	xjac[0,13] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,14] = 	xjac[0,14] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,14] = 	xjac[0,14] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,14] = 	xjac[0,14] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,14] = 	xjac[0,14] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,14] = 	xjac[0,14] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,14] = 	xjac[0,14] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,14] = 	xjac[0,14] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,14] = 	xjac[0,14] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,14] = 	xjac[0,14] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,14] = 	xjac[0,14] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,14] = 	xjac[0,14] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,14] = 	xjac[0,14] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,14] = 	xjac[0,14] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,14] = 	xjac[0,14] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,15] = 	xjac[0,15] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,15] = 	xjac[0,15] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,15] = 	xjac[0,15] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,15] = 	xjac[0,15] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,15] = 	xjac[0,15] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,15] = 	xjac[0,15] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,15] = 	xjac[0,15] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,15] = 	xjac[0,15] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,15] = 	xjac[0,15] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,15] = 	xjac[0,15] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,15] = 	xjac[0,15] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,15] = 	xjac[0,15] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,15] = 	xjac[0,15] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,15] = 	xjac[0,15] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,16] = 	xjac[0,16] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,16] = 	xjac[0,16] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,16] = 	xjac[0,16] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,16] = 	xjac[0,16] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,16] = 	xjac[0,16] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,16] = 	xjac[0,16] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,16] = 	xjac[0,16] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,16] = 	xjac[0,16] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,16] = 	xjac[0,16] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,16] = 	xjac[0,16] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,16] = 	xjac[0,16] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,16] = 	xjac[0,16] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,16] = 	xjac[0,16] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,16] = 	xjac[0,16] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,17] = 	xjac[0,17] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,17] = 	xjac[0,17] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,17] = 	xjac[0,17] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,17] = 	xjac[0,17] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,17] = 	xjac[0,17] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,17] = 	xjac[0,17] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,17] = 	xjac[0,17] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,17] = 	xjac[0,17] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,17] = 	xjac[0,17] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,17] = 	xjac[0,17] - (1e-10*pow(x[0],1.0)*1.0*pow(x[17],0.0))#O_1D + HCl = Cl + OH
	xjac[0,17] = 	xjac[0,17] - (3.6e-11*pow(x[0],1.0)*1.0*pow(x[17],0.0))#O_1D + HCl = ClO + H
	xjac[0,17] = 	xjac[0,17] - (1.35e-11*pow(x[0],1.0)*1.0*pow(x[17],0.0))#O_1D + HCl = O + HCl
	xjac[0,17] = 	xjac[0,17] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,17] = 	xjac[0,17] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,18] = 	xjac[0,18] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,18] = 	xjac[0,18] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,18] = 	xjac[0,18] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,18] = 	xjac[0,18] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,18] = 	xjac[0,18] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,18] = 	xjac[0,18] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,18] = 	xjac[0,18] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,18] = 	xjac[0,18] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,18] = 	xjac[0,18] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,18] = 	xjac[0,18] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,18] = 	xjac[0,18] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,18] = 	xjac[0,18] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,18] = 	xjac[0,18] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,18] = 	xjac[0,18] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,19] = 	xjac[0,19] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,19] = 	xjac[0,19] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,19] = 	xjac[0,19] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,19] = 	xjac[0,19] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,19] = 	xjac[0,19] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,19] = 	xjac[0,19] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,19] = 	xjac[0,19] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,19] = 	xjac[0,19] - (1.55e-10*pow(x[0],1.0)*1.0*pow(x[19],0.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,19] = 	xjac[0,19] - (5.25e-11*pow(x[0],1.0)*1.0*pow(x[19],0.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,19] = 	xjac[0,19] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,19] = 	xjac[0,19] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,19] = 	xjac[0,19] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,19] = 	xjac[0,19] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,19] = 	xjac[0,19] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,20] = 	xjac[0,20] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,20] = 	xjac[0,20] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,20] = 	xjac[0,20] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,20] = 	xjac[0,20] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,20] = 	xjac[0,20] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,20] = 	xjac[0,20] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,20] = 	xjac[0,20] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,20] = 	xjac[0,20] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,20] = 	xjac[0,20] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,20] = 	xjac[0,20] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,20] = 	xjac[0,20] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,20] = 	xjac[0,20] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,20] = 	xjac[0,20] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,20] = 	xjac[0,20] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,21] = 	xjac[0,21] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,21] = 	xjac[0,21] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,21] = 	xjac[0,21] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,21] = 	xjac[0,21] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,21] = 	xjac[0,21] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,21] = 	xjac[0,21] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,21] = 	xjac[0,21] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,21] = 	xjac[0,21] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,21] = 	xjac[0,21] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,21] = 	xjac[0,21] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,21] = 	xjac[0,21] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,21] = 	xjac[0,21] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,21] = 	xjac[0,21] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,21] = 	xjac[0,21] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,22] = 	xjac[0,22] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,22] = 	xjac[0,22] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,22] = 	xjac[0,22] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,22] = 	xjac[0,22] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,22] = 	xjac[0,22] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,22] = 	xjac[0,22] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,22] = 	xjac[0,22] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,22] = 	xjac[0,22] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,22] = 	xjac[0,22] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,22] = 	xjac[0,22] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,22] = 	xjac[0,22] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,22] = 	xjac[0,22] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,22] = 	xjac[0,22] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,22] = 	xjac[0,22] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,23] = 	xjac[0,23] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,23] = 	xjac[0,23] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,23] = 	xjac[0,23] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,23] = 	xjac[0,23] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,23] = 	xjac[0,23] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,23] = 	xjac[0,23] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,23] = 	xjac[0,23] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,23] = 	xjac[0,23] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,23] = 	xjac[0,23] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,23] = 	xjac[0,23] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,23] = 	xjac[0,23] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,23] = 	xjac[0,23] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,23] = 	xjac[0,23] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,23] = 	xjac[0,23] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,24] = 	xjac[0,24] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,24] = 	xjac[0,24] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,24] = 	xjac[0,24] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,24] = 	xjac[0,24] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,24] = 	xjac[0,24] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,24] = 	xjac[0,24] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,24] = 	xjac[0,24] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,24] = 	xjac[0,24] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,24] = 	xjac[0,24] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,24] = 	xjac[0,24] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,24] = 	xjac[0,24] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,24] = 	xjac[0,24] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,24] = 	xjac[0,24] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,24] = 	xjac[0,24] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,25] = 	xjac[0,25] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,25] = 	xjac[0,25] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,25] = 	xjac[0,25] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,25] = 	xjac[0,25] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,25] = 	xjac[0,25] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,25] = 	xjac[0,25] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,25] = 	xjac[0,25] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,25] = 	xjac[0,25] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,25] = 	xjac[0,25] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,25] = 	xjac[0,25] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,25] = 	xjac[0,25] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,25] = 	xjac[0,25] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,25] = 	xjac[0,25] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,25] = 	xjac[0,25] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,26] = 	xjac[0,26] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,26] = 	xjac[0,26] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,26] = 	xjac[0,26] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,26] = 	xjac[0,26] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,26] = 	xjac[0,26] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,26] = 	xjac[0,26] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,26] = 	xjac[0,26] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,26] = 	xjac[0,26] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,26] = 	xjac[0,26] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,26] = 	xjac[0,26] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,26] = 	xjac[0,26] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,26] = 	xjac[0,26] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,26] = 	xjac[0,26] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,26] = 	xjac[0,26] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,27] = 	xjac[0,27] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,27] = 	xjac[0,27] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,27] = 	xjac[0,27] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,27] = 	xjac[0,27] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,27] = 	xjac[0,27] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,27] = 	xjac[0,27] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,27] = 	xjac[0,27] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,27] = 	xjac[0,27] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,27] = 	xjac[0,27] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,27] = 	xjac[0,27] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,27] = 	xjac[0,27] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,27] = 	xjac[0,27] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,27] = 	xjac[0,27] - (3.6e-10*pow(x[0],1.0)*1.0*pow(x[27],0.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,27] = 	xjac[0,27] - (3.6e-10*pow(x[0],1.0)*1.0*pow(x[27],0.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,28] = 	xjac[0,28] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,28] = 	xjac[0,28] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,28] = 	xjac[0,28] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,28] = 	xjac[0,28] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,28] = 	xjac[0,28] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,28] = 	xjac[0,28] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,28] = 	xjac[0,28] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,28] = 	xjac[0,28] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,28] = 	xjac[0,28] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,28] = 	xjac[0,28] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,28] = 	xjac[0,28] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,28] = 	xjac[0,28] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,28] = 	xjac[0,28] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,28] = 	xjac[0,28] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,29] = 	xjac[0,29] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,29] = 	xjac[0,29] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,29] = 	xjac[0,29] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,29] = 	xjac[0,29] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,29] = 	xjac[0,29] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,29] = 	xjac[0,29] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,29] = 	xjac[0,29] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,29] = 	xjac[0,29] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,29] = 	xjac[0,29] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,29] = 	xjac[0,29] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,29] = 	xjac[0,29] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,29] = 	xjac[0,29] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,29] = 	xjac[0,29] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,29] = 	xjac[0,29] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,30] = 	xjac[0,30] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,30] = 	xjac[0,30] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,30] = 	xjac[0,30] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,30] = 	xjac[0,30] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,30] = 	xjac[0,30] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,30] = 	xjac[0,30] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,30] = 	xjac[0,30] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,30] = 	xjac[0,30] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,30] = 	xjac[0,30] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,30] = 	xjac[0,30] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,30] = 	xjac[0,30] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,30] = 	xjac[0,30] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,30] = 	xjac[0,30] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,30] = 	xjac[0,30] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,31] = 	xjac[0,31] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,31] = 	xjac[0,31] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,31] = 	xjac[0,31] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,31] = 	xjac[0,31] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,31] = 	xjac[0,31] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,31] = 	xjac[0,31] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,31] = 	xjac[0,31] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,31] = 	xjac[0,31] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,31] = 	xjac[0,31] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,31] = 	xjac[0,31] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,31] = 	xjac[0,31] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,31] = 	xjac[0,31] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,31] = 	xjac[0,31] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,31] = 	xjac[0,31] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,32] = 	xjac[0,32] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,32] = 	xjac[0,32] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,32] = 	xjac[0,32] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,32] = 	xjac[0,32] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,32] = 	xjac[0,32] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,32] = 	xjac[0,32] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,32] = 	xjac[0,32] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,32] = 	xjac[0,32] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,32] = 	xjac[0,32] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,32] = 	xjac[0,32] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,32] = 	xjac[0,32] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,32] = 	xjac[0,32] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,32] = 	xjac[0,32] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,32] = 	xjac[0,32] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,33] = 	xjac[0,33] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,33] = 	xjac[0,33] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,33] = 	xjac[0,33] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,33] = 	xjac[0,33] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,33] = 	xjac[0,33] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,33] = 	xjac[0,33] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,33] = 	xjac[0,33] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,33] = 	xjac[0,33] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,33] = 	xjac[0,33] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,33] = 	xjac[0,33] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,33] = 	xjac[0,33] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,33] = 	xjac[0,33] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,33] = 	xjac[0,33] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,33] = 	xjac[0,33] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,34] = 	xjac[0,34] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,34] = 	xjac[0,34] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,34] = 	xjac[0,34] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,34] = 	xjac[0,34] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,34] = 	xjac[0,34] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,34] = 	xjac[0,34] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,34] = 	xjac[0,34] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,34] = 	xjac[0,34] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,34] = 	xjac[0,34] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,34] = 	xjac[0,34] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,34] = 	xjac[0,34] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,34] = 	xjac[0,34] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,34] = 	xjac[0,34] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,34] = 	xjac[0,34] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,35] = 	xjac[0,35] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,35] = 	xjac[0,35] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,35] = 	xjac[0,35] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,35] = 	xjac[0,35] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,35] = 	xjac[0,35] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,35] = 	xjac[0,35] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,35] = 	xjac[0,35] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,35] = 	xjac[0,35] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,35] = 	xjac[0,35] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,35] = 	xjac[0,35] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,35] = 	xjac[0,35] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,35] = 	xjac[0,35] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,35] = 	xjac[0,35] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,35] = 	xjac[0,35] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,36] = 	xjac[0,36] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,36] = 	xjac[0,36] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,36] = 	xjac[0,36] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,36] = 	xjac[0,36] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,36] = 	xjac[0,36] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,36] = 	xjac[0,36] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,36] = 	xjac[0,36] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,36] = 	xjac[0,36] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,36] = 	xjac[0,36] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,36] = 	xjac[0,36] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,36] = 	xjac[0,36] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,36] = 	xjac[0,36] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,36] = 	xjac[0,36] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,36] = 	xjac[0,36] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,37] = 	xjac[0,37] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,37] = 	xjac[0,37] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,37] = 	xjac[0,37] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,37] = 	xjac[0,37] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,37] = 	xjac[0,37] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,37] = 	xjac[0,37] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,37] = 	xjac[0,37] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,37] = 	xjac[0,37] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,37] = 	xjac[0,37] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,37] = 	xjac[0,37] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,37] = 	xjac[0,37] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,37] = 	xjac[0,37] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,37] = 	xjac[0,37] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,37] = 	xjac[0,37] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,38] = 	xjac[0,38] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,38] = 	xjac[0,38] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,38] = 	xjac[0,38] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,38] = 	xjac[0,38] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,38] = 	xjac[0,38] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,38] = 	xjac[0,38] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,38] = 	xjac[0,38] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,38] = 	xjac[0,38] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,38] = 	xjac[0,38] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,38] = 	xjac[0,38] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,38] = 	xjac[0,38] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,38] = 	xjac[0,38] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,38] = 	xjac[0,38] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,38] = 	xjac[0,38] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,39] = 	xjac[0,39] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,39] = 	xjac[0,39] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,39] = 	xjac[0,39] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,39] = 	xjac[0,39] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,39] = 	xjac[0,39] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,39] = 	xjac[0,39] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,39] = 	xjac[0,39] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,39] = 	xjac[0,39] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,39] = 	xjac[0,39] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,39] = 	xjac[0,39] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,39] = 	xjac[0,39] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,39] = 	xjac[0,39] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,39] = 	xjac[0,39] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,39] = 	xjac[0,39] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,40] = 	xjac[0,40] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,40] = 	xjac[0,40] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,40] = 	xjac[0,40] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,40] = 	xjac[0,40] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,40] = 	xjac[0,40] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,40] = 	xjac[0,40] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,40] = 	xjac[0,40] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,40] = 	xjac[0,40] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,40] = 	xjac[0,40] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,40] = 	xjac[0,40] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,40] = 	xjac[0,40] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,40] = 	xjac[0,40] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,40] = 	xjac[0,40] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,40] = 	xjac[0,40] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,41] = 	xjac[0,41] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,41] = 	xjac[0,41] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,41] = 	xjac[0,41] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,41] = 	xjac[0,41] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,41] = 	xjac[0,41] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,41] = 	xjac[0,41] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,41] = 	xjac[0,41] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,41] = 	xjac[0,41] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,41] = 	xjac[0,41] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,41] = 	xjac[0,41] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,41] = 	xjac[0,41] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,41] = 	xjac[0,41] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,41] = 	xjac[0,41] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,41] = 	xjac[0,41] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,42] = 	xjac[0,42] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,42] = 	xjac[0,42] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,42] = 	xjac[0,42] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,42] = 	xjac[0,42] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,42] = 	xjac[0,42] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,42] = 	xjac[0,42] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,42] = 	xjac[0,42] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,42] = 	xjac[0,42] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,42] = 	xjac[0,42] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,42] = 	xjac[0,42] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,42] = 	xjac[0,42] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,42] = 	xjac[0,42] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,42] = 	xjac[0,42] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,42] = 	xjac[0,42] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,43] = 	xjac[0,43] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,43] = 	xjac[0,43] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,43] = 	xjac[0,43] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,43] = 	xjac[0,43] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,43] = 	xjac[0,43] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,43] = 	xjac[0,43] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,43] = 	xjac[0,43] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,43] = 	xjac[0,43] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,43] = 	xjac[0,43] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,43] = 	xjac[0,43] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,43] = 	xjac[0,43] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,43] = 	xjac[0,43] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,43] = 	xjac[0,43] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,43] = 	xjac[0,43] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,44] = 	xjac[0,44] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,44] = 	xjac[0,44] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,44] = 	xjac[0,44] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,44] = 	xjac[0,44] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,44] = 	xjac[0,44] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,44] = 	xjac[0,44] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,44] = 	xjac[0,44] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,44] = 	xjac[0,44] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,44] = 	xjac[0,44] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,44] = 	xjac[0,44] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,44] = 	xjac[0,44] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,44] = 	xjac[0,44] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,44] = 	xjac[0,44] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,44] = 	xjac[0,44] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,45] = 	xjac[0,45] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,45] = 	xjac[0,45] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,45] = 	xjac[0,45] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,45] = 	xjac[0,45] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,45] = 	xjac[0,45] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,45] = 	xjac[0,45] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,45] = 	xjac[0,45] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,45] = 	xjac[0,45] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,45] = 	xjac[0,45] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,45] = 	xjac[0,45] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,45] = 	xjac[0,45] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,45] = 	xjac[0,45] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,45] = 	xjac[0,45] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,45] = 	xjac[0,45] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,46] = 	xjac[0,46] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,46] = 	xjac[0,46] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,46] = 	xjac[0,46] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,46] = 	xjac[0,46] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,46] = 	xjac[0,46] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,46] = 	xjac[0,46] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,46] = 	xjac[0,46] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,46] = 	xjac[0,46] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,46] = 	xjac[0,46] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,46] = 	xjac[0,46] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,46] = 	xjac[0,46] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,46] = 	xjac[0,46] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,46] = 	xjac[0,46] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,46] = 	xjac[0,46] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,47] = 	xjac[0,47] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,47] = 	xjac[0,47] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,47] = 	xjac[0,47] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,47] = 	xjac[0,47] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,47] = 	xjac[0,47] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,47] = 	xjac[0,47] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,47] = 	xjac[0,47] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,47] = 	xjac[0,47] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,47] = 	xjac[0,47] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,47] = 	xjac[0,47] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,47] = 	xjac[0,47] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,47] = 	xjac[0,47] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,47] = 	xjac[0,47] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,47] = 	xjac[0,47] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,48] = 	xjac[0,48] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,48] = 	xjac[0,48] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,48] = 	xjac[0,48] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,48] = 	xjac[0,48] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,48] = 	xjac[0,48] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,48] = 	xjac[0,48] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,48] = 	xjac[0,48] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,48] = 	xjac[0,48] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,48] = 	xjac[0,48] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,48] = 	xjac[0,48] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,48] = 	xjac[0,48] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,48] = 	xjac[0,48] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,48] = 	xjac[0,48] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,48] = 	xjac[0,48] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,49] = 	xjac[0,49] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,49] = 	xjac[0,49] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,49] = 	xjac[0,49] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,49] = 	xjac[0,49] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,49] = 	xjac[0,49] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,49] = 	xjac[0,49] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,49] = 	xjac[0,49] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,49] = 	xjac[0,49] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,49] = 	xjac[0,49] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,49] = 	xjac[0,49] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,49] = 	xjac[0,49] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,49] = 	xjac[0,49] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,49] = 	xjac[0,49] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,49] = 	xjac[0,49] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,50] = 	xjac[0,50] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,50] = 	xjac[0,50] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,50] = 	xjac[0,50] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,50] = 	xjac[0,50] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,50] = 	xjac[0,50] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,50] = 	xjac[0,50] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,50] = 	xjac[0,50] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,50] = 	xjac[0,50] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,50] = 	xjac[0,50] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,50] = 	xjac[0,50] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,50] = 	xjac[0,50] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,50] = 	xjac[0,50] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,50] = 	xjac[0,50] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,50] = 	xjac[0,50] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,51] = 	xjac[0,51] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,51] = 	xjac[0,51] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,51] = 	xjac[0,51] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,51] = 	xjac[0,51] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,51] = 	xjac[0,51] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,51] = 	xjac[0,51] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,51] = 	xjac[0,51] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,51] = 	xjac[0,51] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,51] = 	xjac[0,51] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,51] = 	xjac[0,51] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,51] = 	xjac[0,51] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,51] = 	xjac[0,51] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,51] = 	xjac[0,51] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,51] = 	xjac[0,51] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,52] = 	xjac[0,52] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,52] = 	xjac[0,52] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,52] = 	xjac[0,52] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,52] = 	xjac[0,52] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,52] = 	xjac[0,52] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,52] = 	xjac[0,52] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,52] = 	xjac[0,52] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,52] = 	xjac[0,52] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,52] = 	xjac[0,52] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,52] = 	xjac[0,52] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,52] = 	xjac[0,52] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,52] = 	xjac[0,52] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,52] = 	xjac[0,52] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,52] = 	xjac[0,52] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,53] = 	xjac[0,53] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,53] = 	xjac[0,53] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,53] = 	xjac[0,53] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,53] = 	xjac[0,53] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,53] = 	xjac[0,53] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,53] = 	xjac[0,53] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,53] = 	xjac[0,53] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,53] = 	xjac[0,53] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,53] = 	xjac[0,53] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,53] = 	xjac[0,53] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,53] = 	xjac[0,53] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,53] = 	xjac[0,53] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,53] = 	xjac[0,53] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,53] = 	xjac[0,53] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,54] = 	xjac[0,54] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,54] = 	xjac[0,54] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,54] = 	xjac[0,54] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,54] = 	xjac[0,54] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,54] = 	xjac[0,54] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,54] = 	xjac[0,54] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,54] = 	xjac[0,54] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,54] = 	xjac[0,54] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,54] = 	xjac[0,54] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,54] = 	xjac[0,54] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,54] = 	xjac[0,54] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,54] = 	xjac[0,54] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,54] = 	xjac[0,54] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,54] = 	xjac[0,54] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,55] = 	xjac[0,55] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,55] = 	xjac[0,55] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,55] = 	xjac[0,55] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,55] = 	xjac[0,55] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,55] = 	xjac[0,55] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,55] = 	xjac[0,55] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,55] = 	xjac[0,55] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,55] = 	xjac[0,55] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,55] = 	xjac[0,55] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,55] = 	xjac[0,55] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,55] = 	xjac[0,55] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,55] = 	xjac[0,55] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,55] = 	xjac[0,55] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,55] = 	xjac[0,55] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO
	xjac[0,56] = 	xjac[0,56] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[0,56] = 	xjac[0,56] - (1.8e-11*pow(300/T,0)*exp(110/T)*pow(x[0],1.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[0,56] = 	xjac[0,56] - (7.4e-11*pow(300/T,0)*exp(120/T)*pow(x[0],1.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[0,56] = 	xjac[0,56] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[0,56] = 	xjac[0,56] - (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[0,56] = 	xjac[0,56] - (1.1e-10*pow(x[0],1.0)*pow(x[13],1.0))#O_1D + H2 = H + OH
	xjac[0,56] = 	xjac[0,56] - (2.2e-10*pow(x[0],1.0)*pow(x[6],1.0))#O_1D + H2O = 2OH
	xjac[0,56] = 	xjac[0,56] - (1.55e-10*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl + ClO
	xjac[0,56] = 	xjac[0,56] - (5.25e-11*pow(x[0],1.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[0,56] = 	xjac[0,56] - (1e-10*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = Cl + OH
	xjac[0,56] = 	xjac[0,56] - (3.6e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = ClO + H
	xjac[0,56] = 	xjac[0,56] - (1.35e-11*pow(x[0],1.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[0,56] = 	xjac[0,56] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = Cl2 + CO2
	xjac[0,56] = 	xjac[0,56] - (3.6e-10*pow(x[0],1.0)*pow(x[27],1.0))#O_1D + COCl2 = ClCO + ClO

	xjac[1,0] = 	xjac[1,0] - (3.2e-11*pow(300/T,0)*exp(70/T)*1.0*pow(x[0],0.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][0] = 	xjac[1,0] + (3.2e-11*pow(300/T,0)*exp(70/T)*1.0*pow(x[0],0.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][0] = 	xjac[1,0] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,0] = 	xjac[1,0] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][0] = 	xjac[1,0] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][0] = 	xjac[1,0] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][0] = 	xjac[1,0] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][0] = 	xjac[1,0] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][0] = 	xjac[1,0] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,0] = 	xjac[1,0] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,0] = 	xjac[1,0] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][0] = 	xjac[1,0] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,0] = 	xjac[1,0] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,0] = 	xjac[1,0] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,0] = 	xjac[1,0] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,0] = 	xjac[1,0] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][0] = 	xjac[1,0] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][0] = 	xjac[1,0] + (1.2e-10*1.0*pow(x[0],0.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][0] = 	xjac[1,0] + (1.2e-10*1.0*pow(x[0],0.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][0] = 	xjac[1,0] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][0] = 	xjac[1,0] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][0] = 	xjac[1,0] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][0] = 	xjac[1,0] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][0] = 	xjac[1,0] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][0] = 	xjac[1,0] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][0] = 	xjac[1,0] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][0] = 	xjac[1,0] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][0] = 	xjac[1,0] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][0] = 	xjac[1,0] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][0] = 	xjac[1,0] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][0] = 	xjac[1,0] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][0] = 	xjac[1,0] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][0] = 	xjac[1,0] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][0] = 	xjac[1,0] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][0] = 	xjac[1,0] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][0] = 	xjac[1,0] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][0] = 	xjac[1,0] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,0] = 	xjac[1,0] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,0] = 	xjac[1,0] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,0] = 	xjac[1,0] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,0] = 	xjac[1,0] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][0] = 	xjac[1,0] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][0] = 	xjac[1,0] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][0] = 	xjac[1,0] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][0] = 	xjac[1,0] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][0] = 	xjac[1,0] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][0] = 	xjac[1,0] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][0] = 	xjac[1,0] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][0] = 	xjac[1,0] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,0] = 	xjac[1,0] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,0] = 	xjac[1,0] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][0] = 	xjac[1,0] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,0] = 	xjac[1,0] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][0] = 	xjac[1,0] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][0] = 	xjac[1,0] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][0] = 	xjac[1,0] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][0] = 	xjac[1,0] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][0] = 	xjac[1,0] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][0] = 	xjac[1,0] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][0] = 	xjac[1,0] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][0] = 	xjac[1,0] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][0] = 	xjac[1,0] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][0] = 	xjac[1,0] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,1] = 	xjac[1,1] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*1.0*pow(x[1],0.0))#O_1D + O2 = O + O2
	xjac[1][1] = 	xjac[1,1] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*1.0*pow(x[1],0.0))#O_1D + O2 = O + O2
	xjac[1][1] = 	xjac[1,1] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,1] = 	xjac[1,1] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*1.0*pow(x[1],0.0))#O2_1d + O2 = 2O2
	xjac[1][1] = 	xjac[1,1] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*1.0*pow(x[1],0.0))#O2_1d + O2 = 2O2
	xjac[1][1] = 	xjac[1,1] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][1] = 	xjac[1,1] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][1] = 	xjac[1,1] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][1] = 	xjac[1,1] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,1] = 	xjac[1,1] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*1.0*pow(x[1],0.0))#2O + O2 = O3 + O
	xjac[1,1] = 	xjac[1,1] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*2.0*pow(x[1],1.0))#O + 2O2 = O3 + O2
	xjac[1][1] = 	xjac[1,1] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*2.0*pow(x[1],1.0))#O + 2O2 = O3 + O2
	xjac[1,1] = 	xjac[1,1] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*1.0*pow(x[1],0.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,1] = 	xjac[1,1] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*1.0*pow(x[1],0.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,1] = 	xjac[1,1] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*1.0*pow(x[1],0.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,1] = 	xjac[1,1] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*1.0*pow(x[1],0.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][1] = 	xjac[1,1] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][1] = 	xjac[1,1] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][1] = 	xjac[1,1] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][1] = 	xjac[1,1] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][1] = 	xjac[1,1] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][1] = 	xjac[1,1] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][1] = 	xjac[1,1] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][1] = 	xjac[1,1] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][1] = 	xjac[1,1] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][1] = 	xjac[1,1] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][1] = 	xjac[1,1] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][1] = 	xjac[1,1] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][1] = 	xjac[1,1] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][1] = 	xjac[1,1] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][1] = 	xjac[1,1] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][1] = 	xjac[1,1] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][1] = 	xjac[1,1] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][1] = 	xjac[1,1] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][1] = 	xjac[1,1] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][1] = 	xjac[1,1] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][1] = 	xjac[1,1] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,1] = 	xjac[1,1] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*1.0*pow(x[1],0.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,1] = 	xjac[1,1] - (2.3e-12*pow(x[32],1.0)*1.0*pow(x[1],0.0))#S + O2 = SO + O
	xjac[1,1] = 	xjac[1,1] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*1.0*pow(x[1],0.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,1] = 	xjac[1,1] - (2e-15*pow(x[22],1.0)*1.0*pow(x[1],0.0))#ClS + O2 = SO + ClO
	xjac[1][1] = 	xjac[1,1] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][1] = 	xjac[1,1] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][1] = 	xjac[1,1] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][1] = 	xjac[1,1] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][1] = 	xjac[1,1] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][1] = 	xjac[1,1] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][1] = 	xjac[1,1] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][1] = 	xjac[1,1] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,1] = 	xjac[1,1] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*1.0*pow(x[1],0.0))#N + O2 = NO + O
	xjac[1,1] = 	xjac[1,1] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*1.0*pow(x[1],0.0))#HNO + O2 = NO + HO2
	xjac[1][1] = 	xjac[1,1] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,1] = 	xjac[1,1] - (9e-17*1.0*pow(x[1],0.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][1] = 	xjac[1,1] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][1] = 	xjac[1,1] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][1] = 	xjac[1,1] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][1] = 	xjac[1,1] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][1] = 	xjac[1,1] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][1] = 	xjac[1,1] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][1] = 	xjac[1,1] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][1] = 	xjac[1,1] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][1] = 	xjac[1,1] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][1] = 	xjac[1,1] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,2] = 	xjac[1,2] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][2] = 	xjac[1,2] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][2] = 	xjac[1,2] + (2e-16*pow(x[5],1.0)*1.0*pow(x[2],0.0))#O2_1d + O = O2 + O
	xjac[1,2] = 	xjac[1,2] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][2] = 	xjac[1,2] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][2] = 	xjac[1,2] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][2] = 	xjac[1,2] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][2] = 	xjac[1,2] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][2] = 	xjac[1,2] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*2.0*pow(x[2],1.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,2] = 	xjac[1,2] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*2.0*pow(x[2],1.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,2] = 	xjac[1,2] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*1.0*pow(x[2],0.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][2] = 	xjac[1,2] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*1.0*pow(x[2],0.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,2] = 	xjac[1,2] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*1.0*pow(x[2],0.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,2] = 	xjac[1,2] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*1.0*pow(x[2],0.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,2] = 	xjac[1,2] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*1.0*pow(x[2],0.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,2] = 	xjac[1,2] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][2] = 	xjac[1,2] + (8e-12*pow(300/T,0)*exp(-2060/T)*1.0*pow(x[2],0.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][2] = 	xjac[1,2] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][2] = 	xjac[1,2] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][2] = 	xjac[1,2] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][2] = 	xjac[1,2] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][2] = 	xjac[1,2] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][2] = 	xjac[1,2] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][2] = 	xjac[1,2] + (2.2e-11*pow(300/T,0)*exp(120/T)*1.0*pow(x[2],0.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][2] = 	xjac[1,2] + (3e-11*pow(300/T,0)*exp(200/T)*1.0*pow(x[2],0.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][2] = 	xjac[1,2] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][2] = 	xjac[1,2] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][2] = 	xjac[1,2] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][2] = 	xjac[1,2] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][2] = 	xjac[1,2] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][2] = 	xjac[1,2] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][2] = 	xjac[1,2] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*1.0*pow(x[2],0.0))#ClO + O = Cl + O2
	xjac[1][2] = 	xjac[1,2] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][2] = 	xjac[1,2] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][2] = 	xjac[1,2] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][2] = 	xjac[1,2] + (1e-11*1.0*pow(x[2],0.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][2] = 	xjac[1,2] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,2] = 	xjac[1,2] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,2] = 	xjac[1,2] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,2] = 	xjac[1,2] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,2] = 	xjac[1,2] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][2] = 	xjac[1,2] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][2] = 	xjac[1,2] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][2] = 	xjac[1,2] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][2] = 	xjac[1,2] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*1.0*pow(x[2],0.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][2] = 	xjac[1,2] + (3e-14*1.0*pow(x[2],0.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][2] = 	xjac[1,2] + (8e-12*pow(300/T,0)*exp(-9800/T)*1.0*pow(x[2],0.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][2] = 	xjac[1,2] + (1.3e-10*1.0*pow(x[2],0.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][2] = 	xjac[1,2] + (2.32e-16*pow(300/T,0)*exp(-487/T)*1.0*pow(x[2],0.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,2] = 	xjac[1,2] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,2] = 	xjac[1,2] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][2] = 	xjac[1,2] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,2] = 	xjac[1,2] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][2] = 	xjac[1,2] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][2] = 	xjac[1,2] + (5.6e-12*pow(300/T,0)*exp(180/T)*1.0*pow(x[2],0.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][2] = 	xjac[1,2] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][2] = 	xjac[1,2] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][2] = 	xjac[1,2] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][2] = 	xjac[1,2] + (1e-11*1.0*pow(x[2],0.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][2] = 	xjac[1,2] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][2] = 	xjac[1,2] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][2] = 	xjac[1,2] + (4.9e-11*1.0*pow(x[2],0.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][2] = 	xjac[1,2] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,3] = 	xjac[1,3] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][3] = 	xjac[1,3] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][3] = 	xjac[1,3] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,3] = 	xjac[1,3] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][3] = 	xjac[1,3] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][3] = 	xjac[1,3] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][3] = 	xjac[1,3] + (1e-20*pow(x[5],1.0)*1.0*pow(x[3],0.0))#O2_1d + N2 = O2 + N2
	xjac[1][3] = 	xjac[1,3] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][3] = 	xjac[1,3] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,3] = 	xjac[1,3] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,3] = 	xjac[1,3] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][3] = 	xjac[1,3] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,3] = 	xjac[1,3] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*1.0*pow(x[3],0.0))#O + O2 + N2 = O3 + N2
	xjac[1,3] = 	xjac[1,3] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,3] = 	xjac[1,3] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,3] = 	xjac[1,3] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*1.0*pow(x[3],0.0))#H + O2 + N2 = HO2 + N2
	xjac[1][3] = 	xjac[1,3] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][3] = 	xjac[1,3] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][3] = 	xjac[1,3] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][3] = 	xjac[1,3] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][3] = 	xjac[1,3] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][3] = 	xjac[1,3] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][3] = 	xjac[1,3] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][3] = 	xjac[1,3] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][3] = 	xjac[1,3] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][3] = 	xjac[1,3] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][3] = 	xjac[1,3] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][3] = 	xjac[1,3] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][3] = 	xjac[1,3] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][3] = 	xjac[1,3] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][3] = 	xjac[1,3] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][3] = 	xjac[1,3] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][3] = 	xjac[1,3] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][3] = 	xjac[1,3] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][3] = 	xjac[1,3] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][3] = 	xjac[1,3] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][3] = 	xjac[1,3] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,3] = 	xjac[1,3] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,3] = 	xjac[1,3] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,3] = 	xjac[1,3] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,3] = 	xjac[1,3] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][3] = 	xjac[1,3] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][3] = 	xjac[1,3] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][3] = 	xjac[1,3] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][3] = 	xjac[1,3] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][3] = 	xjac[1,3] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][3] = 	xjac[1,3] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][3] = 	xjac[1,3] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][3] = 	xjac[1,3] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,3] = 	xjac[1,3] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,3] = 	xjac[1,3] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][3] = 	xjac[1,3] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,3] = 	xjac[1,3] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][3] = 	xjac[1,3] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][3] = 	xjac[1,3] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][3] = 	xjac[1,3] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][3] = 	xjac[1,3] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][3] = 	xjac[1,3] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][3] = 	xjac[1,3] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][3] = 	xjac[1,3] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][3] = 	xjac[1,3] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][3] = 	xjac[1,3] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][3] = 	xjac[1,3] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,4] = 	xjac[1,4] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][4] = 	xjac[1,4] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][4] = 	xjac[1,4] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,4] = 	xjac[1,4] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][4] = 	xjac[1,4] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][4] = 	xjac[1,4] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][4] = 	xjac[1,4] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][4] = 	xjac[1,4] + (2e-21*pow(x[5],1.0)*1.0*pow(x[4],0.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][4] = 	xjac[1,4] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*1.0*pow(x[4],0.0))#2O + CO2 = O2 + CO2
	xjac[1,4] = 	xjac[1,4] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,4] = 	xjac[1,4] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][4] = 	xjac[1,4] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,4] = 	xjac[1,4] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,4] = 	xjac[1,4] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,4] = 	xjac[1,4] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*1.0*pow(x[4],0.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,4] = 	xjac[1,4] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][4] = 	xjac[1,4] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][4] = 	xjac[1,4] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][4] = 	xjac[1,4] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][4] = 	xjac[1,4] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][4] = 	xjac[1,4] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][4] = 	xjac[1,4] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][4] = 	xjac[1,4] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][4] = 	xjac[1,4] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][4] = 	xjac[1,4] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][4] = 	xjac[1,4] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][4] = 	xjac[1,4] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][4] = 	xjac[1,4] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][4] = 	xjac[1,4] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][4] = 	xjac[1,4] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][4] = 	xjac[1,4] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][4] = 	xjac[1,4] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][4] = 	xjac[1,4] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][4] = 	xjac[1,4] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][4] = 	xjac[1,4] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][4] = 	xjac[1,4] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][4] = 	xjac[1,4] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,4] = 	xjac[1,4] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*1.0*pow(x[4],0.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,4] = 	xjac[1,4] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,4] = 	xjac[1,4] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,4] = 	xjac[1,4] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][4] = 	xjac[1,4] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][4] = 	xjac[1,4] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][4] = 	xjac[1,4] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][4] = 	xjac[1,4] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][4] = 	xjac[1,4] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][4] = 	xjac[1,4] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][4] = 	xjac[1,4] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][4] = 	xjac[1,4] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,4] = 	xjac[1,4] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,4] = 	xjac[1,4] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][4] = 	xjac[1,4] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,4] = 	xjac[1,4] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][4] = 	xjac[1,4] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][4] = 	xjac[1,4] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][4] = 	xjac[1,4] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][4] = 	xjac[1,4] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][4] = 	xjac[1,4] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][4] = 	xjac[1,4] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][4] = 	xjac[1,4] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][4] = 	xjac[1,4] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][4] = 	xjac[1,4] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][4] = 	xjac[1,4] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,5] = 	xjac[1,5] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][5] = 	xjac[1,5] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][5] = 	xjac[1,5] + (2e-16*1.0*pow(x[5],0.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,5] = 	xjac[1,5] - (3.6e-18*pow(300/T,0)*exp(-220/T)*1.0*pow(x[5],0.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][5] = 	xjac[1,5] + (3.6e-18*pow(300/T,0)*exp(-220/T)*1.0*pow(x[5],0.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][5] = 	xjac[1,5] + (4.8e-18*1.0*pow(x[5],0.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][5] = 	xjac[1,5] + (1e-20*1.0*pow(x[5],0.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][5] = 	xjac[1,5] + (2e-21*1.0*pow(x[5],0.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][5] = 	xjac[1,5] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,5] = 	xjac[1,5] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,5] = 	xjac[1,5] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][5] = 	xjac[1,5] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,5] = 	xjac[1,5] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,5] = 	xjac[1,5] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,5] = 	xjac[1,5] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,5] = 	xjac[1,5] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][5] = 	xjac[1,5] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][5] = 	xjac[1,5] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][5] = 	xjac[1,5] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][5] = 	xjac[1,5] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*1.0*pow(x[5],0.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][5] = 	xjac[1,5] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][5] = 	xjac[1,5] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][5] = 	xjac[1,5] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][5] = 	xjac[1,5] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][5] = 	xjac[1,5] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][5] = 	xjac[1,5] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][5] = 	xjac[1,5] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][5] = 	xjac[1,5] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][5] = 	xjac[1,5] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][5] = 	xjac[1,5] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][5] = 	xjac[1,5] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][5] = 	xjac[1,5] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][5] = 	xjac[1,5] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][5] = 	xjac[1,5] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][5] = 	xjac[1,5] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][5] = 	xjac[1,5] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][5] = 	xjac[1,5] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,5] = 	xjac[1,5] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,5] = 	xjac[1,5] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,5] = 	xjac[1,5] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,5] = 	xjac[1,5] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][5] = 	xjac[1,5] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][5] = 	xjac[1,5] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][5] = 	xjac[1,5] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][5] = 	xjac[1,5] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][5] = 	xjac[1,5] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][5] = 	xjac[1,5] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][5] = 	xjac[1,5] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][5] = 	xjac[1,5] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,5] = 	xjac[1,5] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,5] = 	xjac[1,5] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][5] = 	xjac[1,5] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,5] = 	xjac[1,5] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][5] = 	xjac[1,5] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][5] = 	xjac[1,5] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][5] = 	xjac[1,5] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][5] = 	xjac[1,5] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][5] = 	xjac[1,5] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][5] = 	xjac[1,5] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][5] = 	xjac[1,5] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][5] = 	xjac[1,5] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][5] = 	xjac[1,5] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][5] = 	xjac[1,5] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,6] = 	xjac[1,6] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][6] = 	xjac[1,6] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][6] = 	xjac[1,6] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,6] = 	xjac[1,6] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][6] = 	xjac[1,6] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][6] = 	xjac[1,6] + (4.8e-18*pow(x[5],1.0)*1.0*pow(x[6],0.0))#O2_1d + H2O = O2 + H2O
	xjac[1][6] = 	xjac[1,6] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][6] = 	xjac[1,6] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][6] = 	xjac[1,6] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,6] = 	xjac[1,6] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,6] = 	xjac[1,6] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][6] = 	xjac[1,6] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,6] = 	xjac[1,6] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,6] = 	xjac[1,6] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,6] = 	xjac[1,6] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,6] = 	xjac[1,6] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][6] = 	xjac[1,6] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][6] = 	xjac[1,6] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][6] = 	xjac[1,6] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][6] = 	xjac[1,6] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][6] = 	xjac[1,6] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][6] = 	xjac[1,6] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][6] = 	xjac[1,6] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][6] = 	xjac[1,6] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][6] = 	xjac[1,6] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][6] = 	xjac[1,6] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][6] = 	xjac[1,6] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][6] = 	xjac[1,6] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][6] = 	xjac[1,6] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][6] = 	xjac[1,6] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][6] = 	xjac[1,6] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][6] = 	xjac[1,6] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][6] = 	xjac[1,6] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][6] = 	xjac[1,6] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][6] = 	xjac[1,6] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][6] = 	xjac[1,6] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][6] = 	xjac[1,6] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,6] = 	xjac[1,6] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,6] = 	xjac[1,6] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,6] = 	xjac[1,6] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,6] = 	xjac[1,6] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][6] = 	xjac[1,6] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][6] = 	xjac[1,6] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][6] = 	xjac[1,6] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][6] = 	xjac[1,6] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][6] = 	xjac[1,6] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][6] = 	xjac[1,6] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][6] = 	xjac[1,6] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][6] = 	xjac[1,6] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,6] = 	xjac[1,6] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,6] = 	xjac[1,6] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][6] = 	xjac[1,6] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,6] = 	xjac[1,6] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][6] = 	xjac[1,6] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][6] = 	xjac[1,6] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][6] = 	xjac[1,6] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][6] = 	xjac[1,6] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][6] = 	xjac[1,6] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][6] = 	xjac[1,6] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][6] = 	xjac[1,6] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][6] = 	xjac[1,6] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][6] = 	xjac[1,6] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][6] = 	xjac[1,6] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,7] = 	xjac[1,7] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][7] = 	xjac[1,7] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][7] = 	xjac[1,7] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,7] = 	xjac[1,7] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][7] = 	xjac[1,7] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][7] = 	xjac[1,7] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][7] = 	xjac[1,7] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][7] = 	xjac[1,7] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][7] = 	xjac[1,7] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,7] = 	xjac[1,7] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,7] = 	xjac[1,7] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][7] = 	xjac[1,7] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,7] = 	xjac[1,7] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,7] = 	xjac[1,7] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*1.0*pow(x[7],0.0))#O + O2 + CO = O3 + CO
	xjac[1,7] = 	xjac[1,7] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,7] = 	xjac[1,7] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][7] = 	xjac[1,7] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][7] = 	xjac[1,7] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][7] = 	xjac[1,7] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][7] = 	xjac[1,7] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][7] = 	xjac[1,7] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][7] = 	xjac[1,7] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][7] = 	xjac[1,7] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][7] = 	xjac[1,7] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][7] = 	xjac[1,7] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][7] = 	xjac[1,7] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][7] = 	xjac[1,7] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][7] = 	xjac[1,7] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][7] = 	xjac[1,7] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][7] = 	xjac[1,7] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][7] = 	xjac[1,7] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][7] = 	xjac[1,7] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][7] = 	xjac[1,7] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][7] = 	xjac[1,7] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][7] = 	xjac[1,7] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][7] = 	xjac[1,7] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][7] = 	xjac[1,7] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,7] = 	xjac[1,7] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,7] = 	xjac[1,7] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,7] = 	xjac[1,7] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,7] = 	xjac[1,7] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][7] = 	xjac[1,7] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][7] = 	xjac[1,7] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][7] = 	xjac[1,7] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][7] = 	xjac[1,7] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][7] = 	xjac[1,7] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][7] = 	xjac[1,7] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][7] = 	xjac[1,7] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][7] = 	xjac[1,7] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,7] = 	xjac[1,7] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,7] = 	xjac[1,7] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][7] = 	xjac[1,7] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,7] = 	xjac[1,7] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][7] = 	xjac[1,7] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][7] = 	xjac[1,7] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][7] = 	xjac[1,7] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][7] = 	xjac[1,7] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][7] = 	xjac[1,7] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][7] = 	xjac[1,7] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][7] = 	xjac[1,7] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][7] = 	xjac[1,7] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][7] = 	xjac[1,7] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][7] = 	xjac[1,7] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,8] = 	xjac[1,8] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][8] = 	xjac[1,8] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][8] = 	xjac[1,8] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,8] = 	xjac[1,8] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][8] = 	xjac[1,8] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][8] = 	xjac[1,8] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][8] = 	xjac[1,8] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][8] = 	xjac[1,8] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][8] = 	xjac[1,8] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,8] = 	xjac[1,8] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,8] = 	xjac[1,8] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][8] = 	xjac[1,8] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,8] = 	xjac[1,8] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,8] = 	xjac[1,8] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,8] = 	xjac[1,8] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,8] = 	xjac[1,8] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][8] = 	xjac[1,8] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*1.0*pow(x[8],0.0))#O + O3 = 2O2
	xjac[1][8] = 	xjac[1,8] + (1.2e-10*pow(x[0],1.0)*1.0*pow(x[8],0.0))#O_1D + O3 = 2O2
	xjac[1][8] = 	xjac[1,8] + (1.2e-10*pow(x[0],1.0)*1.0*pow(x[8],0.0))#O_1D + O3 = 2O + O2
	xjac[1][8] = 	xjac[1,8] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*1.0*pow(x[8],0.0))#O2_1d + O3 = 2O2 + O
	xjac[1][8] = 	xjac[1,8] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*1.0*pow(x[8],0.0))#H + O3 = OH + O2
	xjac[1][8] = 	xjac[1,8] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*1.0*pow(x[8],0.0))#OH + O3 = HO2 + O2
	xjac[1][8] = 	xjac[1,8] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*1.0*pow(x[8],0.0))#HO2 + O3 = OH + 2O2
	xjac[1][8] = 	xjac[1,8] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][8] = 	xjac[1,8] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][8] = 	xjac[1,8] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][8] = 	xjac[1,8] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][8] = 	xjac[1,8] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][8] = 	xjac[1,8] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][8] = 	xjac[1,8] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*1.0*pow(x[8],0.0))#Cl + O3 = ClO + O2
	xjac[1][8] = 	xjac[1,8] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][8] = 	xjac[1,8] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][8] = 	xjac[1,8] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][8] = 	xjac[1,8] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][8] = 	xjac[1,8] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][8] = 	xjac[1,8] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][8] = 	xjac[1,8] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,8] = 	xjac[1,8] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,8] = 	xjac[1,8] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,8] = 	xjac[1,8] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,8] = 	xjac[1,8] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][8] = 	xjac[1,8] + (1.2e-11*pow(x[32],1.0)*1.0*pow(x[8],0.0))#S + O3 = SO + O2
	xjac[1][8] = 	xjac[1,8] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*1.0*pow(x[8],0.0))#SO + O3 = SO2 + O2
	xjac[1][8] = 	xjac[1,8] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*1.0*pow(x[8],0.0))#SO2 + O3 = SO3 + O2
	xjac[1][8] = 	xjac[1,8] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][8] = 	xjac[1,8] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][8] = 	xjac[1,8] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][8] = 	xjac[1,8] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][8] = 	xjac[1,8] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,8] = 	xjac[1,8] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,8] = 	xjac[1,8] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][8] = 	xjac[1,8] + (2e-16*pow(x[48],1.0)*1.0*pow(x[8],0.0))#N + O3 = NO + O2
	xjac[1,8] = 	xjac[1,8] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][8] = 	xjac[1,8] + (3e-12*pow(300/T,0)*exp(-1500/T)*1.0*pow(x[8],0.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][8] = 	xjac[1,8] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][8] = 	xjac[1,8] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*1.0*pow(x[8],0.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][8] = 	xjac[1,8] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][8] = 	xjac[1,8] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][8] = 	xjac[1,8] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][8] = 	xjac[1,8] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][8] = 	xjac[1,8] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][8] = 	xjac[1,8] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][8] = 	xjac[1,8] + (5e-19*1.0*pow(x[8],0.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,9] = 	xjac[1,9] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][9] = 	xjac[1,9] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][9] = 	xjac[1,9] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,9] = 	xjac[1,9] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][9] = 	xjac[1,9] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][9] = 	xjac[1,9] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][9] = 	xjac[1,9] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][9] = 	xjac[1,9] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][9] = 	xjac[1,9] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,9] = 	xjac[1,9] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,9] = 	xjac[1,9] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][9] = 	xjac[1,9] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,9] = 	xjac[1,9] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,9] = 	xjac[1,9] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,9] = 	xjac[1,9] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,9] = 	xjac[1,9] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*1.0*pow(x[9],0.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][9] = 	xjac[1,9] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][9] = 	xjac[1,9] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][9] = 	xjac[1,9] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][9] = 	xjac[1,9] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][9] = 	xjac[1,9] + (1.4e-10*pow(300/T,0)*exp(-470/T)*1.0*pow(x[9],0.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][9] = 	xjac[1,9] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][9] = 	xjac[1,9] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][9] = 	xjac[1,9] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][9] = 	xjac[1,9] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][9] = 	xjac[1,9] + (7.29e-12*1.0*pow(x[9],0.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][9] = 	xjac[1,9] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][9] = 	xjac[1,9] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][9] = 	xjac[1,9] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][9] = 	xjac[1,9] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][9] = 	xjac[1,9] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][9] = 	xjac[1,9] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][9] = 	xjac[1,9] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][9] = 	xjac[1,9] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][9] = 	xjac[1,9] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][9] = 	xjac[1,9] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][9] = 	xjac[1,9] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,9] = 	xjac[1,9] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*1.0*pow(x[9],0.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,9] = 	xjac[1,9] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,9] = 	xjac[1,9] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,9] = 	xjac[1,9] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][9] = 	xjac[1,9] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][9] = 	xjac[1,9] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][9] = 	xjac[1,9] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][9] = 	xjac[1,9] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][9] = 	xjac[1,9] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][9] = 	xjac[1,9] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][9] = 	xjac[1,9] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][9] = 	xjac[1,9] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,9] = 	xjac[1,9] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,9] = 	xjac[1,9] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][9] = 	xjac[1,9] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,9] = 	xjac[1,9] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][9] = 	xjac[1,9] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][9] = 	xjac[1,9] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][9] = 	xjac[1,9] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][9] = 	xjac[1,9] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][9] = 	xjac[1,9] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][9] = 	xjac[1,9] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][9] = 	xjac[1,9] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][9] = 	xjac[1,9] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][9] = 	xjac[1,9] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][9] = 	xjac[1,9] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,10] = 	xjac[1,10] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][10] = 	xjac[1,10] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][10] = 	xjac[1,10] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,10] = 	xjac[1,10] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][10] = 	xjac[1,10] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][10] = 	xjac[1,10] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][10] = 	xjac[1,10] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][10] = 	xjac[1,10] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][10] = 	xjac[1,10] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,10] = 	xjac[1,10] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,10] = 	xjac[1,10] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][10] = 	xjac[1,10] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,10] = 	xjac[1,10] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,10] = 	xjac[1,10] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,10] = 	xjac[1,10] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,10] = 	xjac[1,10] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][10] = 	xjac[1,10] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][10] = 	xjac[1,10] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][10] = 	xjac[1,10] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][10] = 	xjac[1,10] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][10] = 	xjac[1,10] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][10] = 	xjac[1,10] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][10] = 	xjac[1,10] + (1e-14*pow(300/T,0)*exp(-490/T)*1.0*pow(x[10],0.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][10] = 	xjac[1,10] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][10] = 	xjac[1,10] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*1.0*pow(x[10],0.0))#O + HO2 = OH + O2
	xjac[1][10] = 	xjac[1,10] + (7.29e-12*pow(x[9],1.0)*1.0*pow(x[10],0.0))#H + HO2 = H2 + O2
	xjac[1][10] = 	xjac[1,10] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*1.0*pow(x[10],0.0))#OH + HO2 = H2O + O2
	xjac[1][10] = 	xjac[1,10] + (2.3e-13*pow(300/T,0)*exp(600/T)*2.0*pow(x[10],1.0))#2HO2 = H2O2 + O2
	xjac[1][10] = 	xjac[1,10] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*2.0*pow(x[10],1.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][10] = 	xjac[1,10] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][10] = 	xjac[1,10] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*1.0*pow(x[10],0.0))#Cl + HO2 = HCl + O2
	xjac[1][10] = 	xjac[1,10] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][10] = 	xjac[1,10] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][10] = 	xjac[1,10] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*1.0*pow(x[10],0.0))#ClO + HO2 = HOCl + O2
	xjac[1][10] = 	xjac[1,10] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][10] = 	xjac[1,10] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][10] = 	xjac[1,10] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,10] = 	xjac[1,10] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,10] = 	xjac[1,10] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,10] = 	xjac[1,10] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,10] = 	xjac[1,10] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][10] = 	xjac[1,10] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][10] = 	xjac[1,10] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][10] = 	xjac[1,10] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][10] = 	xjac[1,10] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][10] = 	xjac[1,10] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][10] = 	xjac[1,10] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][10] = 	xjac[1,10] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][10] = 	xjac[1,10] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,10] = 	xjac[1,10] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,10] = 	xjac[1,10] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][10] = 	xjac[1,10] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,10] = 	xjac[1,10] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][10] = 	xjac[1,10] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][10] = 	xjac[1,10] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][10] = 	xjac[1,10] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][10] = 	xjac[1,10] + (5e-16*1.0*pow(x[10],0.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][10] = 	xjac[1,10] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][10] = 	xjac[1,10] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][10] = 	xjac[1,10] + (3.5e-12*1.0*pow(x[10],0.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][10] = 	xjac[1,10] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][10] = 	xjac[1,10] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][10] = 	xjac[1,10] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,11] = 	xjac[1,11] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][11] = 	xjac[1,11] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][11] = 	xjac[1,11] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,11] = 	xjac[1,11] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][11] = 	xjac[1,11] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][11] = 	xjac[1,11] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][11] = 	xjac[1,11] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][11] = 	xjac[1,11] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][11] = 	xjac[1,11] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,11] = 	xjac[1,11] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,11] = 	xjac[1,11] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][11] = 	xjac[1,11] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,11] = 	xjac[1,11] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,11] = 	xjac[1,11] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,11] = 	xjac[1,11] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,11] = 	xjac[1,11] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][11] = 	xjac[1,11] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][11] = 	xjac[1,11] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][11] = 	xjac[1,11] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][11] = 	xjac[1,11] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][11] = 	xjac[1,11] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][11] = 	xjac[1,11] + (1.7e-12*pow(300/T,0)*exp(-940/T)*1.0*pow(x[11],0.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][11] = 	xjac[1,11] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][11] = 	xjac[1,11] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*1.0*pow(x[11],0.0))#O + OH = O2 + H
	xjac[1][11] = 	xjac[1,11] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][11] = 	xjac[1,11] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][11] = 	xjac[1,11] + (4.8e-11*pow(300/T,0)*exp(250/T)*1.0*pow(x[11],0.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][11] = 	xjac[1,11] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][11] = 	xjac[1,11] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][11] = 	xjac[1,11] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][11] = 	xjac[1,11] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][11] = 	xjac[1,11] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][11] = 	xjac[1,11] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*1.0*pow(x[11],0.0))#ClO + OH = HCl + O2
	xjac[1][11] = 	xjac[1,11] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][11] = 	xjac[1,11] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][11] = 	xjac[1,11] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][11] = 	xjac[1,11] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,11] = 	xjac[1,11] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,11] = 	xjac[1,11] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,11] = 	xjac[1,11] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,11] = 	xjac[1,11] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][11] = 	xjac[1,11] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][11] = 	xjac[1,11] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][11] = 	xjac[1,11] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][11] = 	xjac[1,11] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][11] = 	xjac[1,11] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][11] = 	xjac[1,11] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][11] = 	xjac[1,11] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][11] = 	xjac[1,11] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,11] = 	xjac[1,11] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,11] = 	xjac[1,11] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][11] = 	xjac[1,11] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,11] = 	xjac[1,11] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][11] = 	xjac[1,11] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][11] = 	xjac[1,11] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][11] = 	xjac[1,11] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][11] = 	xjac[1,11] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][11] = 	xjac[1,11] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][11] = 	xjac[1,11] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][11] = 	xjac[1,11] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][11] = 	xjac[1,11] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][11] = 	xjac[1,11] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][11] = 	xjac[1,11] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,12] = 	xjac[1,12] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][12] = 	xjac[1,12] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][12] = 	xjac[1,12] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,12] = 	xjac[1,12] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][12] = 	xjac[1,12] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][12] = 	xjac[1,12] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][12] = 	xjac[1,12] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][12] = 	xjac[1,12] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][12] = 	xjac[1,12] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,12] = 	xjac[1,12] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,12] = 	xjac[1,12] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][12] = 	xjac[1,12] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,12] = 	xjac[1,12] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,12] = 	xjac[1,12] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,12] = 	xjac[1,12] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,12] = 	xjac[1,12] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][12] = 	xjac[1,12] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][12] = 	xjac[1,12] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][12] = 	xjac[1,12] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][12] = 	xjac[1,12] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][12] = 	xjac[1,12] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][12] = 	xjac[1,12] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][12] = 	xjac[1,12] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][12] = 	xjac[1,12] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][12] = 	xjac[1,12] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][12] = 	xjac[1,12] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][12] = 	xjac[1,12] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][12] = 	xjac[1,12] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][12] = 	xjac[1,12] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*1.0*pow(x[12],0.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][12] = 	xjac[1,12] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][12] = 	xjac[1,12] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][12] = 	xjac[1,12] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][12] = 	xjac[1,12] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][12] = 	xjac[1,12] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][12] = 	xjac[1,12] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][12] = 	xjac[1,12] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][12] = 	xjac[1,12] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,12] = 	xjac[1,12] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,12] = 	xjac[1,12] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,12] = 	xjac[1,12] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,12] = 	xjac[1,12] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][12] = 	xjac[1,12] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][12] = 	xjac[1,12] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][12] = 	xjac[1,12] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][12] = 	xjac[1,12] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][12] = 	xjac[1,12] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][12] = 	xjac[1,12] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][12] = 	xjac[1,12] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][12] = 	xjac[1,12] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,12] = 	xjac[1,12] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,12] = 	xjac[1,12] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][12] = 	xjac[1,12] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,12] = 	xjac[1,12] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][12] = 	xjac[1,12] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][12] = 	xjac[1,12] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][12] = 	xjac[1,12] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][12] = 	xjac[1,12] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][12] = 	xjac[1,12] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][12] = 	xjac[1,12] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][12] = 	xjac[1,12] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][12] = 	xjac[1,12] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][12] = 	xjac[1,12] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][12] = 	xjac[1,12] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,13] = 	xjac[1,13] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][13] = 	xjac[1,13] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][13] = 	xjac[1,13] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,13] = 	xjac[1,13] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][13] = 	xjac[1,13] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][13] = 	xjac[1,13] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][13] = 	xjac[1,13] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][13] = 	xjac[1,13] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][13] = 	xjac[1,13] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,13] = 	xjac[1,13] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,13] = 	xjac[1,13] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][13] = 	xjac[1,13] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,13] = 	xjac[1,13] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,13] = 	xjac[1,13] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,13] = 	xjac[1,13] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,13] = 	xjac[1,13] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][13] = 	xjac[1,13] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][13] = 	xjac[1,13] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][13] = 	xjac[1,13] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][13] = 	xjac[1,13] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][13] = 	xjac[1,13] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][13] = 	xjac[1,13] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][13] = 	xjac[1,13] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][13] = 	xjac[1,13] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][13] = 	xjac[1,13] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][13] = 	xjac[1,13] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][13] = 	xjac[1,13] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][13] = 	xjac[1,13] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][13] = 	xjac[1,13] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][13] = 	xjac[1,13] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][13] = 	xjac[1,13] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][13] = 	xjac[1,13] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][13] = 	xjac[1,13] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][13] = 	xjac[1,13] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][13] = 	xjac[1,13] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][13] = 	xjac[1,13] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][13] = 	xjac[1,13] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,13] = 	xjac[1,13] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,13] = 	xjac[1,13] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,13] = 	xjac[1,13] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,13] = 	xjac[1,13] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][13] = 	xjac[1,13] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][13] = 	xjac[1,13] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][13] = 	xjac[1,13] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][13] = 	xjac[1,13] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][13] = 	xjac[1,13] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][13] = 	xjac[1,13] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][13] = 	xjac[1,13] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][13] = 	xjac[1,13] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,13] = 	xjac[1,13] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,13] = 	xjac[1,13] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][13] = 	xjac[1,13] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,13] = 	xjac[1,13] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][13] = 	xjac[1,13] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][13] = 	xjac[1,13] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][13] = 	xjac[1,13] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][13] = 	xjac[1,13] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][13] = 	xjac[1,13] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][13] = 	xjac[1,13] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][13] = 	xjac[1,13] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][13] = 	xjac[1,13] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][13] = 	xjac[1,13] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][13] = 	xjac[1,13] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,14] = 	xjac[1,14] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][14] = 	xjac[1,14] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][14] = 	xjac[1,14] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,14] = 	xjac[1,14] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][14] = 	xjac[1,14] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][14] = 	xjac[1,14] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][14] = 	xjac[1,14] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][14] = 	xjac[1,14] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][14] = 	xjac[1,14] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,14] = 	xjac[1,14] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,14] = 	xjac[1,14] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][14] = 	xjac[1,14] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,14] = 	xjac[1,14] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,14] = 	xjac[1,14] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,14] = 	xjac[1,14] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,14] = 	xjac[1,14] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][14] = 	xjac[1,14] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][14] = 	xjac[1,14] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][14] = 	xjac[1,14] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][14] = 	xjac[1,14] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][14] = 	xjac[1,14] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][14] = 	xjac[1,14] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][14] = 	xjac[1,14] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][14] = 	xjac[1,14] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][14] = 	xjac[1,14] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][14] = 	xjac[1,14] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][14] = 	xjac[1,14] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][14] = 	xjac[1,14] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][14] = 	xjac[1,14] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][14] = 	xjac[1,14] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][14] = 	xjac[1,14] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][14] = 	xjac[1,14] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][14] = 	xjac[1,14] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][14] = 	xjac[1,14] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][14] = 	xjac[1,14] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][14] = 	xjac[1,14] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][14] = 	xjac[1,14] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,14] = 	xjac[1,14] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,14] = 	xjac[1,14] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,14] = 	xjac[1,14] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,14] = 	xjac[1,14] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][14] = 	xjac[1,14] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][14] = 	xjac[1,14] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][14] = 	xjac[1,14] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][14] = 	xjac[1,14] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][14] = 	xjac[1,14] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][14] = 	xjac[1,14] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][14] = 	xjac[1,14] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][14] = 	xjac[1,14] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,14] = 	xjac[1,14] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,14] = 	xjac[1,14] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][14] = 	xjac[1,14] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,14] = 	xjac[1,14] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][14] = 	xjac[1,14] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][14] = 	xjac[1,14] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][14] = 	xjac[1,14] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][14] = 	xjac[1,14] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][14] = 	xjac[1,14] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][14] = 	xjac[1,14] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][14] = 	xjac[1,14] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][14] = 	xjac[1,14] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][14] = 	xjac[1,14] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][14] = 	xjac[1,14] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,15] = 	xjac[1,15] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][15] = 	xjac[1,15] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][15] = 	xjac[1,15] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,15] = 	xjac[1,15] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][15] = 	xjac[1,15] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][15] = 	xjac[1,15] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][15] = 	xjac[1,15] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][15] = 	xjac[1,15] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][15] = 	xjac[1,15] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,15] = 	xjac[1,15] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,15] = 	xjac[1,15] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][15] = 	xjac[1,15] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,15] = 	xjac[1,15] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,15] = 	xjac[1,15] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,15] = 	xjac[1,15] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,15] = 	xjac[1,15] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][15] = 	xjac[1,15] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][15] = 	xjac[1,15] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][15] = 	xjac[1,15] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][15] = 	xjac[1,15] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][15] = 	xjac[1,15] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][15] = 	xjac[1,15] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][15] = 	xjac[1,15] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][15] = 	xjac[1,15] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][15] = 	xjac[1,15] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][15] = 	xjac[1,15] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][15] = 	xjac[1,15] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][15] = 	xjac[1,15] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][15] = 	xjac[1,15] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][15] = 	xjac[1,15] + (2.3e-11*pow(300/T,0)*exp(-200/T)*1.0*pow(x[15],0.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][15] = 	xjac[1,15] + (1.8e-11*pow(300/T,0)*exp(170/T)*1.0*pow(x[15],0.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][15] = 	xjac[1,15] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][15] = 	xjac[1,15] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][15] = 	xjac[1,15] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][15] = 	xjac[1,15] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][15] = 	xjac[1,15] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][15] = 	xjac[1,15] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,15] = 	xjac[1,15] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,15] = 	xjac[1,15] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,15] = 	xjac[1,15] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,15] = 	xjac[1,15] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][15] = 	xjac[1,15] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][15] = 	xjac[1,15] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][15] = 	xjac[1,15] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][15] = 	xjac[1,15] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][15] = 	xjac[1,15] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][15] = 	xjac[1,15] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][15] = 	xjac[1,15] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][15] = 	xjac[1,15] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,15] = 	xjac[1,15] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,15] = 	xjac[1,15] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][15] = 	xjac[1,15] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,15] = 	xjac[1,15] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][15] = 	xjac[1,15] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][15] = 	xjac[1,15] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][15] = 	xjac[1,15] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][15] = 	xjac[1,15] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][15] = 	xjac[1,15] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][15] = 	xjac[1,15] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][15] = 	xjac[1,15] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][15] = 	xjac[1,15] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][15] = 	xjac[1,15] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][15] = 	xjac[1,15] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,16] = 	xjac[1,16] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][16] = 	xjac[1,16] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][16] = 	xjac[1,16] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,16] = 	xjac[1,16] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][16] = 	xjac[1,16] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][16] = 	xjac[1,16] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][16] = 	xjac[1,16] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][16] = 	xjac[1,16] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][16] = 	xjac[1,16] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,16] = 	xjac[1,16] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,16] = 	xjac[1,16] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][16] = 	xjac[1,16] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,16] = 	xjac[1,16] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,16] = 	xjac[1,16] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,16] = 	xjac[1,16] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,16] = 	xjac[1,16] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][16] = 	xjac[1,16] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][16] = 	xjac[1,16] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][16] = 	xjac[1,16] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][16] = 	xjac[1,16] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][16] = 	xjac[1,16] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][16] = 	xjac[1,16] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][16] = 	xjac[1,16] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][16] = 	xjac[1,16] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][16] = 	xjac[1,16] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][16] = 	xjac[1,16] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][16] = 	xjac[1,16] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][16] = 	xjac[1,16] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][16] = 	xjac[1,16] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][16] = 	xjac[1,16] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][16] = 	xjac[1,16] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][16] = 	xjac[1,16] + (3e-11*pow(300/T,0)*exp(70/T)*1.0*pow(x[16],0.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][16] = 	xjac[1,16] + (6e-13*pow(300/T,0)*exp(230/T)*1.0*pow(x[16],0.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][16] = 	xjac[1,16] + (2.7e-12*pow(300/T,0)*exp(220/T)*1.0*pow(x[16],0.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][16] = 	xjac[1,16] + (1e-12*pow(300/T,0)*exp(-1590/T)*2.0*pow(x[16],1.0))#2ClO = Cl2 + O2
	xjac[1][16] = 	xjac[1,16] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][16] = 	xjac[1,16] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,16] = 	xjac[1,16] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,16] = 	xjac[1,16] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,16] = 	xjac[1,16] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,16] = 	xjac[1,16] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][16] = 	xjac[1,16] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][16] = 	xjac[1,16] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][16] = 	xjac[1,16] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][16] = 	xjac[1,16] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][16] = 	xjac[1,16] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][16] = 	xjac[1,16] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][16] = 	xjac[1,16] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][16] = 	xjac[1,16] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,16] = 	xjac[1,16] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,16] = 	xjac[1,16] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][16] = 	xjac[1,16] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,16] = 	xjac[1,16] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][16] = 	xjac[1,16] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][16] = 	xjac[1,16] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][16] = 	xjac[1,16] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][16] = 	xjac[1,16] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][16] = 	xjac[1,16] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][16] = 	xjac[1,16] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][16] = 	xjac[1,16] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][16] = 	xjac[1,16] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][16] = 	xjac[1,16] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][16] = 	xjac[1,16] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,17] = 	xjac[1,17] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][17] = 	xjac[1,17] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][17] = 	xjac[1,17] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,17] = 	xjac[1,17] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][17] = 	xjac[1,17] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][17] = 	xjac[1,17] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][17] = 	xjac[1,17] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][17] = 	xjac[1,17] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][17] = 	xjac[1,17] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,17] = 	xjac[1,17] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,17] = 	xjac[1,17] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][17] = 	xjac[1,17] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,17] = 	xjac[1,17] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,17] = 	xjac[1,17] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,17] = 	xjac[1,17] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,17] = 	xjac[1,17] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][17] = 	xjac[1,17] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][17] = 	xjac[1,17] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][17] = 	xjac[1,17] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][17] = 	xjac[1,17] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][17] = 	xjac[1,17] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][17] = 	xjac[1,17] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][17] = 	xjac[1,17] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][17] = 	xjac[1,17] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][17] = 	xjac[1,17] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][17] = 	xjac[1,17] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][17] = 	xjac[1,17] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][17] = 	xjac[1,17] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][17] = 	xjac[1,17] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][17] = 	xjac[1,17] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][17] = 	xjac[1,17] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][17] = 	xjac[1,17] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][17] = 	xjac[1,17] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][17] = 	xjac[1,17] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][17] = 	xjac[1,17] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][17] = 	xjac[1,17] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][17] = 	xjac[1,17] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,17] = 	xjac[1,17] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,17] = 	xjac[1,17] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,17] = 	xjac[1,17] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,17] = 	xjac[1,17] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][17] = 	xjac[1,17] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][17] = 	xjac[1,17] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][17] = 	xjac[1,17] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][17] = 	xjac[1,17] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][17] = 	xjac[1,17] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][17] = 	xjac[1,17] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][17] = 	xjac[1,17] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][17] = 	xjac[1,17] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,17] = 	xjac[1,17] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,17] = 	xjac[1,17] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][17] = 	xjac[1,17] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,17] = 	xjac[1,17] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][17] = 	xjac[1,17] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][17] = 	xjac[1,17] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][17] = 	xjac[1,17] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][17] = 	xjac[1,17] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][17] = 	xjac[1,17] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][17] = 	xjac[1,17] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][17] = 	xjac[1,17] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][17] = 	xjac[1,17] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][17] = 	xjac[1,17] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][17] = 	xjac[1,17] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,18] = 	xjac[1,18] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][18] = 	xjac[1,18] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][18] = 	xjac[1,18] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,18] = 	xjac[1,18] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][18] = 	xjac[1,18] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][18] = 	xjac[1,18] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][18] = 	xjac[1,18] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][18] = 	xjac[1,18] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][18] = 	xjac[1,18] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,18] = 	xjac[1,18] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,18] = 	xjac[1,18] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][18] = 	xjac[1,18] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,18] = 	xjac[1,18] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,18] = 	xjac[1,18] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,18] = 	xjac[1,18] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,18] = 	xjac[1,18] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][18] = 	xjac[1,18] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][18] = 	xjac[1,18] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][18] = 	xjac[1,18] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][18] = 	xjac[1,18] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][18] = 	xjac[1,18] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][18] = 	xjac[1,18] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][18] = 	xjac[1,18] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][18] = 	xjac[1,18] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][18] = 	xjac[1,18] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][18] = 	xjac[1,18] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][18] = 	xjac[1,18] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][18] = 	xjac[1,18] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][18] = 	xjac[1,18] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][18] = 	xjac[1,18] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][18] = 	xjac[1,18] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][18] = 	xjac[1,18] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][18] = 	xjac[1,18] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][18] = 	xjac[1,18] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][18] = 	xjac[1,18] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][18] = 	xjac[1,18] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][18] = 	xjac[1,18] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,18] = 	xjac[1,18] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,18] = 	xjac[1,18] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,18] = 	xjac[1,18] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,18] = 	xjac[1,18] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][18] = 	xjac[1,18] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][18] = 	xjac[1,18] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][18] = 	xjac[1,18] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][18] = 	xjac[1,18] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][18] = 	xjac[1,18] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][18] = 	xjac[1,18] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][18] = 	xjac[1,18] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][18] = 	xjac[1,18] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,18] = 	xjac[1,18] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,18] = 	xjac[1,18] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][18] = 	xjac[1,18] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,18] = 	xjac[1,18] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][18] = 	xjac[1,18] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][18] = 	xjac[1,18] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][18] = 	xjac[1,18] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][18] = 	xjac[1,18] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][18] = 	xjac[1,18] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][18] = 	xjac[1,18] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][18] = 	xjac[1,18] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][18] = 	xjac[1,18] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][18] = 	xjac[1,18] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][18] = 	xjac[1,18] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,19] = 	xjac[1,19] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][19] = 	xjac[1,19] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][19] = 	xjac[1,19] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,19] = 	xjac[1,19] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][19] = 	xjac[1,19] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][19] = 	xjac[1,19] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][19] = 	xjac[1,19] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][19] = 	xjac[1,19] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][19] = 	xjac[1,19] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,19] = 	xjac[1,19] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,19] = 	xjac[1,19] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][19] = 	xjac[1,19] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,19] = 	xjac[1,19] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,19] = 	xjac[1,19] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,19] = 	xjac[1,19] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,19] = 	xjac[1,19] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][19] = 	xjac[1,19] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][19] = 	xjac[1,19] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][19] = 	xjac[1,19] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][19] = 	xjac[1,19] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][19] = 	xjac[1,19] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][19] = 	xjac[1,19] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][19] = 	xjac[1,19] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][19] = 	xjac[1,19] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][19] = 	xjac[1,19] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][19] = 	xjac[1,19] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][19] = 	xjac[1,19] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][19] = 	xjac[1,19] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][19] = 	xjac[1,19] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][19] = 	xjac[1,19] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][19] = 	xjac[1,19] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][19] = 	xjac[1,19] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][19] = 	xjac[1,19] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][19] = 	xjac[1,19] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][19] = 	xjac[1,19] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][19] = 	xjac[1,19] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][19] = 	xjac[1,19] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,19] = 	xjac[1,19] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,19] = 	xjac[1,19] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,19] = 	xjac[1,19] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,19] = 	xjac[1,19] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][19] = 	xjac[1,19] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][19] = 	xjac[1,19] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][19] = 	xjac[1,19] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][19] = 	xjac[1,19] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][19] = 	xjac[1,19] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][19] = 	xjac[1,19] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][19] = 	xjac[1,19] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][19] = 	xjac[1,19] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,19] = 	xjac[1,19] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,19] = 	xjac[1,19] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][19] = 	xjac[1,19] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,19] = 	xjac[1,19] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][19] = 	xjac[1,19] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][19] = 	xjac[1,19] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][19] = 	xjac[1,19] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][19] = 	xjac[1,19] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][19] = 	xjac[1,19] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][19] = 	xjac[1,19] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][19] = 	xjac[1,19] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][19] = 	xjac[1,19] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][19] = 	xjac[1,19] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][19] = 	xjac[1,19] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,20] = 	xjac[1,20] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][20] = 	xjac[1,20] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][20] = 	xjac[1,20] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,20] = 	xjac[1,20] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][20] = 	xjac[1,20] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][20] = 	xjac[1,20] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][20] = 	xjac[1,20] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][20] = 	xjac[1,20] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][20] = 	xjac[1,20] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,20] = 	xjac[1,20] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,20] = 	xjac[1,20] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][20] = 	xjac[1,20] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,20] = 	xjac[1,20] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,20] = 	xjac[1,20] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,20] = 	xjac[1,20] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,20] = 	xjac[1,20] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][20] = 	xjac[1,20] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][20] = 	xjac[1,20] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][20] = 	xjac[1,20] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][20] = 	xjac[1,20] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][20] = 	xjac[1,20] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][20] = 	xjac[1,20] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][20] = 	xjac[1,20] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][20] = 	xjac[1,20] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][20] = 	xjac[1,20] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][20] = 	xjac[1,20] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][20] = 	xjac[1,20] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][20] = 	xjac[1,20] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][20] = 	xjac[1,20] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][20] = 	xjac[1,20] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][20] = 	xjac[1,20] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][20] = 	xjac[1,20] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][20] = 	xjac[1,20] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][20] = 	xjac[1,20] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][20] = 	xjac[1,20] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][20] = 	xjac[1,20] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][20] = 	xjac[1,20] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,20] = 	xjac[1,20] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,20] = 	xjac[1,20] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,20] = 	xjac[1,20] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,20] = 	xjac[1,20] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][20] = 	xjac[1,20] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][20] = 	xjac[1,20] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][20] = 	xjac[1,20] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][20] = 	xjac[1,20] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][20] = 	xjac[1,20] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][20] = 	xjac[1,20] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][20] = 	xjac[1,20] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][20] = 	xjac[1,20] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,20] = 	xjac[1,20] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,20] = 	xjac[1,20] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][20] = 	xjac[1,20] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,20] = 	xjac[1,20] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][20] = 	xjac[1,20] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][20] = 	xjac[1,20] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][20] = 	xjac[1,20] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][20] = 	xjac[1,20] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][20] = 	xjac[1,20] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][20] = 	xjac[1,20] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][20] = 	xjac[1,20] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][20] = 	xjac[1,20] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][20] = 	xjac[1,20] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][20] = 	xjac[1,20] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,21] = 	xjac[1,21] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][21] = 	xjac[1,21] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][21] = 	xjac[1,21] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,21] = 	xjac[1,21] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][21] = 	xjac[1,21] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][21] = 	xjac[1,21] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][21] = 	xjac[1,21] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][21] = 	xjac[1,21] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][21] = 	xjac[1,21] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,21] = 	xjac[1,21] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,21] = 	xjac[1,21] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][21] = 	xjac[1,21] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,21] = 	xjac[1,21] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,21] = 	xjac[1,21] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,21] = 	xjac[1,21] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,21] = 	xjac[1,21] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][21] = 	xjac[1,21] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][21] = 	xjac[1,21] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][21] = 	xjac[1,21] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][21] = 	xjac[1,21] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][21] = 	xjac[1,21] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][21] = 	xjac[1,21] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][21] = 	xjac[1,21] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][21] = 	xjac[1,21] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][21] = 	xjac[1,21] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][21] = 	xjac[1,21] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][21] = 	xjac[1,21] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][21] = 	xjac[1,21] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][21] = 	xjac[1,21] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][21] = 	xjac[1,21] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][21] = 	xjac[1,21] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][21] = 	xjac[1,21] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][21] = 	xjac[1,21] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][21] = 	xjac[1,21] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][21] = 	xjac[1,21] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][21] = 	xjac[1,21] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][21] = 	xjac[1,21] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,21] = 	xjac[1,21] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,21] = 	xjac[1,21] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,21] = 	xjac[1,21] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,21] = 	xjac[1,21] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][21] = 	xjac[1,21] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][21] = 	xjac[1,21] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][21] = 	xjac[1,21] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][21] = 	xjac[1,21] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][21] = 	xjac[1,21] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][21] = 	xjac[1,21] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][21] = 	xjac[1,21] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][21] = 	xjac[1,21] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,21] = 	xjac[1,21] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,21] = 	xjac[1,21] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][21] = 	xjac[1,21] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,21] = 	xjac[1,21] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][21] = 	xjac[1,21] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][21] = 	xjac[1,21] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][21] = 	xjac[1,21] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][21] = 	xjac[1,21] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][21] = 	xjac[1,21] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][21] = 	xjac[1,21] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][21] = 	xjac[1,21] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][21] = 	xjac[1,21] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][21] = 	xjac[1,21] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][21] = 	xjac[1,21] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,22] = 	xjac[1,22] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][22] = 	xjac[1,22] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][22] = 	xjac[1,22] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,22] = 	xjac[1,22] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][22] = 	xjac[1,22] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][22] = 	xjac[1,22] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][22] = 	xjac[1,22] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][22] = 	xjac[1,22] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][22] = 	xjac[1,22] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,22] = 	xjac[1,22] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,22] = 	xjac[1,22] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][22] = 	xjac[1,22] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,22] = 	xjac[1,22] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,22] = 	xjac[1,22] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,22] = 	xjac[1,22] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,22] = 	xjac[1,22] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][22] = 	xjac[1,22] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][22] = 	xjac[1,22] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][22] = 	xjac[1,22] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][22] = 	xjac[1,22] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][22] = 	xjac[1,22] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][22] = 	xjac[1,22] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][22] = 	xjac[1,22] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][22] = 	xjac[1,22] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][22] = 	xjac[1,22] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][22] = 	xjac[1,22] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][22] = 	xjac[1,22] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][22] = 	xjac[1,22] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][22] = 	xjac[1,22] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][22] = 	xjac[1,22] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][22] = 	xjac[1,22] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][22] = 	xjac[1,22] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][22] = 	xjac[1,22] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][22] = 	xjac[1,22] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][22] = 	xjac[1,22] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][22] = 	xjac[1,22] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][22] = 	xjac[1,22] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,22] = 	xjac[1,22] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,22] = 	xjac[1,22] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,22] = 	xjac[1,22] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,22] = 	xjac[1,22] - (2e-15*1.0*pow(x[22],0.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][22] = 	xjac[1,22] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][22] = 	xjac[1,22] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][22] = 	xjac[1,22] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][22] = 	xjac[1,22] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][22] = 	xjac[1,22] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][22] = 	xjac[1,22] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][22] = 	xjac[1,22] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][22] = 	xjac[1,22] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,22] = 	xjac[1,22] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,22] = 	xjac[1,22] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][22] = 	xjac[1,22] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,22] = 	xjac[1,22] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][22] = 	xjac[1,22] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][22] = 	xjac[1,22] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][22] = 	xjac[1,22] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][22] = 	xjac[1,22] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][22] = 	xjac[1,22] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][22] = 	xjac[1,22] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][22] = 	xjac[1,22] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][22] = 	xjac[1,22] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][22] = 	xjac[1,22] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][22] = 	xjac[1,22] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,23] = 	xjac[1,23] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][23] = 	xjac[1,23] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][23] = 	xjac[1,23] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,23] = 	xjac[1,23] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][23] = 	xjac[1,23] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][23] = 	xjac[1,23] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][23] = 	xjac[1,23] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][23] = 	xjac[1,23] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][23] = 	xjac[1,23] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,23] = 	xjac[1,23] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,23] = 	xjac[1,23] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][23] = 	xjac[1,23] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,23] = 	xjac[1,23] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,23] = 	xjac[1,23] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,23] = 	xjac[1,23] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,23] = 	xjac[1,23] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][23] = 	xjac[1,23] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][23] = 	xjac[1,23] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][23] = 	xjac[1,23] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][23] = 	xjac[1,23] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][23] = 	xjac[1,23] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][23] = 	xjac[1,23] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][23] = 	xjac[1,23] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][23] = 	xjac[1,23] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][23] = 	xjac[1,23] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][23] = 	xjac[1,23] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][23] = 	xjac[1,23] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][23] = 	xjac[1,23] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][23] = 	xjac[1,23] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][23] = 	xjac[1,23] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][23] = 	xjac[1,23] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][23] = 	xjac[1,23] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][23] = 	xjac[1,23] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][23] = 	xjac[1,23] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][23] = 	xjac[1,23] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][23] = 	xjac[1,23] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][23] = 	xjac[1,23] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,23] = 	xjac[1,23] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,23] = 	xjac[1,23] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,23] = 	xjac[1,23] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,23] = 	xjac[1,23] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][23] = 	xjac[1,23] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][23] = 	xjac[1,23] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][23] = 	xjac[1,23] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][23] = 	xjac[1,23] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][23] = 	xjac[1,23] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][23] = 	xjac[1,23] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][23] = 	xjac[1,23] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][23] = 	xjac[1,23] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,23] = 	xjac[1,23] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,23] = 	xjac[1,23] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][23] = 	xjac[1,23] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,23] = 	xjac[1,23] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][23] = 	xjac[1,23] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][23] = 	xjac[1,23] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][23] = 	xjac[1,23] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][23] = 	xjac[1,23] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][23] = 	xjac[1,23] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][23] = 	xjac[1,23] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][23] = 	xjac[1,23] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][23] = 	xjac[1,23] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][23] = 	xjac[1,23] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][23] = 	xjac[1,23] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,24] = 	xjac[1,24] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][24] = 	xjac[1,24] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][24] = 	xjac[1,24] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,24] = 	xjac[1,24] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][24] = 	xjac[1,24] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][24] = 	xjac[1,24] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][24] = 	xjac[1,24] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][24] = 	xjac[1,24] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][24] = 	xjac[1,24] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,24] = 	xjac[1,24] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,24] = 	xjac[1,24] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][24] = 	xjac[1,24] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,24] = 	xjac[1,24] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,24] = 	xjac[1,24] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,24] = 	xjac[1,24] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,24] = 	xjac[1,24] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][24] = 	xjac[1,24] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][24] = 	xjac[1,24] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][24] = 	xjac[1,24] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][24] = 	xjac[1,24] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][24] = 	xjac[1,24] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][24] = 	xjac[1,24] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][24] = 	xjac[1,24] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][24] = 	xjac[1,24] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][24] = 	xjac[1,24] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][24] = 	xjac[1,24] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][24] = 	xjac[1,24] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][24] = 	xjac[1,24] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][24] = 	xjac[1,24] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][24] = 	xjac[1,24] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][24] = 	xjac[1,24] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][24] = 	xjac[1,24] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][24] = 	xjac[1,24] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][24] = 	xjac[1,24] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][24] = 	xjac[1,24] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][24] = 	xjac[1,24] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][24] = 	xjac[1,24] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,24] = 	xjac[1,24] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,24] = 	xjac[1,24] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,24] = 	xjac[1,24] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,24] = 	xjac[1,24] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][24] = 	xjac[1,24] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][24] = 	xjac[1,24] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][24] = 	xjac[1,24] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][24] = 	xjac[1,24] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][24] = 	xjac[1,24] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][24] = 	xjac[1,24] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][24] = 	xjac[1,24] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][24] = 	xjac[1,24] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,24] = 	xjac[1,24] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,24] = 	xjac[1,24] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][24] = 	xjac[1,24] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,24] = 	xjac[1,24] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][24] = 	xjac[1,24] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][24] = 	xjac[1,24] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][24] = 	xjac[1,24] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][24] = 	xjac[1,24] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][24] = 	xjac[1,24] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][24] = 	xjac[1,24] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][24] = 	xjac[1,24] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][24] = 	xjac[1,24] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][24] = 	xjac[1,24] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][24] = 	xjac[1,24] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,25] = 	xjac[1,25] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][25] = 	xjac[1,25] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][25] = 	xjac[1,25] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,25] = 	xjac[1,25] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][25] = 	xjac[1,25] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][25] = 	xjac[1,25] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][25] = 	xjac[1,25] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][25] = 	xjac[1,25] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][25] = 	xjac[1,25] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,25] = 	xjac[1,25] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,25] = 	xjac[1,25] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][25] = 	xjac[1,25] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,25] = 	xjac[1,25] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,25] = 	xjac[1,25] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,25] = 	xjac[1,25] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,25] = 	xjac[1,25] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][25] = 	xjac[1,25] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][25] = 	xjac[1,25] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][25] = 	xjac[1,25] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][25] = 	xjac[1,25] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][25] = 	xjac[1,25] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][25] = 	xjac[1,25] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][25] = 	xjac[1,25] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][25] = 	xjac[1,25] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][25] = 	xjac[1,25] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][25] = 	xjac[1,25] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][25] = 	xjac[1,25] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][25] = 	xjac[1,25] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][25] = 	xjac[1,25] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][25] = 	xjac[1,25] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][25] = 	xjac[1,25] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][25] = 	xjac[1,25] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][25] = 	xjac[1,25] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][25] = 	xjac[1,25] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][25] = 	xjac[1,25] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][25] = 	xjac[1,25] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][25] = 	xjac[1,25] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,25] = 	xjac[1,25] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,25] = 	xjac[1,25] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,25] = 	xjac[1,25] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,25] = 	xjac[1,25] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][25] = 	xjac[1,25] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][25] = 	xjac[1,25] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][25] = 	xjac[1,25] + (3e-12*pow(300/T,0)*exp(-7000/T)*1.0*pow(x[25],0.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][25] = 	xjac[1,25] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][25] = 	xjac[1,25] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][25] = 	xjac[1,25] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*1.0*pow(x[25],0.0))#O + SO2 = SO + O2
	xjac[1][25] = 	xjac[1,25] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*1.0*pow(x[25],0.0))#O(1D) + SO2 = SO + O2
	xjac[1][25] = 	xjac[1,25] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,25] = 	xjac[1,25] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,25] = 	xjac[1,25] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][25] = 	xjac[1,25] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,25] = 	xjac[1,25] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][25] = 	xjac[1,25] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][25] = 	xjac[1,25] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][25] = 	xjac[1,25] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][25] = 	xjac[1,25] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][25] = 	xjac[1,25] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][25] = 	xjac[1,25] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][25] = 	xjac[1,25] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][25] = 	xjac[1,25] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][25] = 	xjac[1,25] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][25] = 	xjac[1,25] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,26] = 	xjac[1,26] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][26] = 	xjac[1,26] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][26] = 	xjac[1,26] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,26] = 	xjac[1,26] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][26] = 	xjac[1,26] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][26] = 	xjac[1,26] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][26] = 	xjac[1,26] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][26] = 	xjac[1,26] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][26] = 	xjac[1,26] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,26] = 	xjac[1,26] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,26] = 	xjac[1,26] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][26] = 	xjac[1,26] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,26] = 	xjac[1,26] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,26] = 	xjac[1,26] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,26] = 	xjac[1,26] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,26] = 	xjac[1,26] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][26] = 	xjac[1,26] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][26] = 	xjac[1,26] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][26] = 	xjac[1,26] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][26] = 	xjac[1,26] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][26] = 	xjac[1,26] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][26] = 	xjac[1,26] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][26] = 	xjac[1,26] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][26] = 	xjac[1,26] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][26] = 	xjac[1,26] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][26] = 	xjac[1,26] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][26] = 	xjac[1,26] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][26] = 	xjac[1,26] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][26] = 	xjac[1,26] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][26] = 	xjac[1,26] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][26] = 	xjac[1,26] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][26] = 	xjac[1,26] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][26] = 	xjac[1,26] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][26] = 	xjac[1,26] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][26] = 	xjac[1,26] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][26] = 	xjac[1,26] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][26] = 	xjac[1,26] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,26] = 	xjac[1,26] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,26] = 	xjac[1,26] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,26] = 	xjac[1,26] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,26] = 	xjac[1,26] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][26] = 	xjac[1,26] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][26] = 	xjac[1,26] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][26] = 	xjac[1,26] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][26] = 	xjac[1,26] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][26] = 	xjac[1,26] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][26] = 	xjac[1,26] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][26] = 	xjac[1,26] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][26] = 	xjac[1,26] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,26] = 	xjac[1,26] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,26] = 	xjac[1,26] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][26] = 	xjac[1,26] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,26] = 	xjac[1,26] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][26] = 	xjac[1,26] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][26] = 	xjac[1,26] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][26] = 	xjac[1,26] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][26] = 	xjac[1,26] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][26] = 	xjac[1,26] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][26] = 	xjac[1,26] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][26] = 	xjac[1,26] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][26] = 	xjac[1,26] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][26] = 	xjac[1,26] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][26] = 	xjac[1,26] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,27] = 	xjac[1,27] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][27] = 	xjac[1,27] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][27] = 	xjac[1,27] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,27] = 	xjac[1,27] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][27] = 	xjac[1,27] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][27] = 	xjac[1,27] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][27] = 	xjac[1,27] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][27] = 	xjac[1,27] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][27] = 	xjac[1,27] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,27] = 	xjac[1,27] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,27] = 	xjac[1,27] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][27] = 	xjac[1,27] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,27] = 	xjac[1,27] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,27] = 	xjac[1,27] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,27] = 	xjac[1,27] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,27] = 	xjac[1,27] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][27] = 	xjac[1,27] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][27] = 	xjac[1,27] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][27] = 	xjac[1,27] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][27] = 	xjac[1,27] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][27] = 	xjac[1,27] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][27] = 	xjac[1,27] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][27] = 	xjac[1,27] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][27] = 	xjac[1,27] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][27] = 	xjac[1,27] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][27] = 	xjac[1,27] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][27] = 	xjac[1,27] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][27] = 	xjac[1,27] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][27] = 	xjac[1,27] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][27] = 	xjac[1,27] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][27] = 	xjac[1,27] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][27] = 	xjac[1,27] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][27] = 	xjac[1,27] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][27] = 	xjac[1,27] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][27] = 	xjac[1,27] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][27] = 	xjac[1,27] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][27] = 	xjac[1,27] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,27] = 	xjac[1,27] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,27] = 	xjac[1,27] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,27] = 	xjac[1,27] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,27] = 	xjac[1,27] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][27] = 	xjac[1,27] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][27] = 	xjac[1,27] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][27] = 	xjac[1,27] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][27] = 	xjac[1,27] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][27] = 	xjac[1,27] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][27] = 	xjac[1,27] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][27] = 	xjac[1,27] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][27] = 	xjac[1,27] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,27] = 	xjac[1,27] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,27] = 	xjac[1,27] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][27] = 	xjac[1,27] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,27] = 	xjac[1,27] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][27] = 	xjac[1,27] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][27] = 	xjac[1,27] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][27] = 	xjac[1,27] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][27] = 	xjac[1,27] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][27] = 	xjac[1,27] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][27] = 	xjac[1,27] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][27] = 	xjac[1,27] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][27] = 	xjac[1,27] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][27] = 	xjac[1,27] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][27] = 	xjac[1,27] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,28] = 	xjac[1,28] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][28] = 	xjac[1,28] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][28] = 	xjac[1,28] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,28] = 	xjac[1,28] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][28] = 	xjac[1,28] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][28] = 	xjac[1,28] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][28] = 	xjac[1,28] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][28] = 	xjac[1,28] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][28] = 	xjac[1,28] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,28] = 	xjac[1,28] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,28] = 	xjac[1,28] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][28] = 	xjac[1,28] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,28] = 	xjac[1,28] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,28] = 	xjac[1,28] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,28] = 	xjac[1,28] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,28] = 	xjac[1,28] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][28] = 	xjac[1,28] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][28] = 	xjac[1,28] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][28] = 	xjac[1,28] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][28] = 	xjac[1,28] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][28] = 	xjac[1,28] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][28] = 	xjac[1,28] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][28] = 	xjac[1,28] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][28] = 	xjac[1,28] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][28] = 	xjac[1,28] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][28] = 	xjac[1,28] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][28] = 	xjac[1,28] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][28] = 	xjac[1,28] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][28] = 	xjac[1,28] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][28] = 	xjac[1,28] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][28] = 	xjac[1,28] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][28] = 	xjac[1,28] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][28] = 	xjac[1,28] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][28] = 	xjac[1,28] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][28] = 	xjac[1,28] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][28] = 	xjac[1,28] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][28] = 	xjac[1,28] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,28] = 	xjac[1,28] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,28] = 	xjac[1,28] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,28] = 	xjac[1,28] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,28] = 	xjac[1,28] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][28] = 	xjac[1,28] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][28] = 	xjac[1,28] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][28] = 	xjac[1,28] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][28] = 	xjac[1,28] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][28] = 	xjac[1,28] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][28] = 	xjac[1,28] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][28] = 	xjac[1,28] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][28] = 	xjac[1,28] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,28] = 	xjac[1,28] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,28] = 	xjac[1,28] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][28] = 	xjac[1,28] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,28] = 	xjac[1,28] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][28] = 	xjac[1,28] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][28] = 	xjac[1,28] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][28] = 	xjac[1,28] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][28] = 	xjac[1,28] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][28] = 	xjac[1,28] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][28] = 	xjac[1,28] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][28] = 	xjac[1,28] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][28] = 	xjac[1,28] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][28] = 	xjac[1,28] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][28] = 	xjac[1,28] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,29] = 	xjac[1,29] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][29] = 	xjac[1,29] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][29] = 	xjac[1,29] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,29] = 	xjac[1,29] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][29] = 	xjac[1,29] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][29] = 	xjac[1,29] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][29] = 	xjac[1,29] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][29] = 	xjac[1,29] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][29] = 	xjac[1,29] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,29] = 	xjac[1,29] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,29] = 	xjac[1,29] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][29] = 	xjac[1,29] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,29] = 	xjac[1,29] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,29] = 	xjac[1,29] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,29] = 	xjac[1,29] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,29] = 	xjac[1,29] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][29] = 	xjac[1,29] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][29] = 	xjac[1,29] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][29] = 	xjac[1,29] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][29] = 	xjac[1,29] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][29] = 	xjac[1,29] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][29] = 	xjac[1,29] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][29] = 	xjac[1,29] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][29] = 	xjac[1,29] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][29] = 	xjac[1,29] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][29] = 	xjac[1,29] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][29] = 	xjac[1,29] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][29] = 	xjac[1,29] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][29] = 	xjac[1,29] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][29] = 	xjac[1,29] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][29] = 	xjac[1,29] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][29] = 	xjac[1,29] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][29] = 	xjac[1,29] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][29] = 	xjac[1,29] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][29] = 	xjac[1,29] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][29] = 	xjac[1,29] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][29] = 	xjac[1,29] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,29] = 	xjac[1,29] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,29] = 	xjac[1,29] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,29] = 	xjac[1,29] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,29] = 	xjac[1,29] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][29] = 	xjac[1,29] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][29] = 	xjac[1,29] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*1.0*pow(x[29],0.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][29] = 	xjac[1,29] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][29] = 	xjac[1,29] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*1.0*pow(x[29],0.0))#O + SO = S + O2
	xjac[1][29] = 	xjac[1,29] + (3e-14*pow(x[2],1.0)*1.0*pow(x[29],0.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][29] = 	xjac[1,29] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][29] = 	xjac[1,29] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][29] = 	xjac[1,29] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,29] = 	xjac[1,29] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,29] = 	xjac[1,29] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][29] = 	xjac[1,29] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,29] = 	xjac[1,29] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][29] = 	xjac[1,29] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][29] = 	xjac[1,29] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][29] = 	xjac[1,29] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][29] = 	xjac[1,29] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][29] = 	xjac[1,29] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][29] = 	xjac[1,29] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][29] = 	xjac[1,29] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][29] = 	xjac[1,29] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][29] = 	xjac[1,29] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][29] = 	xjac[1,29] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,30] = 	xjac[1,30] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][30] = 	xjac[1,30] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][30] = 	xjac[1,30] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,30] = 	xjac[1,30] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][30] = 	xjac[1,30] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][30] = 	xjac[1,30] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][30] = 	xjac[1,30] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][30] = 	xjac[1,30] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][30] = 	xjac[1,30] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,30] = 	xjac[1,30] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,30] = 	xjac[1,30] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][30] = 	xjac[1,30] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,30] = 	xjac[1,30] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,30] = 	xjac[1,30] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,30] = 	xjac[1,30] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,30] = 	xjac[1,30] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][30] = 	xjac[1,30] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][30] = 	xjac[1,30] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][30] = 	xjac[1,30] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][30] = 	xjac[1,30] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][30] = 	xjac[1,30] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][30] = 	xjac[1,30] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][30] = 	xjac[1,30] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][30] = 	xjac[1,30] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][30] = 	xjac[1,30] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][30] = 	xjac[1,30] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][30] = 	xjac[1,30] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][30] = 	xjac[1,30] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][30] = 	xjac[1,30] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][30] = 	xjac[1,30] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][30] = 	xjac[1,30] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][30] = 	xjac[1,30] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][30] = 	xjac[1,30] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][30] = 	xjac[1,30] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][30] = 	xjac[1,30] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][30] = 	xjac[1,30] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][30] = 	xjac[1,30] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,30] = 	xjac[1,30] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,30] = 	xjac[1,30] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,30] = 	xjac[1,30] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,30] = 	xjac[1,30] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][30] = 	xjac[1,30] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][30] = 	xjac[1,30] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][30] = 	xjac[1,30] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][30] = 	xjac[1,30] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][30] = 	xjac[1,30] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][30] = 	xjac[1,30] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][30] = 	xjac[1,30] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][30] = 	xjac[1,30] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*1.0*pow(x[30],0.0))#O + SO3 = SO2 + O2
	xjac[1,30] = 	xjac[1,30] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,30] = 	xjac[1,30] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][30] = 	xjac[1,30] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,30] = 	xjac[1,30] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][30] = 	xjac[1,30] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][30] = 	xjac[1,30] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][30] = 	xjac[1,30] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][30] = 	xjac[1,30] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][30] = 	xjac[1,30] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][30] = 	xjac[1,30] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][30] = 	xjac[1,30] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][30] = 	xjac[1,30] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][30] = 	xjac[1,30] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][30] = 	xjac[1,30] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,31] = 	xjac[1,31] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][31] = 	xjac[1,31] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][31] = 	xjac[1,31] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,31] = 	xjac[1,31] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][31] = 	xjac[1,31] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][31] = 	xjac[1,31] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][31] = 	xjac[1,31] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][31] = 	xjac[1,31] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][31] = 	xjac[1,31] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,31] = 	xjac[1,31] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,31] = 	xjac[1,31] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][31] = 	xjac[1,31] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,31] = 	xjac[1,31] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,31] = 	xjac[1,31] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,31] = 	xjac[1,31] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,31] = 	xjac[1,31] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][31] = 	xjac[1,31] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][31] = 	xjac[1,31] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][31] = 	xjac[1,31] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][31] = 	xjac[1,31] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][31] = 	xjac[1,31] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][31] = 	xjac[1,31] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][31] = 	xjac[1,31] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][31] = 	xjac[1,31] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][31] = 	xjac[1,31] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][31] = 	xjac[1,31] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][31] = 	xjac[1,31] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][31] = 	xjac[1,31] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][31] = 	xjac[1,31] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][31] = 	xjac[1,31] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][31] = 	xjac[1,31] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][31] = 	xjac[1,31] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][31] = 	xjac[1,31] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][31] = 	xjac[1,31] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][31] = 	xjac[1,31] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][31] = 	xjac[1,31] + (1e-11*pow(x[2],1.0)*1.0*pow(x[31],0.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][31] = 	xjac[1,31] + (5e-12*2.0*pow(x[31],1.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,31] = 	xjac[1,31] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,31] = 	xjac[1,31] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,31] = 	xjac[1,31] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,31] = 	xjac[1,31] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][31] = 	xjac[1,31] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][31] = 	xjac[1,31] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][31] = 	xjac[1,31] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][31] = 	xjac[1,31] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][31] = 	xjac[1,31] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][31] = 	xjac[1,31] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][31] = 	xjac[1,31] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][31] = 	xjac[1,31] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,31] = 	xjac[1,31] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,31] = 	xjac[1,31] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][31] = 	xjac[1,31] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,31] = 	xjac[1,31] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][31] = 	xjac[1,31] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][31] = 	xjac[1,31] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][31] = 	xjac[1,31] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][31] = 	xjac[1,31] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][31] = 	xjac[1,31] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][31] = 	xjac[1,31] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][31] = 	xjac[1,31] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][31] = 	xjac[1,31] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][31] = 	xjac[1,31] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][31] = 	xjac[1,31] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,32] = 	xjac[1,32] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][32] = 	xjac[1,32] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][32] = 	xjac[1,32] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,32] = 	xjac[1,32] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][32] = 	xjac[1,32] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][32] = 	xjac[1,32] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][32] = 	xjac[1,32] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][32] = 	xjac[1,32] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][32] = 	xjac[1,32] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,32] = 	xjac[1,32] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,32] = 	xjac[1,32] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][32] = 	xjac[1,32] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,32] = 	xjac[1,32] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,32] = 	xjac[1,32] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,32] = 	xjac[1,32] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,32] = 	xjac[1,32] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][32] = 	xjac[1,32] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][32] = 	xjac[1,32] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][32] = 	xjac[1,32] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][32] = 	xjac[1,32] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][32] = 	xjac[1,32] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][32] = 	xjac[1,32] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][32] = 	xjac[1,32] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][32] = 	xjac[1,32] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][32] = 	xjac[1,32] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][32] = 	xjac[1,32] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][32] = 	xjac[1,32] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][32] = 	xjac[1,32] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][32] = 	xjac[1,32] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][32] = 	xjac[1,32] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][32] = 	xjac[1,32] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][32] = 	xjac[1,32] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][32] = 	xjac[1,32] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][32] = 	xjac[1,32] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][32] = 	xjac[1,32] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][32] = 	xjac[1,32] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][32] = 	xjac[1,32] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,32] = 	xjac[1,32] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,32] = 	xjac[1,32] - (2.3e-12*1.0*pow(x[32],0.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,32] = 	xjac[1,32] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,32] = 	xjac[1,32] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][32] = 	xjac[1,32] + (1.2e-11*1.0*pow(x[32],0.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][32] = 	xjac[1,32] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][32] = 	xjac[1,32] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][32] = 	xjac[1,32] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][32] = 	xjac[1,32] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][32] = 	xjac[1,32] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][32] = 	xjac[1,32] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][32] = 	xjac[1,32] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,32] = 	xjac[1,32] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,32] = 	xjac[1,32] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][32] = 	xjac[1,32] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,32] = 	xjac[1,32] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][32] = 	xjac[1,32] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][32] = 	xjac[1,32] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][32] = 	xjac[1,32] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][32] = 	xjac[1,32] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][32] = 	xjac[1,32] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][32] = 	xjac[1,32] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][32] = 	xjac[1,32] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][32] = 	xjac[1,32] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][32] = 	xjac[1,32] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][32] = 	xjac[1,32] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,33] = 	xjac[1,33] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][33] = 	xjac[1,33] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][33] = 	xjac[1,33] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,33] = 	xjac[1,33] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][33] = 	xjac[1,33] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][33] = 	xjac[1,33] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][33] = 	xjac[1,33] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][33] = 	xjac[1,33] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][33] = 	xjac[1,33] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,33] = 	xjac[1,33] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,33] = 	xjac[1,33] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][33] = 	xjac[1,33] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,33] = 	xjac[1,33] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,33] = 	xjac[1,33] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,33] = 	xjac[1,33] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,33] = 	xjac[1,33] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][33] = 	xjac[1,33] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][33] = 	xjac[1,33] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][33] = 	xjac[1,33] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][33] = 	xjac[1,33] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][33] = 	xjac[1,33] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][33] = 	xjac[1,33] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][33] = 	xjac[1,33] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][33] = 	xjac[1,33] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][33] = 	xjac[1,33] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][33] = 	xjac[1,33] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][33] = 	xjac[1,33] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][33] = 	xjac[1,33] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][33] = 	xjac[1,33] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][33] = 	xjac[1,33] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][33] = 	xjac[1,33] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][33] = 	xjac[1,33] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][33] = 	xjac[1,33] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][33] = 	xjac[1,33] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][33] = 	xjac[1,33] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][33] = 	xjac[1,33] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][33] = 	xjac[1,33] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,33] = 	xjac[1,33] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,33] = 	xjac[1,33] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,33] = 	xjac[1,33] - (1.3e-12*pow(300/T,0)*exp(-330/T)*1.0*pow(x[33],0.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,33] = 	xjac[1,33] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][33] = 	xjac[1,33] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][33] = 	xjac[1,33] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][33] = 	xjac[1,33] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][33] = 	xjac[1,33] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][33] = 	xjac[1,33] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][33] = 	xjac[1,33] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][33] = 	xjac[1,33] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][33] = 	xjac[1,33] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,33] = 	xjac[1,33] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,33] = 	xjac[1,33] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][33] = 	xjac[1,33] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,33] = 	xjac[1,33] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][33] = 	xjac[1,33] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][33] = 	xjac[1,33] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][33] = 	xjac[1,33] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][33] = 	xjac[1,33] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][33] = 	xjac[1,33] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][33] = 	xjac[1,33] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][33] = 	xjac[1,33] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][33] = 	xjac[1,33] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][33] = 	xjac[1,33] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][33] = 	xjac[1,33] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,34] = 	xjac[1,34] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][34] = 	xjac[1,34] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][34] = 	xjac[1,34] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,34] = 	xjac[1,34] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][34] = 	xjac[1,34] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][34] = 	xjac[1,34] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][34] = 	xjac[1,34] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][34] = 	xjac[1,34] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][34] = 	xjac[1,34] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,34] = 	xjac[1,34] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,34] = 	xjac[1,34] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][34] = 	xjac[1,34] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,34] = 	xjac[1,34] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,34] = 	xjac[1,34] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,34] = 	xjac[1,34] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,34] = 	xjac[1,34] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][34] = 	xjac[1,34] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][34] = 	xjac[1,34] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][34] = 	xjac[1,34] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][34] = 	xjac[1,34] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][34] = 	xjac[1,34] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][34] = 	xjac[1,34] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][34] = 	xjac[1,34] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][34] = 	xjac[1,34] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][34] = 	xjac[1,34] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][34] = 	xjac[1,34] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][34] = 	xjac[1,34] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][34] = 	xjac[1,34] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][34] = 	xjac[1,34] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][34] = 	xjac[1,34] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][34] = 	xjac[1,34] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][34] = 	xjac[1,34] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][34] = 	xjac[1,34] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][34] = 	xjac[1,34] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][34] = 	xjac[1,34] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][34] = 	xjac[1,34] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][34] = 	xjac[1,34] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,34] = 	xjac[1,34] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,34] = 	xjac[1,34] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,34] = 	xjac[1,34] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,34] = 	xjac[1,34] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][34] = 	xjac[1,34] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][34] = 	xjac[1,34] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][34] = 	xjac[1,34] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][34] = 	xjac[1,34] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][34] = 	xjac[1,34] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][34] = 	xjac[1,34] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][34] = 	xjac[1,34] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][34] = 	xjac[1,34] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,34] = 	xjac[1,34] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,34] = 	xjac[1,34] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][34] = 	xjac[1,34] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,34] = 	xjac[1,34] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][34] = 	xjac[1,34] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][34] = 	xjac[1,34] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][34] = 	xjac[1,34] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][34] = 	xjac[1,34] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][34] = 	xjac[1,34] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][34] = 	xjac[1,34] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][34] = 	xjac[1,34] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][34] = 	xjac[1,34] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][34] = 	xjac[1,34] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][34] = 	xjac[1,34] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,35] = 	xjac[1,35] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][35] = 	xjac[1,35] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][35] = 	xjac[1,35] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,35] = 	xjac[1,35] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][35] = 	xjac[1,35] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][35] = 	xjac[1,35] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][35] = 	xjac[1,35] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][35] = 	xjac[1,35] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][35] = 	xjac[1,35] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,35] = 	xjac[1,35] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,35] = 	xjac[1,35] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][35] = 	xjac[1,35] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,35] = 	xjac[1,35] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,35] = 	xjac[1,35] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,35] = 	xjac[1,35] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,35] = 	xjac[1,35] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][35] = 	xjac[1,35] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][35] = 	xjac[1,35] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][35] = 	xjac[1,35] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][35] = 	xjac[1,35] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][35] = 	xjac[1,35] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][35] = 	xjac[1,35] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][35] = 	xjac[1,35] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][35] = 	xjac[1,35] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][35] = 	xjac[1,35] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][35] = 	xjac[1,35] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][35] = 	xjac[1,35] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][35] = 	xjac[1,35] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][35] = 	xjac[1,35] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][35] = 	xjac[1,35] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][35] = 	xjac[1,35] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][35] = 	xjac[1,35] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][35] = 	xjac[1,35] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][35] = 	xjac[1,35] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][35] = 	xjac[1,35] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][35] = 	xjac[1,35] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][35] = 	xjac[1,35] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,35] = 	xjac[1,35] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,35] = 	xjac[1,35] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,35] = 	xjac[1,35] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,35] = 	xjac[1,35] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][35] = 	xjac[1,35] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][35] = 	xjac[1,35] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][35] = 	xjac[1,35] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][35] = 	xjac[1,35] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][35] = 	xjac[1,35] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][35] = 	xjac[1,35] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][35] = 	xjac[1,35] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][35] = 	xjac[1,35] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,35] = 	xjac[1,35] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,35] = 	xjac[1,35] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][35] = 	xjac[1,35] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,35] = 	xjac[1,35] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][35] = 	xjac[1,35] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][35] = 	xjac[1,35] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][35] = 	xjac[1,35] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][35] = 	xjac[1,35] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][35] = 	xjac[1,35] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][35] = 	xjac[1,35] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][35] = 	xjac[1,35] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][35] = 	xjac[1,35] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][35] = 	xjac[1,35] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][35] = 	xjac[1,35] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,36] = 	xjac[1,36] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][36] = 	xjac[1,36] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][36] = 	xjac[1,36] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,36] = 	xjac[1,36] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][36] = 	xjac[1,36] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][36] = 	xjac[1,36] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][36] = 	xjac[1,36] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][36] = 	xjac[1,36] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][36] = 	xjac[1,36] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,36] = 	xjac[1,36] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,36] = 	xjac[1,36] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][36] = 	xjac[1,36] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,36] = 	xjac[1,36] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,36] = 	xjac[1,36] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,36] = 	xjac[1,36] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,36] = 	xjac[1,36] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][36] = 	xjac[1,36] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][36] = 	xjac[1,36] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][36] = 	xjac[1,36] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][36] = 	xjac[1,36] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][36] = 	xjac[1,36] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][36] = 	xjac[1,36] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][36] = 	xjac[1,36] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][36] = 	xjac[1,36] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][36] = 	xjac[1,36] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][36] = 	xjac[1,36] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][36] = 	xjac[1,36] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][36] = 	xjac[1,36] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][36] = 	xjac[1,36] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][36] = 	xjac[1,36] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][36] = 	xjac[1,36] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][36] = 	xjac[1,36] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][36] = 	xjac[1,36] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][36] = 	xjac[1,36] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][36] = 	xjac[1,36] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][36] = 	xjac[1,36] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][36] = 	xjac[1,36] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,36] = 	xjac[1,36] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,36] = 	xjac[1,36] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,36] = 	xjac[1,36] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,36] = 	xjac[1,36] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][36] = 	xjac[1,36] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][36] = 	xjac[1,36] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][36] = 	xjac[1,36] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][36] = 	xjac[1,36] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][36] = 	xjac[1,36] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][36] = 	xjac[1,36] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][36] = 	xjac[1,36] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][36] = 	xjac[1,36] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,36] = 	xjac[1,36] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,36] = 	xjac[1,36] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][36] = 	xjac[1,36] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,36] = 	xjac[1,36] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][36] = 	xjac[1,36] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][36] = 	xjac[1,36] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][36] = 	xjac[1,36] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][36] = 	xjac[1,36] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][36] = 	xjac[1,36] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][36] = 	xjac[1,36] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][36] = 	xjac[1,36] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][36] = 	xjac[1,36] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][36] = 	xjac[1,36] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][36] = 	xjac[1,36] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,37] = 	xjac[1,37] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][37] = 	xjac[1,37] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][37] = 	xjac[1,37] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,37] = 	xjac[1,37] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][37] = 	xjac[1,37] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][37] = 	xjac[1,37] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][37] = 	xjac[1,37] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][37] = 	xjac[1,37] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][37] = 	xjac[1,37] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,37] = 	xjac[1,37] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,37] = 	xjac[1,37] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][37] = 	xjac[1,37] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,37] = 	xjac[1,37] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,37] = 	xjac[1,37] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,37] = 	xjac[1,37] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,37] = 	xjac[1,37] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][37] = 	xjac[1,37] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][37] = 	xjac[1,37] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][37] = 	xjac[1,37] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][37] = 	xjac[1,37] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][37] = 	xjac[1,37] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][37] = 	xjac[1,37] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][37] = 	xjac[1,37] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][37] = 	xjac[1,37] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][37] = 	xjac[1,37] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][37] = 	xjac[1,37] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][37] = 	xjac[1,37] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][37] = 	xjac[1,37] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][37] = 	xjac[1,37] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][37] = 	xjac[1,37] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][37] = 	xjac[1,37] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][37] = 	xjac[1,37] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][37] = 	xjac[1,37] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][37] = 	xjac[1,37] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][37] = 	xjac[1,37] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][37] = 	xjac[1,37] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][37] = 	xjac[1,37] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,37] = 	xjac[1,37] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,37] = 	xjac[1,37] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,37] = 	xjac[1,37] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,37] = 	xjac[1,37] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][37] = 	xjac[1,37] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][37] = 	xjac[1,37] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][37] = 	xjac[1,37] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][37] = 	xjac[1,37] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][37] = 	xjac[1,37] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][37] = 	xjac[1,37] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][37] = 	xjac[1,37] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][37] = 	xjac[1,37] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,37] = 	xjac[1,37] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,37] = 	xjac[1,37] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][37] = 	xjac[1,37] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,37] = 	xjac[1,37] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][37] = 	xjac[1,37] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][37] = 	xjac[1,37] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][37] = 	xjac[1,37] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][37] = 	xjac[1,37] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][37] = 	xjac[1,37] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][37] = 	xjac[1,37] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][37] = 	xjac[1,37] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][37] = 	xjac[1,37] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][37] = 	xjac[1,37] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][37] = 	xjac[1,37] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,38] = 	xjac[1,38] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][38] = 	xjac[1,38] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][38] = 	xjac[1,38] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,38] = 	xjac[1,38] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][38] = 	xjac[1,38] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][38] = 	xjac[1,38] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][38] = 	xjac[1,38] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][38] = 	xjac[1,38] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][38] = 	xjac[1,38] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,38] = 	xjac[1,38] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,38] = 	xjac[1,38] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][38] = 	xjac[1,38] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,38] = 	xjac[1,38] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,38] = 	xjac[1,38] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,38] = 	xjac[1,38] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,38] = 	xjac[1,38] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][38] = 	xjac[1,38] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][38] = 	xjac[1,38] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][38] = 	xjac[1,38] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][38] = 	xjac[1,38] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][38] = 	xjac[1,38] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][38] = 	xjac[1,38] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][38] = 	xjac[1,38] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][38] = 	xjac[1,38] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][38] = 	xjac[1,38] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][38] = 	xjac[1,38] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][38] = 	xjac[1,38] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][38] = 	xjac[1,38] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][38] = 	xjac[1,38] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][38] = 	xjac[1,38] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][38] = 	xjac[1,38] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][38] = 	xjac[1,38] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][38] = 	xjac[1,38] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][38] = 	xjac[1,38] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][38] = 	xjac[1,38] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][38] = 	xjac[1,38] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][38] = 	xjac[1,38] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,38] = 	xjac[1,38] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,38] = 	xjac[1,38] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,38] = 	xjac[1,38] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,38] = 	xjac[1,38] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][38] = 	xjac[1,38] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][38] = 	xjac[1,38] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][38] = 	xjac[1,38] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][38] = 	xjac[1,38] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][38] = 	xjac[1,38] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][38] = 	xjac[1,38] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][38] = 	xjac[1,38] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][38] = 	xjac[1,38] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,38] = 	xjac[1,38] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,38] = 	xjac[1,38] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][38] = 	xjac[1,38] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,38] = 	xjac[1,38] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][38] = 	xjac[1,38] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][38] = 	xjac[1,38] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][38] = 	xjac[1,38] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][38] = 	xjac[1,38] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][38] = 	xjac[1,38] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][38] = 	xjac[1,38] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][38] = 	xjac[1,38] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][38] = 	xjac[1,38] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][38] = 	xjac[1,38] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][38] = 	xjac[1,38] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,39] = 	xjac[1,39] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][39] = 	xjac[1,39] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][39] = 	xjac[1,39] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,39] = 	xjac[1,39] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][39] = 	xjac[1,39] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][39] = 	xjac[1,39] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][39] = 	xjac[1,39] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][39] = 	xjac[1,39] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][39] = 	xjac[1,39] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,39] = 	xjac[1,39] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,39] = 	xjac[1,39] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][39] = 	xjac[1,39] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,39] = 	xjac[1,39] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,39] = 	xjac[1,39] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,39] = 	xjac[1,39] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,39] = 	xjac[1,39] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][39] = 	xjac[1,39] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][39] = 	xjac[1,39] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][39] = 	xjac[1,39] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][39] = 	xjac[1,39] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][39] = 	xjac[1,39] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][39] = 	xjac[1,39] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][39] = 	xjac[1,39] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][39] = 	xjac[1,39] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][39] = 	xjac[1,39] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][39] = 	xjac[1,39] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][39] = 	xjac[1,39] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][39] = 	xjac[1,39] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][39] = 	xjac[1,39] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][39] = 	xjac[1,39] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][39] = 	xjac[1,39] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][39] = 	xjac[1,39] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][39] = 	xjac[1,39] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][39] = 	xjac[1,39] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][39] = 	xjac[1,39] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][39] = 	xjac[1,39] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][39] = 	xjac[1,39] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,39] = 	xjac[1,39] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,39] = 	xjac[1,39] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,39] = 	xjac[1,39] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,39] = 	xjac[1,39] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][39] = 	xjac[1,39] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][39] = 	xjac[1,39] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][39] = 	xjac[1,39] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][39] = 	xjac[1,39] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][39] = 	xjac[1,39] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][39] = 	xjac[1,39] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][39] = 	xjac[1,39] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][39] = 	xjac[1,39] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,39] = 	xjac[1,39] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,39] = 	xjac[1,39] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][39] = 	xjac[1,39] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,39] = 	xjac[1,39] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][39] = 	xjac[1,39] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][39] = 	xjac[1,39] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][39] = 	xjac[1,39] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][39] = 	xjac[1,39] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][39] = 	xjac[1,39] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][39] = 	xjac[1,39] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][39] = 	xjac[1,39] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][39] = 	xjac[1,39] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][39] = 	xjac[1,39] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][39] = 	xjac[1,39] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,40] = 	xjac[1,40] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][40] = 	xjac[1,40] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][40] = 	xjac[1,40] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,40] = 	xjac[1,40] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][40] = 	xjac[1,40] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][40] = 	xjac[1,40] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][40] = 	xjac[1,40] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][40] = 	xjac[1,40] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][40] = 	xjac[1,40] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,40] = 	xjac[1,40] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,40] = 	xjac[1,40] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][40] = 	xjac[1,40] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,40] = 	xjac[1,40] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,40] = 	xjac[1,40] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,40] = 	xjac[1,40] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,40] = 	xjac[1,40] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][40] = 	xjac[1,40] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][40] = 	xjac[1,40] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][40] = 	xjac[1,40] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][40] = 	xjac[1,40] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][40] = 	xjac[1,40] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][40] = 	xjac[1,40] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][40] = 	xjac[1,40] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][40] = 	xjac[1,40] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][40] = 	xjac[1,40] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][40] = 	xjac[1,40] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][40] = 	xjac[1,40] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][40] = 	xjac[1,40] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][40] = 	xjac[1,40] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][40] = 	xjac[1,40] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][40] = 	xjac[1,40] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][40] = 	xjac[1,40] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][40] = 	xjac[1,40] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][40] = 	xjac[1,40] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][40] = 	xjac[1,40] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][40] = 	xjac[1,40] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][40] = 	xjac[1,40] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,40] = 	xjac[1,40] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,40] = 	xjac[1,40] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,40] = 	xjac[1,40] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,40] = 	xjac[1,40] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][40] = 	xjac[1,40] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][40] = 	xjac[1,40] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][40] = 	xjac[1,40] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][40] = 	xjac[1,40] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][40] = 	xjac[1,40] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][40] = 	xjac[1,40] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][40] = 	xjac[1,40] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][40] = 	xjac[1,40] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,40] = 	xjac[1,40] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,40] = 	xjac[1,40] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][40] = 	xjac[1,40] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,40] = 	xjac[1,40] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][40] = 	xjac[1,40] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][40] = 	xjac[1,40] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][40] = 	xjac[1,40] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][40] = 	xjac[1,40] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][40] = 	xjac[1,40] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][40] = 	xjac[1,40] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][40] = 	xjac[1,40] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][40] = 	xjac[1,40] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][40] = 	xjac[1,40] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][40] = 	xjac[1,40] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,41] = 	xjac[1,41] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][41] = 	xjac[1,41] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][41] = 	xjac[1,41] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,41] = 	xjac[1,41] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][41] = 	xjac[1,41] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][41] = 	xjac[1,41] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][41] = 	xjac[1,41] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][41] = 	xjac[1,41] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][41] = 	xjac[1,41] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,41] = 	xjac[1,41] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,41] = 	xjac[1,41] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][41] = 	xjac[1,41] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,41] = 	xjac[1,41] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,41] = 	xjac[1,41] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,41] = 	xjac[1,41] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,41] = 	xjac[1,41] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][41] = 	xjac[1,41] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][41] = 	xjac[1,41] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][41] = 	xjac[1,41] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][41] = 	xjac[1,41] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][41] = 	xjac[1,41] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][41] = 	xjac[1,41] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][41] = 	xjac[1,41] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][41] = 	xjac[1,41] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][41] = 	xjac[1,41] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][41] = 	xjac[1,41] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][41] = 	xjac[1,41] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][41] = 	xjac[1,41] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][41] = 	xjac[1,41] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][41] = 	xjac[1,41] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][41] = 	xjac[1,41] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][41] = 	xjac[1,41] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][41] = 	xjac[1,41] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][41] = 	xjac[1,41] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][41] = 	xjac[1,41] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][41] = 	xjac[1,41] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][41] = 	xjac[1,41] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,41] = 	xjac[1,41] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,41] = 	xjac[1,41] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,41] = 	xjac[1,41] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,41] = 	xjac[1,41] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][41] = 	xjac[1,41] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][41] = 	xjac[1,41] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][41] = 	xjac[1,41] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][41] = 	xjac[1,41] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][41] = 	xjac[1,41] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][41] = 	xjac[1,41] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][41] = 	xjac[1,41] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][41] = 	xjac[1,41] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,41] = 	xjac[1,41] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,41] = 	xjac[1,41] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][41] = 	xjac[1,41] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,41] = 	xjac[1,41] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][41] = 	xjac[1,41] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][41] = 	xjac[1,41] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][41] = 	xjac[1,41] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][41] = 	xjac[1,41] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][41] = 	xjac[1,41] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][41] = 	xjac[1,41] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][41] = 	xjac[1,41] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][41] = 	xjac[1,41] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][41] = 	xjac[1,41] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][41] = 	xjac[1,41] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,42] = 	xjac[1,42] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][42] = 	xjac[1,42] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][42] = 	xjac[1,42] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,42] = 	xjac[1,42] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][42] = 	xjac[1,42] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][42] = 	xjac[1,42] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][42] = 	xjac[1,42] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][42] = 	xjac[1,42] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][42] = 	xjac[1,42] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,42] = 	xjac[1,42] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,42] = 	xjac[1,42] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][42] = 	xjac[1,42] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,42] = 	xjac[1,42] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,42] = 	xjac[1,42] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,42] = 	xjac[1,42] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,42] = 	xjac[1,42] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][42] = 	xjac[1,42] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][42] = 	xjac[1,42] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][42] = 	xjac[1,42] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][42] = 	xjac[1,42] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][42] = 	xjac[1,42] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][42] = 	xjac[1,42] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][42] = 	xjac[1,42] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][42] = 	xjac[1,42] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][42] = 	xjac[1,42] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][42] = 	xjac[1,42] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][42] = 	xjac[1,42] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][42] = 	xjac[1,42] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][42] = 	xjac[1,42] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][42] = 	xjac[1,42] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][42] = 	xjac[1,42] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][42] = 	xjac[1,42] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][42] = 	xjac[1,42] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][42] = 	xjac[1,42] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][42] = 	xjac[1,42] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][42] = 	xjac[1,42] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][42] = 	xjac[1,42] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,42] = 	xjac[1,42] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,42] = 	xjac[1,42] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,42] = 	xjac[1,42] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,42] = 	xjac[1,42] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][42] = 	xjac[1,42] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][42] = 	xjac[1,42] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][42] = 	xjac[1,42] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][42] = 	xjac[1,42] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][42] = 	xjac[1,42] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][42] = 	xjac[1,42] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][42] = 	xjac[1,42] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][42] = 	xjac[1,42] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,42] = 	xjac[1,42] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,42] = 	xjac[1,42] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][42] = 	xjac[1,42] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,42] = 	xjac[1,42] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][42] = 	xjac[1,42] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][42] = 	xjac[1,42] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][42] = 	xjac[1,42] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][42] = 	xjac[1,42] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][42] = 	xjac[1,42] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][42] = 	xjac[1,42] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][42] = 	xjac[1,42] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][42] = 	xjac[1,42] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][42] = 	xjac[1,42] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][42] = 	xjac[1,42] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,43] = 	xjac[1,43] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][43] = 	xjac[1,43] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][43] = 	xjac[1,43] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,43] = 	xjac[1,43] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][43] = 	xjac[1,43] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][43] = 	xjac[1,43] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][43] = 	xjac[1,43] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][43] = 	xjac[1,43] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][43] = 	xjac[1,43] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,43] = 	xjac[1,43] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,43] = 	xjac[1,43] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][43] = 	xjac[1,43] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,43] = 	xjac[1,43] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,43] = 	xjac[1,43] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,43] = 	xjac[1,43] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,43] = 	xjac[1,43] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][43] = 	xjac[1,43] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][43] = 	xjac[1,43] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][43] = 	xjac[1,43] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][43] = 	xjac[1,43] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][43] = 	xjac[1,43] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][43] = 	xjac[1,43] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][43] = 	xjac[1,43] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][43] = 	xjac[1,43] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][43] = 	xjac[1,43] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][43] = 	xjac[1,43] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][43] = 	xjac[1,43] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][43] = 	xjac[1,43] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][43] = 	xjac[1,43] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][43] = 	xjac[1,43] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][43] = 	xjac[1,43] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][43] = 	xjac[1,43] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][43] = 	xjac[1,43] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][43] = 	xjac[1,43] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][43] = 	xjac[1,43] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][43] = 	xjac[1,43] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][43] = 	xjac[1,43] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,43] = 	xjac[1,43] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,43] = 	xjac[1,43] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,43] = 	xjac[1,43] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,43] = 	xjac[1,43] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][43] = 	xjac[1,43] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][43] = 	xjac[1,43] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][43] = 	xjac[1,43] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][43] = 	xjac[1,43] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][43] = 	xjac[1,43] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][43] = 	xjac[1,43] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][43] = 	xjac[1,43] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][43] = 	xjac[1,43] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,43] = 	xjac[1,43] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,43] = 	xjac[1,43] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][43] = 	xjac[1,43] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,43] = 	xjac[1,43] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][43] = 	xjac[1,43] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][43] = 	xjac[1,43] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][43] = 	xjac[1,43] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][43] = 	xjac[1,43] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][43] = 	xjac[1,43] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][43] = 	xjac[1,43] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][43] = 	xjac[1,43] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][43] = 	xjac[1,43] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][43] = 	xjac[1,43] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][43] = 	xjac[1,43] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,44] = 	xjac[1,44] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][44] = 	xjac[1,44] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][44] = 	xjac[1,44] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,44] = 	xjac[1,44] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][44] = 	xjac[1,44] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][44] = 	xjac[1,44] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][44] = 	xjac[1,44] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][44] = 	xjac[1,44] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][44] = 	xjac[1,44] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,44] = 	xjac[1,44] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,44] = 	xjac[1,44] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][44] = 	xjac[1,44] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,44] = 	xjac[1,44] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,44] = 	xjac[1,44] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,44] = 	xjac[1,44] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,44] = 	xjac[1,44] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][44] = 	xjac[1,44] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][44] = 	xjac[1,44] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][44] = 	xjac[1,44] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][44] = 	xjac[1,44] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][44] = 	xjac[1,44] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][44] = 	xjac[1,44] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][44] = 	xjac[1,44] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][44] = 	xjac[1,44] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][44] = 	xjac[1,44] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][44] = 	xjac[1,44] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][44] = 	xjac[1,44] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][44] = 	xjac[1,44] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][44] = 	xjac[1,44] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][44] = 	xjac[1,44] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][44] = 	xjac[1,44] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][44] = 	xjac[1,44] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][44] = 	xjac[1,44] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][44] = 	xjac[1,44] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][44] = 	xjac[1,44] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][44] = 	xjac[1,44] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][44] = 	xjac[1,44] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,44] = 	xjac[1,44] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,44] = 	xjac[1,44] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,44] = 	xjac[1,44] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,44] = 	xjac[1,44] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][44] = 	xjac[1,44] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][44] = 	xjac[1,44] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][44] = 	xjac[1,44] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][44] = 	xjac[1,44] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][44] = 	xjac[1,44] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][44] = 	xjac[1,44] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][44] = 	xjac[1,44] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][44] = 	xjac[1,44] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,44] = 	xjac[1,44] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,44] = 	xjac[1,44] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][44] = 	xjac[1,44] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,44] = 	xjac[1,44] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][44] = 	xjac[1,44] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][44] = 	xjac[1,44] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][44] = 	xjac[1,44] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][44] = 	xjac[1,44] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][44] = 	xjac[1,44] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][44] = 	xjac[1,44] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][44] = 	xjac[1,44] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][44] = 	xjac[1,44] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][44] = 	xjac[1,44] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][44] = 	xjac[1,44] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,45] = 	xjac[1,45] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][45] = 	xjac[1,45] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][45] = 	xjac[1,45] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,45] = 	xjac[1,45] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][45] = 	xjac[1,45] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][45] = 	xjac[1,45] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][45] = 	xjac[1,45] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][45] = 	xjac[1,45] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][45] = 	xjac[1,45] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,45] = 	xjac[1,45] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,45] = 	xjac[1,45] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][45] = 	xjac[1,45] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,45] = 	xjac[1,45] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,45] = 	xjac[1,45] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,45] = 	xjac[1,45] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,45] = 	xjac[1,45] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][45] = 	xjac[1,45] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][45] = 	xjac[1,45] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][45] = 	xjac[1,45] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][45] = 	xjac[1,45] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][45] = 	xjac[1,45] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][45] = 	xjac[1,45] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][45] = 	xjac[1,45] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][45] = 	xjac[1,45] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][45] = 	xjac[1,45] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][45] = 	xjac[1,45] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][45] = 	xjac[1,45] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][45] = 	xjac[1,45] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][45] = 	xjac[1,45] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][45] = 	xjac[1,45] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][45] = 	xjac[1,45] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][45] = 	xjac[1,45] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][45] = 	xjac[1,45] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][45] = 	xjac[1,45] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][45] = 	xjac[1,45] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][45] = 	xjac[1,45] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][45] = 	xjac[1,45] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,45] = 	xjac[1,45] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,45] = 	xjac[1,45] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,45] = 	xjac[1,45] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,45] = 	xjac[1,45] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][45] = 	xjac[1,45] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][45] = 	xjac[1,45] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][45] = 	xjac[1,45] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][45] = 	xjac[1,45] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][45] = 	xjac[1,45] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*2.0*pow(x[45],1.0))#O + (SO)2 = S2O + O2
	xjac[1][45] = 	xjac[1,45] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][45] = 	xjac[1,45] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][45] = 	xjac[1,45] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,45] = 	xjac[1,45] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,45] = 	xjac[1,45] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][45] = 	xjac[1,45] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,45] = 	xjac[1,45] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][45] = 	xjac[1,45] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][45] = 	xjac[1,45] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][45] = 	xjac[1,45] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][45] = 	xjac[1,45] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][45] = 	xjac[1,45] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][45] = 	xjac[1,45] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][45] = 	xjac[1,45] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][45] = 	xjac[1,45] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][45] = 	xjac[1,45] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][45] = 	xjac[1,45] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,46] = 	xjac[1,46] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][46] = 	xjac[1,46] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][46] = 	xjac[1,46] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,46] = 	xjac[1,46] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][46] = 	xjac[1,46] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][46] = 	xjac[1,46] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][46] = 	xjac[1,46] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][46] = 	xjac[1,46] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][46] = 	xjac[1,46] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,46] = 	xjac[1,46] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,46] = 	xjac[1,46] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][46] = 	xjac[1,46] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,46] = 	xjac[1,46] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,46] = 	xjac[1,46] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,46] = 	xjac[1,46] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,46] = 	xjac[1,46] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][46] = 	xjac[1,46] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][46] = 	xjac[1,46] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][46] = 	xjac[1,46] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][46] = 	xjac[1,46] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][46] = 	xjac[1,46] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][46] = 	xjac[1,46] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][46] = 	xjac[1,46] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][46] = 	xjac[1,46] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][46] = 	xjac[1,46] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][46] = 	xjac[1,46] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][46] = 	xjac[1,46] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][46] = 	xjac[1,46] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][46] = 	xjac[1,46] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][46] = 	xjac[1,46] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][46] = 	xjac[1,46] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][46] = 	xjac[1,46] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][46] = 	xjac[1,46] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][46] = 	xjac[1,46] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][46] = 	xjac[1,46] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][46] = 	xjac[1,46] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][46] = 	xjac[1,46] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,46] = 	xjac[1,46] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,46] = 	xjac[1,46] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,46] = 	xjac[1,46] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,46] = 	xjac[1,46] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][46] = 	xjac[1,46] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][46] = 	xjac[1,46] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][46] = 	xjac[1,46] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][46] = 	xjac[1,46] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][46] = 	xjac[1,46] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][46] = 	xjac[1,46] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][46] = 	xjac[1,46] + (1.3e-10*pow(x[2],1.0)*1.0*pow(x[46],0.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][46] = 	xjac[1,46] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,46] = 	xjac[1,46] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,46] = 	xjac[1,46] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][46] = 	xjac[1,46] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,46] = 	xjac[1,46] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][46] = 	xjac[1,46] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][46] = 	xjac[1,46] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][46] = 	xjac[1,46] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][46] = 	xjac[1,46] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][46] = 	xjac[1,46] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][46] = 	xjac[1,46] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][46] = 	xjac[1,46] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][46] = 	xjac[1,46] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][46] = 	xjac[1,46] + (4.9e-11*pow(x[2],1.0)*1.0*pow(x[46],0.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][46] = 	xjac[1,46] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,47] = 	xjac[1,47] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][47] = 	xjac[1,47] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][47] = 	xjac[1,47] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,47] = 	xjac[1,47] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][47] = 	xjac[1,47] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][47] = 	xjac[1,47] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][47] = 	xjac[1,47] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][47] = 	xjac[1,47] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][47] = 	xjac[1,47] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,47] = 	xjac[1,47] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,47] = 	xjac[1,47] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][47] = 	xjac[1,47] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,47] = 	xjac[1,47] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,47] = 	xjac[1,47] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,47] = 	xjac[1,47] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,47] = 	xjac[1,47] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][47] = 	xjac[1,47] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][47] = 	xjac[1,47] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][47] = 	xjac[1,47] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][47] = 	xjac[1,47] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][47] = 	xjac[1,47] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][47] = 	xjac[1,47] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][47] = 	xjac[1,47] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][47] = 	xjac[1,47] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][47] = 	xjac[1,47] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][47] = 	xjac[1,47] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][47] = 	xjac[1,47] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][47] = 	xjac[1,47] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][47] = 	xjac[1,47] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][47] = 	xjac[1,47] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][47] = 	xjac[1,47] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][47] = 	xjac[1,47] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][47] = 	xjac[1,47] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][47] = 	xjac[1,47] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][47] = 	xjac[1,47] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][47] = 	xjac[1,47] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][47] = 	xjac[1,47] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,47] = 	xjac[1,47] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,47] = 	xjac[1,47] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,47] = 	xjac[1,47] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,47] = 	xjac[1,47] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][47] = 	xjac[1,47] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][47] = 	xjac[1,47] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][47] = 	xjac[1,47] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][47] = 	xjac[1,47] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][47] = 	xjac[1,47] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][47] = 	xjac[1,47] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][47] = 	xjac[1,47] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][47] = 	xjac[1,47] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,47] = 	xjac[1,47] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,47] = 	xjac[1,47] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][47] = 	xjac[1,47] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,47] = 	xjac[1,47] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][47] = 	xjac[1,47] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][47] = 	xjac[1,47] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][47] = 	xjac[1,47] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][47] = 	xjac[1,47] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][47] = 	xjac[1,47] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][47] = 	xjac[1,47] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][47] = 	xjac[1,47] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][47] = 	xjac[1,47] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][47] = 	xjac[1,47] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][47] = 	xjac[1,47] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,48] = 	xjac[1,48] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][48] = 	xjac[1,48] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][48] = 	xjac[1,48] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,48] = 	xjac[1,48] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][48] = 	xjac[1,48] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][48] = 	xjac[1,48] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][48] = 	xjac[1,48] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][48] = 	xjac[1,48] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][48] = 	xjac[1,48] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,48] = 	xjac[1,48] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,48] = 	xjac[1,48] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][48] = 	xjac[1,48] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,48] = 	xjac[1,48] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,48] = 	xjac[1,48] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,48] = 	xjac[1,48] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,48] = 	xjac[1,48] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][48] = 	xjac[1,48] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][48] = 	xjac[1,48] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][48] = 	xjac[1,48] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][48] = 	xjac[1,48] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][48] = 	xjac[1,48] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][48] = 	xjac[1,48] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][48] = 	xjac[1,48] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][48] = 	xjac[1,48] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][48] = 	xjac[1,48] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][48] = 	xjac[1,48] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][48] = 	xjac[1,48] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][48] = 	xjac[1,48] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][48] = 	xjac[1,48] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][48] = 	xjac[1,48] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][48] = 	xjac[1,48] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][48] = 	xjac[1,48] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][48] = 	xjac[1,48] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][48] = 	xjac[1,48] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][48] = 	xjac[1,48] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][48] = 	xjac[1,48] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][48] = 	xjac[1,48] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,48] = 	xjac[1,48] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,48] = 	xjac[1,48] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,48] = 	xjac[1,48] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,48] = 	xjac[1,48] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][48] = 	xjac[1,48] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][48] = 	xjac[1,48] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][48] = 	xjac[1,48] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][48] = 	xjac[1,48] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][48] = 	xjac[1,48] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][48] = 	xjac[1,48] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][48] = 	xjac[1,48] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][48] = 	xjac[1,48] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,48] = 	xjac[1,48] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*1.0*pow(x[48],0.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,48] = 	xjac[1,48] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][48] = 	xjac[1,48] + (2e-16*1.0*pow(x[48],0.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,48] = 	xjac[1,48] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*1.0*pow(x[48],0.0))#O2(1d) + N = NO + O
	xjac[1][48] = 	xjac[1,48] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][48] = 	xjac[1,48] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][48] = 	xjac[1,48] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][48] = 	xjac[1,48] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][48] = 	xjac[1,48] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][48] = 	xjac[1,48] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][48] = 	xjac[1,48] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][48] = 	xjac[1,48] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][48] = 	xjac[1,48] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][48] = 	xjac[1,48] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,49] = 	xjac[1,49] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][49] = 	xjac[1,49] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][49] = 	xjac[1,49] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,49] = 	xjac[1,49] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][49] = 	xjac[1,49] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][49] = 	xjac[1,49] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][49] = 	xjac[1,49] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][49] = 	xjac[1,49] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][49] = 	xjac[1,49] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,49] = 	xjac[1,49] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,49] = 	xjac[1,49] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][49] = 	xjac[1,49] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,49] = 	xjac[1,49] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,49] = 	xjac[1,49] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,49] = 	xjac[1,49] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,49] = 	xjac[1,49] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][49] = 	xjac[1,49] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][49] = 	xjac[1,49] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][49] = 	xjac[1,49] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][49] = 	xjac[1,49] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][49] = 	xjac[1,49] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][49] = 	xjac[1,49] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][49] = 	xjac[1,49] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][49] = 	xjac[1,49] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][49] = 	xjac[1,49] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][49] = 	xjac[1,49] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][49] = 	xjac[1,49] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][49] = 	xjac[1,49] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][49] = 	xjac[1,49] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][49] = 	xjac[1,49] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][49] = 	xjac[1,49] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][49] = 	xjac[1,49] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][49] = 	xjac[1,49] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][49] = 	xjac[1,49] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][49] = 	xjac[1,49] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][49] = 	xjac[1,49] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][49] = 	xjac[1,49] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,49] = 	xjac[1,49] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,49] = 	xjac[1,49] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,49] = 	xjac[1,49] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,49] = 	xjac[1,49] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][49] = 	xjac[1,49] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][49] = 	xjac[1,49] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][49] = 	xjac[1,49] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][49] = 	xjac[1,49] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][49] = 	xjac[1,49] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][49] = 	xjac[1,49] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][49] = 	xjac[1,49] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][49] = 	xjac[1,49] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,49] = 	xjac[1,49] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,49] = 	xjac[1,49] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][49] = 	xjac[1,49] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,49] = 	xjac[1,49] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][49] = 	xjac[1,49] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*1.0*pow(x[49],0.0))#O3 + NO = NO2 + O2
	xjac[1][49] = 	xjac[1,49] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][49] = 	xjac[1,49] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][49] = 	xjac[1,49] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][49] = 	xjac[1,49] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][49] = 	xjac[1,49] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][49] = 	xjac[1,49] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][49] = 	xjac[1,49] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][49] = 	xjac[1,49] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][49] = 	xjac[1,49] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,50] = 	xjac[1,50] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][50] = 	xjac[1,50] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][50] = 	xjac[1,50] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,50] = 	xjac[1,50] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][50] = 	xjac[1,50] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][50] = 	xjac[1,50] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][50] = 	xjac[1,50] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][50] = 	xjac[1,50] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][50] = 	xjac[1,50] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,50] = 	xjac[1,50] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,50] = 	xjac[1,50] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][50] = 	xjac[1,50] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,50] = 	xjac[1,50] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,50] = 	xjac[1,50] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,50] = 	xjac[1,50] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,50] = 	xjac[1,50] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][50] = 	xjac[1,50] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][50] = 	xjac[1,50] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][50] = 	xjac[1,50] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][50] = 	xjac[1,50] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][50] = 	xjac[1,50] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][50] = 	xjac[1,50] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][50] = 	xjac[1,50] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][50] = 	xjac[1,50] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][50] = 	xjac[1,50] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][50] = 	xjac[1,50] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][50] = 	xjac[1,50] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][50] = 	xjac[1,50] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][50] = 	xjac[1,50] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][50] = 	xjac[1,50] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][50] = 	xjac[1,50] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][50] = 	xjac[1,50] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][50] = 	xjac[1,50] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][50] = 	xjac[1,50] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][50] = 	xjac[1,50] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][50] = 	xjac[1,50] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][50] = 	xjac[1,50] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,50] = 	xjac[1,50] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,50] = 	xjac[1,50] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,50] = 	xjac[1,50] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,50] = 	xjac[1,50] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][50] = 	xjac[1,50] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][50] = 	xjac[1,50] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][50] = 	xjac[1,50] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][50] = 	xjac[1,50] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][50] = 	xjac[1,50] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][50] = 	xjac[1,50] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][50] = 	xjac[1,50] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][50] = 	xjac[1,50] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,50] = 	xjac[1,50] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,50] = 	xjac[1,50] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*1.0*pow(x[50],0.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][50] = 	xjac[1,50] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,50] = 	xjac[1,50] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][50] = 	xjac[1,50] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][50] = 	xjac[1,50] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][50] = 	xjac[1,50] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][50] = 	xjac[1,50] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][50] = 	xjac[1,50] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][50] = 	xjac[1,50] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][50] = 	xjac[1,50] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][50] = 	xjac[1,50] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][50] = 	xjac[1,50] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][50] = 	xjac[1,50] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,51] = 	xjac[1,51] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][51] = 	xjac[1,51] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][51] = 	xjac[1,51] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,51] = 	xjac[1,51] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][51] = 	xjac[1,51] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][51] = 	xjac[1,51] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][51] = 	xjac[1,51] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][51] = 	xjac[1,51] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][51] = 	xjac[1,51] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,51] = 	xjac[1,51] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,51] = 	xjac[1,51] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][51] = 	xjac[1,51] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,51] = 	xjac[1,51] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,51] = 	xjac[1,51] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,51] = 	xjac[1,51] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,51] = 	xjac[1,51] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][51] = 	xjac[1,51] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][51] = 	xjac[1,51] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][51] = 	xjac[1,51] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][51] = 	xjac[1,51] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][51] = 	xjac[1,51] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][51] = 	xjac[1,51] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][51] = 	xjac[1,51] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][51] = 	xjac[1,51] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][51] = 	xjac[1,51] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][51] = 	xjac[1,51] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][51] = 	xjac[1,51] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][51] = 	xjac[1,51] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][51] = 	xjac[1,51] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][51] = 	xjac[1,51] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][51] = 	xjac[1,51] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][51] = 	xjac[1,51] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][51] = 	xjac[1,51] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][51] = 	xjac[1,51] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][51] = 	xjac[1,51] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][51] = 	xjac[1,51] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][51] = 	xjac[1,51] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,51] = 	xjac[1,51] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,51] = 	xjac[1,51] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,51] = 	xjac[1,51] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,51] = 	xjac[1,51] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][51] = 	xjac[1,51] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][51] = 	xjac[1,51] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][51] = 	xjac[1,51] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][51] = 	xjac[1,51] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][51] = 	xjac[1,51] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][51] = 	xjac[1,51] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][51] = 	xjac[1,51] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][51] = 	xjac[1,51] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,51] = 	xjac[1,51] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,51] = 	xjac[1,51] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][51] = 	xjac[1,51] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,51] = 	xjac[1,51] - (9e-17*pow(x[1],1.0)*1.0*pow(x[51],0.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][51] = 	xjac[1,51] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][51] = 	xjac[1,51] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][51] = 	xjac[1,51] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][51] = 	xjac[1,51] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][51] = 	xjac[1,51] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][51] = 	xjac[1,51] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][51] = 	xjac[1,51] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][51] = 	xjac[1,51] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][51] = 	xjac[1,51] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][51] = 	xjac[1,51] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,52] = 	xjac[1,52] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][52] = 	xjac[1,52] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][52] = 	xjac[1,52] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,52] = 	xjac[1,52] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][52] = 	xjac[1,52] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][52] = 	xjac[1,52] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][52] = 	xjac[1,52] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][52] = 	xjac[1,52] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][52] = 	xjac[1,52] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,52] = 	xjac[1,52] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,52] = 	xjac[1,52] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][52] = 	xjac[1,52] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,52] = 	xjac[1,52] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,52] = 	xjac[1,52] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,52] = 	xjac[1,52] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,52] = 	xjac[1,52] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][52] = 	xjac[1,52] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][52] = 	xjac[1,52] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][52] = 	xjac[1,52] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][52] = 	xjac[1,52] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][52] = 	xjac[1,52] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][52] = 	xjac[1,52] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][52] = 	xjac[1,52] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][52] = 	xjac[1,52] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][52] = 	xjac[1,52] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][52] = 	xjac[1,52] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][52] = 	xjac[1,52] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][52] = 	xjac[1,52] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][52] = 	xjac[1,52] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][52] = 	xjac[1,52] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][52] = 	xjac[1,52] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][52] = 	xjac[1,52] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][52] = 	xjac[1,52] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][52] = 	xjac[1,52] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][52] = 	xjac[1,52] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][52] = 	xjac[1,52] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][52] = 	xjac[1,52] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,52] = 	xjac[1,52] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,52] = 	xjac[1,52] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,52] = 	xjac[1,52] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,52] = 	xjac[1,52] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][52] = 	xjac[1,52] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][52] = 	xjac[1,52] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][52] = 	xjac[1,52] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][52] = 	xjac[1,52] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][52] = 	xjac[1,52] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][52] = 	xjac[1,52] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][52] = 	xjac[1,52] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][52] = 	xjac[1,52] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,52] = 	xjac[1,52] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,52] = 	xjac[1,52] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][52] = 	xjac[1,52] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,52] = 	xjac[1,52] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][52] = 	xjac[1,52] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][52] = 	xjac[1,52] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][52] = 	xjac[1,52] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][52] = 	xjac[1,52] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][52] = 	xjac[1,52] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][52] = 	xjac[1,52] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][52] = 	xjac[1,52] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][52] = 	xjac[1,52] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][52] = 	xjac[1,52] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*1.0*pow(x[52],0.0))#O(1D) + N2O = N2 + O2
	xjac[1][52] = 	xjac[1,52] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,53] = 	xjac[1,53] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][53] = 	xjac[1,53] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][53] = 	xjac[1,53] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,53] = 	xjac[1,53] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][53] = 	xjac[1,53] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][53] = 	xjac[1,53] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][53] = 	xjac[1,53] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][53] = 	xjac[1,53] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][53] = 	xjac[1,53] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,53] = 	xjac[1,53] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,53] = 	xjac[1,53] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][53] = 	xjac[1,53] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,53] = 	xjac[1,53] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,53] = 	xjac[1,53] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,53] = 	xjac[1,53] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,53] = 	xjac[1,53] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][53] = 	xjac[1,53] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][53] = 	xjac[1,53] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][53] = 	xjac[1,53] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][53] = 	xjac[1,53] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][53] = 	xjac[1,53] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][53] = 	xjac[1,53] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][53] = 	xjac[1,53] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][53] = 	xjac[1,53] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][53] = 	xjac[1,53] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][53] = 	xjac[1,53] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][53] = 	xjac[1,53] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][53] = 	xjac[1,53] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][53] = 	xjac[1,53] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][53] = 	xjac[1,53] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][53] = 	xjac[1,53] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][53] = 	xjac[1,53] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][53] = 	xjac[1,53] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][53] = 	xjac[1,53] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][53] = 	xjac[1,53] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][53] = 	xjac[1,53] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][53] = 	xjac[1,53] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,53] = 	xjac[1,53] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,53] = 	xjac[1,53] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,53] = 	xjac[1,53] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,53] = 	xjac[1,53] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][53] = 	xjac[1,53] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][53] = 	xjac[1,53] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][53] = 	xjac[1,53] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][53] = 	xjac[1,53] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][53] = 	xjac[1,53] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][53] = 	xjac[1,53] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][53] = 	xjac[1,53] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][53] = 	xjac[1,53] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,53] = 	xjac[1,53] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,53] = 	xjac[1,53] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][53] = 	xjac[1,53] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,53] = 	xjac[1,53] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][53] = 	xjac[1,53] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][53] = 	xjac[1,53] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*1.0*pow(x[53],0.0))#O + NO2 = NO + O2
	xjac[1][53] = 	xjac[1,53] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*1.0*pow(x[53],0.0))#O3 + NO2 = NO3 + O2
	xjac[1][53] = 	xjac[1,53] + (5e-16*pow(x[10],1.0)*1.0*pow(x[53],0.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][53] = 	xjac[1,53] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*1.0*pow(x[53],0.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][53] = 	xjac[1,53] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][53] = 	xjac[1,53] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][53] = 	xjac[1,53] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][53] = 	xjac[1,53] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][53] = 	xjac[1,53] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,54] = 	xjac[1,54] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][54] = 	xjac[1,54] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][54] = 	xjac[1,54] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,54] = 	xjac[1,54] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][54] = 	xjac[1,54] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][54] = 	xjac[1,54] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][54] = 	xjac[1,54] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][54] = 	xjac[1,54] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][54] = 	xjac[1,54] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,54] = 	xjac[1,54] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,54] = 	xjac[1,54] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][54] = 	xjac[1,54] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,54] = 	xjac[1,54] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,54] = 	xjac[1,54] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,54] = 	xjac[1,54] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,54] = 	xjac[1,54] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][54] = 	xjac[1,54] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][54] = 	xjac[1,54] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][54] = 	xjac[1,54] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][54] = 	xjac[1,54] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][54] = 	xjac[1,54] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][54] = 	xjac[1,54] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][54] = 	xjac[1,54] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][54] = 	xjac[1,54] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][54] = 	xjac[1,54] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][54] = 	xjac[1,54] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][54] = 	xjac[1,54] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][54] = 	xjac[1,54] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][54] = 	xjac[1,54] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][54] = 	xjac[1,54] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][54] = 	xjac[1,54] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][54] = 	xjac[1,54] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][54] = 	xjac[1,54] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][54] = 	xjac[1,54] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][54] = 	xjac[1,54] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][54] = 	xjac[1,54] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][54] = 	xjac[1,54] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,54] = 	xjac[1,54] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,54] = 	xjac[1,54] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,54] = 	xjac[1,54] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,54] = 	xjac[1,54] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][54] = 	xjac[1,54] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][54] = 	xjac[1,54] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][54] = 	xjac[1,54] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][54] = 	xjac[1,54] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][54] = 	xjac[1,54] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][54] = 	xjac[1,54] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][54] = 	xjac[1,54] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][54] = 	xjac[1,54] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,54] = 	xjac[1,54] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,54] = 	xjac[1,54] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][54] = 	xjac[1,54] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,54] = 	xjac[1,54] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][54] = 	xjac[1,54] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][54] = 	xjac[1,54] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][54] = 	xjac[1,54] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][54] = 	xjac[1,54] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][54] = 	xjac[1,54] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][54] = 	xjac[1,54] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][54] = 	xjac[1,54] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][54] = 	xjac[1,54] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][54] = 	xjac[1,54] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][54] = 	xjac[1,54] + (5e-19*pow(x[8],1.0)*1.0*pow(x[54],0.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,55] = 	xjac[1,55] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][55] = 	xjac[1,55] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][55] = 	xjac[1,55] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,55] = 	xjac[1,55] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][55] = 	xjac[1,55] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][55] = 	xjac[1,55] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][55] = 	xjac[1,55] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][55] = 	xjac[1,55] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][55] = 	xjac[1,55] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,55] = 	xjac[1,55] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,55] = 	xjac[1,55] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][55] = 	xjac[1,55] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,55] = 	xjac[1,55] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,55] = 	xjac[1,55] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,55] = 	xjac[1,55] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,55] = 	xjac[1,55] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][55] = 	xjac[1,55] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][55] = 	xjac[1,55] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][55] = 	xjac[1,55] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][55] = 	xjac[1,55] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][55] = 	xjac[1,55] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][55] = 	xjac[1,55] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][55] = 	xjac[1,55] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][55] = 	xjac[1,55] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][55] = 	xjac[1,55] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][55] = 	xjac[1,55] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][55] = 	xjac[1,55] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][55] = 	xjac[1,55] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][55] = 	xjac[1,55] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][55] = 	xjac[1,55] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][55] = 	xjac[1,55] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][55] = 	xjac[1,55] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][55] = 	xjac[1,55] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][55] = 	xjac[1,55] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][55] = 	xjac[1,55] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][55] = 	xjac[1,55] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][55] = 	xjac[1,55] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,55] = 	xjac[1,55] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,55] = 	xjac[1,55] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,55] = 	xjac[1,55] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,55] = 	xjac[1,55] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][55] = 	xjac[1,55] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][55] = 	xjac[1,55] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][55] = 	xjac[1,55] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][55] = 	xjac[1,55] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][55] = 	xjac[1,55] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][55] = 	xjac[1,55] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][55] = 	xjac[1,55] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][55] = 	xjac[1,55] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,55] = 	xjac[1,55] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,55] = 	xjac[1,55] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][55] = 	xjac[1,55] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,55] = 	xjac[1,55] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][55] = 	xjac[1,55] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][55] = 	xjac[1,55] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][55] = 	xjac[1,55] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][55] = 	xjac[1,55] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][55] = 	xjac[1,55] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*pow(x[56],1.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][55] = 	xjac[1,55] + (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[1][55] = 	xjac[1,55] + (3.5e-12*pow(x[10],1.0)*pow(x[56],1.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][55] = 	xjac[1,55] + (8.5e-13*pow(300/T,0)*exp(2450/T)*pow(x[56],2.0))#2NO3 = 2NO2 + O2
	xjac[1][55] = 	xjac[1,55] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][55] = 	xjac[1,55] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3
	xjac[1,56] = 	xjac[1,56] - (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][56] = 	xjac[1,56] + (3.2e-11*pow(300/T,0)*exp(70/T)*pow(x[0],1.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[1][56] = 	xjac[1,56] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[1,56] = 	xjac[1,56] - (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][56] = 	xjac[1,56] + (3.6e-18*pow(300/T,0)*exp(-220/T)*pow(x[5],1.0)*pow(x[1],1.0))#O2_1d + O2 = 2O2
	xjac[1][56] = 	xjac[1,56] + (4.8e-18*pow(x[5],1.0)*pow(x[6],1.0))#O2_1d + H2O = O2 + H2O
	xjac[1][56] = 	xjac[1,56] + (1e-20*pow(x[5],1.0)*pow(x[3],1.0))#O2_1d + N2 = O2 + N2
	xjac[1][56] = 	xjac[1,56] + (2e-21*pow(x[5],1.0)*pow(x[4],1.0))#O2_1d + CO2 = O2 + CO2
	xjac[1][56] = 	xjac[1,56] + (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[1,56] = 	xjac[1,56] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[1,56] = 	xjac[1,56] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1][56] = 	xjac[1,56] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[1,56] = 	xjac[1,56] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[1,56] = 	xjac[1,56] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[1,56] = 	xjac[1,56] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[1,56] = 	xjac[1,56] - (((5.7e-32*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((5.7e-32*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#H + O2 + N2 = HO2 + N2
	xjac[1][56] = 	xjac[1,56] + (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[1][56] = 	xjac[1,56] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O2
	xjac[1][56] = 	xjac[1,56] + (1.2e-10*pow(x[0],1.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[1][56] = 	xjac[1,56] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[1][56] = 	xjac[1,56] + (1.4e-10*pow(300/T,0)*exp(-470/T)*pow(x[9],1.0)*pow(x[8],1.0))#H + O3 = OH + O2
	xjac[1][56] = 	xjac[1,56] + (1.7e-12*pow(300/T,0)*exp(-940/T)*pow(x[11],1.0)*pow(x[8],1.0))#OH + O3 = HO2 + O2
	xjac[1][56] = 	xjac[1,56] + (1e-14*pow(300/T,0)*exp(-490/T)*pow(x[10],1.0)*pow(x[8],1.0))#HO2 + O3 = OH + 2O2
	xjac[1][56] = 	xjac[1,56] + (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[1][56] = 	xjac[1,56] + (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[1][56] = 	xjac[1,56] + (7.29e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2 + O2
	xjac[1][56] = 	xjac[1,56] + (4.8e-11*pow(300/T,0)*exp(250/T)*pow(x[11],1.0)*pow(x[10],1.0))#OH + HO2 = H2O + O2
	xjac[1][56] = 	xjac[1,56] + (2.3e-13*pow(300/T,0)*exp(600/T)*pow(x[10],2.0))#2HO2 = H2O2 + O2
	xjac[1][56] = 	xjac[1,56] + (1.7e-33*pow(300/T,0)*exp(1000/T)*N*pow(x[10],2.0)*pow(x[12],1.0))#2HO2 + M = H2O2 + O2 + M
	xjac[1][56] = 	xjac[1,56] + (2.3e-11*pow(300/T,0)*exp(-200/T)*pow(x[15],1.0)*pow(x[8],1.0))#Cl + O3 = ClO + O2
	xjac[1][56] = 	xjac[1,56] + (1.8e-11*pow(300/T,0)*exp(170/T)*pow(x[15],1.0)*pow(x[10],1.0))#Cl + HO2 = HCl + O2
	xjac[1][56] = 	xjac[1,56] + (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[1][56] = 	xjac[1,56] + (6e-13*pow(300/T,0)*exp(230/T)*pow(x[16],1.0)*pow(x[11],1.0))#ClO + OH = HCl + O2
	xjac[1][56] = 	xjac[1,56] + (2.7e-12*pow(300/T,0)*exp(220/T)*pow(x[16],1.0)*pow(x[10],1.0))#ClO + HO2 = HOCl + O2
	xjac[1][56] = 	xjac[1,56] + (1e-12*pow(300/T,0)*exp(-1590/T)*pow(x[16],2.0))#2ClO = Cl2 + O2
	xjac[1][56] = 	xjac[1,56] + (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[1][56] = 	xjac[1,56] + (5e-12*pow(x[31],2.0))#2ClCO3 = 2Cl + 2CO2 + O2
	xjac[1,56] = 	xjac[1,56] - (((2e-31*pow(300/T,1.6)*exp(0/T))*N*(7.5e-11*pow(300/T,0)*exp(0/T)))/(((2e-31*pow(300/T,1.6)*exp(0/T))*N) +(7.5e-11*pow(300/T,0)*exp(0/T)))*pow(x[9],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#H + O2 + CO2 = HO2 + CO2
	xjac[1,56] = 	xjac[1,56] - (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[1,56] = 	xjac[1,56] - (1.3e-12*pow(300/T,0)*exp(-330/T)*pow(x[33],1.0)*pow(x[1],1.0))#HSO3 + O2 = = HO2 + SO3
	xjac[1,56] = 	xjac[1,56] - (2e-15*pow(x[22],1.0)*pow(x[1],1.0))#ClS + O2 = SO + ClO
	xjac[1][56] = 	xjac[1,56] + (1.2e-11*pow(x[32],1.0)*pow(x[8],1.0))#S + O3 = SO + O2
	xjac[1][56] = 	xjac[1,56] + (4.5e-12*pow(300/T,0)*exp(-1170/T)*pow(x[29],1.0)*pow(x[8],1.0))#SO + O3 = SO2 + O2
	xjac[1][56] = 	xjac[1,56] + (3e-12*pow(300/T,0)*exp(-7000/T)*pow(x[25],1.0)*pow(x[8],1.0))#SO2 + O3 = SO3 + O2
	xjac[1][56] = 	xjac[1,56] + (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[1][56] = 	xjac[1,56] + (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[1][56] = 	xjac[1,56] + (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[1][56] = 	xjac[1,56] + (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[1][56] = 	xjac[1,56] + (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[1,56] = 	xjac[1,56] - (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[1,56] = 	xjac[1,56] - (5.25e-12*pow(300/T,0)*exp(-1510/T)*pow(x[50],1.0)*pow(x[1],1.0))#HNO + O2 = NO + HO2
	xjac[1][56] = 	xjac[1,56] + (2e-16*pow(x[48],1.0)*pow(x[8],1.0))#N + O3 = NO + O2
	xjac[1,56] = 	xjac[1,56] - (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[1][56] = 	xjac[1,56] + (3e-12*pow(300/T,0)*exp(-1500/T)*pow(x[8],1.0)*pow(x[49],1.0))#O3 + NO = NO2 + O2
	xjac[1][56] = 	xjac[1,56] + (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[1][56] = 	xjac[1,56] + (1.2e-13*pow(300/T,0)*exp(-2450/T)*pow(x[8],1.0)*pow(x[53],1.0))#O3 + NO2 = NO3 + O2
	xjac[1][56] = 	xjac[1,56] + (5e-16*pow(x[10],1.0)*pow(x[53],1.0))#HO2 + NO2 = HNO2 + O2
	xjac[1][56] = 	xjac[1,56] + (4.5e-14*pow(300/T,0)*exp(-1260/T)*1.0*pow(x[56],0.0)*pow(x[53],1.0))#NO3 + NO2 = NO + NO2 + O2
	xjac[1][56] = 	xjac[1,56] + (1e-11*pow(x[2],1.0)*1.0*pow(x[56],0.0))#O + NO3 = O2 + NO2
	xjac[1][56] = 	xjac[1,56] + (3.5e-12*pow(x[10],1.0)*1.0*pow(x[56],0.0))#HO2 + NO3 = HNO3 + O2
	xjac[1][56] = 	xjac[1,56] + (8.5e-13*pow(300/T,0)*exp(2450/T)*2.0*pow(x[56],1.0))#2NO3 = 2NO2 + O2
	xjac[1][56] = 	xjac[1,56] + (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
	xjac[1][56] = 	xjac[1,56] + (5e-19*pow(x[8],1.0)*pow(x[54],1.0))#O3 + HNO2 = O2 + HNO3

	xjac[2][0] = 	xjac[2,0] + (3.2e-11*pow(300/T,0)*exp(70/T)*1.0*pow(x[0],0.0)*pow(x[1],1.0))#O_1D + O2 = O + O2
	xjac[2][0] = 	xjac[2,0] + (1.8e-11*pow(300/T,0)*exp(110/T)*1.0*pow(x[0],0.0)*pow(x[3],1.0))#O_1D + N2 = O + N2
	xjac[2][0] = 	xjac[2,0] + (7.4e-11*pow(300/T,0)*exp(120/T)*1.0*pow(x[0],0.0)*pow(x[4],1.0))#O_1D + CO2 = O + CO2
	xjac[2,0] = 	xjac[2,0] - (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[2][0] = 	xjac[2,0] + (2e-16*pow(x[5],1.0)*pow(x[2],1.0))#O2_1d + O = O2 + O
	xjac[2,0] = 	xjac[2,0] - (2.898e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2 + CO2
	xjac[2,0] = 	xjac[2,0] - (8.712e-23*pow(1/T,2)*exp(0/T)*N*pow(x[2],2.0)*pow(x[4],1.0))#2O + CO2 = O2_1d + CO2
	xjac[2,0] = 	xjac[2,0] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[2][0] = 	xjac[2,0] + (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],2.0)*pow(x[1],1.0))#2O + O2 = O3 + O
	xjac[2,0] = 	xjac[2,0] - (((5.9e-34*pow(300/T,2.4)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.9e-34*pow(300/T,2.4)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],2.0))#O + 2O2 = O3 + O2
	xjac[2,0] = 	xjac[2,0] - (((5.95e-34*pow(300/T,2.3)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((5.95e-34*pow(300/T,2.3)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[3],1.0))#O + O2 + N2 = O3 + N2
	xjac[2,0] = 	xjac[2,0] - (((6.7e-34*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((6.7e-34*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[7],1.0))#O + O2 + CO = O3 + CO
	xjac[2,0] = 	xjac[2,0] - (((1.4e-33*pow(300/T,2.5)*exp(0/T))*N*(2.8e-12*pow(300/T,0)*exp(0/T)))/(((1.4e-33*pow(300/T,2.5)*exp(0/T))*N) +(2.8e-12*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[1],1.0)*pow(x[4],1.0))#O + O2 + CO2 = O3 + CO2
	xjac[2,0] = 	xjac[2,0] - (8e-12*pow(300/T,0)*exp(-2060/T)*pow(x[2],1.0)*pow(x[8],1.0))#O + O3 = 2O2
	xjac[2][0] = 	xjac[2,0] + (1.2e-10*1.0*pow(x[0],0.0)*pow(x[8],1.0))#O_1D + O3 = 2O + O2
	xjac[2][0] = 	xjac[2,0] + (5.2e-11*pow(300/T,0)*exp(-2840/T)*pow(x[5],1.0)*pow(x[8],1.0))#O2_1d + O3 = 2O2 + O
	xjac[2,0] = 	xjac[2,0] - (1.3e-29*pow(1/T,1)*exp(0/T)*N*pow(x[2],1.0)*pow(x[9],1.0)*pow(x[12],1.0))#O + H + M = OH + M
	xjac[2,0] = 	xjac[2,0] - (8.5e-20*pow(300/T,-2.7)*exp(-3160/T)*pow(x[2],1.0)*pow(x[13],1.0))#O + H2 = OH + H
	xjac[2,0] = 	xjac[2,0] - (2.2e-11*pow(300/T,0)*exp(120/T)*pow(x[2],1.0)*pow(x[11],1.0))#O + OH = O2 + H
	xjac[2][0] = 	xjac[2,0] + (4.2e-12*pow(300/T,0)*exp(-240/T)*pow(x[11],2.0))#2OH = H2O + O
	xjac[2,0] = 	xjac[2,0] - (3e-11*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2
	xjac[2,0] = 	xjac[2,0] - (6e-13*pow(300/T,0)*exp(200/T)*pow(x[2],1.0)*pow(x[10],1.0))#O + HO2 = OH + O2_1d
	xjac[2][0] = 	xjac[2,0] + (1.62e-12*pow(x[9],1.0)*pow(x[10],1.0))#H + HO2 = H2O + O
	xjac[2,0] = 	xjac[2,0] - (1.4e-12*pow(300/T,0)*exp(-2000/T)*pow(x[2],1.0)*pow(x[14],1.0))#O + H2O2 = OH + HO2
	xjac[2,0] = 	xjac[2,0] - (5e-32*N*pow(x[15],1.0)*pow(x[2],1.0)*pow(x[12],1.0))#Cl + O + M = ClO + M
	xjac[2][0] = 	xjac[2,0] + (8.33e-12*pow(300/T,0)*exp(-2790/T)*pow(x[15],1.0)*pow(x[11],1.0))#Cl + OH = HCl + O
	xjac[2,0] = 	xjac[2,0] - (7.4e-12*pow(300/T,0)*exp(-1650/T)*pow(x[2],1.0)*pow(x[19],1.0))#O + Cl2 = ClO + Cl
	xjac[2][0] = 	xjac[2,0] + (5.25e-11*1.0*pow(x[0],0.0)*pow(x[19],1.0))#O_1D + Cl2 = Cl2 + O
	xjac[2,0] = 	xjac[2,0] - (3e-11*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2
	xjac[2,0] = 	xjac[2,0] - (6e-13*pow(300/T,0)*exp(70/T)*pow(x[16],1.0)*pow(x[2],1.0))#ClO + O = Cl + O2_1d
	xjac[2,0] = 	xjac[2,0] - (1e-11*pow(300/T,0)*exp(-3300/T)*pow(x[2],1.0)*pow(x[17],1.0))#O + HCl = OH + Cl
	xjac[2][0] = 	xjac[2,0] + (1.35e-11*1.0*pow(x[0],0.0)*pow(x[17],1.0))#O_1D + HCl = O + HCl
	xjac[2,0] = 	xjac[2,0] - (1.7e-13*pow(x[2],1.0)*pow(x[18],1.0))#O + HOCl = OH + ClO
	xjac[2,0] = 	xjac[2,0] - (3e-11*pow(x[2],1.0)*pow(x[20],1.0))#O + ClCO = Cl + CO2
	xjac[2,0] = 	xjac[2,0] - (3e-12*pow(x[2],1.0)*pow(x[20],1.0))#O + ClCO = CO + ClO
	xjac[2,0] = 	xjac[2,0] - (1e-11*pow(x[2],1.0)*pow(x[31],1.0))#O + ClCO3 = Cl + O2 + CO2
	xjac[2,0] = 	xjac[2,0] - (((1.7e-33*pow(300/T,0)*exp(-1510/T))*N*(2.66e-14*pow(300/T,0)*exp(-1459/T)))/(((1.7e-33*pow(300/T,0)*exp(-1510/T))*N) +(2.66e-14*pow(300/T,0)*exp(-1459/T)))*pow(x[2],1.0)*pow(x[7],1.0)*pow(x[12],1.0))#O + CO + M = CO2 + M
	xjac[2,0] = 	xjac[2,0] - (6.5e-33*pow(300/T,0)*exp(-2180/T)*pow(x[2],1.0)*pow(x[7],2.0))#O + 2CO = CO2 + CO
	xjac[2,0] = 	xjac[2,0] - (3.4e-33*pow(300/T,0)*exp(-2180/T)*pow(x[2],2.0)*pow(x[7],1.0))#2O + CO = CO2 + O
	xjac[2][0] = 	xjac[2,0] + (3.4e-33*pow(300/T,0)*exp(-2180/T)*pow(x[2],2.0)*pow(x[7],1.0))#2O + CO = CO2 + O
	xjac[2,0] = 	xjac[2,0] - (1.5e-34*pow(300/T,0)*exp(900/T)*N*pow(x[32],1.0)*pow(x[2],1.0)*pow(x[12],1.0))#S + O + M = SO + M
	xjac[2][0] = 	xjac[2,0] + (2.3e-12*pow(x[32],1.0)*pow(x[1],1.0))#S + O2 = SO + O
	xjac[2,0] = 	xjac[2,0] - (1.6e-13*pow(300/T,0)*exp(-2280/T)*pow(x[29],1.0)*pow(x[2],1.0))#SO + O = SO2 + O
	xjac[2][0] = 	xjac[2,0] + (1.6e-13*pow(300/T,0)*exp(-2280/T)*pow(x[29],1.0)*pow(x[2],1.0))#SO + O = SO2 + O
	xjac[2,0] = 	xjac[2,0] - (2.2e-11*pow(300/T,0)*exp(-84/T)*pow(x[2],1.0)*pow(x[24],1.0))#O + S2 = SO + S
	xjac[2,0] = 	xjac[2,0] - (8e-11*pow(x[2],1.0)*pow(x[36],1.0))#O + S3 = SO + S2
	xjac[2,0] = 	xjac[2,0] - (8e-11*pow(x[2],1.0)*pow(x[38],1.0))#O + S4 = SO + S3
	xjac[2,0] = 	xjac[2,0] - (8e-11*pow(300/T,0)*exp(-200/T)*pow(x[2],1.0)*pow(x[39],1.0))#O + S5 = S4 + SO
	xjac[2,0] = 	xjac[2,0] - (8e-11*pow(300/T,0)*exp(-300/T)*pow(x[2],1.0)*pow(x[40],1.0))#O + S6 = S5 + SO
	xjac[2,0] = 	xjac[2,0] - (8e-11*pow(300/T,0)*exp(-200/T)*pow(x[2],1.0)*pow(x[41],1.0))#O + S7 = S6 + SO
	xjac[2,0] = 	xjac[2,0] - (8e-11*pow(300/T,0)*exp(-400/T)*pow(x[2],1.0)*pow(x[42],1.0))#O + S8 = S7 + SO
	xjac[2,0] = 	xjac[2,0] - (6.6e-13*pow(300/T,0)*exp(-2760/T)*pow(x[2],1.0)*pow(x[29],1.0))#O + SO = S + O2
	xjac[2,0] = 	xjac[2,0] - (((5.1e-31*pow(300/T,0)*exp(0/T))*N*(5.3e-11*pow(300/T,0)*exp(0/T)))/(((5.1e-31*pow(300/T,0)*exp(0/T))*N) +(5.3e-11*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[12],1.0))#O + SO + M = SO2 + M
	xjac[2,0] = 	xjac[2,0] - (1e-11*pow(x[43],1.0)*pow(x[2],1.0)*pow(x[44],1.0)*pow(x[29],1.0))#ClC(O)OO + SO = Cl + SO2 + CO2
	xjac[2,0] = 	xjac[2,0] - (3e-14*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = S2O + O2
	xjac[2,0] = 	xjac[2,0] - (3e-15*pow(x[2],1.0)*pow(x[29],1.0)*pow(x[45],2.0))#O + (SO)2 = SO + SO2
	xjac[2,0] = 	xjac[2,0] - (8e-12*pow(300/T,0)*exp(-9800/T)*pow(x[2],1.0)*pow(x[25],1.0))#O + SO2 = SO + O2
	xjac[2,0] = 	xjac[2,0] - (7e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO2 + O
	xjac[2][0] = 	xjac[2,0] + (7e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO2 + O
	xjac[2,0] = 	xjac[2,0] - (1.3e-10*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[25],1.0))#O(1D) + SO2 = SO + O2
	xjac[2,0] = 	xjac[2,0] - (1.32e-31*pow(300/T,0)*exp(-1000/T)*N*pow(x[2],1.0)*pow(x[25],1.0)*pow(x[12],1.0))#O + SO2 + M = SO3 + M
	xjac[2,0] = 	xjac[2,0] - (1e-15*pow(x[43],1.0)*pow(x[2],1.0)*pow(x[44],1.0)*pow(x[25],1.0))#ClC(O)OO + SO2 = Cl + SO3 + CO2
	xjac[2,0] = 	xjac[2,0] - (2.32e-16*pow(300/T,0)*exp(-487/T)*pow(x[2],1.0)*pow(x[30],1.0))#O + SO3 = SO2 + O2
	xjac[2,0] = 	xjac[2,0] - (1.7e-12*pow(x[2],1.0)*pow(x[37],1.0))#O + S2O = 2SO
	xjac[2,0] = 	xjac[2,0] - (1.2e-10*pow(x[2],1.0)*pow(x[22],1.0))#O + ClS = SO + Cl
	xjac[2,0] = 	xjac[2,0] - (1e-13*pow(x[2],1.0)*pow(x[23],1.0))#O + ClS2 = SO + ClS
	xjac[2,0] = 	xjac[2,0] - (5e-11*pow(300/T,0)*exp(-600/T)*pow(x[2],1.0)*pow(x[28],1.0))#O + OSCl = SO2 + Cl
	xjac[2,0] = 	xjac[2,0] - (2e-11*pow(300/T,0)*exp(-600/T)*pow(x[2],1.0)*pow(x[28],1.0))#O + OSCl = SO + ClO
	xjac[2,0] = 	xjac[2,0] - (1e-11*pow(x[2],1.0)*pow(x[26],1.0))#O + ClSO2 = SO2 + ClO
	xjac[2,0] = 	xjac[2,0] - (1.6e-11*pow(300/T,0)*exp(-2150/T)*pow(x[2],1.0)*pow(x[21],1.0))#O + OCS = SO + CO
	xjac[2][0] = 	xjac[2,0] + (1.5e-11*pow(300/T,0)*exp(-3600/T)*pow(x[48],1.0)*pow(x[1],1.0))#N + O2 = NO + O
	xjac[2][0] = 	xjac[2,0] + (9e-17*pow(x[1],1.0)*pow(x[51],1.0)*pow(x[48],1.0))#O2(1d) + N = NO + O
	xjac[2,0] = 	xjac[2,0] - (3.5e-37*pow(300/T,0.6)*exp(0/T)*N*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[3],1.0)*pow(x[12],1.0))#O(1D) + N2 + M = N2O + M
	xjac[2,0] = 	xjac[2,0] - (((9e-31*pow(300/T,1.5)*exp(0/T))*N*(3e-11*pow(300/T,0)*exp(0/T)))/(((9e-31*pow(300/T,1.5)*exp(0/T))*N) +(3e-11*pow(300/T,0)*exp(0/T)))*pow(x[2],1.0)*pow(x[49],1.0)*pow(x[12],1.0))#O + NO + M = NO2 + M
	xjac[2][0] = 	xjac[2,0] + (2.1e-11*pow(300/T,0)*exp(100/T)*pow(x[48],1.0)*pow(x[49],1.0))#N + NO = N2 + O
	xjac[2,0] = 	xjac[2,0] - (5.6e-12*pow(300/T,0)*exp(180/T)*pow(x[2],1.0)*pow(x[53],1.0))#O + NO2 = NO + O2
	xjac[2,0] = 	xjac[2,0] - (((2.5e-31*pow(300/T,1.8)*exp(0/T))*N*(2.2e-11*pow(300/T,0.7)*exp(0/T)))/(((2.5e-31*pow(300/T,1.8)*exp(0/T))*N) +(2.2e-11*pow(300/T,0.7)*exp(0/T)))*pow(x[2],1.0)*pow(x[53],1.0)*pow(x[12],1.0))#O + NO2 + M = HNO3 + M
	xjac[2][0] = 	xjac[2,0] + (5.8e-12*pow(300/T,0)*exp(220/T)*pow(x[48],1.0)*pow(x[53],1.0))#N + NO2 = N2O + O
	xjac[2,0] = 	xjac[2,0] - (1e-11*pow(x[2],1.0)*pow(x[56],1.0))#O + NO3 = O2 + NO2
	xjac[2,0] = 	xjac[2,0] - (6.7e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = 2NO
	xjac[2,0] = 	xjac[2,0] - (4.9e-11*pow(x[2],1.0)*pow(x[46],1.0)*pow(x[52],1.0))#O(1D) + N2O = N2 + O2
