from pylab import *
from scipy.interpolate import interp1d
import pandas as pd

#Read stellar irradiance, and wavelength
wav = []
irr = []
with open('../data/stellar/sun.ir', 'r') as file:
    # Skip the first line
     next(file)
     for line in file:
        # Split the line into two columns of numbers
        column1, column2 = map(float, line.strip().split())
        
        # Do something with the numbers, for example, print them
        wav.append(double(column1)*1e-10)
        irr.append(double(column2)*1e10)

#Read O2 photolysis cross section, and branching ratio
O2_xcross = "../data/VULCAN/O2/O2_cross.csv"

O2 = pd.read_csv(
    O2_xcross, skiprows=2, usecols=[0, 1, 2, 3], names=["wav", "abs", "diss", "ion"]
)

w = O2["wav"][:]*1e-9
diss_x = O2["diss"][:] * 1e-4
f_diss = interp1d(w, diss_x, fill_value=(0, 0), bounds_error=False)

O2_branch = "../data/VULCAN/O2/O2_branch.csv"

O2_br  = pd.read_csv(
    O2_branch, skiprows=2, usecols=[0, 1, 2], names=["wav", "p1", "p2"]
)

w_br = O2_br["wav"][:]*1e-9
br1 = O2_br["p1"][:]
br2 = O2_br["p2"][:]
f1 = interp1d(w_br, br1,fill_value=(0, 0), bounds_error=False)
f2 = interp1d(w_br, br2,fill_value=(0, 0), bounds_error=False)
#f1 = interp1d(w_br, br1,fill_value=(br1[0], array(br1)[-1]), bounds_error=False)
#f2 = interp1d(w_br, br2,fill_value=(br2[0], array(br2)[-1]), bounds_error=False)
O2_x1 = f1(wav)*f_diss(wav)
O2_x2 = f2(wav)*f_diss(wav)

'''
O2_x1 = f1(wav)*f_diss(wav)
O2_x2 = f2(wav)*f_diss(wav)
plot(w*1e9, diss_x, 'b-')
plot(array(wav)*1e9, O2_x1, 'r-')
plot(array(wav)*1e9, O2_x2, 'r--')
yscale('log')
xlim([0, 500])
show()
'''
#Read O3 photolysis cross section, and branching ratio
O3_xcross = "../data/VULCAN/O3/O3_cross.csv"

O3 = pd.read_csv(
    O3_xcross, skiprows=2, usecols=[0, 1, 2, 3], names=["wav", "abs", "diss", "ion"]
)

w = O3["wav"][:]*1e-9
diss_x = O3["diss"][:] * 1e-4

f_diss = interp1d(w, diss_x,fill_value=(0, 0), bounds_error=False)

O3_branch = "../data/VULCAN/O3/O3_branch.csv"

O3_br  = pd.read_csv(
    O3_branch, skiprows=2, usecols=[0, 1, 2], names=["wav", "p1", "p2"]
)

w_br = O3_br["wav"][:]*1e-9
br3 = O3_br["p1"][:]
br4 = O3_br["p2"][:]
#f3 = interp1d(w_br, br3,fill_value=(br3[0], array(br3)[-1]), bounds_error=False)
#f4 = interp1d(w_br, br4,fill_value=(br4[0], array(br4)[-1]), bounds_error=False)
f3 = interp1d(w_br, br3,fill_value=(0, 0), bounds_error=False)
f4 = interp1d(w_br, br4,fill_value=(0, 0), bounds_error=False)

O3_x1 = f3(wav)*f_diss(wav)
O3_x2 = f4(wav)*f_diss(wav)


plot(w*1e9, diss_x, 'b-')
plot(array(wav)*1e9, O3_x1, 'r-')
plot(array(wav)*1e9, O3_x2, 'r--')
yscale('log')
xlim([0, 500])
show()



#Compute photolysis cross section
h = 6.626E-34
c = 3E8

aF = array(irr)*array(wav)/(h*c)
R1 = 0
R2 = 0
R3 = 0
R4 = 0

for i in range(0, len(aF) -1):
    R1 = R1 + 0.5*((aF[i]*O2_x1[i]) + (aF[i+1]*O2_x1[i+1]))*(wav[i+1] - wav[i])
    R2 = R2 + 0.5*((aF[i]*O2_x2[i]) + (aF[i+1]*O2_x2[i+1]))*(wav[i+1] - wav[i])
    R3 = R3 + 0.5*((aF[i]*O3_x1[i]) + (aF[i+1]*O3_x1[i+1]))*(wav[i+1] - wav[i])
    R4 = R4 + 0.5*((aF[i]*O3_x2[i]) + (aF[i+1]*O3_x2[i+1]))*(wav[i+1] - wav[i])


print("Photolysis rates")
print(R1)
print(R2)
print(R3)
print(R4)


