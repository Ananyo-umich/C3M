#The program initiates the reaction network solver for 1-d along z direction
#Author: Ananyo Bhattacharya
#Email: ananyo@umich.edu
#Affiliation: University of Michigan, Ann Arbor
from pylab import *
from glob import glob
import multiprocessing as mp
import concurrent.futures
import re
from equation import rate_equation

#Reaction structure
class reaction:
  def __init__(self, idNum):
    self.Data = []
    self.Reactants = []
    self.Products = []
    self.Rstoic = []
    self.Pstoic = []
    self.molecules = []
    self.idNum = idNum
    self.RateFlag = []
    self.ReactantID = []
    self.ProductID = []
    self.Equation = []
    
    
#Functions
def stoich(word, pattern):
  coeff = re.findall(pattern, word)
  if(coeff == []):
    ans = ['1']
  else:
    ans = coeff
  return ans

#Reading input file
inpFile = "scheme.inp"
reactionFileSearch = "reaction_database"
rateFileSearch = "rate_database"
reactionFile = []
rateFile = []
Temp = []
with open(inpFile) as f:
  InpLines = f.read().splitlines()
  f.close()
  
for i in range(len(InpLines)):
  if(re.findall(reactionFileSearch, InpLines[i]) != []):
    reactionFile = re.findall("[A-Za-z0-9.\/_]+",InpLines[i].split("=")[1])
  if(re.findall(rateFileSearch, InpLines[i]) != []):
    rateFile = re.findall("[A-Za-z0-9.\/_]+",InpLines[i].split("=")[1])
  if(re.findall("temp", InpLines[i]) != []):  
    Temp = float(re.findall("[0-9.]+",InpLines[i].split("=")[1])[0])

reaction_array = []
with open(reactionFile[0]) as f:
    lines = f.read().splitlines()
    f.close()

num = len(lines)
#Regular expression to identify reactants and products
chem_pattern = '[A-Za-z0-9()]+' #Chemistry pattern to read reactants and products
reactant_pattern = '[A-Za-z0-9()]+' 
product_pattern = '[A-Za-z0-9()]+'
stoich_pattern = re.compile("\d")
#Varibales in the chemical reaction network
chemicals = []
for i in range(num):
  line = (lines[i]).split(",")
  reactionID = line[0]
  rxnClass = reaction(i)
  rxnClass.Data = line
  rxnClass.RateFlag = line[2]
  #Regex expression to search for reactants and products
  rxnClass.Equation = line[1]
  rxnClass.Reactants = re.findall(reactant_pattern, (line[1].split("="))[0]) 
  rxnClass.Products = re.findall(product_pattern, (line[1].split("="))[1]) 
  rxnClass.ReactantID = zeros(len(rxnClass.Reactants))
  rxnClass.ProductID = zeros(len(rxnClass.Products))
  rxnClass.Rstoic = zeros(len(rxnClass.Reactants))
  rxnClass.Pstoic = zeros(len(rxnClass.Products))
  for rts in range(len(rxnClass.Reactants)):
    var = stoich_pattern.match(rxnClass.Reactants[rts]) 
    if(var == None):
      rxnClass.Rstoic[rts] = 1
    if(var != None):
      rxnClass.Rstoic[rts] = float(var.group(0))
      rxnClass.Reactants[rts] = rxnClass.Reactants[rts].replace(var.group(0), "",1)
  for pts in range(len(rxnClass.Products)):
    var = stoich_pattern.match(rxnClass.Products[pts])
    if(var == None):
      rxnClass.Pstoic[pts] = 1
    if(var != None):
      rxnClass.Pstoic[pts] = float(var.group(0))
      rxnClass.Products[pts] = rxnClass.Products[pts].replace(var.group(0), "",1)
  rxnClass.molecules = concatenate((rxnClass.Reactants, rxnClass.Products))
  reaction_array.append(rxnClass)
  #Regex to find stoichiometric coefficients
  for molecule in rxnClass.molecules:
        if molecule not in chemicals:
            chemicals.append(molecule)
  stoi_reac = re.findall(stoich_pattern, line[1]) 
  


variables = array(chemicals)  
print(variables)
for varIndex in range(len(variables)):
  for rxn in reaction_array:
    for rts in range(len(rxn.Reactants)):
      if(variables[varIndex] == rxn.Reactants[rts]):
        rxn.ReactantID[rts] = varIndex 
    for pts in range(len(rxn.Products)):
      if(variables[varIndex] == rxn.Products[pts]): 
        rxn.ProductID[pts] = varIndex 
 

#accessing a netcdf file for atmosphere
#InitData = nc.Dataset(AtmosFile, 'r' ,format = 'NETCDF4')
#Pres = maxEUV['lat'][:]
#Temp = maxEUV['lon'][:]
#for i in variables:
#vardata = InitData[str(variable)
#Kzz
#


#Writing python code for chemical reaction network solver
#Searching variables and reactions to create ODEs
with open('chemnet.py', 'w') as code:
#Import packages
  code.write('from pylab import * \nimport netCDF4 as nc \nfrom glob import glob\n\n' )
  code.write('\nT = '+ str(Temp) +' \n\n')
  #for varIndex in range(len(variables)):
  #  code.write('i' +variables[varIndex].replace(["(", ")"],"_")+ '= '+ str(varIndex) + '\n')
  code.write('def Jacobian(x):\n')
  code.write("\txjac = zeros(["+str(len(variables))+","+str(len(variables))+"])\n")
  for varIndex in range(len(variables)):
    for varIndex2 in range(len(variables)):
      for rn in range(len(reaction_array)):
        if variables[varIndex] in reaction_array[rn].Reactants:
          rate = rate_equation.ChemicalReactionRateCoefficient(reaction_array[rn].Data) 
          jac = rate          
          for rts in range(len(reaction_array[rn].Reactants)): 
            if(variables[varIndex2] == reaction_array[rn].Reactants[rts]):
              jac = jac + "*" +str(reaction_array[rn].Rstoic[rts])+ "*pow(x[" + str(int(reaction_array[rn].ReactantID[rts])) + "]," + str(reaction_array[rn].Rstoic[rts] - 1) + ")"
            if(variables[varIndex2] != reaction_array[rn].Reactants[rts]):
              jac = jac + "*pow(x[" + str(int(reaction_array[rn].ReactantID[rts])) + "]," + str(reaction_array[rn].Rstoic[rts]) + ")"          
          loss =  "(" + jac + ")"        
          code.write("\txjac[" +str(varIndex)+","+str(varIndex2)+ "] = " + "\txjac[" +str(varIndex)+","+str(varIndex2)+ "] - " + loss  + '#' + str(reaction_array[rn].Equation) + '\n' )
        
        if variables[varIndex] in reaction_array[rn].Products:
          rate = rate_equation.ChemicalReactionRateCoefficient(reaction_array[rn].Data) 
          jac = rate 
          for pts in range(len(reaction_array[rn].Products)):
            if(variables[varIndex2] == reaction_array[rn].Products[pts]): 
              jac = jac + "*" +str(reaction_array[rn].Pstoic[pts])+ "*pow(x[" + str(int(reaction_array[rn].ProductID[pts])) + "]," + str(reaction_array[rn].Pstoic[pts] - 1) + ")"
            if(variables[varIndex2] != reaction_array[rn].Products[pts]): 
              jac = jac + "*pow(x[" + str(int(reaction_array[rn].ProductID[pts])) + "]," + str(reaction_array[rn].Pstoic[pts]) + ")"
          production =  "(" + jac + ")"
          code.write("\txjac[" +str(varIndex)+"]["+str(varIndex2)+ "] = " + "\txjac[" +str(varIndex)+","+str(varIndex2)+ "] - " + production  + '#' + str(reaction_array[rn].Equation) + '\n' )
    code.write('\n' )
  code.write('\treturn xjac\n\n' )
  code.write('def ReactionRate(x):\n')
  code.write("\tdx = zeros("+str(len(variables))+")\n")
  for varIndex in range(len(variables)):
    for rn in range(len(reaction_array)):
      if variables[varIndex] in reaction_array[rn].Reactants:
        rate = rate_equation.ChemicalReactionRateCoefficient(reaction_array[rn].Data) 
        for rts in range(len(reaction_array[rn].Reactants)): 
          rate = rate + "*pow(x[" + str(int(reaction_array[rn].ReactantID[rts])) + "]," + str(reaction_array[rn].Rstoic[rts]) + ")"            
        loss =  "(" + rate + ")"        
        code.write("\tdx[" +str(varIndex)+ "] = " + "\tdx[" +str(varIndex)+ "] - " + loss  + '#' + str(reaction_array[rn].Equation) + '\n' )
        
      if variables[varIndex] in reaction_array[rn].Products:
        rate = rate_equation.ChemicalReactionRateCoefficient(reaction_array[rn].Data) 
        for pts in range(len(reaction_array[rn].Products)): 
          rate =  rate + "*pow(x[" + str(int(reaction_array[rn].ProductID[pts])) + "]," + str(reaction_array[rn].Pstoic[pts]) + ")" #concentration of products
        production =  "(" + rate + ")"
        code.write("\tdx[" +str(varIndex)+ "] = " + "\tdx[" +str(varIndex)+ "] + " + production  + '#' + str(reaction_array[rn].Equation) + '\n' )
    code.write('\n' )
  code.write('\treturn dx\n\n' )
#Initial conditions
  code.write("t_init = 0\n")
  code.write("t_final = 100\n")
  code.write("n = 10000\n")
  code.write("x = abs(randn("+ str(len(variables))+",n))\nrate = zeros("+ str(len(variables)) +")\njacob = zeros(["+str(len(variables))+","+str(len(variables))+"])\n")
  for i in range(len(InpLines)):
    for varIndex in range(len(variables)):
      if(re.findall("#", InpLines[i]) == []):
        if(re.findall("[A-Za-z0-9-]+", InpLines[i].split("=")[0]) == [variables[varIndex]]):
          print(InpLines[i])
          val = re.findall("[0-9e.-]+", InpLines[i].split("=")[1])
          code.write("x[" + str(varIndex)+"] ="+ val[0] +"\n")
  
  
  code.write("dt = (t_final - t_init)/n\n")
  code.write("for i in range(n-1):\n")
  code.write("\ta_matrix = Jacobian(x[:,i])\n")
  code.write("\tprint(i)\n")
  code.write("\tprint(a_matrix)\n")
  code.write("\tx_pred = x[:,i]\n")
  code.write("\tg = x_pred-(ReactionRate(x[:,i])*dt)-x[:,i]\n")
  code.write("\tb_matrix = matmul(linalg.inv(Jacobian(x[:,i])),g)\n")
  code.write("\tx[:,i+1] = x[:,i] - (b_matrix)\n")
  code.close()
 


      
        

#Template code for Python
'''
from pylab import *
import netcdf
from glob import glob
import netCDF4 as nc

#Initiate concentration variable (maybe a separate file for monitoring variables)
x = zeros(len(variable))
dx = zeros(len(variable)) #Growth or loss rate
dt = 0.1 #Time step (V. Imp.)
nt = 1000 #Time scale

initial conditions needed
for time in range(nt):
#List of reactions
 writing the Jacobian of the reactions
 making the matrix (if there are n reactions the dimensions are n x n) 
 solving it



solver
'''