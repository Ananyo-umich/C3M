#The program initiates the reaction network solver for box model
#Author: Ananyo Bhattacharya
#Email: ananyo@umich.edu
#Affiliation: University of Michigan, Ann Arbor

from pylab import *
from glob import glob
import multiprocessing as mp
import concurrent.futures
import re
from equation import rate_equation
from config import *
from glob import glob
import os

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
    
class atmosphere:
  def _init_(self):
    self.Temp = []
    self.Alt = []
    self.Press = []
    self.NDensity = []    
    self.Eddy = []


#Functions
def stoich(word, pattern):
  coeff = re.findall(pattern, word)
  if(coeff == []):
    ans = ['1']
  else:
    ans = coeff
  return ans

#Reading input file
inpFile = str(args['i'])
reactionFileSearch = "reaction_database"
rateFileSearch = "rate_database"
atmosFileSearch = "planet_file"
reactionFile = []
rateFile = []
Temp = []
Press = []
Tstart = []
Tend = []
nTime = []
with open(inpFile) as f:
  InpLines = f.read().splitlines()
  f.close()
  
for i in range(len(InpLines)):
  if(re.findall(reactionFileSearch, InpLines[i]) != []):
    reactionFile = re.findall("[A-Za-z0-9.\/_]+",InpLines[i].split("=")[1])
  if(re.findall(atmosFileSearch, InpLines[i]) != []):
    atmosFile = re.findall("[A-Za-z0-9.\/_]+",InpLines[i].split("=")[1])
  if(re.findall(rateFileSearch, InpLines[i]) != []):
    rateFile = re.findall("[A-Za-z0-9.\/_]+",InpLines[i].split("=")[1])
  if(re.findall("temp", InpLines[i]) != []):  
    Temp = float(re.findall("[+0-9.Ee-]+",InpLines[i].split("=")[1])[0])
  if(re.findall("press", InpLines[i]) != []):  
    Press = float(re.findall("[+0-9.Ee-]+",InpLines[i].split("=")[1])[0])
  if(re.findall("tstart", InpLines[i]) != []):  
    Tstart = float(re.findall("[+0-9.Ee-]+",InpLines[i].split("=")[1])[0])
  if(re.findall("tend", InpLines[i]) != []):  
    Tend = float(re.findall("[+0-9.Ee-]+",InpLines[i].split("=")[1])[0])
  if(re.findall("nTime", InpLines[i]) != []):  
    nTime = int(re.findall("[+0-9.Ee-]+",InpLines[i].split("=")[1])[0])

reaction_array = []
with open(reactionFile[0]) as f:
    lines = f.read().splitlines()
    f.close()

num = len(lines)
#Regular expression to identify reactants and products
chem_pattern = '[A-Za-z0-9()]+' #Chemistry pattern to read reactants and products
reactant_pattern = '[A-Za-z0-9_]+' 
product_pattern = '[A-Za-z0-9_]+'
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
  
#Physical properties of the atmosphere
with open(atmosFile[0]) as f:
    lines = f.read().splitlines()
    f.close()
size = len(lines)
planet = atmosphere(0)
for i in range(size):
  data = lines[i].split(",")
  planet.Temp.append(data[0])
  planet.Alt.append(data[1])
  planet.Press.append(data[2])
  planet.NDensity.append(data[3])
  planet.Eddy.append(data[4])
  

#Reacrion ID for substrates
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
 

#Writing python code for chemical reaction network solver
#Searching variables and reactions to create ODEs

with open(args['ex'], 'w') as code:
#Import packages
  code.write('from pylab import * \nimport netCDF4 as nc \nfrom glob import glob\n\n' )
  code.write('AtmosData = genfromtxt"Venus.txt"\n')
  code.write('Temp = AtmosData[:,0]\n')
  code.write('Alt = AtmosData[:,1]\n')
  code.write('Press = AtmosData[:,2]\n')
  code.write('NDensity = AtmosData[:,3]\n')
  code.write('Eddy = AtmosData[:,4]\n')
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
          for rts in range(len(reaction_array[rn].Reactants)):
            if(variables[varIndex2] == reaction_array[rn].Reactants[rts]): 
              jac = jac + "*" +str(reaction_array[rn].Rstoic[rts])+ "*pow(x[" + str(int(reaction_array[rn].ReactantID[rts])) + "]," + str(reaction_array[rn].Rstoic[rts] - 1) + ")"
            if(variables[varIndex2] != reaction_array[rn].Reactants[rts]): 
              jac = jac + "*pow(x[" + str(int(reaction_array[rn].ReactantID[rts])) + "]," + str(reaction_array[rn].Rstoic[rts]) + ")"
          production =  "(" + jac + ")"
          code.write("\txjac[" +str(varIndex)+"]["+str(varIndex2)+ "] = " + "\txjac[" +str(varIndex)+","+str(varIndex2)+ "] + " + production  + '#' + str(reaction_array[rn].Equation) + '\n' )
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
        for rts in range(len(reaction_array[rn].Reactants)): 
          rate =  rate + "*pow(x[" + str(int(reaction_array[rn].ReactantID[rts])) + "]," + str(reaction_array[rn].Rstoic[rts]) + ")" #concentration of products
        production =  "(" + rate + ")"
        code.write("\tdx[" +str(varIndex)+ "] = " + "\tdx[" +str(varIndex)+ "] + " + production  + '#' + str(reaction_array[rn].Equation) + '\n' )
    code.write('\n' )
  code.write('\treturn dx\n\n' )
#Initial conditions
  code.write("t_init = "+str(Tstart)+"\n")
  code.write("t_final = "+str(Tend)+"\n")
  code.write("n = "+str(nTime)+"\n")
  code.write("x = 1e-6*N*abs(randn("+ str(len(variables))+",n))\nrate = zeros("+ str(len(variables)) +")\njacob = zeros(["+str(len(variables))+","+str(len(variables))+"])\n")
  '''
  for i in range(len(InpLines)):
    for varIndex in range(len(variables)):
      if(re.findall("#", InpLines[i]) == []):
        if(re.findall("[A-Za-z0-9-]+", InpLines[i].split("=")[0]) == [variables[varIndex]]):
          print(InpLines[i])
          val = re.findall("[0-9e.-]+", InpLines[i].split("=")[1])
          code.write("x[" + str(varIndex)+"] ="+ val[0] +"\n")
  '''
  
  
  code.write("dt = "+str((Tend-Tstart)/nTime)+"\n")
  code.write("for i in range(n-1):\n")
  code.write("\tj_matrix = Jacobian(x[:,i])\n")
  code.write("\tx_pred = x[:,i]\n")
  code.write("\trate = ReactionRate(x[:,i])\n")
  code.write("\ti_matrix = identity(len(x[:,i]))\n")
  code.write("\tx_prev = x_pred\n")
  code.write("\tx_next = x_pred\n")
  code.write("#Backward Euler Integration Scheme (Li and Chen, 2020)\n")
  code.write("\twhile True:\n")
  code.write("\t\tb_matrix = (rate - ((x_prev - x[:,i])/dt))\n")
  code.write("\t\tx_int = x_next\n")
  code.write("\t\tx_next = x_prev + matmul(linalg.inv((i_matrix/dt) - j_matrix),b_matrix)\n")
  code.write("\t\tx_prev = x_int\n")
  code.write("\t\tprint(str(abs(x_next - x_prev)/x_prev))\n")
  code.write("\t\tif (abs(x_next - x_prev)/x_prev <= 1e-10).all():\n")
  code.write("\t\t\tbreak\n")
  code.write("\tx[:,i+1] = x_next\n")
  code.write("\t\n")
  #code.write("\tg = x_pred-(ReactionRate(x[:,i])*dt)-x[:,i]\n")
  #code.write("\tb_matrix = matmul(linalg.inv(Jacobian(x[:,i])),g)\n")
  #code.write("\tx[:,i+1] = x[:,i] - (b_matrix)\n")
  code.write("\tprint('Chemical Reaction Network solved for '+str(i)+' Time Step')\n")
  code.write("\n#Writing the Output into NETCDF format\n")
  code.write("OutFile = '"+str(args['o'])+"' \n")
  code.write("data = nc.Dataset(OutFile, 'w' ,format = 'NETCDF4')\n")
  code.write("pres = data.createDimension('press',1)\n")
  code.write("temp = data.createDimension('temp',1)\n")
  code.write("P = data.createVariable('press', 'f8', ('press',))\n")
  code.write("T = data.createVariable('temp', 'f8', ('temp',))\n")
  code.write("time = data.createDimension('time',n)\n")
  code.write("t = data.createVariable('time', 'f8', ('time',))\n")
  for varIndex in range(len(variables)):
    code.write(str(variables[varIndex]) + "= data.createVariable('"+ str(variables[varIndex]) + "', 'f8', ('time',))\n")
    code.write(str(variables[varIndex]) + "[:] = x["+str(varIndex)+",:]\n")
  code.write("data.close()")
  code.close()
