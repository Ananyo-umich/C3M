#The program will generate and analyze graphical network of chemical reactions
from pylab import *
from glob import glob
import multiprocessing as mp
import concurrent.futures
import re
import networkx as nx

#Reaction structure
class reaction:
  def __init__(self, idNum):
    self.Reactants = []
    self.Products = []
    self.Rstoic = []
    self.Pstoic = []
    self.molecules = []
    self.idNum = idNum
    self.RateFlag = []
    self.acoeff = []
    self.bcoeff = []
    self.ccoeff = []
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
    reactionFile = re.findall("[a-z.\/]+",InpLines[i].split("=")[1])
  if(re.findall(rateFileSearch, InpLines[i]) != []):
    rateFile = re.findall("[a-z.\/]+",InpLines[i].split("=")[1])
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
  reactionID = int(line[0])
  rxnClass = reaction(reactionID)
  rxnClass.RateFlag = line[2]
  #Regex expression to search for reactants and products
  rxnClass.Equation = line[1]
  rxnClass.acoeff = float(line[3])
  rxnClass.bcoeff = float(line[4])
  rxnClass.ccoeff = float(line[5])
  print(rxnClass.idNum)
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
      rxnClass.Reactants[rts] = rxnClass.Reactants[rts].replace(var.group(0), "")
  for pts in range(len(rxnClass.Products)):
    var = stoich_pattern.match(rxnClass.Products[pts])
    if(var == None):
      rxnClass.Pstoic[pts] = 1
    if(var != None):
      rxnClass.Pstoic[pts] = float(var.group(0))
      rxnClass.Products[pts] = rxnClass.Products[pts].replace(var.group(0), "")
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
        
        
chemnet = nx.DiGraph()
#Directed graph between reactants and products
#Edge comment based on reactants
#Combining subgraphs