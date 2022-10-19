#The program generates am executable python program for a chemistry network
#Notes:
#Input - inp file, spectrum file, chem network file
#Author: Ananyo Bhattacharya
#Email: ananyo@umich.edu
#Affiliation: University of Michigan, Ann Arbor

from pylab import *
from glob import glob
import multiprocessing as mp
import concurrent.futures
import re
from config import *
from glob import glob
import os
import cantera as ct

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
inpFile = str(args['i'])
reactionFileSearch = "reaction_database"
rateFileSearch = "rate_database"
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
  rxnClass.Reactants = re.findall(reactant_pattern, (line[1].split("<=>"))[0]) 
  rxnClass.Products = re.findall(product_pattern, (line[1].split("<=>"))[1]) 
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
            if molecule[0] != 'M':
              chemicals.append(molecule)
  stoi_reac = re.findall(stoich_pattern, line[1]) 
  
  


variables = array(chemicals)  
for varIndex in range(len(variables)):
  for rxn in reaction_array:
    for rts in range(len(rxn.Reactants)):
      if(variables[varIndex] == rxn.Reactants[rts]):
        rxn.ReactantID[rts] = varIndex 
    for pts in range(len(rxn.Products)):
      if(variables[varIndex] == rxn.Products[pts]): 
        rxn.ProductID[pts] = varIndex 

elements = ['C','H','O']
elementList = '[C, H, O]'

#Template for YAML files
'''
phases:
- name: gas1
  thermo: ideal-gas
  elements: [O, H, N, Ar]
  species: [H2, H, O, O2, OH, H2O, HO2, H2O2, N2, AR]
  kinetics: gas
  reactions: all
  transport: mixture-averaged
  state:
    T: 300.0
    P: 1.01325e+05
- name: gas2
  thermo: ideal-gas
  kinetics: gas
  transport: mixture-averaged
  state: {T: 300.0, 1 atm}

species:
- H2: ...
- H: ...
...
- AR: ...

reactions:
...

'''
#The stoichiometric coefficients have to be kept apart from the molecule name

#Writing yaml file for substrates and elements
with open(args['o'], 'w') as chem:
    ReactionSpecies = ''
    for cix in range(len(variables)):
        if(cix !=0):
            ReactionSpecies = ReactionSpecies + ", "
        ReactionSpecies = ReactionSpecies + variables[cix]
    chem.write('phases:')
    chem.write('\n- name: gas')
    chem.write('\n  thermo: ideal-gas')
    chem.write('\n  elements: ' + str(elementList))
    chem.write('\n  species: [' + str(ReactionSpecies) + ']')
    chem.write('\n  kinetics: gas')
    chem.write('\n  reactions: all')
    chem.write('\n  skip-undeclared-elements: true')
    chem.write('\n\nspecies:')
    for cix in range(len(variables)):
      composition = ''
      for eix in range(len(elements)):
        search_param = ''+ str(elements[eix]) + '[0-9]+|'+ str(elements[eix]) + ''
        VarIntr = re.findall(search_param, variables[cix])
        if(VarIntr != []):
          Atoms = zeros(len(VarIntr))
          for inx in range(len(VarIntr)):
            Atom = re.findall('[0-9]+', VarIntr[inx])
            if(Atom == []):
              Atoms[inx] = 1
            if(Atom != []):
              Atoms[inx] = int(Atom[0])
          if(composition != ''):
            composition = composition + ","
          composition = composition + str(elements[eix]) + ": " + str(sum(Atoms))
      chem.write('\n- name: ' + str(variables[cix]))
      chem.write('\n  composition: {' + str(composition) + '}')
      chem.write('\n  enthalpy: 0')
      chem.write('\n  thermo:')
      chem.write('\n    model: NASA7')
      chem.write('\n    temperature-ranges: [200.0, 1000.0, 3500.0]')
      chem.write('\n    data:')
      chem.write('\n    - [2.34433112, 7.98052075e-03, -1.9478151e-05, 2.01572094e-08, -7.37611761e-12,')
      chem.write('\n      -917.935173, 0.683010238]')
      chem.write('\n    - [3.3372792, -4.94024731e-05, 4.99456778e-07, -1.79566394e-10, 2.00255376e-14,')
      chem.write('\n      -950.158922, -3.20502331]')
    chem.write('\n\nreactions:')
    for rxn in reaction_array:
      reactant_string = ''
      product_string = ''
      for inx in range(len(rxn.Reactants)):
        if(rxn.Reactants[inx] != 'M'):
          if(reactant_string != ''):
            reactant_string = reactant_string + ' + ' + str(rxn.Rstoic[inx]) + ' ' + str(rxn.Reactants[inx])
          if(reactant_string == ''):
            reactant_string =  str(rxn.Rstoic[inx]) + ' ' + str(rxn.Reactants[inx])
      for inx in range(len(rxn.Products)):
        if(rxn.Products[inx] != 'M'):
          if(product_string != ''):
            product_string = product_string + ' + ' + str(rxn.Pstoic[inx]) + ' ' + str(rxn.Products[inx]) 
          if(product_string == ''):
            product_string =  str(rxn.Pstoic[inx]) + ' ' + str(rxn.Products[inx])
#Types of Rate Flags currently incorporated in the code:
# p_JPL, p_SWRI, p_VULCAN, t, t0, m0, m1, m2, m3
      if(rxn.RateFlag == 't'):
        chem.write("\n- equation: '" + reactant_string + " <=> " + product_string + " '")
        RateCoeff = 'A: ' + str(double(rxn.Data[3])) + ',T0: ' + str(double(rxn.Data[4])) + ',b: ' + str(double(rxn.Data[5])) + ',Ea: ' + str(double(rxn.Data[6]))
        chem.write('\n  rate-constant: {' + RateCoeff + ' }')
#rate-constant: {A: 6.E-34, T0: 300, b: -2.3}
      if(rxn.RateFlag == 't0'):
        chem.write("\n- equation: '" + reactant_string + " <=> " + product_string + " '")
        RateCoeff = 'A: ' + str(double(rxn.Data[3])) + ',T0: 300.0 ,b: 0.0 ,Ea: 0.0'
        chem.write('\n  rate-constant: {' + RateCoeff + ' }')
#Three-body reaction
      if(rxn.RateFlag == 'm0'):
        chem.write("\n- equation: '" + reactant_string + " + M <=> " + product_string + " + M '")
        chem.write('\n  type: three-body')
        RateCoeff = 'A: ' + str(double(rxn.Data[3])) + ',T0: ' + str(double(rxn.Data[4])) + ',b: ' + str(double(rxn.Data[5])) + ',Ea: ' + str(double(rxn.Data[6]))
        chem.write('\n  rate-constant: {' + RateCoeff + ' }')
        #chem.write('\n  efficiencies: {O2: 0.4, CO: 0.75, CO2: 1.5, H2O: 6.5, CH4: 0.48, H2: 1.0,')
        #chem.write('\n    C2H6: 1.0, C2H4: 1.6}')
#Three-body reaction
      if(rxn.RateFlag == 'm1'):
        chem.write("\n- equation: '" + reactant_string + " + M <=> " + product_string + " + M'")
        chem.write('\n  type: three-body')
        RateCoeff = 'A: ' + str(double(rxn.Data[3])) + ',T0: ' + str(double(rxn.Data[4])) + ',b: ' + str(double(rxn.Data[5])) + ',Ea: ' + str(double(rxn.Data[6]))

        chem.write('\n  rate-constant: {' + RateCoeff + ' }')
        #chem.write('\n  efficiencies: {O2: 0.4, CO: 0.75, CO2: 1.5, H2O: 6.5, CH4: 0.48, H2: 1.0,')
        #chem.write('\n    C2H6: 1.0, C2H4: 1.6}')
#Three-body reaction
      if(rxn.RateFlag == 'm2'):
        chem.write("\n- equation: '" + reactant_string + " (+ M) <=> " + product_string + " (+ M)'")
        chem.write('\n  type: falloff')
        RateCoeff_lowP = 'A: ' + str(double(rxn.Data[3])) + ',T0: ' + str(double(rxn.Data[4])) + ',b: ' + str(double(rxn.Data[5])) + ',Ea: ' + str(double(rxn.Data[6]))
        RateCoeff_highP = 'A: ' + str(double(rxn.Data[7])) + ',T0: ' + str(double(rxn.Data[8])) + ',b: ' + str(double(rxn.Data[9])) + ',Ea: ' + str(double(rxn.Data[10]))
        chem.write('\n  low-P-rate-constant: {'+ RateCoeff_lowP +'}')
        chem.write('\n  high-P-rate-constant: {'+ RateCoeff_highP +'}')
        chem.write('\n  Troe: {A: 1, T3: 94.0, T1: 10000.0}')
        #chem.write('\n  efficiencies: {O2: 0.4, CO: 0.75, CO2: 1.5, H2O: 6.5, CH4: 0.48, H2: 1.0,')
        #chem.write('\n    C2H6: 1.0, C2H4: 1.6}')

#Writing the reactions from the reaction database
#Using the Arrhenius term for reaction rate
    chem.close()
 
#Writing Executable Python File for Cantera
with open(args['ex'], 'w') as code:
    code.write('from pylab import *\nimport netCDF4 as nc\nfrom glob import glob\nimport cantera as ct\n\n')
    code.write('\nT = '+ str(Temp))
    code.write('\nP = '+ str(Press))
    code.write('\n#Creatng chemistry network')
    ReactionList = ''
    ReactionSpecies = ''
    for cix in range(len(variables)):
        if(cix !=0):
            ReactionSpecies = ReactionSpecies + ", "
        ReactionSpecies = ReactionSpecies + variables[cix]
    for rxn in reaction_array:
        if(reaction_array.index(rxn) !=0 ):
            ReactionList = ReactionList + " , "
        #code.write('#Creating reaction ' + str(rxn.Equation)+'\n' )
#Reactant string for Cantera elementary reaction
        Reactant_string = ''
        for rix in range(len(rxn.Reactants)):
            if(rix != 0):
                Reactant_string = Reactant_string + " , "
            Reactant_string = Reactant_string + "'"+str(rxn.Reactants[rix])+ "'"+" : "+str(rxn.Rstoic[rix])
        Product_string = ''
        for pix in range(len(rxn.Products)):
            if(pix != 0):
                Product_string = Product_string + " , "
            Product_string = Product_string + "'"+str(rxn.Products[pix])+ "'"+" : "+str(rxn.Pstoic[pix])
        ReactionIX = 'r'+str(reaction_array.index(rxn)+1)
        ReactionList = ReactionList + ReactionIX
        #code.write('\nr'+str(reaction_array.index(rxn)+1)+' = ct.ElementaryReaction({'+Reactant_string+'}, {'+Product_string+'}) \n')
    code.write('\ngas = ct.Solution("'+args['o']+'")')
    #code.write('\nrun("dot {0} -Tpng -o{1} -Gdpi=200".format(dot_file, img_file).split())')
    code.close()


    


