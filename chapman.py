from pylab import *
import netCDF4 as nc
from glob import glob
import cantera as ct


T = 200.0
P = 100000.0
#Creatng chemistry network#Creating reaction O2 <=> 2O

r1 = ct.ElementaryReaction({'O2' : 1.0}, {'O' : 2.0}) 
#Creating reaction O + O2 <=> O3

r2 = ct.ElementaryReaction({'O' : 1.0 , 'O2' : 1.0}, {'O3' : 1.0}) 
#Creating reaction O3 <=> O + O2

r3 = ct.ElementaryReaction({'O3' : 1.0}, {'O' : 1.0 , 'O2' : 1.0}) 
#Creating reaction O3 <=> O_1D + O2

r4 = ct.ElementaryReaction({'O3' : 1.0}, {'O_1D' : 1.0 , 'O2' : 1.0}) 
#Creating reaction O + O3 <=> 2O2

r5 = ct.ElementaryReaction({'O' : 1.0 , 'O3' : 1.0}, {'O2' : 2.0}) 
#Creating reaction O_1D <=> O

r6 = ct.ElementaryReaction({'O_1D' : 1.0}, {'O' : 1.0}) 

gas = ct.Solution("example.yaml")
r = ct.IdealGasReactor(gas)
net = ct.ReactorNet([r])
T = r.T
gas.equilibrate("HP")
element = "S"
diagram = ct.ReactionPathDiagram(gas,element)
diagram.title = "Reaction path diagram following {0}".format(element)
diagram.label_threshold = 0.01
img_file = "rxnpath.png"
dot_file = "rxnpath.dot"
diagram.write_dot(dot_file)