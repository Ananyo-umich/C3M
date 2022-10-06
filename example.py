from pylab import *
import netCDF4 as nc
from glob import glob
import cantera as ct
import scipy.integrate

class ReactorOde:
    def __init__(self, gas):
        # Parameters of the ODE system and auxiliary data are stored in the
        # ReactorOde object.
        self.gas = gas
        self.P = gas.P
        self.T = gas.T

    def __call__(self, t, y):
        """the ODE function, y' = f(t,y) """

        # State vector is [T, Y_1, Y_2, ... Y_K]
        self.gas.set_unnormalized_mass_fractions(y[1:])
        #self.gas.TP = self.T, self.P
        rho = self.gas.density

        wdot = self.gas.net_production_rates
        dYdt = wdot * self.gas.molecular_weights / rho
        return hstack((0,dYdt))



#Creatng chemistry network#Creating reaction H2SO4 + H2O <=> SO3 + 2H2O

r1 = ct.ElementaryReaction({'H2SO4' : 1.0 , 'H2O' : 1.0}, {'SO3' : 1.0 , 'H2O' : 2.0}) 
#Creating reaction SO3 + CO <=> SO2 + CO2

r2 = ct.ElementaryReaction({'SO3' : 1.0 , 'CO' : 1.0}, {'SO2' : 1.0 , 'CO2' : 1.0}) 
#Creating reaction SO3 + OCS <=> S2O2 + CO2

r3 = ct.ElementaryReaction({'SO3' : 1.0 , 'OCS' : 1.0}, {'S2O2' : 1.0 , 'CO2' : 1.0}) 

gas = ct.Solution("example.yaml")
gas.TPY = 200, 10, 'SO2:1e-5,CO2:0.96,H2O:1e-4,S2O2:1e-7,H2SO4:1e-5,CO:1e-6,SO3:1e-12,OCS:1e-9'
r = ct.IdealGasReactor(gas)
#net = ct.ReactorNet([r])
#gas.equilibrate("HP")
y0 = hstack((gas.T,gas.X))
#Setting up the ODE solver
ode = ReactorOde(gas)
solver = scipy.integrate.ode(ode)
solver.set_integrator('vode', method='bdf', with_jacobian=True)
solver.set_initial_value(y0, 0.0)

#print(gas.net_production_rates)
# Integrate the equations, keeping T(t) and Y(k,t)
t_end = 1e-1
states = ct.SolutionArray(gas, 1, extra={'t': [0.0]})
dt = 1e-5
while solver.successful() and solver.t < t_end:
    solver.integrate(solver.t + dt)
    gas.X = solver.y[1:]
    states.append(gas.state, t=solver.t)
    
 


L1 = plt.plot(states.t, states('H2O').Y, color='b', label='T', lw=2)
L2 = plt.plot(states.t, states('OCS').Y, color='r', label='T', lw=2)
L3 = plt.plot(states.t, states('SO2').Y, color='C1', label='T', lw=2)
print(states('H2O').Y)
xlabel('time(s)')
xscale('log')
yscale('log')
ylabel('mole fraction')
savefig('sample_result.png')


