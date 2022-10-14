from pylab import * 
import netCDF4 as nc 
from glob import glob


T = 227.0
P = 1197.0
Kb = 1.38e-23
N = P/(Kb*T) 

def Jacobian(x):
	xjac = zeros([4,4])
	xjac[0,0] = 	xjac[0,0] - (6e-11*1.0*pow(x[0],0.0))#O2 = 2O
	xjac[0,0] = 	xjac[0,0] - (6e-34*pow(300/T,2.4)*exp(0/T)*N*pow(x[1],1.0)*1.0*pow(x[0],0.0))#O + O2 = O3
	xjac[0][0] = 	xjac[0,0] + (7e-4*pow(x[2],1.0))#O3 = O + O2
	xjac[0][0] = 	xjac[0,0] + (5e-4*pow(x[2],1.0))#O3 = O_1D + O2
	xjac[0][0] = 	xjac[0,0] + (8e-12*pow(300/T,0)*exp(-2060/T)*N*pow(x[1],1.0)*pow(x[2],1.0))#O + O3 = 2O2
	xjac[0,1] = 	xjac[0,1] - (6e-11*pow(x[0],1.0))#O2 = 2O
	xjac[0,1] = 	xjac[0,1] - (6e-34*pow(300/T,2.4)*exp(0/T)*N*1.0*pow(x[1],0.0)*pow(x[0],1.0))#O + O2 = O3
	xjac[0][1] = 	xjac[0,1] + (7e-4*pow(x[2],1.0))#O3 = O + O2
	xjac[0][1] = 	xjac[0,1] + (5e-4*pow(x[2],1.0))#O3 = O_1D + O2
	xjac[0][1] = 	xjac[0,1] + (8e-12*pow(300/T,0)*exp(-2060/T)*N*1.0*pow(x[1],0.0)*pow(x[2],1.0))#O + O3 = 2O2
	xjac[0,2] = 	xjac[0,2] - (6e-11*pow(x[0],1.0))#O2 = 2O
	xjac[0,2] = 	xjac[0,2] - (6e-34*pow(300/T,2.4)*exp(0/T)*N*pow(x[1],1.0)*pow(x[0],1.0))#O + O2 = O3
	xjac[0][2] = 	xjac[0,2] + (7e-4*1.0*pow(x[2],0.0))#O3 = O + O2
	xjac[0][2] = 	xjac[0,2] + (5e-4*1.0*pow(x[2],0.0))#O3 = O_1D + O2
	xjac[0][2] = 	xjac[0,2] + (8e-12*pow(300/T,0)*exp(-2060/T)*N*pow(x[1],1.0)*1.0*pow(x[2],0.0))#O + O3 = 2O2
	xjac[0,3] = 	xjac[0,3] - (6e-11*pow(x[0],1.0))#O2 = 2O
	xjac[0,3] = 	xjac[0,3] - (6e-34*pow(300/T,2.4)*exp(0/T)*N*pow(x[1],1.0)*pow(x[0],1.0))#O + O2 = O3
	xjac[0][3] = 	xjac[0,3] + (7e-4*pow(x[2],1.0))#O3 = O + O2
	xjac[0][3] = 	xjac[0,3] + (5e-4*pow(x[2],1.0))#O3 = O_1D + O2
	xjac[0][3] = 	xjac[0,3] + (8e-12*pow(300/T,0)*exp(-2060/T)*N*pow(x[1],1.0)*pow(x[2],1.0))#O + O3 = 2O2

	xjac[1][0] = 	xjac[1,0] + (6e-11*1.0*pow(x[0],0.0))#O2 = 2O
	xjac[1,0] = 	xjac[1,0] - (6e-34*pow(300/T,2.4)*exp(0/T)*N*pow(x[1],1.0)*1.0*pow(x[0],0.0))#O + O2 = O3
	xjac[1][0] = 	xjac[1,0] + (7e-4*pow(x[2],1.0))#O3 = O + O2
	xjac[1,0] = 	xjac[1,0] - (8e-12*pow(300/T,0)*exp(-2060/T)*N*pow(x[1],1.0)*pow(x[2],1.0))#O + O3 = 2O2
	xjac[1][0] = 	xjac[1,0] + (3.2e-11*pow(300/T,0)*exp(0/T)*N*pow(x[3],1.0))#O_1D = O
	xjac[1][1] = 	xjac[1,1] + (6e-11*pow(x[0],1.0))#O2 = 2O
	xjac[1,1] = 	xjac[1,1] - (6e-34*pow(300/T,2.4)*exp(0/T)*N*1.0*pow(x[1],0.0)*pow(x[0],1.0))#O + O2 = O3
	xjac[1][1] = 	xjac[1,1] + (7e-4*pow(x[2],1.0))#O3 = O + O2
	xjac[1,1] = 	xjac[1,1] - (8e-12*pow(300/T,0)*exp(-2060/T)*N*1.0*pow(x[1],0.0)*pow(x[2],1.0))#O + O3 = 2O2
	xjac[1][1] = 	xjac[1,1] + (3.2e-11*pow(300/T,0)*exp(0/T)*N*pow(x[3],1.0))#O_1D = O
	xjac[1][2] = 	xjac[1,2] + (6e-11*pow(x[0],1.0))#O2 = 2O
	xjac[1,2] = 	xjac[1,2] - (6e-34*pow(300/T,2.4)*exp(0/T)*N*pow(x[1],1.0)*pow(x[0],1.0))#O + O2 = O3
	xjac[1][2] = 	xjac[1,2] + (7e-4*1.0*pow(x[2],0.0))#O3 = O + O2
	xjac[1,2] = 	xjac[1,2] - (8e-12*pow(300/T,0)*exp(-2060/T)*N*pow(x[1],1.0)*1.0*pow(x[2],0.0))#O + O3 = 2O2
	xjac[1][2] = 	xjac[1,2] + (3.2e-11*pow(300/T,0)*exp(0/T)*N*pow(x[3],1.0))#O_1D = O
	xjac[1][3] = 	xjac[1,3] + (6e-11*pow(x[0],1.0))#O2 = 2O
	xjac[1,3] = 	xjac[1,3] - (6e-34*pow(300/T,2.4)*exp(0/T)*N*pow(x[1],1.0)*pow(x[0],1.0))#O + O2 = O3
	xjac[1][3] = 	xjac[1,3] + (7e-4*pow(x[2],1.0))#O3 = O + O2
	xjac[1,3] = 	xjac[1,3] - (8e-12*pow(300/T,0)*exp(-2060/T)*N*pow(x[1],1.0)*pow(x[2],1.0))#O + O3 = 2O2
	xjac[1][3] = 	xjac[1,3] + (3.2e-11*pow(300/T,0)*exp(0/T)*N*1.0*pow(x[3],0.0))#O_1D = O

	xjac[2][0] = 	xjac[2,0] + (6e-34*pow(300/T,2.4)*exp(0/T)*N*pow(x[1],1.0)*1.0*pow(x[0],0.0))#O + O2 = O3
	xjac[2,0] = 	xjac[2,0] - (7e-4*pow(x[2],1.0))#O3 = O + O2
	xjac[2,0] = 	xjac[2,0] - (5e-4*pow(x[2],1.0))#O3 = O_1D + O2
	xjac[2,0] = 	xjac[2,0] - (8e-12*pow(300/T,0)*exp(-2060/T)*N*pow(x[1],1.0)*pow(x[2],1.0))#O + O3 = 2O2
	xjac[2][1] = 	xjac[2,1] + (6e-34*pow(300/T,2.4)*exp(0/T)*N*1.0*pow(x[1],0.0)*pow(x[0],1.0))#O + O2 = O3
	xjac[2,1] = 	xjac[2,1] - (7e-4*pow(x[2],1.0))#O3 = O + O2
	xjac[2,1] = 	xjac[2,1] - (5e-4*pow(x[2],1.0))#O3 = O_1D + O2
	xjac[2,1] = 	xjac[2,1] - (8e-12*pow(300/T,0)*exp(-2060/T)*N*1.0*pow(x[1],0.0)*pow(x[2],1.0))#O + O3 = 2O2
	xjac[2][2] = 	xjac[2,2] + (6e-34*pow(300/T,2.4)*exp(0/T)*N*pow(x[1],1.0)*pow(x[0],1.0))#O + O2 = O3
	xjac[2,2] = 	xjac[2,2] - (7e-4*1.0*pow(x[2],0.0))#O3 = O + O2
	xjac[2,2] = 	xjac[2,2] - (5e-4*1.0*pow(x[2],0.0))#O3 = O_1D + O2
	xjac[2,2] = 	xjac[2,2] - (8e-12*pow(300/T,0)*exp(-2060/T)*N*pow(x[1],1.0)*1.0*pow(x[2],0.0))#O + O3 = 2O2
	xjac[2][3] = 	xjac[2,3] + (6e-34*pow(300/T,2.4)*exp(0/T)*N*pow(x[1],1.0)*pow(x[0],1.0))#O + O2 = O3
	xjac[2,3] = 	xjac[2,3] - (7e-4*pow(x[2],1.0))#O3 = O + O2
	xjac[2,3] = 	xjac[2,3] - (5e-4*pow(x[2],1.0))#O3 = O_1D + O2
	xjac[2,3] = 	xjac[2,3] - (8e-12*pow(300/T,0)*exp(-2060/T)*N*pow(x[1],1.0)*pow(x[2],1.0))#O + O3 = 2O2

	xjac[3][0] = 	xjac[3,0] + (5e-4*pow(x[2],1.0))#O3 = O_1D + O2
	xjac[3,0] = 	xjac[3,0] - (3.2e-11*pow(300/T,0)*exp(0/T)*N*pow(x[3],1.0))#O_1D = O
	xjac[3][1] = 	xjac[3,1] + (5e-4*pow(x[2],1.0))#O3 = O_1D + O2
	xjac[3,1] = 	xjac[3,1] - (3.2e-11*pow(300/T,0)*exp(0/T)*N*pow(x[3],1.0))#O_1D = O
	xjac[3][2] = 	xjac[3,2] + (5e-4*1.0*pow(x[2],0.0))#O3 = O_1D + O2
	xjac[3,2] = 	xjac[3,2] - (3.2e-11*pow(300/T,0)*exp(0/T)*N*pow(x[3],1.0))#O_1D = O
	xjac[3][3] = 	xjac[3,3] + (5e-4*pow(x[2],1.0))#O3 = O_1D + O2
	xjac[3,3] = 	xjac[3,3] - (3.2e-11*pow(300/T,0)*exp(0/T)*N*1.0*pow(x[3],0.0))#O_1D = O

	return xjac

def ReactionRate(x):
	dx = zeros(4)
	dx[0] = 	dx[0] - (6e-11*pow(x[0],1.0))#O2 = 2O
	dx[0] = 	dx[0] - (6e-34*pow(300/T,2.4)*exp(0/T)*N*pow(x[1],1.0)*pow(x[0],1.0))#O + O2 = O3
	dx[0] = 	dx[0] + (7e-4*pow(x[2],1.0))#O3 = O + O2
	dx[0] = 	dx[0] + (5e-4*pow(x[2],1.0))#O3 = O_1D + O2
	dx[0] = 	dx[0] + (8e-12*pow(300/T,0)*exp(-2060/T)*N*pow(x[1],1.0)*pow(x[2],1.0))#O + O3 = 2O2

	dx[1] = 	dx[1] + (6e-11*pow(x[0],1.0))#O2 = 2O
	dx[1] = 	dx[1] - (6e-34*pow(300/T,2.4)*exp(0/T)*N*pow(x[1],1.0)*pow(x[0],1.0))#O + O2 = O3
	dx[1] = 	dx[1] + (7e-4*pow(x[2],1.0))#O3 = O + O2
	dx[1] = 	dx[1] - (8e-12*pow(300/T,0)*exp(-2060/T)*N*pow(x[1],1.0)*pow(x[2],1.0))#O + O3 = 2O2
	dx[1] = 	dx[1] + (3.2e-11*pow(300/T,0)*exp(0/T)*N*pow(x[3],1.0))#O_1D = O

	dx[2] = 	dx[2] + (6e-34*pow(300/T,2.4)*exp(0/T)*N*pow(x[1],1.0)*pow(x[0],1.0))#O + O2 = O3
	dx[2] = 	dx[2] - (7e-4*pow(x[2],1.0))#O3 = O + O2
	dx[2] = 	dx[2] - (5e-4*pow(x[2],1.0))#O3 = O_1D + O2
	dx[2] = 	dx[2] - (8e-12*pow(300/T,0)*exp(-2060/T)*N*pow(x[1],1.0)*pow(x[2],1.0))#O + O3 = 2O2

	dx[3] = 	dx[3] + (5e-4*pow(x[2],1.0))#O3 = O_1D + O2
	dx[3] = 	dx[3] - (3.2e-11*pow(300/T,0)*exp(0/T)*N*pow(x[3],1.0))#O_1D = O

	return dx

t_init = 0.0
t_final = 1000.0
n = 100000
x = 1e-10*N*abs(randn(4,n))
rate = zeros(4)
jacob = zeros([4,4])
x[1] =3e7
x[2] =1e12
x[0] =7e16
dt = 0.01
for i in range(n-1):
	j_matrix = Jacobian(x[:,i])
	x_pred = x[:,i]
	rate = ReactionRate(x[:,i])
	i_matrix = identity(len(x[:,i]))
	x_prev = x_pred
	x_next = x_pred
#Backward Euler Integration Scheme (Li and Chen, 2020)
	while True:
		b_matrix = (rate - ((x_prev - x[:,i])/dt))
		x_int = x_next
		x_next = x_prev + matmul(linalg.inv((i_matrix/dt) - j_matrix),b_matrix)
		x_prev = x_int
		print(str(abs(x_next - x_prev)/x_prev))
		if (abs(x_next - x_prev)/x_prev <= 1e-10).all():
			break
	x[:,i+1] = x_next
	
	print('Chemical Reaction Network solved for '+str(i)+' Time Step')

#Writing the Output into NETCDF format
OutFile = 'chapman.nc' 
data = nc.Dataset(OutFile, 'w' ,format = 'NETCDF4')
pres = data.createDimension('press',1)
temp = data.createDimension('temp',1)
P = data.createVariable('press', 'f8', ('press',))
T = data.createVariable('temp', 'f8', ('temp',))
time = data.createDimension('time',n)
t = data.createVariable('time', 'f8', ('time',))
O2= data.createVariable('O2', 'f8', ('time',))
O2[:] = x[0,:]
O= data.createVariable('O', 'f8', ('time',))
O[:] = x[1,:]
O3= data.createVariable('O3', 'f8', ('time',))
O3[:] = x[2,:]
O_1D= data.createVariable('O_1D', 'f8', ('time',))
O_1D[:] = x[3,:]
data.close()