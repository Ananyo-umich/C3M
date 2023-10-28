from pylab import *

data_1 = genfromtxt('O2diss_linear.txt')
data_2 = genfromtxt('O2diss_log.txt')


time1 = data_1[:,0]
time2 = data_2[:,0]

#O2
O2_d1 = data_1[:,1]
O2_d2 = data_2[:,1]
print(time)


#O
O_d1 = data_1[:,2]
O_d2 = data_2[:,2]

plot(time1, O2_d1, 'r-')
plot(time2, O2_d2, 'b-')



plot(time1, O_d1, 'r--')
plot(time2, O_d2, 'b--')
legend(['0.7 AU (O$_{2}$) linear','0.7 AU (O$_{2}$) log','0.7 AU (O) linear','0.7 AU (O) log'], loc = 'best')
yscale('log')
xscale('log')
xlabel('Time (s)')
ylabel('X')
savefig('O2diss_integ_comp.png')