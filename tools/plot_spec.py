from pylab import *


data = genfromtxt("sun_spec.inp")
wavelength = data[:, 0]
irradiance = data[:, 1]

print(irradiance)
plot(wavelength, irradiance)
# yscale('log')
# xlim([250,25000])
xlabel("wavelength")
ylabel("irradiance")
savefig("solar.png")

sum = 0
for i in range(len(irradiance) - 1):
    sum = sum + (irradiance[i] * (wavelength[i + 1] - wavelength[i]))

print(sum)
