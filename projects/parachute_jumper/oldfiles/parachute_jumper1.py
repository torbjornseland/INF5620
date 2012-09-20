#/usr/bin/env python
import vertical_motion as vm
from numpy import linspace, zeros
import	matplotlib.pyplot as plt

if __name__ == '__main__':
	rho = 0.79 #density of air 20 degrees C [kg/m^3]
	rho_b = 1003 #density of human body [kg/m^3]
	v0 = 0; #starting velosity [m/s]
	my = 18.27*10**(-6) # viscosity of air [Pa*s]
	m = 80 #mass of avrage human [kg]
	d = 0.5 #diameter of human [m]
	C_D = 1. #drag coefficient, human [-]
	V = 0.08 #volume of human [m^3]
	A = 0.9 #Cross section of human [m^2]

	theta = 0.5

	T = 50 #final time[s]
	N = 100000 #nr of datapoints
	t_list = linspace(0,T,N);
	dt = T/float(N-1);

	epsilon = 0.1 #reynolds number stokes drag limit

	v_list = zeros(N);
	v_list[0] = v0;

	flag = True
	for i in range(N-1):
		v_list[i+1] = vm.crank_nicholsen(rho,rho_b,my,d,C_D,V,A,epsilon,dt, v_list[i], g = -9.81)

    #Plotting
	plt.plot(t_list,v_list)
	plt.show()
	raw_input();
