#/usr/bin/env python

import nose.tools as nt
import vertical_motion as vm
from numpy import linspace, zeros

def test_crank_nicholsen():
	#results = [list[0], v_list[100], v_list[221], v_list[-100], v_list[-1]]
	results = [0.0000000000000000, -1.1755187262322282, -2.5955339812878586, \
			-44.3344607484978610, -44.3344607484978610]

	rho = 1.2041 #density of air 20 degrees C [kg/m^3]
	rho_b = 1000 #density of human body [kg/m^3]
	v0 = 0; #starting velosity [m/s]
	my = 18.27*10**(-6) # viscosity of air [Pa*s]
	m = 75 #mass of avrage human [kg]
	d = 0.5 #diameter of human [m]
	C_D = 1.15 #long sylinder[-]
	V = 0.075 #volume of human [m^3]
	A = 0.3*1.8 #0.0314 #Cross section of human [m^2]

	theta = 0.5

	T_final = 120 #final time[s]
	N = 100000 #nr of datapoints
	t_list = linspace(0,T_final,N);
	dt = T_final/float(N-1);

	epsilon = 0.1 #reynolds number stokes drag limit

	v_list = zeros(N);
	v_list[0] = v0;

	for i in range(N-1):
		v_list[i+1] = vm.crank_nicholsen(rho,rho_b,my,d,C_D,V,A,epsilon,dt,v_list[i], g = -9.81)

	nt.assert_almost_equal(v_list[0], results[0], delta=1E-10)
	nt.assert_almost_equal(v_list[100], results[1], delta=1E-10)
	nt.assert_almost_equal(v_list[221], results[2], delta=1E-10)
	nt.assert_almost_equal(v_list[-100], results[3], delta=1E-10)
	nt.assert_almost_equal(v_list[-1], results[4], delta=1E-10)
