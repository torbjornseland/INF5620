#/usr/bin/env python

import parachute_jumper2 as pj
import nose.tools as nt

def test_foward_euler_parachute_jumper():
	#res_z = [z[0], z[20], z[100], z[-1]]
	#res_v = [v[0], v[20], v[100], v[-1]]
	res_z = [5000.0000000000000000, 4999.8138034722796874, 4995.1816863121885035, -0.0886208033146170]
	res_v = [0.0000000000000000, -1.9593960265106347, -9.6643865086380423, -10.5020300807767963]

	#init values
	par = {}
	par['dt'] = 0.01 #Timestep
	par['v0'] = 0 #Starting velocity
	par['z0'] = 5000 #Starting altitude
	par['rho0'] = 0.79 #density of air [kg/m^3]
	par['my'] = 18.27*10**(-6) #Viscosity of air [Pa*s]
	par['d'] = 0.5 #diameter of human [m]
	par['C_D'] = 1 #Drag coefficient
	par['A'] = 0.9 #Area of human [m^2]
	par['V'] = 0.08 #Volume of human [m^3]
	par['m'] = 80 #Mass of human [kg]
	par['T'] = 288 #Temperatur [K]
	par['L'] = -0.0065 #
	par['Rs'] = 8.314 #
	par['M'] = 0.029 #
	par['Ts'] = 288 #
	par['R'] = 8.314 #
	par['Tp'] = 20 # Time when parachute opens [s]
	par['Ap'] = 12 # Area of parachute [m^2]
	par['C_Dp'] = 1.5 # Drag coefficient parachute 

	z,v = pj.Forward_Euler_parachute_jumper(par).solve()

	nt.assert_almost_equal(z[0], res_z[0], delta=1E-10)
	nt.assert_almost_equal(z[20], res_z[1], delta=1E-10)
	nt.assert_almost_equal(z[100], res_z[2], delta=1E-10)
	nt.assert_almost_equal(z[-1], res_z[3], delta=1E-10)
	nt.assert_almost_equal(v[0], res_v[0], delta=1E-10)
	nt.assert_almost_equal(v[20], res_v[1], delta=1E-10)
	nt.assert_almost_equal(v[100], res_v[2], delta=1E-10)
	nt.assert_almost_equal(v[-1], res_v[3], delta=1E-10)

