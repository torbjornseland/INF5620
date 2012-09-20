# /usr/bin/env python

import nose.tools as nt
import vertical_motion as vm
from numpy import zeros, sqrt

def analytical(v, g, dt):
	g = -9.81
	for i in range(1, len(v)):
		v[i] = v[0] + g*dt*i
	return v


def test_init():
	dic = {}
	dic['rho'] = 0 
	dic['my'] = 18.27*10**(-6)
	dic['T0'] =  288
	dic['Rs'] = 8.34
	dic['L'] = 0.0065
	dic['M'] =  0.029
	dic['rho_b'] =  1003
	dic['m'] = 80
	dic['d'] = 0.0
	dic['C_D'] =  1.0 
	dic['V'] = 0.08
	dic['A'] =  0.9
	dic['C_Dp'] =  1.5
	dic['Tp'] =  20
	dic['Ap'] = 12

	v0 = 0
	T = 20
	dt = 0.001
	g = -9.81
	N = int(T/float(dt))
	time, v = vm.init(v0, T, dt, N, dic)
	v1 = zeros(N)
	v1[0] = v0
	v1 = analytical(v1,g, dt)

	# Testing for linear solution
	for i in range(len(v)):
		nt.assert_almost_equal(v[i], v1[i], delta=1E-10)

	# Testing terminal velocity
	dic['rho'] = 0.79
	dic['d'] = 0.5
	T = 50
	N = int(T/float(dt))
	time, v = vm.init(v0, T, dt, N, dic)

	a = .5*dic['C_D']*(dic['rho']*dic['A'])/(dic['rho_b']*dic['V'])
	b = g*(1-(dic['rho']/dic['rho_b']))
	t_m = -sqrt(abs(b)/a)

	nt.assert_almost_equal(v[-1], t_m, delta=1E-5)

	
	




