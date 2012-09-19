#/usr/bin/env python
from numpy import array, pi

def force_gravity(m, g = -9.81):
	return m*g;

def force_buoyancy(rho,V,g = -9.81):
	return -rho*V*g;

def force_stokes_drag(d, my, v):
	return -3*pi*d*my*v

def force_quadratic_drag(C_D, rho, A, v):
	return -.5*C_D*rho*A*abs(v)*v

def re(v,d,rho,my):
    return (d*abs(v)*rho)/my;

