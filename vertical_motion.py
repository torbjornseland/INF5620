from numpy import  array, linspace, zeros, pi;
import matplotlib.pyplot as plt
from parachutemodule import re


def crank_nicholsen_stokes_drag(a,b,dt,v_1):
    return ((-a*dt*0.5 + 1)*v_1 + b*dt)/float(1 + a*dt*0.5)

def crank_nicholsen_quadratic_drag(a,b,dt,v_1):
    return (v_1 + dt*b)/float(1+dt*a*abs(v_1)) 

def crank_nicholsen(rho,rho_b,my,d,C_D,V,A,epsilon,dt,v_1, g = -9.81):
    """
    Compute one velosity step of a vertical movment in fluid
    :param dt: Timestep
    :param v_1: Initial velosity
    :param rho: Density of fluid
    :param rho_b: Density of object
    :param my: Dynamic viscosity
    :param d: Diameter of object perpendicular to flow
    :param C_D: Drag coefficient
    :param V: Volume of object
    :param A: Cross-sectional area of object
    """
    b = g*(1-rho/float(rho_b))
    if(re(v_1,d,rho,my)< epsilon):
        a = (3*pi*d*my)/(V*rho_b);
        return crank_nicholsen_stokes_drag(a,b,dt,v_1);
    else:
        a = 0.5*C_D*(rho*A)/(rho_b*V);
        return crank_nicholsen_quadratic_drag(a,b,dt,v_1);

if __name__ == '__main__':
	rho = 1.2041 #density of air 20 degrees C [kg/m^3]
	rho_b = 1000 #density of human body [kg/m^3]
	v0 = 0; #starting velosity [m/s]
	my = 18.27*10**(-6) # viscosity of air [Pa*s]
	m = 75 #mass of avrage human [kg]
	d = 0.5 #diameter of human [m]
	C_D = 1.15 #drag coefficient, long sylinder[-]
	V = 0.075 #volume of human [m^3]
	A = 0.3*1.8 #Cross section of human [m^2]

	T_final = 120 #final time[s]
	N = 100000 #nr of datapoints
	t_list = linspace(0,T_final,N);
	dt = T_final/float(N-1);

	epsilon = 0.1 #reynolds number stokes drag limit

	v_list = zeros(N);
	v_list[0] = v0;

	for i in range(N-1):
		v_list[i+1] = crank_nicholsen(rho,rho_b,my,d,C_D,V,A,epsilon,dt,v_list[i], g = -9.81)

    #Plotting
	plt.plot(t_list,v_list)
	plt.show()
	print v_list
	raw_input();


