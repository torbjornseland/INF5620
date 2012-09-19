#/usr/bin/env python
import matplotlib.pyplot as plt
from parachutemodule import *
from numpy import linspace, zeros
from vertical_motion import crank_nicholsen

def get_drag_force(v):
    if(re(v,d,rho,my)< epsilon):
        #use stokes drag model
        return force_stokes_drag(d,my,v);
    else:
        #use quadratic drag model
        return force_quadratic_drag(C_D,rho,A,v);



rho = 1.2041 #density of air 20 degrees C [kg/m^3]
rho_b = 1000 #density of human body [kg/m^3]
v0 = 0; #starting velosity [m/s]
my = 18.27*10**(-6) # viscosity of air [Pa*s]
m = 75 #mass of avrage human [kg]
d = 0.5 #diameter of human [m]
C_D = 1.15 #long sylinder[-]
V = 0.075 #volume of human [m^3]
A = 0.3*1.8 #Cross section of human [m^2]

T_final = 120 #final time[s]
N = 100000 #nr of datapoints
t_list = linspace(0,T_final,N);
dt = T_final/float(N-1); #timestep[s]

epsilon = 0.1 #reynolds number stokes drag limit

v_list = zeros(N);
v_list[0] = v0;
Fd_list = zeros(N);
Fd_list[0] = get_drag_force(v_list[0]);

#make list of the constant forces
gravity = zeros(N)+force_gravity(m);
buoyancy = zeros(N)+force_buoyancy(rho,V);

for i in range(N-1):
    v_list[i+1] = crank_nicholsen(rho,rho_b,my,d,C_D,V,A,epsilon,dt,v_list[i], g = -9.81)
    Fd_list[i+1] = get_drag_force(v_list[i+1]);

p1, = plt.plot(t_list,Fd_list,)
p2, = plt.plot(t_list,gravity)
p3, = plt.plot(t_list,buoyancy)
plt.legend([p1,p2,p3],["drag","gravity","buoyancy"]);
plt.show()
raw_input();


