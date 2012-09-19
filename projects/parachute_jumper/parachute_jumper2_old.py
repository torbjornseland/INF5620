#/usr/bin/env python
import vertical_motion_force_plot as vmfp
import vertical_motion as vm
from scitools.all import *



class Forwad_Euler_parachute_jumper:
    """
    Simulate a parachute jumper jumping from a sepsified height 
    until he/she hits the ground with the forward euler method
    """
    def __init__(self,parameters):
        """
        :param parameters: All parameters needed to simulate a parachute 
                    jump in a dictionary. 
        Dictionary should include
                dt - timetep
                rho0 - desity of air at z0 height
                my - dynamic viscosity of air
                d - diameter of parachuter
                C_D - drag coeffisient
                A - cross-section area of parachuter
                V - volume of parachuter
                m - mass of parachuter
                T - temprature
                L - TODO
                Rs - TODO
                M - TODO
                Ts - TODO
                R - TODO
        Dictionary can include
                Tg - time when parachute opens
                Ag - area of parachute
                C_Dg - drag coefisient of parachute
        """

        par = parameters
        self.dt = float(par['dt'])
        self.rho0 = float(par['rho0'])
        self.my = float(par['my'])
        self.d = float(par['d'])
        self.C_D = float(par['C_D'])
        self.A = float(par['A'])
        self.V = float(par['V'])
        self.m = float(par['m'])
        self.T = float(par['T'])
        self.L = float(par['L'])
        self.Rs = float(par['Rs'])
        self.M = float(par['M'])
        self.Ts = float(par['Ts'])
        self.R = float(par['R'])
        self.z = [float(par['z0'])]
        self.v = [float(par['v0'])]

        self.p = (self.rho0*self.R*self.Ts)/self.M
        self.rho = rho0
        self.epsilon = 0.1

        if('Tp' in par.keys()):
            self.Tp = float(par['Tp'])
            self.Ap = float(par['Ap'])
            self.C_Dp = float(par['C_Dp'])
		
		else:
		    self.Tp = None;
	
    def iteration(self):
        """Solve for one iteration"""
        self.z.append(altitude(self.z[-1], self.v[-1], self.dt))
        #p has to be updated before rho is updated, 
        #use z[-2] to take the correct timestep
        self.p = decay_of_pressure(self.p, self.z[-2], self.M, self.Rs, self.T, self.L, self.dt)
        self.v.append(velocity(v[-1], self.rho, self.d, self.my, self.C_D, self.A, self.dt, self.V, self.m, self.epsilon))
        self.rho = density(self.p, self.M, self.R, self.Ts)

    def solve(self):
        """Slove until patashuter hits the ground"""
        #Tg is given, if test for optimization
        if(not self.Tg == None):
            flag = True;
            count = 0;
            while(self.z[-1] > 0):
                #parachute open
                if(dt*count >= Tg and flag):
                    self.A = self.Ap;
                    self.C_D = self.C_Dp;
                    flag = False;
                self.iteration();
                count += 1;
        
        #Tg not given, iteration as normal
        else:
            while(self.z[-1] > 0):
                self.iteration();

        return self.z, self.v



def altitude(z_1, v_1, dt):
    """One step of alitude with forward euler"""
	return z_1 + dt*v_1

def decay_of_pressure(p_1, z_1, M, Rs, T, L, dt, g=-9.81):
    """One step of pressure with forward euler"""
	return p_1 + (dt*M*g)/(Rs*(T+L*z_1))

def density(p, M, R, Ts):
    """Convertion of pressure to density"""
	return p*M/float(R*Ts)

def velocity(v_1, rho_1, d, my, C_D, A, dt, V, m, epsilon):
    """One step of velosity(vertical motion) with forward euler"""
	re = vm.re(v_1,d,rho_1,my)
	f_b = vmfp.force_buoyancy(rho_1,V)
	f_g = vmfp.force_gravity(m)
	if re<epsilon:
	    #Use stokes drag model
		f_d = vmfp.force_stokes_drag(d, my, v_1)
	else:
	    #Use quadratic drag model
		f_d = vmfp.force_quadratic_drag(C_D, rho_1, A, v_1)
	v = v_1 + (dt/float(m))*(f_d + f_b + f_g)
	return v

def forward_euler_parachute_jumper(parameters):
	par = parameters
	dt = float(par['dt'])
	v0 = float(par['v0'])
	z0 = float(par['z0'])
	rho0 = float(par['rho0'])
	my = float(par['my'])
	d = float(par['d'])
	C_D = float(par['C_D'])
	A = float(par['A'])
	V = float(par['V'])
	m = float(par['m'])
	T = float(par['T'])
	L = float(par['L'])
	Rs = float(par['Rs'])
	M = float(par['M'])
	Ts = float(par['Ts'])
	R = float(par['R'])

	v = [v0]
	z = [z0]
	p = (rho0*R*Ts)/M
	rho = rho0
	epsilon = 0.1

	if 'Tp' in par.keys():
		count = 0
		flag = True
		Ap = par['Ap']
		Tp = par['Tp']
		C_Dp = par['C_Dp']
		while z[-1]>0:
			if dt*count >= Tp and flag:
				A = Ap
				C_D = C_Dp
				flag = False
			#p has to be updated before rho is updated
			z.append(altitude(z[-1], v[-1], dt))
			p = decay_of_pressure(p, z[-2], M, Rs, T, L, dt)
			v.append(velocity(v[-1], rho, d, my, C_D, A, dt, V, m, epsilon))
			rho = density(p, M, R, Ts)
			count += 1
	else:	
		while z[-1]>0:
			#p has to be updated before rho is updated
			z.append(altitude(z[-1], v[-1], dt))
			p = decay_of_pressure(p, z[-2], M, Rs, T, L, dt)
			v.append(velocity(v[-1], rho, d, my, C_D, A, dt, V, m, epsilon))
			rho = density(p, M, R, Ts)
	return z, v

if __name__ == '__main__':
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

	z,v  = forward_euler_parachute_jumper(par)
	t = linspace(0, par['dt']*len(v), len(v))
	file = open('results.txt', 'w')
	file.write('z\n')
	print z[0]
	print z[20]
	print z[100]
	print z[-1]

	figure()
	plot(t,z)
	figure()
	plot(t, v)

	fe = Forwad_Euler_parachute_jumper(par)
	z,v = fe.solve()
	figure()
	plot(t,z)
	figure()
	plot(t, v)
	raw_input()


