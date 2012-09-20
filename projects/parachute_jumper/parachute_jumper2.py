
class Problem:
    def __init__(self,parameters):
        self.parameters = parameters
        self.g = -9.81;

    def initial_condition(self, init_dict):
        par = self.parameters
        init_dict['p0'] = init_dict['rho0']*par['R']*self.T(init_dict['z0'])/float(par['M'])
        self.init_dict = init_dict

    def get_initial_condition(self):
        return init_dict;

    def dz(self, v):
        return v

    def dp(self,z):
        par = self.parameters
        return par['M']*self.g/float(par['Rs']*(par['T0']+par['L']*z))

    def density(self,p,z):
        par = self.parameters
        return p*par['M']/float(par['R']*self.T(z))

    def dv(self, v, rho):
        par = self.parameters
        # Update rho in parameters such that it can be used to 
        # calculate drag force
        par['rho'] = rho;
        f_d = self.get_drag_force(v)
        return f_d/par['m'] + self.g - (rho/par['rho_b'])*self.g

    def use_parascute(self):
        return 'Tp' in self.parameters.keys()

    def get_parachute(self):
        return self.parameters['Tp']

    def update_parameters(self):
        self.parameters['A'] = self.parameters['Ap']
        self.parameters['C_D'] = self.parameters['C_Dp']

    def T(self,z):
        par = self.parameters
        return par['T0'] - par['L']*z


class Solver:
    def __init__(self, problem):
        self.problem = problem
        init_cond = problem.get_initial_condition()
        self.z = [init_cond['v0']]
        self.p = init_cond['p0']
        self.rho = init_cond['rho0']
        self.v = [init_cond['v0']]
        self.dt = 0.1

    def set_timestep(self, dt):
        self.dt = dt

    def solve(self):

        problem = self.problem

        if problem.use_paraschute():
            counter = 0
            flag = True
            # Run as until parachuter hits the ground
            while(z[-1] > 0):
                if(counter*dt >= problem.get_parachute() and flag):
                    problem.update_parameters()
                    flag = False
                self.iteration()
                counter += 1

        else:
            while(z[-1] > 0):
                self.iteration();

        return self.z,self.v


    def iteration(self):
            dt = self.dt
            problem = self.problem
            z = self.z
            v = self.v
            p = self.p
            rho = self.rho

            # Use forward euler to solve system of equations
            p = p + dt*self.problem.dp(z[-1])
            z.append(z[-1] + dt*problem.dz(v[-1]))
            v.append(v[-1] + dt*problem.dv(v[-1],rho))
            rho = problem.density(p)



if __name__ == '__main__':
    """
	t = linspace(0,par['dt']*len(v),len(v))
	plot(t,z)
	figure()
	plot(t, v)"""
