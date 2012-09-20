

class Solver:
    """
    Class to solve vertical motion with a given problem
    Example of usage:

    c = Solver(problem) # Initiate solver with the Problem class
    c.set_timestep(0.001) # Set the time step of the problem
    v = c.solve(10) # Computed velocity
    """
    def __init__(self, problem):
        self.problem = problem
        self.dt = 0.1
        self.v = [problem.v0]

    def set_timestep(self,dt):
        self.dt = dt

    def solve(self,T):
        """Solves vertical motion problem to timestep T"""
        problem = self.problem
        v = self.v
        N = int(T/float(dt))
        for i in range(N-1):
            # Use Stokes drag model
            if(problem.Re(v[i]) < 1):
                v[i+1] = ((-problem.a_s*dt*0.5 + 1)*v[i] + problem.b_s*dt)/float(1 + problem.a_s*dt*0.5)
            # Use quadratic drag model
            else:
                v[i+1] = (v[i] + dt*problem.b_q)/float(1+dt*problem.a_q*abs(v[i])) 

        return v

