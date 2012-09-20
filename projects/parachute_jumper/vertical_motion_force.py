from Problem import Problem
from solver_vertical_motion import Solver
from read_file import read_file
from scitools.std import plot
from numpy import linspace, array, zeros
import sys

v0 = 0
T = 100
dt = 0.01
N = int(T/float(dt))

parameters = read_file(sys.argv[1:])
problem = Problem(parameters)
problem.set_initial_condition(v0);

c = Solver(problem)
c.set_timestep(dt)
v = c.solve(T)

v = array(v)
F_d = zeros(N)
for i in range(N):
        F_d[i] = problem.get_drag_force(v[i])

F_g = zeros(N) + problem.get_gravity_force()
F_b = zeros(N) + problem.get_buoyancy_force()

time = linspace(0,dt*N,N);

plot(time,F_d,time,F_g,time,F_b, legend=['drag','gravity','buoyancy'])
raw_input()

