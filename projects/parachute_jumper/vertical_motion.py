from Problem import Problem
from solver_vertical_motion import Solver
from read_file import read_file
from scitools.std import plot
from numpy import linspace
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

time = linspace(0,dt*N,N);

plot(time,v)
raw_input()



