from Problem import Problem
from solver_vertical_motion import Solver
from read_file import read_file
from scitools.std import plot
from numpy import linspace
import sys


def init(v0, T, dt, N, parameters):
	problem = Problem(parameters)
	problem.set_initial_condition(v0);

	c = Solver(problem)
	c.set_timestep(dt)
	v = c.solve(T)

	time = linspace(0,dt*N,N);
	return time, v


if __name__ == '__main__':
	parameters = read_file(sys.argv[1:])
	v0 = 0
	T = 100
	dt = 0.01
	N = int(T/float(dt))

	time, v = init(v0, T, dt, N, parameters)

	plot(time,v)
	raw_input()

