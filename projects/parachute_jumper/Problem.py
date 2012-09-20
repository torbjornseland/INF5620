from numpy import pi

class Problem:
	"""This class takes in a dictionary with all the parameters.
	It gives back the drag force after checking the re.
	""" 
	def __init__(self, param):
		self.param = param
		self.g = -9.81
		self.a_s = (3*pi*param['d']*param['my'])/(param['V']*param['rho_b'])
		self.a_q =  0.5*param['C_D']*(param['rho']*param['A'])/(param['rho_b']*param['V'])
		self.b = self.g*(1-param['rho']/float(param['rho_b']))

	def set_initial_condition(self,v0):
		self.param['v0']= v0

	def get_initial_condition(self):
		return self.param['v0']

	def re(self,v):
		param = self.param
		self.Re = (param['d']*abs(v)*param['rho'])/param['my']
		return self.Re
    
	def get_gravity_force(self):
		return self.param['m']*self.g

	def get_buoyancy_force(self):
	        param = self.param
	        return -param['V']*param['rho']*self.g

	def get_drag_force_stokes(self,v):
		self.v = v
		self.stokes = -3*pi*self.param['d']*self.param['my']*self.v
		return self.stokes

	def get_drag_force_quadratic(self,v):
		self.v = v
		self.quadratic =-.5*self.param['C_D']*self.param['rho']*self.param['A']*abs(self.v)*self.v
		return self.quadratic

	def get_drag_force(self,v):
		self.v = v
		if self.re(self.v) < 1:
			return self.get_drag_force_stokes(self.v)
		else:
			return self.get_drag_force_quadratic(self.v)
		

if __name__ == '__main__':
	param = {'rho':1.2041, 'rho_b':1000, 'v0':0, 'my':18.27*10**(-6), 'm':75, 'd':0.5, 'C_D':1.15, 'V':0.075, 'A':0.3*1.8}
	test = problem(param)
	print test.re(1)
	print test.force()
	print test.get_drag_force_stokes(2)
	print test.get_drag_force_quadratic(2)

