from numpy import pi

class problem:
	def __init__(self, param):
		self.param = param
		self.g = -9.81
		self.a_s = (3*pi*param['d']*param['my'])/(param['V']*param['rho_b'])
		self.a_q =  0.5*param['C_D']*(param['rho']*param['A'])/(param['rho_b']*param['V'])
		self.b = self.g*(1-param['rho']/float(param['rho_b']))

	def set_initial_condition(self,v0):
		self.param['v0']= v0

	def re(self,v):
		self.v = v
		self.Re = (param['d']*abs(self.v)*param['rho'])/param['my']
		return self.Re

	def force(self):
		self.g = -9.81
		self.gravity = self.param['m']*self.g
		self.buoyancy = -self.param['rho']*self.param['V']*self.g
		return self.gravity+self.buoyancy

	def get_drag_force_stokes(self,v):
		self.v = v
		self.stokes = -3*pi*self.param['d']*self.param['my']*self.v
		return self.stokes+ self.force()

	def get_drag_force_quadratic(self,v):
		self.v = v
		self.quadratic =-.5*self.param['C_D']*self.param['rho']*self.param['A']*abs(self.v)*self.v
		return self.quadratic+ self.force()
		
if __name__ == '__main__':
	param = {'rho':1.2041, 'rho_b':1000, 'v0':0, 'my':18.27*10**(-6), 'm':75, 'd':0.5, 'C_D':1.15, 'V':0.075, 'A':0.3*1.8}
	test = problem(param)
	print test.re(1)
	print test.force()
	print test.get_drag_force_stokes(2)
	print test.get_drag_force_quadratic(2)

