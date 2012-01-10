import scipy as sp
import numpy as np

# Common transformations 
def noop(x):
	return x

def podds(p):
	return p/(1.0-p)

def oddsp(o):
	return o/(1.0+o)

def plogodds(p):
	return sp.log(podds(p))

def logoddsp(lo):
	return oddsp(sp.exp(lo))

transform_inverses = {sp.exp:sp.log, noop:noop, podds:oddsp, plogodds:logoddsp}
for (k,v) in transform_inverses.items():
	transform_inverses[v] = k

class Parameter(object):
	"""A parameter in an optimization scheme. Given a transform function and an inverse, can have the transformed value be optimized."""
	def __init__(self, name, value=None, transform=None, deriv=None, hess=None):
		self._name = name
		self._value = value
		self._transform = (transform,)[0]
		# DAD: catch error here if there's no known inverse? For now, allow error to propagate normally.
		if len((transform,)) == 1:
			self._inv_transform = transform_inverses[self._transform]
		else:
			self._inv_transform = (transform,)[1]
		self._deriv_fxn = deriv
		self._hess_fxn = hess
		self._optderiv = None
		self._opthess = None
	
	@property
	def value(self):
		"""Get the value of the parameter"""
		return self._value

	@value.setter
	def value(self, v):
		"""Set the value of the parameter"""
		self._value = v
	
	@property
	def name(self):
		"""Get the parameter's name"""
		return self._name
	
	@property
	def transform(self):
		"""Get the transfom function"""
		return self._transform
	
	@property
	def optvalue(self):
		"""Get the optimizer-prepared value"""
		return self._transform(self._value)
	
	@optvalue.setter
	def optvalue(self, ov):
		"""Store the untransformed value given the optimizer-prepared value"""
		self._value = self._inv_transform(ov)
	
	@property
	def optderivative(self):
		return self._optderiv #self._deriv_fxn(self.optvalue, data)
	
	@optderivative.setter
	def optderivative(self, od):
		"""Store the specified derivative of the optimizer-ready value"""
		self._optderiv = od
	
	@property
	def opthessian(self):
		return self._opthess #self._hess_fxn(self.optvalue)
	
	def __str__(self):
		return "{me.name} {me.value} ({me.optvalue}), f'={me.optderivative}".format(me=self)
	
class ParameterTable(object):
	"""A collection of parameters"""
	def __init__(self, parameter_list=[]):
		self._parameter_dict = dict([(p.name,p) for p in parameter_list])
		self._parameter_names = sorted(self._parameter_dict.keys())
	
	def add(self, parameter):
		self._parameter_dict[parameter.name] = parameter
		self._parameter_names = sorted(self._parameter_dict.keys())
	
	def get(self, param_name, default=None):
		return self._parameter_dict.get(param_name,default)
	
	def __getitem__(self, param_name):
		return self._parameter_dict[param_name].value
	
	def toOptimizer(self):
		#opt_params = [self._parameter_dict[id].optvalue for id in self._parameter_names]
		return self.optvalues
	
	def fromOptimizer(self, param_list):
		assert len(param_list) == len(self._parameter_names)
		for (id, v) in zip(self._parameter_names, param_list):
			self._parameter_dict[id].optvalue = v
	
	@property
	def optderivatives(self):
		return = sp.r_[[self._parameter_dict[id].optderivative for id in self.names]]
	
	@property
	def opthess(self):
		return = sp.r_[[self._parameter_dict[id].hessian for id in self.names]]
	
	@property
	def names(self):
		return self._parameter_names
	
	@property
	def keys(self):
		return self._parameter_names
	
	@property
	def values(self):
		return [self._parameter_dict[k].value for k in self.names]

	@property
	def optvalues(self):
		return [self._parameter_dict[k].optvalue for k in self.names]
	
	def __str__(self):
		"""Write out the parameter name, value, and optimizer value in tabular form"""
		s = ''
		for k in self._parameter_names:
			p = self.get(k)
			s += str(p) + '\n'
		return s[:-1]

class LogLikelihoodManager(object):
	"""Generic base object for LogLikehood wrapping. Contains parameters, and has methods for fitting and calculating log-likelihood."""
	def __init__(self):
		self._parameters = ParameterTable()
	
	def getParameter(self, id, default=None):
		return self._parameters.get(id,default)
	
	def __getitem__(self, id):
		return self.getParameter(id).value

	def updateParameters(self, opt_param_list):
		"""Pull parameters from the optimizer-ready list"""
		self._parameters.fromOptimizer(opt_param_list)
			
	@property
	def parameters(self):
		"""Get the parameters"""
		return self._parameters
		
	def logLikelihood(self, opt_param_list):
		"""Compute the log-likelihood given a list of parameters"""
		# First, pull parameters from optimizer.
		self.updateParameters(opt_param_list)
		# Then, compute log likelihood, however that is done.
		raise Exception, "The logLikelihood() method must be overridden!"
	
	def derivative(self, opt_param_list):
		raise Exception, "The derivative() method must be overridden!"
	
	def hessian(self, opt_param_list):
		raise Exception, "The hessian() method must be overridden!"
	
	def fit(self):
		"""Fit the specified model"""
		starting_params = self.parameters.toOptimizer()
		fitted_params = sp.optimize.fmin(self.logLikelihood, starting_params, disp=False)
		self.updateParameters(fitted_params)
		return self.parameters
