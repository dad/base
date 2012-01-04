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
	def __init__(self, name, value=None, transform=None):
		self._name = name
		self._value = value
		self._transform = transform
		# DAD: catch error here?
		self._inv_transform = transform_inverses[self._transform]
	
	@property
	def value(self):
		return self._value
	
	@property
	def name(self):
		return self._name
	
	@property
	def transform(self):
		return self._transform
	
	@property
	def optvalue(self):
		return self._transform(self._value)
	
	@optvalue.setter
	def optvalue(self, ov):
		"""Retrieve the untransformed value from the optimizer-prepared value"""
		self._value = self._inv_transform(ov)
	
class ParameterTable(object):
	def __init__(self, parameter_list=[]):
		self._parameter_dict = dict([(p.name,p) for p in parameter_list])
		self._parameter_names = sorted(self._parameter_dict.keys())
	
	def add(self, parameter):
		self._parameter_dict[parameter.name] = parameter
		self._parameter_names = sorted(self._parameter_dict.keys())
	
	def get(self, param_name, default=None):
		return self._parameter_dict.get(param_name,default)
	
	def __getitem__(self, param_name):
		return self._parameter_dict[param_name]
	
	def toOptimizer(self):
		opt_params = [self._parameter_dict[id].optvalue for id in self._parameter_names]
		return opt_params
	
	def fromOptimizer(self, param_list):
		assert len(param_list) == len(self._parameter_names)
		for (id, v) in zip(self._parameter_names, param_list):
			self._parameter_dict[id].optvalue = v
	
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
