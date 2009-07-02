import time

def printTiming(func):
	def wrapper(*arg):
		t1 = time.time()
		res = func(*arg)
		t2 = time.time()
		print '%s took %0.3f ms' % (func.func_name, (t2-t1)*1000.0)
		return res
	return wrapper

class OutStreams:
	streams = []
	def __init__(self, stream_list):
		self.streams = stream_list

	def write(self, line):
		for outs in self.streams:
			outs.write(line)
	
