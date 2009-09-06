import time, os

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
	def __init__(self, stream_list=[]):
		if not isinstance(stream_list, list):
			stream_list = [stream_list]
		self.streams = stream_list

	def write(self, line):
		for outs in self.streams:
			outs.write(line)

	def addStream(self, stream):
		self.streams.append(stream)

	def removeStream(self, stream):
		self.streams.remove(stream)

class ReaderError(Exception):
	"""Exception class for Readers"""

class DelimitedLineReader:
	handler_dict = {"s":str, "f":float, "d":int}
	
	def __init__(self, in_file, header=True, field_defs=None, delim="\t", strip=True, comment_str="#", custom_handler_dict=None):
		self.infile = in_file
		self.delim = delim
		self.strip = strip
		self.comment_str = comment_str
		self.cur_flds = None
		self.cur_line = self.infile.readline()
		self.n_lines_read = 1
		self.has_header = header
		self.headers = None
		self.handlers = []
		
		if field_defs is None:
			# Attempt to infer the handlers
			if self.has_header:
				self.getHeader()
				line = self.next(process=False)
			else:
				line = self.getRawLine()
			while self.isComment():
				line = self.next(process=False)
			assert self.isValid()
			self.inferHandlers(line)
		else:
			for h in field_defs:
				try:
					self.handlers.append(self.handler_dict[h])
				except KeyError: # unrecognized field definition
					if not custom_handler_dict is None:
						try:
							self.handlers.append(custom_handler_dict[h])
						except KeyError, ke:
							raise FileReaderException("No custom handler provided for field-type %s" % h)

	def next(self, process=True):
		self.cur_line = self.infile.readline()
		if self.isValid():
			self.n_lines_read += 1
		# Read line until we find something 
		while self.isComment():
			self.cur_line = self.file.readline()
			self.n_lines_read += 1
		res = None
		if self.isValid():
			if process:
				self.cur_flds = self.process()
				res = self.cur_flds
			else:
				res = self.cur_line
		#print "$%d-%s^" % (self.getNumRead(), self.getRawLine())
		return res

	def process(self, apply_handlers=True):
		res = None
		if self.isValid():
			res = self.processLine(self.cur_line, apply_handlers)
		return res
			
	def processLine(self, line, apply_handlers=True):
		if self.strip:
			line = line.strip()
		flds = line.split(self.delim)
		if apply_handlers:
			if not self.handlers:
				# Infer handlers are strings
				self.handlers = [str for i in range(len(flds))]
			res = [self.handlers[i](flds[i]) for i in range(len(flds))]
		else:
			res = flds
		return res

	def isValid(self):
		return self.cur_line
	
	def getRawLine(self):
		return self.cur_line

	def get(self):
		return self.cur_flds

	def getField(self, field):
		if self.headers:
			res = self.cur_flds[self.headers.index(field)]
		else:
			res = self.cur_flds[field]
		return res

	def getNumRead(self):
		return self.n_lines_read

	def isComment(self):
		res = False
		if self.isValid():
			line = self.cur_line.strip()
			if len(line)>0:
				res = line[0] == self.comment_str
		return res

	def getHeader(self):
		if not self.headers:
			while self.isComment():
				s = self.next(process=False)
			if self.isValid() and not self.isComment():
				self.headers = self.process(apply_handlers=False)
			else:
				raise ReaderError("Never encountered a header line")
		return self.headers

	def inferHandlers(self, line):
		if self.strip:
			line = line.strip()
		flds = line.split(self.delim)
		self.handlers = [None]*len(flds)
		inferred_string = ""
		for fi in range(len(flds)):
			fld = flds[fi]
			# int -> float -> string
			for handler_key in "dfs":
				handler = self.handler_dict[handler_key]
				try:
					res = handler(fld)
					# If we get here without an exception, add the handler
					self.handlers[fi] = handler
					inferred_string += handler_key
					break
				except ValueError:
					continue
		#print inferred_string
		return inferred_string

def test001():
	inf = file(os.path.expanduser("~/research/data/scerevisiae/scer-trna-anticodons.txt"),'r')
	# Infer the types
	fp = DelimitedLineReader(inf)
	vals = []
	v = fp.process()
	while fp.isValid():
		vals.append(v[2])
		v = fp.next()
	assert sum(vals) == 273
	inf.close()
	print "** infer header types"
		
def test002():
	inf = file(os.path.expanduser("~/research/data/scerevisiae/scer-trna-anticodons.txt"),'r')
	fp = DelimitedLineReader(inf, "ssdss")
	fp.next()
	while fp.isValid():
		flds = fp.get()
		fp.next()
	inf.close()
	print "** read through"

def test003():
	inf = file(os.path.expanduser("~/research/data/scerevisiae/scer-trna-anticodons.txt"),'r')
	fp = DelimitedLineReader(inf, "ssdss")
	header = fp.getHeader()
	assert '-'.join(header) == 'aa-dna.anticodon-count-anticodon-codons'
	inf.close()
	print "** header parsing"

def test004():
	inf = file(os.path.expanduser("~/research/data/scerevisiae/scer-trna-anticodons.txt"),'r')
	inf.readline()
	# Infer the types -- no header
	fp = DelimitedLineReader(inf, header=False)
	vals = []
	v = fp.process()
	while fp.isValid():
		vals.append(v[2])
		v = fp.next()
	#print sum(vals)
	assert sum(vals) == 273
	inf.close()
	print "** infer header types 2"

if __name__=="__main__":
	# Tests
	test001()
	test002()
	test003()
	test004()
	print "** All tests passed **"
