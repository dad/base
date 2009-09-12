import time, os, random, string

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

def looseIntParser(x):
	v = None
	try:
		v = int(x)
	except ValueError:
		if x =='NA' or x =='':
			v = None
		else:
			v = looseFloatParser(x)
	return v

def looseFloatParser(x):
	v = None
	try:
		v = float(x)
	except ValueError:
		if x =='NA' or x =='':
			v = None
		else:
			v = x
	return v

class LineCache:
	
	"""Class for caching lines read by DelimitedLineReader."""
	def __init__(self, instream):
		self.cache = []
		self.instream = instream
		self.cache_size = 100
		self.refill()
	
	def add(self, line):
		self.cache.append(line)
	
	def pop(self):
		# Use standard Python queue pattern
		res = self.cache.pop(0)
		if len(self.cache) == 0:
			self.refill()
		return res

	def refill(self):
		# Refill the cache
		for li in range(self.cache_size):
			line = self.instream.readline()
			# If end of file, bail.
			if not line:
				break
			else:
				self.add(line)
	
	def getLine(self, index):
		while index >= len(self.cache):
			self.refill()
		return self.cache[index]
	
	def isEmpty(self):
		return len(self.cache) == 0
	
	def len(self):
		return len(self.cache)

class ReaderError(Exception):
	"""Exception class for Readers"""


class DelimitedLineReader:
	"""
	Class which parses delimited files line-by-line.
	Typical usage:
	dlr = DelimitedLineReader(file(...,'r'), delim='\t')
	headers = dlr.getHeader()
	while dlr.isValid():
	    flds = dlr.get()
	"""
	handler_dict = {"s":str, "f":looseFloatParser, "d":looseIntParser}

	def __init__(self, in_file, header=True, field_defs=None, sep="\t", strip=True, comment_str="#", custom_handler_dict=None):
		self.infile = in_file
		self.delim = sep
		self.strip = strip
		self.comment_str = comment_str
		# Data
		self.cache = LineCache(self.infile)
		self.cur_flds = None
		self.cur_line = None
		self.n_lines_read = 0
		self.has_header = header
		self.headers = None
		self.field_defs = None
		self.handlers = []

		if self.has_header:
			self.getHeader()
			line = self.next(process=False)
		else:
			line = self.getRawLine()
		if field_defs is None:
			# Attempt to infer the handlers
			while self.isComment():
				line = self.next(process=False)
			assert self.isValid()
			self.field_defs = self.inferHandlers(line)
		else:
			self.field_defs = field_defs
			for h in self.field_defs:
				try:
					self.handlers.append(self.handler_dict[h])
				except KeyError: # unrecognized field definition
					if not custom_handler_dict is None:
						try:
							self.handlers.append(custom_handler_dict[h])
						except KeyError, ke:
							raise ReaderError("No custom handler provided for field-type %s" % h)

	def next(self, process=True):
		if not self.atEnd():
			self.cur_line = self.cache.pop()
			self.n_lines_read += 1
		else:
			raise ReaderError("Attempt to read past end of stream")
		# Read line until we find something
		while self.isComment() and not self.atEnd():
			self.cur_line = self.cache.pop() #self.file.readline()
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
		if line != "":
			flds = line.split(self.delim)
			if apply_handlers:
				assert len(flds) <= len(self.handlers)
				if not self.handlers:
					# Infer handlers are strings
					self.handlers = [str for i in range(len(flds))]
				res = [self.handlers[i](flds[i]) for i in range(len(flds))]
			else:
				res = flds
		else:
			res = None
		return res

	def isValid(self):
		return not self.cur_line is None

	def atEnd(self):
		return self.cache.isEmpty()

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
			#print self.cur_line, self.cache.len()
			line = self.cur_line.strip()
			if len(line)>0:
				res = line[0] == self.comment_str
		return res

	def getHeader(self):
		if not self.headers:
			li = 0
			self.cur_line = self.cache.getLine(li)
			while self.isComment():
				li += 1
				self.cur_line = self.cache.getLine(li)
			if self.isValid() and not self.isComment():
				self.headers = self.process(apply_handlers=False)
			else:
				raise ReaderError("Never encountered a header line")
		return self.headers

	def inferHandlers(self, line):
		# DAD: run through fields until we've seen at least one non-NA for each.
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
	# Normal
	n_lines = 100
	header_list = ["str","float","int","str"]
	fname = "tmp_normal.txt"
	makeFile(fname, header_list, "sfds", n_lines, '\t', 0.0)
	inf = file(fname, 'r')
	# Infer the types
	fp = DelimitedLineReader(inf)
	header = fp.getHeader()
	vals = []
	while not fp.atEnd():
		#print fp.cache.cache[-1],
		v = fp.next()
		#print fp.isValid(), v
		vals.append(v[2])
	assert sum(vals) > n_lines
	inf.close()
	os.remove(fname)
	print "** infer header types"

def test002():
	# Normal
	n_lines = 100
	header_list = ["str","float","int","str"]
	fname = "tmp_normal.txt"
	makeFile(fname, header_list, "sfds", n_lines, '\t', 0.0)
	inf = file(fname, 'r')
	# Infer the types
	fp = DelimitedLineReader(inf)
	while not fp.atEnd():
		fp.next()
	inf.close()
	os.remove(fname)
	print "** read through"

def test003():
	# Normal
	n_lines = 100
	header_list = ["str","float","int","str"]
	fname = "tmp_normal.txt"
	makeFile(fname, header_list, "sfds", n_lines, '\t', 0.0)
	inf = file(fname, 'r')
	# Infer the types
	fp = DelimitedLineReader(inf)
	header = fp.getHeader()
	assert header == header_list
	inf.close()
	os.remove(fname)
	print "** header parsing"

def test004():
	# Normal
	n_lines = 100
	header_list = ["str","float","int","str"]
	fname = "tmp_normal.txt"
	makeFile(fname, None, "sfds", n_lines, '\t', 0.0)
	inf = file(fname, 'r')
	# Infer the types
	fp = DelimitedLineReader(inf)
	header = fp.getHeader()
	vals = []
	while not fp.atEnd():
		#print fp.cache.cache[-1],
		v = fp.next()
		#print fp.isValid(), v
		vals.append(v[2])
	assert sum(vals) > n_lines
	inf.close()
	os.remove(fname)
	print "** infer header types 2"

def test005():
	# Normal
	n_lines = random.randint(10,300)
	header_list = ["str","float","int","str"]
	fname = "tmp_normal.txt"
	makeFile(fname, None, "sfds", n_lines, '\t', 0.0)
	inf = file(fname, 'r')
	# Infer the types
	fp = DelimitedLineReader(inf)
	header = fp.getHeader()
	while not fp.atEnd():
		v = fp.next()
	assert fp.getNumRead() == n_lines
	inf.close()
	os.remove(fname)
	print "** lines read"

def test006():
	# NA's at some frequency -- infer types.
	n_lines = 100
	header_list = ["str","float","int","str"]
	fname = "tmp_normal.txt"
	makeFile(fname, header_list, "sfds", n_lines, '\t', 0.1)
	inf = file(fname, 'r')
	# Infer the types
	fp = DelimitedLineReader(inf)
	header = fp.getHeader()
	vals = []
	while not fp.atEnd():
		#print fp.cache.cache[-1],
		v = fp.next()
		#print v
		if not v[2] is None:
			vals.append(v[2])
	assert sum(vals) > n_lines
	inf.close()
	os.remove(fname)
	print "** infer header types with NA's"

def randString():
	return ''.join(random.sample(string.letters, 10))
def randInt():
	return random.randint(1,10000)
def randFloat():
	return random.random()
def makeFile(fname, headers, fld_types, n_lines, sep, na_prob):
	dtype = {}
	dtype['f'] = randFloat
	dtype['d'] = randInt
	dtype['s'] = randString
	outf = file(fname,'w')
	if headers:
		assert len(headers) == len(fld_types)
		outf.write("%s\n" % sep.join(headers))
	for li in range(n_lines):
		line = ''
		if random.random() < na_prob:
			line = "NA"
		else:
			line = str(dtype[fld_types[0]]())
		if len(fld_types) > 1:
			for f in fld_types[1:]:
				if random.random() < na_prob:
					line = sep.join([line, "NA"])
				else:
					line = sep.join([line, str(dtype[f]())])
		outf.write("%s\n" % line)
	return outf	
				
				

if __name__=="__main__":
	# Make test files
	# Normal
	makeFile("tmp_normal.txt", ["str","float","int","str"], "sfds", 100, '\t', 0.0)
	# No header
	makeFile("tmp_noheader.txt", None, "sfds", 100, '\t', 0.0)
	# Ends with comments
	f = makeFile("tmp_endcomments.txt", ["str","float","int","str"], "sfds", 100, '\t', 0.0)
	f.write("# wakka\n")
	f.write("# foo\n")
	f.close()
	# All comments
	outf = file("tmp_allcomments.txt",'w')
	for i in range(100):
		outf.write("# Comment!\n")
	outf.close()
	# Randomly placed NA's in first several lines
	
	# Tests
	test001()
	test002()
	test003()
	test004()
	test005()
	test006()
	print "** All tests passed **"
