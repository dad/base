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
	def __init__(self, stream_list=[]):
		if not isinstance(stream_list, list):
			stream_list = [stream_list]
		self.streams = stream_list[:]

	def write(self, line):
		for outs in self.streams:
			outs.write(line)

	def addStream(self, stream):
		self.streams.append(stream)

	def removeStream(self, stream):
		self.streams.remove(stream)

def isBlank(s):
	return s.strip() == ''

def isComment(s, comment_char='#'):
	st = s.strip()
	res = False
	if len(st) > 0:
		res = st[0] == comment_char
	return res

def isNA(x):
	return x is None or x == 'NA' or x == ''

def strNA(x):
	if isNA(x):
		return "NA"
	else:
		return str(x)

def formatNA(x, format, sep=None):
	if isinstance(x, list):
		flds = format.split(sep)
		if sep is None:
			sep = '\t'
		assert len(flds) == len(x)
		res = []
		for i in range(len(x)):
			if isNA(x[i]):
				res.append("NA")
			else:
				res.append(flds[i] % x[i])
		return sep.join(res)
	else:
		if isNA(x):
			return "NA"
		else:
			return format % x

def looseIntParser(x):
	v = None
	try:
		v = int(x)
	except ValueError:
		if isNA(x):
			v = None
		else:
			v = naFloatParser(x)
	return v

def naIntParser(x):
	v = None
	try:
		v = int(x)
	except ValueError, ve:
		if isNA(x):
			v = None
		else:
			raise ve
	return v

def naFloatParser(x):
	v = None
	try:
		v = float(x)
	except ValueError, ve:
		if isNA(x):
			v = None
		else:
			raise ve
	return v

def naStringParser(x):
	"""A parser that respects NA's."""
	v = None
	if not isNA(x):
		v = str(x)
	return v

def dictParser(x, entry_sep=';', key_sep='='):
	entries = x.split(entry_sep)
	return dict([entry.split(key_sep) for entry in entries])

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
		eof = False
		for li in range(self.cache_size):
			line = self.instream.readline()
			# If end of file, bail.
			if not line:
				eof = True
				break
			else:
				self.add(line)
		return eof

	def getLine(self, index):
		eof = False
		while not eof and (index >= len(self.cache)):
			eof = self.refill()
		if index >= len(self.cache):
			raise ReaderEOFError("Attempt to read past end of stream (%d, %d)" % (index, len(self.cache)))
		return self.cache[index]

	def isEmpty(self):
		return len(self.cache) == 0

	def len(self):
		return len(self.cache)

class ReaderError(Exception):
	"""Exception class for Readers"""

class ReaderEOFError(ReaderError):
	"""Error when Reader goes past end of file"""

def basicHeaderFixer(flds):
	# Make headers unique
	if len(flds) > len(set(flds)):
		new_flds = ['']*len(flds)
		for fi in range(len(flds)):
			f = flds[fi]
			f_count = flds.count(f)
			if f_count == 1 or fi == 0:
				new_flds[fi] = f
			else:  # more than one instance
				i = 1
				new_fld_name = '%s.%d' % (f,i)
				# only look at field names up to this one.
				# if this is the first instance of a repeated field name,
				# it will not be changed.
				while new_fld_name in set(new_flds + flds[0:fi-1]):
					i += 1
					new_fld_name = '%s.%d' % (f,i)
				new_flds[fi] = new_fld_name
		# Ensure we've actually fixed the problem
		assert len(new_flds) == len(set(new_flds))
		res = new_flds
	else:
		res = flds
	return res

class DelimitedLineReader:
	"""
	Class which parses delimited files line-by-line.
	Typical usage:
	dlr = DelimitedLineReader(file(...,'r'), delim='\t')
	headers = dlr.getHeader()
	while not dlr.atEnd():
	    flds = dlr.next()
	"""
	handler_dict = {"s":str, "f":naFloatParser, "d":naIntParser}

	def __init__(self, in_file, header=True, field_defs=None, sep="\t", strip=False, comment_str="#", custom_handler_dict=None, header_name_processor=basicHeaderFixer):
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
		self.field_defs = field_defs
		self.handlers = []
		self.header_name_processor = header_name_processor

		if self.has_header:
			# Position reader right after header
			self.getHeader(move_to_data=True)

		if field_defs is None:
			# Attempt to infer the handlers
			self.field_defs = self.inferHandlers()
		else:
			for h in self.field_defs:
				try:
					self.handlers.append(self.handler_dict[h])
				except KeyError: # unrecognized field definition
					if not custom_handler_dict is None:
						try:
							self.handlers.append(custom_handler_dict[h])
						except KeyError, ke:
							raise ReaderError("No custom handler provided for field-type %s" % h)

	def next(self, process=True, apply_handlers=True):
		if not self.atEnd():
			self.cur_line = self.cache.pop()
			self.n_lines_read += 1
		else:
			raise ReaderEOFError("Attempt to read past end of stream")
		# Read line until we find something
		while self.isComment() and not self.atEnd():
			self.cur_line = self.cache.pop() #self.file.readline()
			self.n_lines_read += 1
		res = None
		if self.isValid():
			if process:
				self.cur_flds = self.process(apply_handlers)
				res = self.cur_flds
			else:
				res = self.cur_line
		#print "$%d-%s^" % (self.getNumRead(), self.getRawLine())
		return res

	def nextDict(self, apply_handlers=True):
		## Return fields as a dictionary, keyed by the header names
		headers = self.getHeader()
		res = None
		if not headers is None:
			flds = self.next(True, apply_handlers)
			res = dict(zip(headers, flds))
		else:
			raise ReaderError, "Attempt to return dictionary of fields, but header is empty"
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
			flds[-1] = flds[-1].strip() # Get rid of \n
			if apply_handlers:
				assert len(flds) <= len(self.handlers)
				if not self.handlers:
					# Infer handlers are strings
					self.handlers = [str for i in range(len(flds))]
				done_processing = False
				while not done_processing:
					try:
						# Apply the new handlers to get the data.
						res = [self.handlers[hi](flds[hi]) for hi in range(len(flds))]
						done_processing = True
					except ValueError, ve:
						#print "updating handler %d" % hi
						# Adaptively update handlers if there was a value error.
						handler_key = self.inferHandlerKey(flds[hi])
						self.handlers[hi] = self.handler_dict[handler_key]
						# Reapply the new handlers to get the data.
						#res = [self.handlers[i](flds[i]) for i in range(len(flds))]

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

	def getHeader(self, move_to_data=True):
		if self.has_header:
			if not self.headers:
				li = 0
				self.cur_line = self.cache.getLine(li)
				while self.isComment():
					li += 1
					self.cur_line = self.cache.getLine(li)
				if self.isValid() and not self.isComment():
					self.headers = self.process(apply_handlers=False)
					if move_to_data:
						# Position reader at line following header
						for pop_li in range(li+1):
							self.cache.pop()
							self.n_lines_read += 1
				else:
					raise ReaderError("Never encountered a header line")
			# Now process the headers -- e.g., de-dupe column names.
			self.headers = self.header_name_processor(self.headers)
			return self.headers
		else:
			return None # No header

	def inferHandlerKey(self, fld):
		# int -> float -> string
		found_key = None
		for handler_key in "dfs":  # DAD: should include custom handlers up front
			handler = self.handler_dict[handler_key]
			try:
				res = handler(fld)
				# If we get here without an exception, save the key
				found_key = handler_key
				break
			except ValueError:
				continue
		return found_key


	def inferHandlers(self):
		# DAD: run through fields until we've seen at least one non-NA for each.
		handlers_identified = False
		li = 0
		self.cur_line = self.cache.getLine(li)
		self.handlers = None
		inferred_string = []
		while not handlers_identified and self.isValid():
			if not self.isComment():
				# Not a comment line -- parse it.
				if self.strip:
					self.cur_line = self.cur_line.strip()
				flds = self.cur_line.split(self.delim)
				flds[-1] = flds[-1].strip() # Get rid of \n
				# Initialize empty handler list if we haven't done so already
				if self.handlers is None:
					self.handlers = [None]*len(flds)
					inferred_string = ['X']*len(flds)
				if len(flds) != len(self.handlers):
					print flds
				assert len(flds) == len(self.handlers), "Number of fields %d not equal to number of handlers %d" % (len(flds), len(self.handlers))
				for hi in range(len(self.handlers)):
					fld = flds[hi]
					if self.handlers[hi] is None:
						if not isNA(fld):
							handler_key = self.inferHandlerKey(fld)
							inferred_string[hi] = handler_key
							self.handlers[hi] = self.handler_dict[handler_key]
					else: # handler has already been found; just confirm, and upgrade if necessary
						try:
							val = self.handlers[hi](fld)
						except ValueError:
							#print "upgrading handler", inferred_string[hi],
							handler_key = self.inferHandlerKey(fld)
							inferred_string[hi] = handler_key
							self.handlers[hi] = self.handler_dict[handler_key]
							#print "to", handler_key

				# We're finished when all handlers are not None.
				handlers_identified = len([h for h in self.handlers if h is None]) == 0
			if not handlers_identified:
				#print "cache:", self.cache.cache
				li += 1
				try:
					self.cur_line = self.cache.getLine(li)
				except ReaderEOFError:
					# We've reached the end of the file with an inconclusive result -- some fields
					# still can't have types inferred.
					# Just assume everything's a string.
					for hi in range(len(self.handlers)):
						if self.handlers[hi] is None:
							self.handlers[hi] = self.handler_dict["s"]
					handlers_identified = True
		#print inferred_string
		inferred_string = ''.join(inferred_string)
		return inferred_string

	def setHandlerType(self, handler_index, type_string):
		try:
			self.handlers[handler_index] = self.handler_dict[type_string]
		except KeyError, ke:
			raise ReaderError, "Unknown handler type %s" % type_string
		except IndexError:
			raise ReaderError, "Bad handler index %d" % handler_index


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
	print "** 001 infer header types"

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
	print "** 002 read through"

def test003():
	# Normal
	n_lines = 100
	header_list = ["str","float","int","anotherStr"]
	fname = "tmp_normal.txt"
	makeFile(fname, header_list, "sfds", n_lines, '\t', 0.0)
	inf = file(fname, 'r')
	# Infer the types
	fp = DelimitedLineReader(inf)
	header = fp.getHeader()
	assert header == header_list
	inf.close()
	os.remove(fname)
	print "** 003 header parsing"

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
	print "** 004 infer header types 2"

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
	print "** 005 lines read"

def test006():
	# NA's at some frequency -- infer types.
	n_lines = 100
	header_list = ["str","float","int","str"]
	fname = "tmp_normal.txt"
	# Put an NA in the first line -- test the ability to look past to infer types
	first_lines = [randString() + "\t1.00\tNA\t" + randString() + "\n"]
	makeFile(fname, None, "sfds", n_lines, '\t', 0.1, first_lines=first_lines, last_lines=None)
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
	print "** 006 infer header types with NA's"

def test007():
	# NA's in every line -- infer types.
	n_lines = 0
	header_list = ["str","float","int","str"]
	fname = "tmp_normal.txt"
	# Put an NA in the first line -- test the ability to look past to infer types
	first_lines = [randString() + "\t59.3\tNA\t" + randString() + "\n",
		randString() + "\tNA\t12\t" + randString() + "\n"]
	makeFile(fname, None, "sfds", n_lines, '\t', 0.1, first_lines=first_lines, last_lines=None)
	inf = file(fname,'r')
	# Infer the types
	fp = DelimitedLineReader(inf, header=False)
	header = fp.getHeader()
	vals = []
	while not fp.atEnd():
		#print fp.cache.cache[-1],
		v = fp.next()
		#print v
		if not v[2] is None:
			vals.append(v[2])
	assert sum(vals) == 12
	inf.close()
	os.remove(fname)
	print "** 007 infer header types with NA's in every line"

def test008():
	# NA's in every line -- infer types.
	n_lines = 0
	header_list = ["str","float","int","str"]
	fname = "tmp_normal.txt"
	# Put an NA in the first line -- test the ability to look past to infer types
	first_lines = [randString() + "\t59.3\tNA\t" + randString() + "\n"]*100
	makeFile(fname, None, "sfds", n_lines, '\t', 0.1, first_lines=first_lines, last_lines=None)
	inf = file(fname,'r')
	# Infer the types
	fp = DelimitedLineReader(inf, header=False)
	header = fp.getHeader()
	inf.close()
	os.remove(fname)
	print "** 008 infer header types with NA's in one full column"

def maxQuantHeader(header_flds):
	# Header line
	header_line = '\t'.join(header_flds)
	header_line = header_line.lower()
	header_line = header_line.replace("h/l", "hl")
	header_line = header_line.replace(" ", ".")
	header_line = header_line.replace("(", ".")
	header_line = header_line.replace(".[%]", "")
	header_line = header_line.replace("[", ".")
	header_line = header_line.replace("]", "")
	header_line = header_line.replace(")", "")
	header_line = header_line.replace("/", '.')
	header_line = header_line.replace("+", "plus")
	header_line = header_line.replace(".\t", '\t')
	header_line = header_line.replace(".\n", '\n')
	header_line = header_line.replace("..", '.')
	header_flds = header_line.strip().split('\t')
	return header_flds

def test009():
	# Header processing
	n_lines = 100
	header_list = ["int","int.1","int","int"]
	fname = "tmp_normal.txt"
	makeFile(fname, header_list, "ffds", n_lines, '\t', 0.1)
	inf = file(fname,'r')
	# Infer the types
	fp = DelimitedLineReader(inf)
	header = fp.getHeader()
	assert header[0] == 'int'
	assert header[1] == 'int.1'
	assert header[2] == 'int.2'
	assert header[3] == 'int.3'
	inf.close()
	os.remove(fname)
	print "** 009 header processing"

def test010():
	# Header processing
	n_lines = 100
	header_list = ["Ratio (H/L)","float","int","This + That"]
	fname = "tmp_normal.txt"
	makeFile(fname, header_list, "ffds", n_lines, '\t', 0.1)
	inf = file(fname,'r')
	# Infer the types
	fp = DelimitedLineReader(inf, header=True, header_name_processor=maxQuantHeader)
	header = fp.getHeader()
	assert header[0] == 'ratio.hl'
	assert header[-1] == 'this.plus.that'
	inf.close()
	os.remove(fname)
	print "** 010 header processing: custom header processor"

def test011():
	# Adaptive field redefinition
	n_lines = 100
	header_list = ["float","float","int","str"]
	fname = "tmp_normal.txt"
	last_lines = ["0.2\t0.3\tnot.an.int\twakka\n"]
	makeFile(fname, header_list, "ffds", n_lines, '\t', 0.01, last_lines=last_lines)
	inf = file(fname,'r')
	# Infer the types
	fp = DelimitedLineReader(inf, header=True, header_name_processor=maxQuantHeader)
	header = fp.getHeader()
	while not fp.atEnd():
		flds = fp.next()
	inf.close()
	os.remove(fname)
	print "** 011 adaptive handler updating"


def randString():
	return ''.join(random.sample(string.letters, 10))
def randInt():
	return random.randint(1,10000)
def randFloat():
	return random.random()
def makeFile(fname, headers, fld_types, n_lines, sep, na_prob, first_lines=None, last_lines=None):
	dtype = {}
	dtype['f'] = randFloat
	dtype['d'] = randInt
	dtype['s'] = randString
	outf = file(fname,'w')
	if first_lines:
		for line in first_lines:
			outf.write(line)
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
	if last_lines:
		for line in last_lines:
			outf.write(line)
	outf.close()
	return fname



if __name__=="__main__":
	# Tests
	test001()
	test002()
	test003()
	test004()
	test005()
	test006()
	test007()
	test008()
	test009()
	test010()
	test011()
	print "** All tests passed **"
