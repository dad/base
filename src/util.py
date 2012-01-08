import time, os, random, string, sys, math
import na

# Float equality.  On my system (WinXP, Python 2.6), smallest distinguishable float difference is 7.45e-9.
def feq(f1,f2,eps=1e-8):
	return abs(f1-f2)<eps

def printTiming(func):
	"""Use as follows. Given a function foo(arg1,arg2), define as:
	@util.printTiming
	def foo(arg1,arg2):
		<some stuff>
	When foo(a,b) is called, printTiming will execute.
	"""
	# Define a wrapper function for func which stores the time before and after calling func.
	def wrapper(*arg):
		# Store the current time, in microseconds
		t1 = time.time()
		# Call the function
		res = func(*arg)
		# Store the current time again, in microseconds
		t2 = time.time()
		print '{0} took {1:0.3f} ms'.format(func.func_name, (t2-t1)*1000.0)
		return res
	# Return the wrapper function
	return wrapper

@printTiming
def test_printTiming(s1):
	time.sleep(0.5)
	print 'The next message should show at least 500ms of wait time:\n\t',

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

	def flush(self):
		for outs in self.streams:
			outs.flush()

def isNA(x):
	print "util.isNA() should be updated to na.isNA()"
	return na.isNA(x)


""" Class for applying formatting automatically.
"""
class FieldFormatter:
	def __init__(self, var, format, transform=lambda x: x):
		self.var = var
		self.format = format
		self.transform = transform

	def __str__(self):
		res = None
		if not na.isNA(self.var):
			try:
				trans_var = self.transform(self.var)
				res = self.format.format(trans_var)
			except ValueError:
				pass
			except TypeError:
				pass
		else:
			res = na.NA
		return res


def dictParser(x, entry_sep=';', key_sep='='):
	entries = x.split(entry_sep)
	return dict([entry.split(key_sep) for entry in entries])

##########################
# Simple testing framework
#
# Example usage:
# if __name__=='__main__':
#	harness = util.TestHarness()
#	harness.add(AnExampleTestCase())
#	harness.add(TestWrapper(anExampleTestFunction, "hello", "hello"))
#	harness.add(TestWrapper(anExampleTestFunction, "hello", "world"))
#	harness.run()
##########################

class AnExampleTestCase(object):
	"""Example test case."""
	def run(self):
		# Because this function returns False, it will declare that it failed.
		# Returning True would suppress the failure message and produce success.
		this_test_passed = False
		self.message = "Goodness, I failed!"
		return this_test_passed

def anExampleTestFunction(arg1,arg2):
	"""An example test function."""
	return arg1==arg2

class TestWrapper(object):
	"""A wrapper class which allows any function returning a bool to be used as a test."""
	def __init__(self, fun, *args):
		self.fun = fun
		self.args = args
		self.__doc__ = fun.__doc__
	
	@property
	def name(self):
		return self.fun.__name__

	def run(self):
		return self.fun(*(self.args))

class TestHarness(object):
	"""A simple framework for registering and running tests."""
	def __init__(self):
		# Maintain a list of tests.
		self.testcases = []
	
	def add(self, test):
		# Add a test to the list.
		self.testcases.append(test)
	
	def run(self, stream=sys.stdout):
		"""Run the registered tests. Defaults to sending information to sys.stdout; pass in stream=file(fname,'w') for file output.
		@return a tuple containing the number of tests, the number that passed, and the total time elapsed.
		"""
		# Run the tests, time them, and keep track of the results.
		n_tests = 0
		n_passed = 0
		total_time = 0.0
		if len(self.testcases) > 0:
			stream.write("Running {} tests\n".format(len(self.testcases)) + "-"*40 + "\n")
			for test in self.testcases:
				# Report on the test
				test_name = test.__class__.__name__
				if test_name == 'TestWrapper':
					test_name = test.name
				stream.write("{} ({}): {}...".format(n_tests+1, test_name, test.__doc__))
				stream.flush()
				# Store the current time again, in seconds
				t1 = time.clock()
				# Run the test and store the result.
				test_result = test.run()
				# Store the current time again, in seconds
				t2 = time.clock()
				timing = '[{:.1f} sec]'.format(t2-t1)
				# Update the total time
				total_time += (t2-t1)
				# One more test done...
				n_tests += 1
				if not test_result:
					# If test failed, print timing and its message, if any.
					line = "failed. {}\n".format(timing)
					if not getattr(test, "message", None) is None:
						line += "\t{}\n".format(test.message)
					stream.write(line)
				else:
					# If test passed, update passed test count and print the timing.
					n_passed += 1
					stream.write("passed {}\n".format(timing))
			# Write out a brief summary of the test results.
			stream.write("-"*40 + "\n")
			stream.write("{0} test{1} completed, {2} passed, {3:.1f} seconds.\n".format(n_tests, "" if n_tests==1 else "s", n_passed, total_time))
		return (n_tests, n_passed, total_time)
		
#########################
# DelimitedLineReader
#########################

def isBlank(s):
	return s.strip() == ''

def isComment(s, comment_char='#'):
	st = s.strip()
	res = st.startswith(comment_char)
	return res

def looseIntParser(x):
	v = None
	try:
		v = int(x)
	except ValueError:
		if not na.isNA(x):
			v = naFloatParser(x)
	return v

def naIntParser(x):
	v = None
	try:
		v = int(x)
	except ValueError, ve:
		if not na.isNA(x):
			raise ve
	return v

def naFloatParser(x):
	v = None
	try:
		v = float(x)
	except ValueError, ve:
		if not na.isNA(x):
			raise ve
	return v

def naStringParser(x):
	"""A parser that respects NA's."""
	v = None
	if not na.isNA(x):
		v = str(x)
	return v

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

def defaultHeader(header_flds):
	return header_flds


class LineCache:
	"""Class for caching lines read by DelimitedLineReader."""
	def __init__(self, instream, comment_str='#'):
		self.cache = []
		self.instream = instream
		self.cache_size = 100
		self.comment_str = comment_str
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
				if not line.strip().startswith(self.comment_str):
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
	res = flds
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
		self.cache = LineCache(self.infile, comment_str=self.comment_str)
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
		## If no headers, return dictionary keyed by column numbers
		headers = self.getHeader()
		res = None
		flds = self.next(True, apply_handlers)
		if not headers is None:
			res = dict(zip(headers, flds))
		else:
			res = dict(zip(range(len(flds)),flds))
			#raise ReaderError, "Attempt to return dictionary of fields, but header is empty"
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
				if len(flds) > len(self.handlers):
					print len(flds), len(self.handlers)
					print flds
					sys.exit()
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
		max_lines = 100
		self.cur_line = self.cache.getLine(li)
		self.handlers = None
		inferred_string = []
		while not handlers_identified and li < max_lines and self.isValid():
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
						if not na.isNA(fld):
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
		if not handlers_identified and li >= max_lines:
			# Went past the allowed number of lines to look ahead; set all unset handlers to strings
			for hi in range(len(self.handlers)):
				if self.handlers[hi] is None:
					self.handlers[hi] = self.handler_dict["s"]
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

class test001:
	"""infer header types"""
	def run(self):
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
		res = sum(vals) > n_lines
		inf.close()
		os.remove(fname)
		return res

class test002:
	"""read through"""
	def run(self):
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
		return True

class test003:
	"""header parsing"""
	def run(self):
		# Normal
		n_lines = 100
		header_list = ["str","float","int","anotherStr"]
		fname = "tmp_normal.txt"
		makeFile(fname, header_list, "sfds", n_lines, '\t', 0.0)
		inf = file(fname, 'r')
		# Infer the types
		fp = DelimitedLineReader(inf)
		header = fp.getHeader()
		res = header == header_list
		inf.close()
		os.remove(fname)
		return res

class test004:
	"""infer header types 2"""
	def run(self):
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
		res = sum(vals) > n_lines
		inf.close()
		os.remove(fname)
		return res

class test005():
	"""lines read"""
	def run(self):
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
		res = fp.getNumRead() == n_lines
		inf.close()
		os.remove(fname)
		return res

class test006:
	"""infer header types with NA's"""
	def run(self):
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
		res = sum(vals) > n_lines
		inf.close()
		os.remove(fname)
		return res

class test007:
	"""infer header types with NA's in every line"""
	def run(self):
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
		res = sum(vals) == 12
		inf.close()
		os.remove(fname)
		return res

class test008:
	"""infer header types with NA's in one full column"""
	def run(self):
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
		return True

class test009:
	"""header processing"""
	def run(self):
		# Header processing
		n_lines = 100
		header_list = ["int","int.1","int","int"]
		fname = "tmp_normal.txt"
		makeFile(fname, header_list, "ffds", n_lines, '\t', 0.1)
		inf = file(fname,'r')
		# Infer the types
		fp = DelimitedLineReader(inf)
		header = fp.getHeader()
		res = True
		res = res and header[0] == 'int'
		res = res and header[1] == 'int.1'
		res = res and header[2] == 'int.2'
		res = res and header[3] == 'int.3'
		inf.close()
		os.remove(fname)
		return res

class test010:
	"""header processing: custom header processor"""
	def run(self):
		# Header processing
		n_lines = 100
		header_list = ["Ratio (H/L)","float","int","This + That"]
		fname = "tmp_normal.txt"
		makeFile(fname, header_list, "ffds", n_lines, '\t', 0.1)
		inf = file(fname,'r')
		# Infer the types
		fp = DelimitedLineReader(inf, header=True, header_name_processor=maxQuantHeader)
		header = fp.getHeader()
		res = True
		res = res and header[0] == 'ratio.hl'
		res = res and header[-1] == 'this.plus.that'
		inf.close()
		os.remove(fname)
		return res

class test011:
	"""adaptive handler updating"""
	def run(self):
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
		return True

class test012:
	"""comment as last line"""
	def run(self):
		# Comment as last line
		n_lines = 10
		header_list = ["str","float","int","str"]
		fname = "tmp_comment_last.txt"
		makeFile(fname, header_list, "sfds", n_lines, '\t', 0.0)
		inf = file(fname, 'a')
		inf.write("# comment\n")
		inf.close()
		inf = file(fname, 'r')
		# Infer the types
		fp = DelimitedLineReader(inf)
		while not fp.atEnd():
			flds = fp.nextDict()
		inf.close()
		os.remove(fname)
		return True

class test013:
	"""no headers but call to nextDict"""
	def run(self):
		# No headers but call to nextDict
		n_lines = 10
		fname = "tmp_no_header.txt"
		field_types = "sfds"
		makeFile(fname, None, field_types, n_lines, '\t', 0.0)
		inf = file(fname, 'r')
		# Infer the types
		fp = DelimitedLineReader(inf, header=False)
		while not fp.atEnd():
			flds = fp.nextDict()
			assert set(flds.keys()) == set(range(len(field_types)))
		inf.close()
		os.remove(fname)
		return True

class test014:
	"""comments as first lines"""
	def run(self):
		n_lines = 10
		header_list = ["str","float","int","str"]
		fname = "tmp_comment_first.txt"
		makeFile(fname, header_list, "sfds", n_lines, '\t', 0.0, first_lines="# first line comment\n# second line comment")
		inf = file(fname, 'r')
		# Infer the types
		fp = DelimitedLineReader(inf)
		while not fp.atEnd():
			flds = fp.nextDict()
		inf.close()
		os.remove(fname)
		return True


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

# Read

class LightDataFrame:
	def __init__(self, headers, data):
		self._headers = basicHeaderFixer(headers)
		# Data must be a list of lists, each list is a column.
		self._data = data
		self._header_lookup = dict([(h, self._headers.index(h)) for h in self._headers])

	def getRow(self, row_index):
		return [d[row_index] for d in self._data]

	def getRowDict(self, row_index):
		return dict(zip(self._headers,self.getRow(row_index)))

	def getCol(self, col_index):
		"""Index by column, given by column_key"""
		col_index = self._header_lookup[column_key]
		return [d[col_index] for d in self._data]

	def __getitem__(self, column_key):
		"""Index by column, given by column_key"""
		col_index = self._header_lookup[column_key]
		return self._data[col_index]

	def getNumRows(self):
		return len(self._data[0])

	def getNumCols(self):
		return len(self._data)

	def getHeaders(self):
		return self._headers

	def __str__(self):
		return "** not implemented **"

	nrow = property(getNumRows)
	ncol = property(getNumCols)
	headers = property(getHeaders)


def readTable(fname, header=True, sep='\t', header_name_processor=defaultHeader):
	inf = file(fname,'r')
	dlr = DelimitedLineReader(inf, header=header, sep=sep, header_name_processor=header_name_processor)
	header_flds = dlr.getHeader()
	data_dict = {}
	initialized = False
	while not dlr.atEnd():
		flds = dlr.nextDict()
		if not initialized:
			for h in header_flds:
				data_dict[h] = [flds[h]]
			initialized = True
		else:
			for h in header_flds:
				data_dict[h].append(flds[h])
	inf.close()
	return LightDataFrame(header_flds, [data_dict[h] for h in header_flds])

if __name__=="__main__":
	# Utility function tests
	test_printTiming('test')
	harness = TestHarness()
	# DelimitedLineReader tests
	harness.add(test001())
	harness.add(test002())
	harness.add(test003())
	harness.add(test004())
	harness.add(test005())
	harness.add(test006())
	harness.add(test007())
	harness.add(test008())
	harness.add(test009())
	harness.add(test010())
	harness.add(test011())
	harness.add(test012())
	harness.add(test013())
	harness.add(test014())
	harness.add(TestWrapper(anExampleTestFunction, "hello", "hello"))
	harness.run()
