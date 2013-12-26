import time, os, random, string, sys, math, traceback, unittest
import util, na

@util.printTiming
def test_printTiming():
	time.sleep(0.5)
	print 'The next message should show at least 500ms of wait time:\n\t',

class test001(unittest.TestCase):
	"""infer header types"""
	def test_run(self):
		# Normal
		n_lines = 100
		header_list = ["str","float","int","str"]
		fname = "tmp_normal.txt"
		makeFile(fname, header_list, "sfds", n_lines, '\t', 0.0)
		inf = file(fname, 'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf)
		header = fp.getHeader()
		vals = []
		while not fp.atEnd():
			#print fp.cache.cache[-1],
			v = fp.next()
			#print fp.isValid(), v
			vals.append(v[2])
		self.assertTrue(sum(vals) > n_lines)
		inf.close()
		os.remove(fname)

class test002(unittest.TestCase):
	"""read through"""
	def test_run(self):
		# Normal
		n_lines = 100
		header_list = ["str","float","int","str"]
		fname = "tmp_normal.txt"
		makeFile(fname, header_list, "sfds", n_lines, '\t', 0.0)
		inf = file(fname, 'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf)
		while not fp.atEnd():
			fp.next()
		inf.close()
		os.remove(fname)

class test003(unittest.TestCase):
	"""header parsing"""
	def test_run(self):
		# Normal
		n_lines = 100
		header_list = ["str","float","int","anotherStr"]
		fname = "tmp_normal.txt"
		makeFile(fname, header_list, "sfds", n_lines, '\t', 0.0)
		inf = file(fname, 'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf)
		header = fp.getHeader()
		res = header == header_list
		inf.close()
		os.remove(fname)
		return res

class test004(unittest.TestCase):
	"""infer header types 2"""
	def test_run(self):
		# Normal
		n_lines = 100
		header_list = ["str","float","int","str"]
		fname = "tmp_normal.txt"
		makeFile(fname, None, "sfds", n_lines, '\t', 0.0)
		inf = file(fname, 'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf)
		header = fp.getHeader()
		vals = []
		while not fp.atEnd():
			#print fp.cache.cache[-1],
			v = fp.next()
			#print fp.isValid(), v
			vals.append(v[2])
		self.assertTrue(sum(vals) > n_lines)
		inf.close()
		os.remove(fname)

class test005(unittest.TestCase):
	"""lines read"""
	def test_run(self):
		# Normal
		n_lines = random.randint(10,300)
		header_list = ["str","float","int","str"]
		fname = "tmp_normal.txt"
		makeFile(fname, None, "sfds", n_lines, '\t', 0.0)
		inf = file(fname, 'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf)
		header = fp.getHeader()
		while not fp.atEnd():
			v = fp.next()
		res = fp.getNumRead() == n_lines
		inf.close()
		os.remove(fname)
		return res

class test006(unittest.TestCase):
	"""infer header types with NA's"""
	def test_run(self):
		# NA's at some frequency -- infer types.
		n_lines = 100
		header_list = ["str","float","int","str"]
		fname = "tmp_normal.txt"
		# Put an NA in the first line -- test the ability to look past to infer types
		first_lines = [randString() + "\t1.00\tNA\t" + randString() + "\n"]
		makeFile(fname, None, "sfds", n_lines, '\t', 0.1, first_lines=first_lines, last_lines=None)
		inf = file(fname, 'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf)
		header = fp.getHeader()
		vals = []
		while not fp.atEnd():
			#print fp.cache.cache[-1],
			v = fp.next()
			#print v
			if not v[2] is None:
				vals.append(v[2])
		self.assertTrue(sum(vals) > n_lines)
		inf.close()
		os.remove(fname)

class test007(unittest.TestCase):
	"""infer header types with NA's in every line"""
	def test_run(self):
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
		fp = util.DelimitedLineReader(inf, header=False)
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

class test008(unittest.TestCase):
	"""infer header types with NA's in one full column"""
	def test_run(self):
		# NA's in every line -- infer types.
		n_lines = 0
		header_list = ["str","float","int","str"]
		fname = "tmp_normal.txt"
		# Put an NA in the first line -- test the ability to look past to infer types
		first_lines = [randString() + "\t59.3\tNA\t" + randString() + "\n"]*100
		makeFile(fname, None, "sfds", n_lines, '\t', 0.1, first_lines=first_lines, last_lines=None)
		inf = file(fname,'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf, header=False)
		header = fp.getHeader()
		inf.close()
		os.remove(fname)

class test009(unittest.TestCase):
	"""header processing"""
	def test_run(self):
		# Header processing
		n_lines = 100
		header_list = ["int","int.1","int","int"]
		fname = "tmp_normal.txt"
		makeFile(fname, header_list, "ffds", n_lines, '\t', 0.1)
		inf = file(fname,'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf)
		header = fp.getHeader()
		self.assertTrue(header[0] == 'int')
		self.assertTrue(header[1] == 'int.1')
		self.assertTrue(header[2] == 'int.2')
		self.assertTrue(header[3] == 'int.3')
		inf.close()
		os.remove(fname)

class test010(unittest.TestCase):
	"""header processing: custom header processor"""
	def test_run(self):
		# Header processing
		n_lines = 100
		header_list = ["Ratio (H/L)","float","int","This + That"]
		fname = "tmp_normal.txt"
		makeFile(fname, header_list, "ffds", n_lines, '\t', 0.1)
		inf = file(fname,'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf, header=True, header_name_processor=util.maxQuantHeader)
		header = fp.getHeader()
		self.assertTrue(header[0] == 'ratio.hl')
		self.assertTrue(header[-1] == 'this.plus.that')
		inf.close()
		os.remove(fname)

class test011(unittest.TestCase):
	"""adaptive handler updating"""
	def test_run(self):
		# Adaptive field redefinition
		n_lines = 100
		header_list = ["float","float","int","str"]
		fname = "tmp_normal.txt"
		last_lines = ["0.2\t0.3\tnot.an.int\twakka\n"]
		makeFile(fname, header_list, "ffds", n_lines, '\t', 0.01, last_lines=last_lines)
		inf = file(fname,'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf, header=True, header_name_processor=util.maxQuantHeader)
		header = fp.getHeader()
		while not fp.atEnd():
			flds = fp.next()
		inf.close()
		os.remove(fname)

class test012(unittest.TestCase):
	"""comment as last line"""
	def test_run(self):
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
		fp = util.DelimitedLineReader(inf)
		while not fp.atEnd():
			flds = fp.nextDict()
		inf.close()
		os.remove(fname)

class test013(unittest.TestCase):
	"""no headers but call to nextDict"""
	def test_run(self):
		# No headers but call to nextDict
		n_lines = 10
		fname = "tmp_no_header.txt"
		field_types = "sfds"
		makeFile(fname, None, field_types, n_lines, '\t', 0.0)
		inf = file(fname, 'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf, header=False)
		while not fp.atEnd():
			flds = fp.nextDict()
			self.assertTrue(set(flds.keys()) == set(range(len(field_types))))
		inf.close()
		os.remove(fname)

class test014(unittest.TestCase):
	"""comments as first lines"""
	def test_run(self):
		n_lines = 10
		header_list = ["str","float","int","str"]
		fname = "tmp_comment_first.txt"
		makeFile(fname, header_list, "sfds", n_lines, '\t', 0.0, first_lines="# first line comment\n# second line comment")
		inf = file(fname, 'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf)
		while not fp.atEnd():
			flds = fp.nextDict()
		# if we make it this far, we did not throw an error.
		self.assertTrue(True)
		inf.close()
		os.remove(fname)

class test015(unittest.TestCase):
	"""multiple delimiters -- spaces -- with sep=None"""
	def test_run(self):
		n_lines = 10
		header_list = ["str","float","int","str"]
		fname = "tmp_comment_first.txt"
		inf = file(fname, 'w')
		inf.write("one\ttwo\tthree\n")
		inf.write("  a  b c\n")
		inf.write(" a   b  c\n")
		inf.close()
		inf = file(fname, 'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf, header=True, sep=None)
		while not fp.atEnd():
			flds = fp.nextDict()
			self.assertTrue(flds['one']=='a')
			self.assertTrue(flds['two']=='b')
			self.assertTrue(flds['three']=='c')
			self.assertFalse(flds['three']=='b')
		inf.close()
		os.remove(fname)

class test016(unittest.TestCase):
	"""iterator over entries"""
	def test_run(self):
		n_lines = 10
		header_list = ["str","float","int","str"]
		fname = "tmp_comment_first.txt"
		inf = file(fname, 'w')
		inf.write("one\ttwo\tthree\n")
		inf.write("  a  b c\n")
		inf.write(" a   b  c\n")
		inf.close()
		inf = file(fname, 'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf, header=True, sep=None)
		for flds in fp.entries:
			self.assertTrue(flds[0]=='a')
			self.assertTrue(flds[1]=='b')
			self.assertTrue(flds[2]=='c')
		inf.close()
		os.remove(fname)

class test017(unittest.TestCase):
	"""iterator over dict entries"""
	def test_run(self):
		n_lines = 10
		header_list = ["str","float","int","str"]
		fname = "tmp_comment_first.txt"
		inf = file(fname, 'w')
		inf.write("one\ttwo\tthree\n")
		inf.write("  a  b c\n")
		inf.write(" a   b  c\n")
		inf.close()
		inf = file(fname, 'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf, header=True, sep=None)
		for flds in fp.dictentries:
			self.assertTrue(flds['one']=='a')
			self.assertTrue(flds['two']=='b')
			self.assertTrue(flds['three']=='c')
			self.assertFalse(flds['three']=='b')
		inf.close()
		os.remove(fname)


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
		outf.write("{}\n".format(sep.join(headers)))
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
		outf.write("{}\n".format(line))
	if last_lines:
		for line in last_lines:
			outf.write(line)
	outf.close()
	return fname

if __name__=="__main__":
	test_printTiming()
	unittest.main(verbosity=2)
