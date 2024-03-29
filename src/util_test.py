import time, os, random, string, sys, math, traceback, unittest
import util, na

@util.printTiming
def test_printTiming():
	time.sleep(0.5)
	print('The next message should show at least 500ms of wait time:\n\t', file=sys.stdout)

class test001(unittest.TestCase):
	def test_run(self):
		"""infer header types"""
		# Normal
		n_lines = 100
		header_list = ["str","float","int","str"]
		fname = "tmp_normal.txt"
		makeFile(fname, header_list, "sfds", n_lines, '\t', 0.0)
		inf = open(fname, 'r')
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
	def test_run(self):
		"""read through"""
		# Normal
		n_lines = 100
		header_list = ["str","float","int","str"]
		fname = "tmp_normal.txt"
		makeFile(fname, header_list, "sfds", n_lines, '\t', 0.0)
		inf = open(fname, 'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf)
		while not fp.atEnd():
			fp.next()
		inf.close()
		os.remove(fname)

class test003(unittest.TestCase):
	def test_run(self):
		"""header parsing"""
		# Normal
		n_lines = 100
		header_list = ["str","float","int","anotherStr"]
		fname = "tmp_normal.txt"
		makeFile(fname, header_list, "sfds", n_lines, '\t', 0.0)
		inf = open(fname, 'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf)
		header = fp.getHeader()
		res = header == header_list
		inf.close()
		os.remove(fname)
		return res

class test004(unittest.TestCase):
	def test_run(self):
		"""infer header types 2"""
		# Normal
		n_lines = 100
		header_list = ["str","float","int","str"]
		fname = "tmp_normal.txt"
		makeFile(fname, None, "sfds", n_lines, '\t', 0.0)
		inf = open(fname, 'r')
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
	def test_run(self):
		"""lines read"""
		# Normal
		n_lines = random.randint(10,300)
		header_list = ["str","float","int","str"]
		fname = "tmp_normal.txt"
		makeFile(fname, None, "sfds", n_lines, '\t', 0.0)
		inf = open(fname, 'r')
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
	def test_run(self):
		"""infer header types with NA's"""
		# NA's at some frequency -- infer types.
		n_lines = 100
		header_list = ["str","float","int","str"]
		fname = "tmp_normal.txt"
		# Put an NA in the first line -- test the ability to look past to infer types
		first_lines = [randString() + "\t1.00\tNA\t" + randString() + "\n"]
		makeFile(fname, None, "sfds", n_lines, '\t', 0.1, first_lines=first_lines, last_lines=None)
		inf = open(fname, 'r')
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
	def test_run(self):
		"""infer header types with NA's in every line"""
		# NA's in every line -- infer types.
		n_lines = 0
		header_list = ["str","float","int","str"]
		fname = "tmp_normal.txt"
		# Put an NA in the first line -- test the ability to look past to infer types
		first_lines = [randString() + "\t59.3\tNA\t" + randString() + "\n",
			randString() + "\tNA\t12\t" + randString() + "\n"]
		makeFile(fname, None, "sfds", n_lines, '\t', 0.1, first_lines=first_lines, last_lines=None)
		inf = open(fname,'r')
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
	def test_run(self):
		"""infer header types with NA's in one full column"""
		# NA's in every line -- infer types.
		n_lines = 0
		header_list = ["str","float","int","str"]
		fname = "tmp_normal.txt"
		# Put an NA in the first line -- test the ability to look past to infer types
		first_lines = [randString() + "\t59.3\tNA\t" + randString() + "\n"]*100
		makeFile(fname, None, "sfds", n_lines, '\t', 0.1, first_lines=first_lines, last_lines=None)
		inf = open(fname,'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf, header=False)
		header = fp.getHeader()
		inf.close()
		os.remove(fname)

class test009(unittest.TestCase):
	def test_run(self):
		"""header processing"""
		# Header processing
		n_lines = 100
		header_list = ["int","int.1","int","int"]
		fname = "tmp_normal.txt"
		makeFile(fname, header_list, "ffds", n_lines, '\t', 0.1)
		inf = open(fname,'r')
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
	def test_run(self):
		"""header processing: custom header processor"""
		# Header processing
		n_lines = 100
		header_list = ["Ratio (H/L)","float","int","This + That"]
		fname = "tmp_normal.txt"
		makeFile(fname, header_list, "ffds", n_lines, '\t', 0.1)
		inf = open(fname,'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf, header=True, header_name_processor=util.maxQuantHeader)
		header = fp.getHeader()
		self.assertTrue(header[0] == 'ratio.hl')
		self.assertTrue(header[-1] == 'this.plus.that')
		inf.close()
		os.remove(fname)

class test011(unittest.TestCase):
	def test_run(self):
		"""adaptive handler updating"""
		# Adaptive field redefinition
		n_lines = 100
		header_list = ["float","float","int","str"]
		fname = "tmp_normal.txt"
		last_lines = ["0.2\t0.3\tnot.an.int\twakka\n"]
		makeFile(fname, header_list, "ffds", n_lines, '\t', 0.01, last_lines=last_lines)
		inf = open(fname,'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf, header=True, header_name_processor=util.maxQuantHeader)
		header = fp.getHeader()
		while not fp.atEnd():
			flds = fp.next()
		inf.close()
		os.remove(fname)

class test012(unittest.TestCase):
	def test_run(self):
		"""comment as last line"""
		# Comment as last line
		n_lines = 10
		header_list = ["str","float","int","str"]
		fname = "tmp_comment_last.txt"
		makeFile(fname, header_list, "sfds", n_lines, '\t', 0.0)
		inf = open(fname, 'a')
		inf.write("# comment\n")
		inf.close()
		inf = open(fname, 'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf)
		while not fp.atEnd():
			flds = fp.nextDict()
		inf.close()
		os.remove(fname)

class test013(unittest.TestCase):
	def test_run(self):
		"""no headers but call to nextDict"""
		# No headers but call to nextDict
		n_lines = 10
		fname = "tmp_no_header.txt"
		field_types = "sfds"
		makeFile(fname, None, field_types, n_lines, '\t', 0.0)
		inf = open(fname, 'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf, header=False)
		while not fp.atEnd():
			flds = fp.nextDict()
			self.assertTrue(set(flds.keys()) == set(range(len(field_types))))
		inf.close()
		os.remove(fname)

class test014(unittest.TestCase):
	def test_run(self):
		"""comments as first lines"""
		n_lines = 10
		header_list = ["str","float","int","str"]
		fname = "tmp_comment_first.txt"
		makeFile(fname, header_list, "sfds", n_lines, '\t', 0.0, first_lines="# first line comment\n# second line comment")
		inf = open(fname, 'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf)
		while not fp.atEnd():
			flds = fp.nextDict()
		# if we make it this far, we did not throw an error.
		self.assertTrue(True)
		inf.close()
		os.remove(fname)

class test015(unittest.TestCase):
	def test_run(self):
		"""multiple delimiters -- spaces -- with sep=None"""
		n_lines = 10
		header_list = ["str","float","int","str"]
		fname = "tmp_comment_first.txt"
		inf = open(fname, 'w')
		inf.write("one\ttwo\tthree\n")
		inf.write("  a  b c\n")
		inf.write(" a   b  c\n")
		inf.close()
		inf = open(fname, 'r')
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
	def test_run(self):
		"""iterator over entries"""
		n_lines = 10
		header_list = ["str","float","int","str"]
		fname = "tmp_comment_first.txt"
		inf = open(fname, 'w')
		inf.write("one\ttwo\tthree\n")
		inf.write("  a  b c\n")
		inf.write(" a   b  c\n")
		inf.close()
		inf = open(fname, 'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf, header=True, sep=None)
		for flds in fp.entries:
			self.assertTrue(flds[0]=='a')
			self.assertTrue(flds[1]=='b')
			self.assertTrue(flds[2]=='c')
		inf.close()
		os.remove(fname)

class test017(unittest.TestCase):
	def test_run(self):
		"""iterator over dict entries"""
		n_lines = 10
		header_list = ["str","float","int","str"]
		fname = "tmp_comment_first.txt"
		inf = open(fname, 'w')
		inf.write("one\ttwo\tthree\n")
		inf.write("  a  b c\n")
		inf.write(" a   b  c\n")
		inf.close()
		inf = open(fname, 'r')
		# Infer the types
		fp = util.DelimitedLineReader(inf, header=True, sep=None)
		for flds in fp.dictentries:
			self.assertTrue(flds['one']=='a')
			self.assertTrue(flds['two']=='b')
			self.assertTrue(flds['three']=='c')
			self.assertFalse(flds['three']=='b')
		inf.close()
		os.remove(fname)

class test018(unittest.TestCase):
	def test_run(self):
		"""set handler type"""
		fname = "tmp_types.txt"
		inf = open(fname, 'w')
		inf.write("one\ttwo\tthree\n")
		#inf.write("  NA	NA	NA\n")
		inf.write(" 1.0  two  3\n")
		inf.close()
		inf = open(fname, 'r')
		# Assume all types are string
		fp = util.DelimitedLineReader(inf, header=True, field_defs='sss', sep=None)
		fp.setColumnType("three", "d")
		flds = fp.nextDict()
		self.assertTrue(flds['three']==3)
		inf.close()
		os.remove(fname)


def randString():
	return ''.join(random.sample(string.ascii_letters, 10))
	
def randInt():
	return random.randint(1,10000)

def randFloat():
	return random.random()

def makeFile(fname, headers, fld_types, n_lines, sep, na_prob, first_lines=None, last_lines=None):
	dtype = {}
	dtype['f'] = randFloat
	dtype['d'] = randInt
	dtype['s'] = randString
	outf = open(fname,'w')
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


class test019(unittest.TestCase):
	def test_run(self):
		"""listdict"""
		ld = util.listdict()
		ld['a'].append(1)
		ld['a'].append(2)
		self.assertTrue(ld['a'] == [1,2])

class test020(unittest.TestCase):
	def test_run(self):
		"""listdict initialization"""
		d = {'a':1, 'b':2}
		ld = util.listdict(d)
		ld['a'].append(2)
		self.assertTrue(ld['a'] == [1,2])
		self.assertTrue(ld['b'] == [2])

class test021(unittest.TestCase):
	def test_run(self):
		"""readTable basic"""
		fname = "tmp_lightdataframe.txt"
		inf = open(fname, 'w')
		inf.write("one\ttwo\tthree\n")
		inf.write("a\tb\t3\n")
		inf.write("a\tb\t33\n")
		inf.close()
		inf = open(fname, 'r')
		ldf = util.readTable(inf, header=True)
		self.assertTrue(ldf['three'][0] == 3)
		self.assertTrue(ldf['three'][1] == 33)
		inf.close()
		os.remove(fname)


class test022(unittest.TestCase):
	def test_run(self):
		"""readTable dictrows"""
		fname = "tmp_lightdataframe.txt"
		inf = open(fname, 'w')
		inf.write("one\ttwo\tthree\n")
		inf.write("a\tb\t3\n")
		inf.write("a\tb\t33\n")
		inf.close()
		inf = open(fname, 'r')
		ldf = util.readTable(inf, header=True)
		for (ri, flds) in enumerate(ldf.dictrows):
			if ri == 0:
				self.assertTrue(flds['three'] == 3)
			if ri == 1:
				self.assertTrue(flds['three'] == 33)
		inf.close()
		os.remove(fname)

class test023(unittest.TestCase):
	def test_run(self):
		"""skip, positive case"""
		fname = "tmp_skip.txt"
		inf = open(fname, 'w')
		inf.write("blah\nblah\nblah\n")
		inf.write("one\ttwo\tthree\n")
		inf.write("a\tb\t3\n")
		inf.write("a\tb\t33\n")
		inf.close()
		with open(fname,'r') as inf:
			dlr = util.DelimitedLineReader(inf, header=True, skip=3)
			for (ri, flds) in enumerate(dlr.dictentries):
				if ri == 0:
					self.assertTrue(flds['three'] == 3)
				if ri == 1:
					self.assertTrue(flds['three'] == 33)
		inf.close()
		os.remove(fname)

class test024(unittest.TestCase):
	def test_run(self):
		"""skip, negative case"""
		fname = "tmp_skip.txt"
		inf = open(fname, 'w')
		inf.write("blah\nblah\nblah\n")
		inf.write("one\ttwo\tthree\n")
		inf.write("a\tb\t3\n")
		inf.write("a\tb\t33\n")
		inf.close()
		with open(fname,'r') as inf:
			dlr = util.DelimitedLineReader(inf, header=True, skip=3)
			for (ri, flds) in enumerate(dlr.dictentries):
				try:
					flds['three']
				except KeyError:
					self.assertTrue(True)
		inf.close()
		os.remove(fname)

class test025(unittest.TestCase):
	def test_run(self):
		"""skip, off by one"""
		fname = "tmp_skip.txt"
		inf = open(fname, 'w')
		inf.write("blah\nblah\nblah\n")
		inf.write("one\ttwo\tthree\n")
		inf.write("a\tb\t3\n")
		inf.write("a\tb\t33\n")
		inf.close()
		with open(fname,'r') as inf:
			dlr = util.DelimitedLineReader(inf, header=True, skip=2)
			for (ri, flds) in enumerate(dlr.dictentries):
				try:
					flds['three']
				except KeyError:
					self.assertTrue(True)
		inf.close()
		os.remove(fname)

class test026(unittest.TestCase):
	def test_run(self):
		"""set header names"""
		fname = "tmp_setheader.txt"
		inf = open(fname, 'w')
		inf.write("one\ttwo\tthree\n")
		inf.write("a\tb\t3\n")
		inf.write("a\tb\t33\n")
		inf.close()
		with open(fname,'r') as inf:
			dlr = util.DelimitedLineReader(inf, header=True)
			dlr.setHeaderNames(['a','b','c'])
			for (ri, flds) in enumerate(dlr.dictentries):
				if ri == 0:
					self.assertTrue(flds['c'] == 3)
				if ri == 1:
					self.assertTrue(flds['c'] == 33)
		inf.close()
		os.remove(fname)

class test027(unittest.TestCase):
	def test_run(self):
		"""readTable header"""
		fname = "tmp_lightdataframe.txt"
		inf = open(fname, 'w')
		inf.write("one\ttwo\tthree\n")
		inf.write("a\tb\t3\n")
		inf.write("a\tb\t33\n")
		inf.close()
		with open(fname,'r') as inf:
			ldf = util.readTable(inf, header=True)
			h = ldf.header
			self.assertTrue(h[1]=='two')
		os.remove(fname)

class test028(unittest.TestCase):
	def test_run(self):
		"""formatLine"""
		# DAD: todo
		#self.assertFalse(True)
		dout = util.DelimitedOutput()
		dout.addHeader("one","first",'d')
		dout.addHeader("two","second",'s')
		entry = {'one':2, 'two':'this'}
		line = dout.formatLine(entry)

class test029(unittest.TestCase):
	def test_run(self):
		"""encoding"""
		fname = "tmp_lightdataframe.txt"
		inf = open(fname, 'w')
		inf.write("one\ttwo\tthree\n")
		inf.write("a\tb\t3\n")
		inf.write("a\tb\t“extra%20TFIIIC%29”\n")
		inf.close()
		with open(fname,'r') as inf:
			ldf = util.DelimitedLineReader(inf, header=True)
			#h = ldf.header
			for line in ldf.entries:
				print(line)
			#self.assertTrue(h[1]=='two')
		os.remove(fname)


if __name__=="__main__":
	test_printTiming()
	unittest.main(verbosity=2)
