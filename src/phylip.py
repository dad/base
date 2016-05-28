import sys, os, string, re

class Scanner:
	def feed(self, handle, consumer):
		if isinstance(handle, File.UndoHandle):
			uhandle = handle
		else:
			uhandle = File.UndoHandle(handle)
			self._scan_tree(uhandle, consumer)
			
	def _scan_tree_uhandle(self, uhandle, consumer):
		tree_text = "".join(uhandle.readlines()).replace("\n","")
		tree_text = re.sub('\s', '', tree_text)
		#print "text: ", tree_text
		_scan_tree_text(self, tree_text, consumer)

	def _scan_tree_text(self, tree_text, consumer):
		pos = 0
		rooted = is_rooted(tree_text)
		consumer.start_tree(rooted)
		
		while True:
			c = tree_text[pos]
			if c == '(':
				consumer.begin_node()
				pos += 1
			elif c == ')':
				name = ''
				pos += 1
				c = tree_text[pos]
				while re.match('[_\.\w]',c):
					name += c
					pos += 1
					c = tree_text[pos]
				consumer.end_node(name)
			elif c == ',':
				pos += 1
			elif c == ':':
				# ready to process branch length
				pos += 1
				c = tree_text[pos]
				length = ''
				while re.match('[\.\d]',c):
					length += c
					pos += 1
					c = tree_text[pos]
				consumer.branch_length(float(length))
			elif c == ';':
				consumer.end_tree()
				break
			elif c == "'":
				pos += 1
			else:
				name = ''
				while re.match('[_\.\w]',c):
					name += c
					pos += 1
					c = tree_text[pos]
				consumer.leaf(name)

def is_rooted(tree_string):
	"""Reports whether a tree is rooted by determining whether there is only one comma in the depth-1 node."""
	pos = 0
	c = tree_string[pos]
	depth = 0
	comma = 0
	for pos in range(0,len(tree_string)):
		c = tree_string[pos]
		if c == '(':
			depth += 1
		elif c == ')':
			depth -= 1
		elif c == ',':
			if depth == 1:
				comma += 1
	return comma == 1


class Consumer: # (AbstractConsumer):
	tab = ''
	def __init__(self):
		self.stack = Stack()
		self.data = None

	def printit(self, s):
		#return
		print self.tab,s
		
	def start_tree(self, rooted=False):
		self.rooted = rooted
		self.pos = 0
		self.depth = 0
		self.stack.empty()
		
	def begin_node(self):
		#self.printit("begin")
		self.tab += '\t'
		node = Node()
		self.stack.push(node)
		self.depth = self.depth + 1
		
	def end_node(self, name):
		# 3 cases: depth 1 and unrooted tree (3 nodes)
		#       or internal node             (2 nodes)
		#       or (depth 1 and rooted tree) (2 nodes)
		if self.depth == 1 and not self.rooted:
			right = self.stack.pop()
			middle = self.stack.pop()
			left = self.stack.pop()
			parent = self.stack.top()
			parent.left = left
			parent.middle = middle
			parent.right = right
			parent.name = name
			right.parent = parent
			left.parent = parent
			middle.parent = parent			
		else:
			# internal node or rooted tree
			right = self.stack.pop()
			left = self.stack.pop()
			parent = self.stack.top()
			parent.left = left
			parent.right = right
			parent.name = name
			right.parent = parent
			left.parent = parent
			self.depth -= 1
		if False: # debugging
			self.printit(right.name)
			self.printit(left.name)
			self.tab = self.tab[0:-1]
			self.printit(parent.name)
			#self.printit("end")
			
	def leaf(self, name):
		node = Node()
		node.name = name
		self.stack.push(node)
		
	def branch_length(self, l):
		node = self.stack.top()
		node.length = l
		if False: # debugging
			self.printit("%f"%node.length)

	def end_tree(self):
		self.data = self.stack.pop()
		
class Node:
	""" A simplified representation of a tree (just nodes) """
	def __init__(self, name=None, left=None, middle=None, right=None, length=None, parent=None):
		self.name = name
		self.left = left
		self.right = right
		self.length = length
		self.middle = middle
		self.parent = parent
		#self.children = [x for x in [self.left, self.right, self.middle] if not x is None]
		
	def __str__(self):
		return self.recur_str()+";"

	# Newick tree format:
	# http://evolution.genetics.washington.edu/phylip/newicktree.html
	# Ancestral nodes may be named.
	def recur_str(self):
		if self.middle:
			# Middle child only for unrooted tree at base.
			assert (self.left is not None) and (self.right is not None)
			return "(" + str(self.left.recur_str()) + "," + str(self.middle.recur_str()) + "," + str(self.right.recur_str()) + ")"
		elif self.left:
			# If there's a left child, then there must be a right child
			assert (self.left is not None) and (self.right is not None)			
			if not self.name:
				name = ""
			else:
				name = self.name
			out_str = "(" + str(self.left.recur_str()) + "," + str(self.right.recur_str()) + ')%s' % name
			if self.length is not None:
				out_str += ":%1.6f" % self.length
			return out_str
		else:
			# No children
			out_str = self.name
			if self.length is not None:
				out_str += ":%1.6f" % self.length
			return out_str

	def children(self):
		return [x for x in [self.left, self.right, self.middle] if not x is None]
	
	def tree_nodes(self):
		nodes = []
		self.subtree(nodes)
		return nodes

	def leaves(self):
		return [n for n in self.tree_nodes() if len(n.children())==0]
	
	def subtree(self, nodes):
		nodes.append(self)
		#for c in self.children:
		#	c.subtree(nodes)
		if self.middle:
			self.middle.subtree(nodes)
		if self.left:
			self.left.subtree(nodes)
		if self.right:
			self.right.subtree(nodes)			

	def get_lineage(self):
		lineage = [self]
		top_parent = self
		while top_parent.parent:
			lineage.append(top_parent.parent)
			top_parent = top_parent.parent
		return lineage

	def distance(self, node):
		my_lineage = self.get_lineage()
		her_lineage = node.get_lineage()
		#ps(my_lineage)
		#ps(her_lineage)
		branch_nodes = set(my_lineage).symmetric_difference(set(her_lineage))
		#ps(branch_nodes)
		length = sum([x.length for x in branch_nodes])
		return length

def ps(l):
	for x in l:
		print x.name, x.length
	print "  **"

class Stack:
	def __init__(self):
		self._l = []
		
	def push(self, item):
		self._l.append(item)
			
	def pop(self):
		last = self._l[-1]
		self._l = self._l[:-1]
		return last

	def top(self):
		return self._l[-1]
	
	def empty(self):
		self._l = []
		
				
class Parser: #(AbstractParser):
	def __init__(self):
		self._scanner = Scanner()
		self._consumer = Consumer()

	def read(self, handle):
		self._scanner.feed(handle, self._consumer)
		return self._consumer.data

	def parse(self, tree_text):
		self._scanner._scan_tree_text(tree_text, self._consumer)
		return self._consumer.data

def parse_tree(tree_string):
	p = Parser()
	return p.parse(tree_string.replace(' ',''))

def read_tree(tree_filename):
	tree_string = ''
	for line in file(tree_filename,'r').readlines():
		if line[0] != "#":
			tree_string += line.strip()
	return parse_tree(tree_string)
					
def test():
	tree_strings = ["(((sbay:50.000000,scas:0.000004)8:11.818649,(smik:6.760546,spar:8.815752)9:19.825863)7:0.000000,sklu:0.000004)6;",\
					"(Bovine:0.693950,(Hylobates:0.360790,(Pongo:0.336360,(G._Gorilla:0.171470, (P._paniscus:0.192680,H._sapiens:0.119270):0.083860):0.061240):0.150570):0.549390, Rodent:1.214600);",\
					"A;",\
					"(x1,x2,x3);"]
	for tree in tree_strings:
		#print tree
		t = parse_tree(tree)
		new_tree = "%s" % t
		#print new_tree
		assert tree.replace(" ",'') == new_tree
	print "All phylip tests passed"
	
