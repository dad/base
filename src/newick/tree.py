'''
A Python module for parsing Newick files.

Copyright (C) 2003-2006, Thomas Mailund <mailund@birc.au.dk>

This module contains the representation of trees and a parser for
creating trees from a Newick string or file.

Updated, extended and maintained by D. Allan Drummond <dadrummond@gmail.com>
'''

import lexer
import parser

class TreeError(Exception):
	"""Error manipulating a newick.tree object."""


class Tree(object):
	'''
	Python representation of a tree.
	'''

	def __init__(self):
		self._children = []
		self._leaves_cache = None
		self.name = None
		self.length = None
		self.parent = None

	def addChild(self, c):
		self._children.append(c)
		c.parent = self
		# we need to invalidate this when we add children
		self._leaves_cache = None

	def removeLeaf(self, c):
		# DAD: fix this!
		#print self.name, "removing", c.name
		if not c in self.leaves:
			raise TreeError, "%s is not a leaf of %s" % (c, self)
		par = c.parent
		par._children.remove(c)
		par._leaves_cache = None
		root = self.lineage[-1]
		if len(par._children) == 1:
			lone_kid = par.children[0]
			pp = par.parent
			if not pp is None:
				pp._children.remove(par)
				pp.addChild(lone_kid)
			else:
				# par is the root
				# re-root the tree at lone_kid
				lone_kid.parent = None
				root = lone_kid
		return root

	def unroot(self):
		"""Unroot the tree in the Newick sense."""
		root = self.root
		if root.isRooted() and len(root.leaves) > 2:
			# add one child's children directly to the root
			# find first child that's not childless
			kkids = [x for x in root.children if not x.isLeaf()]
			target_kid = kkids[0]
			for c in target_kid.children:
				root.addChild(c)
			root._children.remove(target_kid)

	def isRooted(self):
		root = self.root
		res = len(root._children) == 2
		if res:
			# only way to not be rooted is to have only 2 leaves.
			res = not (len(root.leaves) == 2)
		return res

	def getRoot(self):
		return self.lineage[-1]

	def getChildren(self):
		'''
		get_children() -- return the list of child leaves/sub-trees.
		'''
		return self._children

	def dfs_traverse(self,visitor):
		'''
		dfs_traverse(visitor) -- do a depth first traversal.
		Part of the Visitor Pattern; performs a depth first traversal,
		calling methods in visitor along the way.
		'''
		visitor.pre_visit_tree(self)
		for c in self._children:
			visitor.pre_visit_child(self, c)
			c.dfs_traverse(visitor)
			visitor.post_visit_child(self, c)
		visitor.post_visit_tree(self)

	def getLeaves(self):
		'''
		getLeaves() --  return list of leaves in this (sub-)tree.
		'''
		if self.isLeaf():
			return [self]
		elif self._leaves_cache is None:
			self._leaves_cache = []
			for c in self._children:
				self._leaves_cache.extend(c.leaves)
		return self._leaves_cache

	def getNodes(self):
		res = [self]
		for c in self._children:
			res.extend(c.getNodes())
		return res

	def getLineage(self):
		lineage = [self]
		top_parent = self
		while top_parent.parent:
			lineage.append(top_parent.parent)
			top_parent = top_parent.parent
		return lineage

	def getMostRecentCommonAncestor(self, node2):
		lineage1 = self.lineage
		lineage2 = set(node2.lineage)
		mrca = None
		for n in lineage1:
			if n in lineage2:
				mrca = n
				break
		return mrca

	def measureFrom(self, node, measureFxn, combineFxn=sum):
		my_lineage = set(self.lineage)
		her_lineage = set(node.lineage)
		branch_nodes = my_lineage.symmetric_difference(her_lineage)
		meas = combineFxn([measureFxn(x) for x in branch_nodes ])
		return meas

	def distanceFrom(self, node):
		return self.measureFrom(node, lambda x: x.length, sum)

	def isLeaf(self):
		return len(self._children)==0

	# special functions and accessors...
	def __repr__(self):
		#print "# kids = ", len(self.children)
		tree_str = ""
		if not self.isLeaf():
			tree_str = '('
			sep = ''
			for c in self.children:
				tree_str += sep+str(c)
				sep = ','
			tree_str += ")"
		#if not self.name is None:
		if self.isLeaf():
			tree_str += "%s" % self.name
		if not self.length is None:
			tree_str += ":%1.6f" % self.length
		if self.parent == None: # I'm the root!
			tree_str += ";"
		return tree_str

	leaves = property(getLeaves, None, None,
					  "List of leaves in this subtree.")
	children = property(getChildren, None, None,
					  "List of child nodes for this subtree.")
	lineage = property(getLineage, doc="Line of descent of this subtree back to root.")

	nodes = property(getNodes, None, None, doc="List of nodes in this subtree.")

	root = property(getRoot, None, None, doc="Root of this tree object (valid even for unrooted trees).")

class TreeVisitor(object):
	'''
	Part of the Visitor Pattern.
	'''

	def __init__(self):
		pass

	def pre_visit_tree(self,t):
		'''
		pre_visit_tree(t) -- callback called before exploring (sub-)tree t.
		'''
		pass

	def post_visit_tree(self,t):
		'''
		post_visit_tree(t) -- callback called after exploring (sub-)tree t.
		'''
		pass

	def pre_visit_child(self,child):
		'''
		pre_visit_edge(src, child)
			-- callback called before exploring an edge.

		'''
		pass

	def post_visit_child(self,child):
		'''
		post_visit_edge(src, child)
			-- callback called before exploring an edge.

		'''
		pass

	def visit_leaf(self,l):
		'''
		visit_leaf(l) -- callback called when exploring leaf l.
		'''
		pass


class _TreeBuilder(parser.AbstractHandler):
	"""
	Builds a tree
	"""
	def __init__(self):
		self.stack = []
		self.root = None

	def newTreeBegin(self):
		#print "tb.newTreeBegin", len(self.stack)
		t = Tree()
		if len(self.stack) == 0:
			# This is the root of the tree
			self.root = t
			t.parent = None
		else:
			# Not the root -- add as a child
			parent = self.stack[-1]
			parent.addChild(t)
		# Top of the stack is always the most recent parent
		self.stack.append(t)

	def newTreeEnd(self, bootstrap, length):
		#print "tb.newTreeEnd", len(self.stack), bootstrap, length
		parent = self.stack.pop()
		if bootstrap is not None:
			parent.name = "%s" % bootstrap
		parent.length = length

	def newLeaf(self,name,length):
		#print "newLeaf", name, length
		child = Tree()
		child.name = name
		child.length = length
		parent = self.stack[-1]
		parent.addChild(child)

	def getResult(self):
		return self.root


def parseTree(input):
	'''Parse input as a Newick description of a tree and return the
	tree in a tree data structure.'''
	return parser.parse(input,_TreeBuilder())

def readTree(file):
	string = ''
	for line in file.readlines():
		if not line.startswith("#"):
			string += line.strip()
	return parseTree(string)


def add_parent_links(tree):
	'''Extend all nodes (except for the root, of course) with a parent link.'''
	class V(TreeVisitor):
		def __init__(self):
			self.cur_parent = [None]
		def pre_visit_edge(self,src,b,l,dst):
			dst.parent = src
		def pre_visit_tree(self,t):
			t.parent = self.cur_parent[-1]
			self.cur_parent.append(t)
		def post_visit_tree(self,t):
			self.cur_parent.pop()
		def visit_leaf(self,l):
			l.parent = self.cur_parent[-1]
	tree.dfs_traverse(V())

def add_distance_from_root(tree):
	'''Extend all nodes with the distance (branch length) from the root'''
	tree.distance_from_root = 0.0	   # 'tree' is the root...
	class V(TreeVisitor):
		def pre_visit_edge(self,src,b,l,dst):
			if l is None: l = 0
			dst.distance_from_root = src.distance_from_root + l
	tree.dfs_traverse(V())




