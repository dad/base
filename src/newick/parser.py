'''
A Python module for parsing Newick files.

Copyright (C) 2003-2006, Thomas Mailund <mailund@birc.au.dk>

This module contains the functionality for grammatical analysis. '''

import tokens

class ParserError:
	'''Exception thrown if the parser encounters an error.'''
	def __init__(self,err):
		self.err = err

	def __repr__(self):
		return "ParserError: "+self.err


class AbstractHandler:
	'''Interface (and NO-OP implementations) of event handlers for
	parsing trees.  A handler can be used for extracting information
	from a tree without explicitly representing the tree in a
	datastructure.'''
   
	def newTreeBegin(self):
		'''Callback called when the parsing begins.'''
		pass

	def newTreeEnd(self):
		'''Callback called when the tree is completely parsed.'''
		pass

	def newLeaf(self, name):
		'''Callback called when a leaf is parsed.  A name is always
		provided, although it can be the empty string if an identifier
		was not explicitly given in the input.'''
		pass

	def newChild(self, name):
		'''Callback called when a new child is parsed.  A name is always
		provided, although it can be the empty string if an identifier
		was not explicitly given in the input.'''
		pass

class _Parser:
	'''State of the parser during parsing.  Should not be used
	directly by users of this package.'''

	def __init__(self, lexer, handler):
		self.lexer = lexer
		self.handler = handler

	def parse(self):
		result = None
		if self.lexer.peek_token(tokens.LParen):
			result = self.parseNode()
		else:
			result = self.parseChild()

		remaining = self.lexer.remaining()
		if remaining != '' and not self.lexer.peek_token(tokens.SemiColon):
			raise ParserError("Unexpected token following tree: " +
							  self.lexer.remaining())
		return result

	def parseNode(self):
		#print "parse_node", self.lexer.remaining()
		''' Parse node on the form ( <edge list> ) '''
		self.handler.newTreeBegin()
		self.lexer.read_token(tokens.LParen)
		self.parseChildren()
		self.lexer.read_token(tokens.RParen)
		# Read bootstrap and length data for this node, if present
		bootstrap = None
		length = None
		if self.lexer.peek_token(tokens.Number):
			bootstrap = int(self.lexer.read_token(tokens.Number).get_number())
		if self.lexer.peek_token(tokens.Colon):
			self.lexer.read_token(tokens.Colon)
			length = self.lexer.read_token(tokens.Number).get_number()
		self.handler.newTreeEnd(bootstrap, length)
		#self.handler.finishNode(bootstrap,length)


	def parseChildren(self):
		#print "parse_children", self.lexer.remaining()
		''' parse a comma-separated list of children. '''
		while True:
			self.parseChild()
			# Break out when the next character isn't a comma
			if self.lexer.peek_token(tokens.Comma):
				self.lexer.read_token(tokens.Comma)
			else:
				break

	def parseChild(self):
		#print "parseChild", self.lexer.remaining()
		# If next character indicates a subtree, parse the subtree
		if self.lexer.peek_token(tokens.LParen):
			self.parseNode()
		else:
			# Otherwise, parse the leaf
			self.parseLeaf()

	def parseLeaf(self):
		#print "parse_leaf", self.lexer.remaining()
		''' Parse a node of the form "id:length" '''
		if self.lexer.peek_token(tokens.Comma) or \
			   self.lexer.peek_token(tokens.RParen):
			# blank name
			self.handler.newLeaf("",None)
			return

		# special case for when the id is just a number
		'''
		if self.lexer.peek_token(tokens.Number):
			id = str(int(self.lexer.read_token(tokens.Number).get_number()))
			self.handler.new_leaf(id)
			return
		'''
		id = self.lexer.read_token(tokens.ID).get_name()

		length = None
		if self.lexer.peek_token(tokens.Colon):
			self.lexer.read_token(tokens.Colon)
			length = self.lexer.read_token(tokens.Number).get_number()

		name = id
		if id == '_':
			# blank name
			name = ''
		self.handler.newLeaf(name,length)

def parse(input, event_handler):
	'''Parse input and invoke callbacks in event_handler.  If
	event_handler implements a getResult() method, parse will return
	the result of calling this after complete parsing, otherwise None
	is returned.'''
	import lexer
	l = lexer.Lexer(input)
	_Parser(l,event_handler).parse()
	if hasattr(event_handler,"getResult"):
		return event_handler.getResult()
	

    
