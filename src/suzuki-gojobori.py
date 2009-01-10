#!/usr/bin/python

import codon

class Node:
	def __init__( self, sequence=None, left=None, right=None ):
		self.sequence = sequence
		self.left = left
		self.right = right

	def printTree( self, level=0 ):
		if self.left != None:
			self.left.printTree( level+1 )
		print "      "*level, self.sequence
		if self.right != None:
			self.right.printTree( level+1 )
	
	def printCodonTree( self, pos, level=0 ):
		if self.left != None:
			self.left.printCodonTree( pos, level+1 )
		print "      "*level, self.sequence[pos:pos+3]
		if self.right != None:
			self.right.printCodonTree( pos, level+1 )
	
	def allCodonsValid( self, pos ):
		if not codon.codonValid( self.sequence[pos:pos+3] ):
			return False
		else:
			valid = True
			if self.left != None:
				valid = valid and self.right.allCodonsValid( pos )
			if valid:
				if self.right != None:
					return valid and self.left.allCodonsValid( pos )
				else:
					return True
			else:
				return False


def printData( c1, c2, Dn, Ds, N, S ):
	"For debugging purposes only"
	if debug:
		print "%s %s Dn=%g, Ds=%g, N=%g, S=%g, Dn/N=%g, Ds/S=%g" % ( c1, c2, Dn, Ds, N, S, Dn/N, Ds/S )

def treeDnDsAtCodon( node, pos ):
	"Calculate Dn and Ds over entire tree, at position 'pos' only. "
	if node.left == None or node.right == None:
		return ( 0., 0. )
	croot = node.sequence[pos:pos+3]
	cleft = node.left.sequence[pos:pos+3]
	cright = node.right.sequence[pos:pos+3]
	( N, S ) = codon.calcNS( croot )
	if S==0:
		return ( -100000., -100000. )
	( Dn1, Ds1 ) = codon.calcDnDs( croot, cleft )
#	printData( croot, cleft, Dn1, Ds1, N, S )
	( Dn2, Ds2 ) = codon.calcDnDs( croot, cright )
#	printData( croot, cright, Dn2, Ds2, N, S )
	( Dnleft, Dsleft ) = treeDnDsAtCodon( node.left, pos )
	( Dnright, Dsright ) = treeDnDsAtCodon( node.right, pos )
	return ( Dnleft + Dnright + (Dn1 + Dn2)/N, Dsleft + Dsright + (Ds1 + Ds2)/S )

#
# All code below this line is for demonstration purposes only
#
########################################

from random import *
from array import *

bases = ( 'A', 'C', 'T', 'G' )


def mutate( s, n ):
	l = len( s )
	snew = array( 'c', s )
	for i in range( n ):
		snew[randint(0,l-1)] = bases[randint(0,3)]
	return snew.tostring()


def main():
	n = 10000
	m = 30000
	s1 = ''
	for i in range( 3*n ):
		s1 += bases[randint(0,3)]
	# create mutated sequences
	s2 = mutate( s1, m )
	s3 = mutate( s1, m )
	s4 = mutate( s2, m )
	s5 = mutate( s2, m )
	s6 = mutate( s3, m )
	s7 = mutate( s3, m )
	s8 = mutate( s6, m )
	s9 = mutate( s6, m )
	# add some gaps to test removal of sites with gaps
	s1 = s1[0:1] + '-' + s1[5:3*n]
	s2 = s2[0:4] + '-' + s2[8:3*n]
	s3 = s3[0:7] + '-' + s3[11:3*n]
	s4 = s4[0:10] + '-' + s4[14:3*n]
	s5 = s5[0:13] + '-' + s5[17:3*n]
	s6 = s6[0:16] + '-' + s6[20:3*n]
	s7 = s7[0:19] + '-' + s7[23:3*n]
	s8 = s8[0:22] + '-' + s8[26:3*n]
	s9 = s9[0:25] + '-' + s9[29:3*n]
	tree = Node( s1, Node( s2, Node( s4 ), Node( s5 ) ), Node( s3, Node( s6, Node( s8 ), Node( s9 ) ), Node( s7 ) ) )
#	tree.printTree()
	print "site\tDn/N\tDs/S"
	for i in range( n ):
#		tree.printCodonTree( 3*i )
		valid = tree.allCodonsValid( 3*i )
#		print valid
		if valid:
			( Dn, Ds ) = treeDnDsAtCodon( tree, 3*i )
			if Dn > 0:
				print "%i\t%g\t%g\t" % ( i, Dn, Ds )

main()

