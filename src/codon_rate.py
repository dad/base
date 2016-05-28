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
	#if S==0:
	#	return ( -100000., -100000. )
	( Dn1, Ds1 ) = codon.calcDnDs( croot, cleft )
	#printData( croot, cleft, Dn1, Ds1, N, S )
	( Dn2, Ds2 ) = codon.calcDnDs( croot, cright )
	#printData( croot, cright, Dn2, Ds2, N, S )
	( Dnleft, Dsleft ) = treeDnDsAtCodon( node.left, pos )
	( Dnright, Dsright ) = treeDnDsAtCodon( node.right, pos )
	codon_ds = 0.0
	if S > 0:
		codon_ds = (Ds1 + Ds2)/S
	return ( Dnleft + Dnright + (Dn1 + Dn2)/N, Dsleft + Dsright + codon_ds )

