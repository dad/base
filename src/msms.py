import sys, os, math, string, random, re
from optparse import OptionParser
import stats, util, biofile, translate, cai
from pysvg.filter import *
from pysvg.gradient import *
from pysvg.linking import *
from pysvg.script import *
from pysvg.shape import *
from pysvg.structure import *
from pysvg.style import *
from pysvg.text import *
from pysvg.builders import *


if __name__=='__main__':
	parser = OptionParser(usage="%prog [options] <evidence filename>")
	parser.add_option("-o", "--out", dest="out_fname", type="string", default=None, help="output filename")
	(options, args) = parser.parse_args()
	
	# Set up some output
	info_outs = util.OutStreams(sys.stdout)
	data_outs = util.OutStreams()
	if not options.out_fname is None:
		outf = file(options.out_fname,'w')
		data_outs.addStream(outf)
	else:
		data_outs.addStream(sys.stdout)

	# Data to take in:
	# masses, intensities, labels, sequence
	
	oh = ShapeBuilder()
	s = svg("test")
  
	s.addElement(oh.createRect(0, 0, 400, 200, 12, 12, strokewidth=2, stroke='navy'))
	s.addElement(oh.createRect(100, 50, 200, 100, strokewidth=2, stroke='navy', fill='yellow'))
	s.addElement(oh.createCircle(700, 500, 50, strokewidth=5, stroke='red'))
	s.addElement(oh.createCircle(810, 500, 50, strokewidth=5, stroke='yellow', fill='#AAAAAA'))
	s.addElement(oh.createEllipse(600, 50, 50, 30, strokewidth=5, stroke='red'))
	s.addElement(oh.createEllipse(700, 50, 50, 30, strokewidth=5, stroke='yellow', fill='#00AABB'))
	s.addElement(oh.createLine(0, 0, 300, 300, strokewidth=2, stroke="black"))
  
	myStyle = StyleBuilder()
	myStyle.setStrokeWidth(2)
	myStyle.setStroke('black')
	l = line(0, 0, 300, 300)
	l.set_style(myStyle.getStyle())
	s.addElement(l)
	#easier method with ShapeBuilder
	oh = ShapeBuilder()
	s.addElement(oh.createLine(10, 10, 300, 300, strokewidth=2, stroke="blue"))
	s.save('./5_Line.svg')
