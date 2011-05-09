#!/usr/bin/python
# -*- coding: iso-8859-1 -*-


from pysvg.filter import *
from pysvg.gradient import *
from pysvg.linking import *
from pysvg.script import *
from pysvg.shape import *
from pysvg.structure import *
from pysvg.style import *
from pysvg.text import *
from pysvg.builders import *


    oh = ShapeBuilder()
    s = svg("test")
  
    s.addElement(oh.createRect(0, 0, 400, 200, 12, 12, strokewidth=2, stroke='navy'))
    s.addElement(oh.createRect(100, 50, 200, 100, strokewidth=2, stroke='navy', fill='yellow'))
    s.addElement(oh.createCircle(700, 500, 50, strokewidth=5, stroke='red'))
    s.addElement(oh.createCircle(810, 500, 50, strokewidth=5, stroke='yellow', fill='#AAAAAA'))
    s.addElement(oh.createEllipse(600, 50, 50, 30, strokewidth=5, stroke='red'))
    s.addElement(oh.createEllipse(700, 50, 50, 30, strokewidth=5, stroke='yellow', fill='#00AABB'))
    s.addElement(oh.createLine(0, 0, 300, 300, strokewidth=2, stroke="black"))
    s.save('./testoutput/4_Shapes.svg')
  
def Line():
    s = svg("test")
    myStyle = StyleBuilder()
    myStyle.setStrokeWidth(2)
    myStyle.setStroke('black')
    l = line(0, 0, 300, 300)
    l.set_style(myStyle.getStyle())
    s.addElement(l)
    #easier method with ShapeBuilder
    oh = ShapeBuilder()
    s.addElement(oh.createLine(10, 0, 300, 300, strokewidth=2, stroke="blue"))
    s.save('./testoutput/5_Line.svg')
