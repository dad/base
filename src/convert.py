#! /usr/bin/python

import sys
import string

f = file(sys.argv[1], 'r')

lines = f.readlines()

for l in lines:
    ls = l.split('\r')
    for i in ls:
        print i

