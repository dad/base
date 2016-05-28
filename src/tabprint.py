#!/usr/bin/python

# tabprint.py -- Print a tab-delimited text file into a fixed width format for easy viewing/printing. 
# Used in conjunction with tabview.py, which provides no printing facilities or ability to
# deal with pipes.

# Example usage
# cat foo.txt | tabprint.py | less -S
# tabprint.py foo.txt bar.txt | less -S
# tabprint.py foo.txt bar.txt > foobar.txt; tabview.py foobar.txt

import os, sys, re, string, tabview, fileinput #we don't currently use the tabview library, but one could imagine better integration here
COLUMN_WIDTH = 20
padstr = "".join([" " for nn in range(0,COLUMN_WIDTH+1)]) #make a blank string for padding


#0) Get the input from either stdin or all input streams
if len(sys.argv) == 1:
    fh = sys.stdin
else:
    fh = fileinput.input() #See http://docs.python.org/lib/module-fileinput.html

#1) Iterate over that input and write it in a fixed width format to stdout
def cropstring(xx):
    xlen = len(xx)
    if(xlen > COLUMN_WIDTH):
        return(xx[0:COLUMN_WIDTH]) #crop
    else:
        return(xx + padstr[0:(COLUMN_WIDTH-xlen)]) #pad

for line in fh:
    toks = line.rstrip("\n").split("\t")
    try:
        sys.stdout.write("\t".join([cropstring(xx) for xx in toks]) + "\n")
    except:
        #If pipe is broken, end gracefully
        sys.exit()
    

        





    

