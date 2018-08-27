import sys
import random
import os
from os import listdir
from os.path import isfile, join


'''
Does various utilitary functionalities, depending on the mode argument.
Here's what this script can do:

> python hyppo_utils.py --mode=relations_reformat --infile=[original relations file] --outfile=[output file]
Takes an original relations file and converts it to a relations file in a format in which each line is a single relation.

'''

mode = ""

infile = ""
outfile = ""

for arg in sys.argv[1:]:
  if arg.startswith("--mode="):
    mode = arg.split("=")[1]
  if arg.startswith("--infile="):
    infile = arg.split("=")[1]
  if arg.startswith("--outfile="):
    outfile = arg.split("=")[1]
  
	
if mode == "":
	print("Please specify a mode.")
	sys.exit()

	
if mode =="relations_reformat":
	outstr = ""
	f = open(infile, 'r')
	for line in f:
		line = line.replace("\n", "").replace("\r", "")
		
		if line.startswith("RELATIONS="):
			outstr = line.replace("RELATIONS=", "")
			outstr = outstr.replace(";;", "\n")
			outstr = outstr.replace("\t", " ")
			
	f.close()
	
	if outfile == "":
		print(outstr)
	else:
		f = open(outfile, 'w')
		f.write(outstr)
		f.close()
	
