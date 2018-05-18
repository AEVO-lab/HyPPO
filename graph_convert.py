import sys
import random
import os
from os import listdir
from os.path import isfile, join


inclusters = ""

for arg in sys.argv[1:]:
  if arg.startswith("--inclusters="):
    pz = arg.split("=")
    inclusters = pz[1]
	
if inclusters == "":
	print("This program converts a .clusters file to a .gexf file.  These are used by various graph visualization software such as gephi.  Please specify as clusters file with the --inclusters= argument.  The gexf file is printed directly to the standard output.")
	sys.exit()
	
f = open(inclusters, 'r')


groups = []

curgroup = None
for line in f:
	line = line.replace("\n", "")
	
	
	if line.startswith("<GROUP>"):
		curgroup = []
	elif line.startswith("</GROUP>"):
		groups.append(curgroup)
		curgroup = None
	else:
		if curgroup != None:
			curgroup.append(line)
	
	
f.close()

strout = '<?xml version="1.0" encoding="UTF-8"?>\n'
strout += '<gexf xmlns="http://www.gexf.net/1.2draft" version="1.2">\n'
'''<meta lastmodifieddate="2009-03-20">
	<creator>Gexf.net</creator>
	<description>A hello world! file</description>
</meta>'''
strout += '<graph mode="static" defaultedgetype="directed">\n'
        
nodestr = ''
edgestr = ''


cptnode = 0
cptedge = 0

for g in groups:
	
	
	first_id_in_group = cptnode
	for c in g:
		nodestr += '<node id="' + str(cptnode) + '" label="' + c + '" />\n'
		
		for i in range(first_id_in_group, cptnode):
			edgestr += '<edge id="' + str(cptedge) + '" source="' + str(i) + '" target="' + str(cptnode) + '" />\n'
			cptedge += 1
	
		cptnode += 1
		
		
strout += "<nodes>\n" + nodestr + "</nodes>\n"
strout += "<edges>\n" + edgestr + "</edges>\n"

strout += "</graph>\n"
strout += "</gexf>"

print(strout)
		