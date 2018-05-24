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



color_list = ["0000FF", "FF0000", "00FF00", "FFFF00", "00FFFF", "FF00FF", "000000", "CD6155", "AF7AC5", "3498DB", "17A589", "28B463", "F1C40F", "E67E22", "2C3E50", "33FF33"]

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
strout += '<gexf xmlns="http://www.gexf.net/1.3draft" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:viz="http://www.gexf.net/1.3draft/viz"    xsi:schemaLocation="http://www.gexf.net/1.3draft http://www.gexf.net/1.3draft/gexf.xsd"      version="1.3">\n'
'''<meta lastmodifieddate="2009-03-20">
	<creator>Gexf.net</creator>
	<description>A hello world! file</description>
</meta>'''
strout += '<graph mode="static" defaultedgetype="directed">\n'
        
nodestr = ''
edgestr = ''


cptnode = 0
cptedge = 0

color_counter = 0

for g in groups:
	
	
	first_id_in_group = cptnode
	for c in g:
		nodestr += '<node id="' + str(cptnode) + '" label="' + c + '">\n'
		nodestr += '<viz:color hex="#' + color_list[color_counter] + '" />\n'
		#nodestr += '<viz:color r="255" g="0" b="0" />\n'
		nodestr += '</node>\n'
		
		for i in range(first_id_in_group, cptnode):
			edgestr += '<edge id="' + str(cptedge) + '" source="' + str(i) + '" target="' + str(cptnode) + '">\n'
			edgestr += '<viz:color hex="#' + color_list[color_counter] + '" />\n'
			edgestr += '</edge>\n'
			cptedge += 1
	
		cptnode += 1
		
	
	color_counter = (color_counter + 1) % len(color_list)
		
		
strout += "<nodes>\n" + nodestr + "</nodes>\n"
strout += "<edges>\n" + edgestr + "</edges>\n"

strout += "</graph>\n"
strout += "</gexf>"

print(strout)
		