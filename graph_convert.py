import sys
import random
import os
from os import listdir
from os.path import isfile, join


inclusters = ""
inrelations = ""

for arg in sys.argv[1:]:
  if arg.startswith("--inclusters="):
    inclusters = arg.split("=")[1]
  if arg.startswith("--inrelations="):
    inrelations = arg.split("=")[1]
    
	
if inclusters == "":
	print("This program converts a .clusters file to a .gexf file.  These are used by various graph visualization software such as gephi.  Please specify as clusters file with the --inclusters= argument.  The gexf file is printed directly to the standard output.")
	sys.exit()
	
f = open(inclusters, 'r')



#color_list = ["0000FF", "FF0000", "00FF00", "FFFF00", "00FFFF", "FF00FF", "000000", "CD6155", "AF7AC5", "3498DB", "17A589", "28B463", "F1C40F", "E67E22", "2C3E50", "33FF33"]

color_list = ['222222', 'F3C300', '875692', 'F38400', 'A1CAF1', 'BE0032', 'C2B280', '848482', '008856', 'E68FAC', '0067A5', 'F99379', '604E97', 'F6A600', 'B3446C', 'DCD300', '882D17', '8DB600', '654522', 'E25822', '2B3D26']


groups = []
gene_group = {}
cptgroup = 0

curgroup = None
for line in f:
	line = line.replace("\n", "")
	
	
	if line.startswith("<GROUP>"):
		curgroup = []
	elif line.startswith("</GROUP>"):
		groups.append(curgroup)
		curgroup = None
		cptgroup += 1
	else:
		if curgroup != None:
			curgroup.append(line)
			gene_group[line] = cptgroup
	
	
f.close()


all_orthologs = []
if inrelations != "":
	
	f = open(inrelations, 'r')
	for line in f:
		line = line.replace("\n", "")
		if line.startswith("RELATIONS="):
			rz = line[10:].split(";;")
			for r in rz:
				pz = r.split("\t")
				if pz[2] == "Orthologs":
					all_orthologs.append( [pz[0], pz[1]] )
	f.close()			



strout = '<?xml version="1.0" encoding="UTF-8"?>\n'
strout += '<gexf xmlns="http://www.gexf.net/1.3draft" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:viz="http://www.gexf.net/1.3draft/viz"    xsi:schemaLocation="http://www.gexf.net/1.3draft http://www.gexf.net/1.3draft/gexf.xsd"      version="1.3">\n'
'''<meta lastmodifieddate="2009-03-20">
	<creator>Gexf.net</creator>
	<description>A hello world! file</description>
</meta>'''
strout += '<graph mode="static" defaultedgetype="undirected">\n'
        
nodestr = ''
edgestr = ''


cptnode = 0
cptedge = 0

color_counter = 0

gene_node_ids = {}

for g in groups:
	
	
	first_id_in_group = cptnode
	for c in g:
		nodestr += '<node id="' + str(cptnode) + '" label="' + c + '">\n'
		nodestr += '<viz:color hex="#' + color_list[color_counter] + '" />\n'
		nodestr += '</node>\n'
		
		for i in range(first_id_in_group, cptnode):
			edgestr += '<edge id="' + str(cptedge) + '" weight="1" source="' + str(i) + '" target="' + str(cptnode) + '">\n'
			edgestr += '<viz:color hex="#' + color_list[color_counter] + '" />\n'
			edgestr += '</edge>\n'
			cptedge += 1
		

		gene_node_ids[c] = cptnode
		
		
	
		cptnode += 1
		
	
	color_counter = (color_counter + 1) % len(color_list)


	
for pair in all_orthologs:
	g1 = gene_group[pair[0]]
	g2 = gene_group[pair[1]]
	

	if g1 != g2:
		
			edgestr += '<edge id="' + str(cptedge) + '" weight="0.1" source="' + str(gene_node_ids[pair[0]]) + '" target="' + str(gene_node_ids[pair[1]]) + '">\n'
			edgestr += '<viz:color hex="#000000" />\n'
			edgestr += '</edge>\n'
			cptedge += 1


		
strout += "<nodes>\n" + nodestr + "</nodes>\n"
strout += "<edges>\n" + edgestr + "</edges>\n"

strout += "</graph>\n"
strout += "</gexf>"

print(strout)
		