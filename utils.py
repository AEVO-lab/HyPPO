import sys

#IMPORTANT NOTE
#This script compares the relations output by hyppo with the relations from swisstree.
#swisstree likes to use lots of different formats across families, and the main 
#job of this script is to attempt to handle this mess.
#This script is only intended for usage in the experiments, and is destined to be thrown away after.
#Therefore, do not expect quality code, let alon understandable code.

#arg 1 must be of type hyppo output
#arg 2 must be of type 1 per line (the true relations)
#eg 
#python3 utils.py ./data/swisstree/out/sequences001.fst.mafft.relations  ./data/swisstree/ST001/gene_relationships.txt 
#python3 utils.py ./data/swisstree/out/sequences003.fst.hyppo_sptree.relations  ./data/swisstree/ST003/gene_relationships.txt
#python3 utils.py ./data/swisstree/out_oma/oma.relations  ./data/swisstree/ST001/gene_relationships.txt --shortformat1


args = sys.argv[1:]

f1 = ""
f2 = ""

waitForRelations = True
shortFormat2 = False
synonymsFile = ""
tryToFixNames = False

for a in args:
	if not a.startswith("--"):
		if f1 == "":
			f1 = a
		elif f2 == "":
			f2 = a
	else:
		if a == "--shortformat1":
			waitForRelations = False
		elif a == "--shortformat2":
			shortFormat2 = True
		elif a == "--tryToFixNames":
			tryToFixNames = True
		elif a.startswith("--synonymsfile="):
			synonymsFile = a.split("=")[1]
		

allSynonyms = {}
if synonymsFile != "":
	f = open(synonymsFile, "r")
	for line in f:
		line = line.replace("\n", "")
		elems = []
		pz = line.split()
		for p in pz:
			px = p.split(",")
			for q in px:
				elems.append(q)
				
		for e in elems:
			allSynonyms[e] = elems
				
	#print(allSynonyms)	
	f.close()
	
	


#loads output from hyppo
rels1 = {}
f1_f = open(f1)
for line in f1_f:
	line = line.replace("\n", "")
	
	if line != "":
		if line.startswith("RELATIONS=") or not waitForRelations:
			if line.startswith("RELATIONS="):
				rstr = line.split("=")[1]
			else:
				rstr = line
			
			rz = rstr.split(";;")
			for r in rz:
				if r != "":
					qz = r.split("\t")
					
					
					wz1 = qz[0].split("_")
					name1 = wz1[0] + "_" + wz1[1]
					wz2 = qz[1].split("_")
					name2 = wz2[0] + "_" + wz2[1]
					
					key1 = name1 + ";;" + name2
					key2 = name2 + ";;" + name1
					
					rels1[key1] = qz[2]
					rels1[key2] = qz[2]
			
f1_f.close()





#loads output from one per line
rels2 = {}
f2_f = open(f2)
for line in f2_f:
	line = line.replace("\n", "")
	
	if not shortFormat2:
		if line != "":
			pz = line.split()
			key1 = pz[0].upper() + ";;" + pz[1].upper()
			key2 = pz[1].upper() + ";;" + pz[0].upper()
			rtype = pz[2].replace("ortholog", "Orthologs").replace("paralog", "Paralogs")
			rels2[key1] = rtype
			rels2[key2] = rtype
	else:
		if line != "":
			rz = line.split(";;")
			
			for r in rz:
				if r != "":
					qz = r.split("\t")
					wz1 = qz[0].split("_")
					name1 = wz1[0] + "_" + wz1[1]
					wz2 = qz[1].split("_")
					name2 = wz2[0] + "_" + wz2[1]
					
					key1 = name1 + ";;" + name2
					key2 = name2 + ";;" + name1
					rels2[key1] = qz[2]
					rels2[key2] = qz[2]
f2_f.close()

'''
if tryToFixNames:
	#print("Trying to fix names")
	for key in rels2:
		pz = key.split(";;")
		for i in [0, 1]:
			qz = pz[i].split("_")
			syns = [pz[i]]
			if len(qz) == 3:
				syns.append(qz[0] + "_" + qz[1])
				syns.append(qz[2] + "_" + qz[1])
				#syns.append(qz[0].replace("Hox", "HX").upper() + "_" + qz[1])
				#syns.append(qz[2].replace("Hox", "HX").upper() + "_" + qz[1])
				

				
				
			else:
				pass
				syns.append(qz[0].replace("Hox", "HX").upper() + "_" + qz[1])	
				
			for e in syns:
				if not e in allSynonyms:
					allSynonyms[e] = []
				else:
					print(e + " already there")
					print(allSynonyms[e])
					print(syns)
				for s in syns:
					if s not in allSynonyms[e]:
						allSynonyms[e].append(s)
						
				print(allSynonyms[e])
#print(allSynonyms)
#sys.exit()
'''

nbrels = 0
nbsame = 0

nbSameOrtho = 0		#TP
nbPara1Ortho2 = 0	#FN
nbOrtho1Para2 = 0	#FP

#NOTE: due to crappy swisstree datasets, some relations will be in rels2 but not in rels1.  We just care about those in rels1
for rkey in rels2:

	#check if rkey is missing because of stupid synonyms (names in the dataset are inconsistent)
	r1key = rkey
	syns1 = [rkey.split(";;")[0]]
	syns2 = [rkey.split(";;")[1]]
	if not r1key in rels1:
		genes = rkey.split(";;")
		for elem in allSynonyms:
			if genes[0] in allSynonyms[elem]:
				syns1 = allSynonyms[elem]
			elif genes[1] in allSynonyms[elem]:
				syns2 = allSynonyms[elem]
	
	for e1 in syns1:
		parts1 = e1.split("_")
		possibles1 = [e1, e1.replace("Hox", "HX").upper()]
		if len(parts1) >= 2:
			possibles1.append(parts1[0] + "_" + parts1[1])
		if len(parts1) == 3:
			possibles1.append(parts1[2] + "_" + parts1[1])
			possibles1.append(parts1[0] + "_" + parts1[1])
		
		for p1 in possibles1:
			for e2 in syns2:
				parts2 = e2.split("_")
				
				possibles2 = [e2, e2.replace("Hox", "HX").upper()]
				if len(parts2) >= 2:
					possibles2.append(parts2[0] + "_" + parts2[1])
				
				if len(parts2) == 3:
					possibles2.append(parts2[2] + "_" + parts2[1])
					possibles2.append(parts2[0] + "_" + parts2[1])
				
				for p2 in possibles2:
					
					if p1 + ";;" + p2 in rels1:
						r1key = p1 + ";;" + p2
					

	
	#remove 3rd thing
	'''if not r1key in rels1 and tryToFixNames:
		genes = rkey.split(";;")
		pz1 = genes[0].split("_")
		pz2 = genes[1].split("_")
		name1 = pz1[0] + "_" + pz1[1]
		name2 = pz2[0] + "_" + pz2[1]
		r1key = name1 + ";;" + name2
	'''	
	
				
				
	if r1key in rels1:
		r1 = rels1[r1key]
		r2 = rels2[rkey]
		
		if r1 == r2:
			nbsame += 1
			if r1 == "Orthologs":
				nbSameOrtho += 1
		else:
			if r1 == "Orthologs":
				nbOrtho1Para2 += 1
			elif r1 == "Paralogs":
				nbPara1Ortho2 += 1
		
		nbrels += 1
	else:
		pass
		#print(rkey + " not found in rels1")
		#print(syns1)
		#print(syns2)
	

precision = float(nbSameOrtho)/(float(nbSameOrtho) + float(nbOrtho1Para2))

recall = float(nbSameOrtho)/(float(nbSameOrtho) + float(nbPara1Ortho2))
	
acc = float(nbsame) / float(nbrels)

print("Precision " + str(precision))
print("Recall " + str(recall))
print("Accuracy " + str(nbsame/2) + "/" + str(nbrels/2) + "    " + str(acc))



