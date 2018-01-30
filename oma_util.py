import sys
import os
from os.path import isfile, join, exists
from os import listdir
from sequences import Sequence, Distances

#eg python3 fasta_splitter.py --seqfile=./simphy/DL0000005/01/mult/alignment_01_1.5000000000.fas --outdir=./testdir

omadir = ""
outfile_relations = ""
outfile_groups = ""
infasta = ""
mode = "parse_relations"
genes_list = ""
speciestree_newick = ""

for arg in sys.argv[1:]:
	if arg.startswith("--omadir="):
		pz = arg.split("=")
		omadir = pz[1]
		
	if arg.startswith("--infasta="):
		pz = arg.split("=")
		infasta = pz[1]
		
	if arg.startswith("--outfile_relations="):
		pz = arg.split("=")
		outfile_relations = pz[1]
		
	if arg.startswith("--mode="):
		pz = arg.split("=")
		mode = pz[1]
	
	if arg.startswith("--speciestree="):
		pz = arg.split("=")
		speciestree_newick = pz[1]

	if arg.startswith("--outfile_groups="):
		pz = arg.split("=")
		outfile_groups = pz[1]

def handleMode(mode):
	######################################################################################
	# parse mode = get pairwise orthologs inferred by oma, put them in a file
	######################################################################################
	if mode == "parse_relations":
		orthodir = join(omadir, "Output/PairwiseOrthologs")
		list = listdir(orthodir)

		orthologs = []
		ortho_set = set()

		for filename in list:
			
			f = open(join(orthodir, filename))
			
			for line in f:
				line = line.replace("\n", "")
				
				if not line.startswith("#"):
					pz = line.split("\t")
					if len(pz) >= 4:
						g1 = pz[2].replace(" ", "").replace("\r", "")
						g2 = pz[3].replace(" ", "").replace("\r", "")
						orthologs.append([g1, g2])
						ortho_set.add(g1 + "-" + g2)
						ortho_set.add(g2 + "-" + g1)
			f.close()

		trout = ""   	#trout instead of the usual, les funny strout
		for pair in orthologs:
			trout += pair[0] + "\t" + pair[1] + "\tOrthologs;;"
			
		
		for i in range(len(genes_list)):
			for j in range(i + 1, len(genes_list)):
				key = genes_list[i] + "-" + genes_list[j]
				if not key in ortho_set:
					trout += genes_list[i] + "\t" + genes_list[j] + "\tParalogs;;"
		
		trout = trout[0:-2]	#extra ;; at end
		
		if outfile_relations == "":
			print(trout)
		else:
			f = open(outfile_relations, 'w')
			f.write(trout)
			f.close()
	######################################################################################
	# parse mode = get ortholog groups inferred by oma, put them in a file
	######################################################################################
	elif mode == "parse_groups":
		strout = ""
		seen_genes = set()
		f = open( join(omadir, "Output/OrthologousGroups.txt") )
		for line in f:
			line = line.replace("\n", "")
			if line.startswith("OMA"):
				pz = line.split("\t")
				
				strout += "<GROUP>\n"
				for i in range(1, len(pz)): #ignore first item
					gene = pz[i].split(":")[1]
					seen_genes.add(gene)
					strout += gene + "\n"
				strout += "</GROUP>\n"
				
		for g in genes_list:
			if not g in seen_genes:
				strout += "<GROUP>\n" + g + "\n</GROUP>\n"
				
		strout = strout[0:-1]	#extra \n
		
		if outfile_groups == "":
			print(trout)
		else:
			f = open(outfile_groups, 'w')
			f.write(strout)
			f.close()
		
	######################################################################################
	# run mode = setup oma directory and run it
	######################################################################################
	elif mode == "run":
		#first, we'll get the modded species tree in a temp file
		tmpfile = infasta + ".tmpspeciestree"
		
		species_list = ""
		species_dict = {}
		for s in sequences:
			if "__" in s.name:
				
				name = s.name.split("__")[0]
				if not name in species_dict:
					if species_list != "":
						species_list += ","
					species_list += name
					species_dict[name] = 1
					
		cmd = "./sgutils/sgutils -m restrict_species_tree -l \"" + species_list + "\" -s \"" + speciestree_newick + "\" -o \"" + tmpfile + "\""
		print("EXEC " + cmd)
		os.system(cmd)
		
		f = open(tmpfile)
		speciestree_newick_restricted = f.readline().replace("\n", "")
		f.close()
		
		#now, get into oma dir and execute it
		cwd = os.getcwd()
		
		if not os.path.exists(omadir):
			os.mkdir(omadir)
		if not os.path.exists(join(omadir, "DB")):
			os.mkdir(join(omadir, "DB"))
		
		cmd = "python3 fasta_splitter.py --seqfile=" + infasta + " --outdir=" + join(omadir, "DB")
		print("EXEC " + cmd)
		os.system(cmd)
		
		
		os.chdir(omadir)
		
		os.system("rm parameters.drw")
		os.system("OMA -p") #this creates a default config file.
		#we must edit the config to put "DNA", and give the species tree restrcted to the genes at hand.
		outcfg = ""
		f = open("parameters.drw")
		for line in f:
			line = line.replace("\n", "")
			if line.startswith("InputDataType"):
				line = "InputDataType := 'DNA';"
			elif line.startswith("SpeciesTree"):
				line = "SpeciesTree := '" + speciestree_newick_restricted + "';"
			
			outcfg += line + "\n"
		f.close()
		f = open("parameters.drw", 'w')
		f.write(outcfg)
		f.close()
		print("Config edited")
		
		cmd = "OMA -n 5"
		print("EXEC " + cmd)
		os.system(cmd)
		
		os.chdir(cwd)
		
		
 
modes = mode.split(",")

genes_list = []
sequences = None
if infasta != "":
	sequences = Sequence.readSequencesFromFastaFile( infasta )
	for s in sequences:
		genes_list.append(s.name.replace(" ", "").replace("\r", ""))
	
#print(genes_list)
for mode in modes:
	handleMode(mode)
		
		
