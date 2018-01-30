from sequences import Sequence, Distances
import sys
import random
import hyppo_classes
import os
import numpy
import math
from os import listdir
from os.path import isfile, join



def output_family_scores_to_file(outfile, format, gene_family, overwrite = True):
		
	if isfile(outfile) and not overwrite:
		return
	
	if format == "matrix":
		matrixstr = Distances.getMatrixString(gene_family.sequences, gene_family.scores)
		with open( outfile, 'w') as file:
			file.write(matrixstr)
			
	elif format == "edgelist":
		outstr = ""
		for i in range(len(gene_family.sequences)):
			for j in range(i + 1, len(gene_family.sequences)):
				#if gene_family.scores[i][j] > 0:
				outstr += gene_family.sequences[i].name + " " + gene_family.sequences[j].name + " " + str(gene_family.scores[i][j]) + "\n"
		
		with open( outfile, 'w') as file:
			file.write(outstr)
		


def get_clusters_from_file(file, name_suffix = ""):
	clusters = []
	f = open(file, 'r')
	curgroup = None
	for line in f:
		line = line.replace("\n", "")
		if line.startswith("<GROUP>"):
			curgroup = set()
		elif line.startswith("</GROUP>"):
			clusters.append(curgroup)
		elif line != "" and curgroup != None:
			curgroup.add(line + name_suffix)
	f.close()
	
	return clusters

	
def write_bbh_distances(outfilename, sequences, scores, species_separator = "__", species_index = 1):
	
	species = set()
	for g in sequences:
		sp = g.name.split(species_separator)[species_index]
		if not sp in species:
			species.add(sp)
	
	favBySpecies = []
	for i in range(len(sequences)):
		isp = sequences[i].name.split(species_separator)[species_index]
		
		favBySpecies[i] = {}
		
		for j in range(len(sequences)):
			jsp = sequences[i].name.split(species_separator)[species_index]
				
			if jsp not in favBySpecies[i]:
				favBySpecies[i][jsp] = j
			else:
				curmax = favBySpecies[i][jsp]
				if scores[i][j] > curmax:
					favBySpecies[i][jsp] = j
	
	
	
	
def write_clusters(file, clusters):
	outstr = ""
	for c in clusters:
		outstr += "<GROUP>\n"
		for g in c:
			outstr += g + "\n"
		outstr += "</GROUP>\n"
	f = open(file, 'w')
	f.write(outstr)
	f.close()

	

	
class GeneFamily:
	

	sequences = []
	scores = None
	clusters_filename = ""
	seq_file = ""
	workdir = ""
	
	def load_from_alignment_file(self, filename):
		self.seq_file = filename
	
		self.sequences = Sequence.readSequences( filename )
		

		
class HyPPO_FullPredictor:

	clusters_filenames = []
	relations_filenames = []
	other_args = {}

	def get_clusters_filenames(self):
		return clusters_filenames
		
	def get_relations_filenames(self):
		return relations_filenames
	
	def predict_orthologs(self, files, workdir, speciestree_file, other_args):
		pass
		
class HyPPO_Scores:
	
	other_args = {}

	def get_scores_filename(self, format, gene_family):
		pass 
	
	def load_scores(self, gene_family):
		pass
		
		
		
		
class HyPPO_ClusterPredictor:
	
	all_gene_families = None
	other_args = {}
		
	def get_clusters_filename(self, gene_family, no_iter, speciestree_filename = None):
		pass
	
	def infer_clusters(self, hyppo_scores, gene_family, no_iter, speciestree_filename = None):
		pass
		
		
		
class HyPPO_SpeciesTreeMaker:
	
	other_args = {}
	
	def get_speciestree_filename(self, gene_families, no_iter):
		pass
		
	def infer_species_tree(self, hyppo_scores, gene_families, no_iter):
		pass
		
		
class HyPPO_InterClusterPredictor:

	other_args = {}

	def get_relations_filename(self, gene_family, speciestree_filename = None):
		pass
		
	def infer_intercluster_relations(self, hyppo_scores, gene_family, speciestree_filename = None):
		pass
		

		
		
#TODO: next 5 classes are just copy pasta...
class MusclePctID(HyPPO_Scores):
	
	sequences = []
	
	def get_scores_filename(self, format, gene_family):
		return join(gene_family.workdir, os.path.basename(gene_family.seq_file) + ".musclepctid." + format)
	
	def load_scores(self, gene_family):
		
		dealignedfile = gene_family.seq_file + ".dealigned.fa"
		Sequence.outputSequencesToFasta(gene_family.sequences, dealignedfile)
		
		
		musclefile = gene_family.seq_file + ".muscle.fa"
		
		cmd = "muscle -in " + dealignedfile + " -out " + musclefile
 
		print("EXEC " + cmd)
		
		if not os.path.isfile(musclefile) or not "use_cache" in self.other_args or os.path.getsize(musclefile) <= 10:
			os.system(cmd)
		else:
			print("Actually, file exists and we'll use it.")
		
		seqs = Sequence.readSequences(musclefile)
		
		
		gene_family.scores = Distances.getPairwisePctID(sequences = seqs, verbose = False, run_nw_algorithm = False)

		
		
class TCoffeePctID(HyPPO_Scores):
	
	sequences = []
	
	def get_scores_filename(self, format, gene_family):
		return join(gene_family.workdir, os.path.basename(gene_family.seq_file) + ".tcoffeepctid." + format)
	
	def load_scores(self, gene_family):
		
		dealignedfile = gene_family.seq_file + ".dealigned.fa"
		Sequence.outputSequencesToFasta(gene_family.sequences, dealignedfile)
		
		
		tcoffeefile = gene_family.seq_file + ".tcoffee.fa"
		
		cmd = "t_coffee -in " + dealignedfile + " -output fasta_aln -outfile " + tcoffeefile
 
		print("EXEC " + cmd)
		
		if not os.path.isfile(tcoffeefile) or not "use_cache" in self.other_args or os.path.getsize(tcoffeefile) <= 10:
			os.system(cmd)
		else:
			print("Actually, file exists and we'll use it.")
		
		seqs = Sequence.readSequences(tcoffeefile)
		
		
		gene_family.scores = Distances.getPairwisePctID(sequences = seqs, verbose = False, run_nw_algorithm = False)
		
		
		
class FSAPctID(HyPPO_Scores):
	
	sequences = []
	
	def get_scores_filename(self, format, gene_family):
		return join(gene_family.workdir, os.path.basename(gene_family.seq_file) + ".fsapctid." + format)
	
	def load_scores(self, gene_family):
		
		dealignedfile = gene_family.seq_file + ".dealigned.fa"
		Sequence.outputSequencesToFasta(gene_family.sequences, dealignedfile)
		
		
		fsafile = gene_family.seq_file + ".fsa.fa"
		
		cmd = "fsa " + dealignedfile  + " > " + fsafile
 
		print("EXEC " + cmd)
		
		if not os.path.isfile(fsafile) or not "use_cache" in self.other_args or os.path.getsize(fsafile) <= 10:
			os.system(cmd)
		else:
			print("Actually, file exists and we'll use it.")
		
		seqs = Sequence.readSequences(fsafile)
		
		
		gene_family.scores = Distances.getPairwisePctID(sequences = seqs, verbose = False, run_nw_algorithm = False)
		
		
class MAFFTPctID(HyPPO_Scores):
	
	sequences = []
	
	def get_scores_filename(self, format, gene_family):
		return join(gene_family.workdir, os.path.basename(gene_family.seq_file) + ".mafftpctid." + format)
	
	def load_scores(self, gene_family):
		
		dealignedfile = gene_family.seq_file + ".dealigned.fa"
		Sequence.outputSequencesToFasta(gene_family.sequences, dealignedfile)
		
		
		mafftfile = gene_family.seq_file + ".mafft.fa"
		
		cmd = "mafft --localpair --maxiterate 5000 " + dealignedfile + " > " + mafftfile
 
		print("EXEC " + cmd)
		
		if not os.path.isfile(mafftfile) or not "use_cache" in self.other_args or os.path.getsize(mafftfile) <= 10:
			os.system(cmd)
		else:
			print("Actually, file exists and we'll use it.")
		
		seqs = Sequence.readSequences(mafftfile)
		
		
		gene_family.scores = Distances.getPairwisePctID(sequences = seqs, verbose = False, run_nw_algorithm = False)
		
		

class ClustalPctID(HyPPO_Scores):
	
	sequences = []
	
	def get_scores_filename(self, format, gene_family):
		return join(gene_family.workdir, os.path.basename(gene_family.seq_file) + ".clupctid." + format)
	
	def load_scores(self, gene_family):
		
		scores = [[0 for col in range(len(gene_family.sequences))] for row in range(len(gene_family.sequences))]
		
		#cluster won't work if nbseqs <= 2.  In this case, we'll just throw out some scores
		if len(gene_family.sequences) == 1:
			scores[0][0] = 1
			gene_family.scores = scores
			return
		elif len(gene_family.sequences) == 2:
			scores[0][0] = 1
			scores[0][1] = 1
			scores[1][0] = 1
			scores[1][1] = 1
			gene_family.scores = scores
			return
			
		cldealignedfile = gene_family.seq_file + ".dealigned.fa"
		Sequence.outputSequencesToFasta(gene_family.sequences, cldealignedfile)
		
		
		clfile = gene_family.seq_file + ".clustal.fa"
		pfile = gene_family.seq_file + ".clustal.pctid"
		cmd = "clustalo -i " + cldealignedfile + " -o " + clfile + " --distmat-out=" + pfile + " --percent-id --full"

		print("EXEC " + cmd)
		if not os.path.isfile(pfile) or not "use_cache" in self.other_args:
			os.system(cmd)
		else:
			print("Actually, file exists and we'll use it.")
		
		d = Distances.readMatrixFile(pfile, "count", " ")
		genes = d["labels"]
		matrix = d["matrix"]
		
		genes_pos = {}
		for i in range(len(genes)):
			genes_pos[genes[i]] = i
		
		
		
		for i in range(len(gene_family.sequences)):
			for j in range(len(gene_family.sequences)):
				scores[i][j] = matrix[ genes_pos[gene_family.sequences[i].name] ] [  genes_pos[gene_family.sequences[j].name] ]
		
		gene_family.scores = scores
		
		
class ScoresPctID(HyPPO_Scores):
	
	sequences = []
	
	def get_scores_filename(self, format, gene_family):
		return join(gene_family.workdir, os.path.basename(gene_family.seq_file) + ".pctid." + format)
	
	def load_scores(self, gene_family):
		
		compute_pctid = False
		if "compute_pctid" in self.other_args and self.other_args["compute_pctid"] == "1":
			print("Computing pct id")
			compute_pctid = True
		
		the_len = len(gene_family.sequences[0].alignedSeq)
		for s in gene_family.sequences:
			if len(s.alignedSeq) != the_len:
				compute_pctid = True
				print("Sequences in " + gene_family.seq_file + " differ in length.  We will use the pct id as our scores.")
				break
		
		use_cache = False
		if "use_cache" in self.other_args:
			use_cache = True
		
		if not compute_pctid:
			gene_family.scores = Distances.getPairwisePctID(sequences = gene_family.sequences, verbose = False, run_nw_algorithm = False)
		else:
			tmpfile = gene_family.seq_file + ".ocredgelist"
			
			if not os.path.isfile(tmpfile) or not use_cache:
				cmd = "./OCR/OCR -m pctid -f " + gene_family.seq_file + " -o " + tmpfile
				print("EXEC " + cmd)
				os.system(cmd)
			else:
				print("Using cached file " + tmpfile)
			
			gene_family.scores = Distances.getScoresFromEdgeList(gene_family.sequences, tmpfile)
			

			
			
class MaxScoreClusterPredictor(HyPPO_ClusterPredictor):
	
	def get_clusters_filename(self, gene_family, no_iter, speciestree_filename = None):
		return join(gene_family.workdir, os.path.basename(gene_family.seq_file) + ".maxscore_" + str(no_iter) + ".clusters")
	
	def infer_clusters(self, hyppo_scores, gene_family, no_iter, speciestree_filename = None):
		
		scores_file = hyppo_scores.get_scores_filename("matrix", gene_family)
		output_family_scores_to_file(scores_file, "matrix", gene_family, True)
		
		outfile = self.get_clusters_filename(gene_family, no_iter, speciestree_filename)
		#call greedy on scores_file
		cmd = "python3 ./greedy_clusterer.py --infile=" + scores_file + " --outfile=" + outfile + " --species_separator=\"" + self.other_args["species_separator"] + "\" --species_index=" + self.other_args["species_index"]
		print("EXEC " + cmd)
		os.system(cmd)
		
		
		#gene_family.clusters_filename = outfile
		
		
		
class GivenClusterPredictor(HyPPO_ClusterPredictor):
	
	
	def get_given_cluster_filename(self, gene_family):
		return other_args[gene_family.seq_file]
	
	def get_clusters_filename(self, gene_family, no_iter, speciestree_filename = None):
		return join(gene_family.workdir, os.path.basename(gene_family.seq_file) + ".given.clusters")
		
	def infer_clusters(self, hyppo_scores, gene_family, no_iter, speciestree_filename = None):
		
		filename = self.get_given_cluster_filename(gene_family, other_args)
		
		suffix = other_args[gene_family.seq_file + "_suffix"]
		clusters = get_clusters_from_file(filename, suffix)
		write_clusters(get_clusters_filename(gene_family, no_iter, speciestree_filename, other_args))
		
		#gene_family.clusters_filename = filename
		
		
class BottomUpClusterPredictor(HyPPO_ClusterPredictor):
	
	def get_clusters_filename(self, gene_family, no_iter, speciestree_filename = None):
		return join(gene_family.workdir, os.path.basename(gene_family.seq_file) + ".bottomup_" + str(no_iter) + ".clusters")
	
	def infer_clusters(self, hyppo_scores, gene_family, no_iter, speciestree_filename = None):
		
		if speciestree_filename == None:
			print("ERROR: to use BottomUpClusterPredictor, you need a species tree.  Program will exit.")
			sys.exit()
		else:
		
			scores_file = hyppo_scores.get_scores_filename("edgelist", gene_family)
			output_family_scores_to_file(scores_file, "edgelist", gene_family, True)
			
			outfile = self.get_clusters_filename(gene_family, no_iter, speciestree_filename)
			#call OCR on scores_file + sp tree
			cmd = "./OCR/OCR -m find_clusters -w " + scores_file + " -s " + speciestree_filename + " -o " + outfile
			cmd += " -spsep \"" + self.other_args["species_separator"] + "\" -spindex " + self.other_args["species_index"]
			
			print("EXEC " + cmd)
			os.system(cmd)
		
			#gene_family.clusters_filename = outfile
		
		
class BottomUpSpeciesTreeMaker(HyPPO_SpeciesTreeMaker):
	
	def get_speciestree_filename(self, gene_families, no_iter):
		strout = join(gene_families[0].workdir, "speciestree_bottomup_" + str(no_iter))
		if "mode_string" in self.other_args:
			strout += "." + self.other_args["mode_string"]
		strout += ".nw"
		return strout
	
	def infer_species_tree(self, hyppo_scores, gene_families, no_iter):
		
		files = ""
		for g in gene_families:
			if files != "":
				files += ";;"
			files += g.clusters_filename
		
		outfile = self.get_speciestree_filename(gene_families, no_iter)
		
		#if os.path.isfile(outfile) and "use_species_tree_cache" in self.other_args:
		#	print("Speciess tree file exists.  Will use.")
		#	return outfile
		
		#call OCR on these clusters
		cmd = "./OCR/OCR -c \"" + files + "\" -m species_tree -o " + outfile	
		cmd += " -spsep \"" + self.other_args["species_separator"] + "\" -spindex " + self.other_args["species_index"]
		#cmd += " -v"  
		print("EXEC " + cmd)
		os.system(cmd)
		
		return outfile



		
		
		
class GreedyInterClusterPredictor(HyPPO_InterClusterPredictor):

	def get_relations_filename(self, gene_family, speciestree_filename = None):
		return join(gene_family.workdir, os.path.basename(gene_family.seq_file) + ".greedy.relations")
		
	def infer_intercluster_relations(self, hyppo_scores, gene_family, speciestree_filename = None):
		clusters_filename = gene_family.clusters_filename
		outfile = self.get_relations_filename(gene_family, speciestree_filename)
		
		scores_file = hyppo_scores.get_scores_filename("edgelist", gene_family)
		output_family_scores_to_file(scores_file, "edgelist", gene_family, True)
		
		#call OCR 					
		cmd = "./OCR/OCR -s " + speciestree_filename + " -c " + clusters_filename
		cmd += " -w " + scores_file
		cmd += " -o " + outfile
		cmd += " -p0"
		cmd += " -spsep \"" + self.other_args["species_separator"] + "\" -spindex " + self.other_args["species_index"]
		#cmd += " -v"
		
		print("EXEC " + cmd)
		os.system(cmd)
		
		return outfile
		
		
		
		
		
class OMAPredictor(HyPPO_FullPredictor):

	clusters_filenames = ""
	relations_filenames = ""

	def predict_orthologs(self, files, workdir, speciestree_file):
		
		runOMA = True 
		
		if os.path.isfile(join(workdir, "Output/OrthologousGroups.txt")) and "use_cache" in self.other_args:
			print("Will use cached output files from oma")
			runOMA = False
		
		dbdir = join(workdir, "DB")
		if not os.path.exists(dbdir):
			os.mkdir(dbdir)
		
		seqs_by_species = Sequence.combineSequenceFilesBySpecies(files, self.other_args["species_separator"], int(self.other_args["species_index"]))
		
		if speciestree_file == "":
			print("OMA is unpredicatable without a species tree.")

		spprefix = ""
		if "species_prefix" in self.other_args:
			spprefix = self.other_args["species_prefix"]
		
		doAA = False
		if "convertToAA" in self.other_args:
			doAA = True
		
		species_list = ""
		genes_list = []
		for key in seqs_by_species:
			outfile = join(dbdir, spprefix + key + ".fa")
			Sequence.outputSequencesToFasta(sequences = seqs_by_species[key], filename = outfile, name_suffix = "", aligned = False, convertToAA = doAA, name_prefix = spprefix)
			
			if species_list != "":
				species_list += ","
			species_list += key
			
			for seq in seqs_by_species[key]:
				genes_list.append(seq.name)
			
			
		#if we use a species tree, OMA requires removing species not appearing in the files
		#sgutils can do that
		if speciestree_file != "":
			oma_sptree_file = join(workdir, "oma_species_tree.nw")
			
			with open(speciestree_file, 'r') as myfile:
				speciestree_newick = myfile.read().replace('\n', '')
			
			cmd = "./sgutils/sgutils -m restrict_species_tree -l \"" + species_list + "\" -s \"" + speciestree_newick + "\" -o \"" + oma_sptree_file + "\""
			print("EXEC " + cmd)
			os.system(cmd)
		
			f = open(oma_sptree_file)
			speciestree_newick_restricted = f.readline().replace("\n", "")
			f.close()
		
			################################################################
			#special case here: if only one species present, oma makes an error.
			#In this case, we make 1 gene = 1 cluster
			if "," not in speciestree_newick_restricted:
				print("Only one species found.  Will make all genes paralogs.")
				self.clusters_filenames = [join(workdir, "oma.clusters")]
				self.relations_filenames = [join(workdir, "oma.relations")]
				
				clusters = []
				for g in genes_list:
					clusters.append([g])
				write_clusters(self.clusters_filenames[0], clusters)
				
				f = open(self.relations_filenames[0], 'w')
				tmp = ""
				for a in range(len(genes_list)):
					for b in range(a + 1, len(genes_list)):
						tmp += genes_list[a] + "\t" + genes_list[b] + "\t" + "Paralogs;;"
				tmp = tmp[0:-2]
				f.write(tmp)
				f.close()
				
				return
			################################################################
		
		
		input_type = "DNA"
		if "seqtype" in self.other_args and self.other_args["seqtype"] == "AA":
			print("Input type set to AA")
			input_type = self.other_args["seqtype"]
		
		#now, get into oma dir and execute it
		cwd = os.getcwd()
		
		
		os.chdir(workdir)
		
		os.system("rm parameters.drw")
		os.system("OMA -p") #this creates a default config file.
		#we must edit the config to put "DNA", and give the species tree restrcted to the genes at hand.
		outcfg = ""
		f = open("parameters.drw")
		for line in f:
			line = line.replace("\n", "")
			if line.startswith("InputDataType"):
				line = "InputDataType := '" + input_type + "';"
			elif line.startswith("SpeciesTree") and speciestree_file != "":
				line = "SpeciesTree := '" + speciestree_newick_restricted + "';"
			
			outcfg += line + "\n"
		f.close()
		f = open("parameters.drw", 'w')
		f.write(outcfg)
		f.close()
		print("Config edited")
		
		cmd = "OMA -n 7"
		print("EXEC " + cmd)
		
		if runOMA:
			os.system(cmd)
		else:
			print("Not really - using cache instead.")
		
		os.chdir(cwd)
		
		
		
		self.parse_groups(genes_list, workdir)
		self.parse_relations(genes_list, workdir)
	
		self.clusters_filenames = [join(workdir, "oma.clusters")]
		self.relations_filenames = [join(workdir, "oma.relations")]
	
	def parse_groups(self, genes_list, workdir):
		strout = ""
		seen_genes = set()
		f = open( join(workdir, "Output/OrthologousGroups.txt") )
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
		
		
		f = open(join(workdir, "oma.clusters"), 'w')
		f.write(strout)
		f.close()

	def parse_relations(self, genes_list, workdir):
		orthodir = join(workdir, "Output/PairwiseOrthologs")
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

		trout = ""   	#trout instead of the usual, less funny strout
		
		for pair in orthologs:
			trout += pair[0] + "\t" + pair[1] + "\tOrthologs;;"
			
		#if we missed some relations, mark them as Paralogs
		for i in range(len(genes_list)):
			for j in range(i + 1, len(genes_list)):
				key = genes_list[i] + "-" + genes_list[j]
				if not key in ortho_set:
					trout += genes_list[i] + "\t" + genes_list[j] + "\tParalogs;;"
					#print(key + " marked as Unknown")
		trout = trout[0:-2]	#extra ;; at end
		
		f = open(join(workdir, "oma.relations"), 'w')
		f.write(trout)
		f.close()
		
		
		
class OMAClusterPredictor(HyPPO_ClusterPredictor):
	
	def get_clusters_filename(self, gene_family, no_iter, speciestree_filename = None):
		return join(gene_family.workdir, os.path.basename(gene_family.seq_file) + ".omaclp_" + str(no_iter) + ".clusters")
	
	def infer_clusters(self, hyppo_scores, gene_family, no_iter, speciestree_filename = None):
		
		oma = OMAPredictor()
		workdir_tmp = join(gene_family.workdir, "oma_tmp")
		oma.predict_orthologs([gene_family.seq_file], workdir_tmp, speciestree_filename, None)
		
		tmp_clusters_file = oma.clusters_filenames
		clusters = get_clusters_from_file(tmp_clusters_file)
		g.clusters = clusters
		
		#outfile = self.get_clusters_filename(gene_family, no_iter, speciestree_filename)

		
		
		
		
#as of now, only works on Manuel's computer.  Please harass me at mlafond2@uOttawa.ca if you want to use this. 
#It does assume that you have everything setup to run orthomcl, and that you are using the orthomcl-pipeline utility.
class OrthoMCLPredictor(HyPPO_FullPredictor):

	clusters_filenames = ""
	relations_filenames = ""
	cached_clusters = []

	def predict_orthologs(self, files, workdir, speciestree_file):
		
		workdir_seqs = join(workdir, "in")
		workdir_out = join(workdir, "out")
		
		allseqs = []
		for f in files:
			seqs = Sequence.readSequences( f )
			for s in seqs:
				allseqs.append(s)
		
		if len(self.cached_clusters) == 0:
			seqs_by_species = Sequence.combineSequenceFilesBySpecies(files, self.other_args["species_separator"], int(self.other_args["species_index"]))
			
			isdna = True
			if "seqtype" in self.other_args and self.other_args["seqtype"] == "AA":
				print("Input type set to AA")
				isdna = False
			
			if not os.path.exists(workdir_seqs):
				os.mkdir(workdir_seqs)
			for key in seqs_by_species:
				outfile = join(workdir_seqs, key + ".fasta")
				Sequence.outputSequencesToFasta(sequences = seqs_by_species[key], filename = outfile, name_suffix = "", aligned = False, convertToAA = isdna, name_prefix = "")
		
			cmd = "/u/lafonman/src/orthomcl-pipeline/bin/orthomcl-pipeline -i " + workdir_seqs + " -o " + workdir_out + " -m /u/lafonman/src/orthomcl-pipeline/orthomcl.conf --nocompliant --yes"
			print("EXEC " + cmd)
			os.system(cmd)
		
			seen_genes = set()
			clusters = []
			clfile = join(workdir_out, "groups/groups.txt")
			f = open(clfile, 'r')
			for line in f:
				line = line.replace("\n", "")
				if line != "":
					gz = line.split(":")[1].split()
					cluster = set()
					for g in gz:
						gname = g.split("|")[1]
						cluster.add(gname)
						seen_genes.add(gname)
					clusters.append(cluster)
			
			f.close()
			
			for s in allseqs:
				name = s.name
				if not name in seen_genes:
					cl = set()
					cl.add(name)
					clusters.append(cl)
					
			self.cached_clusters = clusters
		else:
			print("USING CACHED CLUSTERS")
			clusters = self.cached_clusters
		
		
		#restrict clusters to current family
		#f_set = set()
		#for s in gene_family.sequences:
		#	f_set.add(s.name)
		#f_clusters = []
		#for cl in clusters:
		#	inter = f_set.intersection(cl)
		#	if len(inter) > 0:
		#		f_clusters.append(inter)
		
			
		self.clusters_filenames = [join(workdir, "orthomcl.clusters")]
		self.relations_filenames = [join(workdir, "orthomcl.relations")]
		
		write_clusters(self.clusters_filenames[0], clusters)
		
		#output relations
		relstr = ""
		seen_keys = {}
		for c in clusters:
			for c1 in c:
				for c2 in c:
					if c1 != c2:
					
						key1 = c1 + ";;" + c2
						key2 = c2 + ";;" + c1
						
						if key1 not in seen_keys and key2 not in seen_keys:
							seen_keys[key1] = 1
							seen_keys[key2] = 1
							relstr += c1 + "\t" + c2 + "\t"
							
							sp1 = c1.split(self.other_args["species_separator"])[int(self.other_args["species_index"])]
							sp2 = c2.split(self.other_args["species_separator"])[int(self.other_args["species_index"])]
							
							if sp1 != sp2:
								relstr += "Orthologs"
							else:
								relstr += "Paralogs"
							
							relstr += ";;"
		relstr = relstr[0:-2]	#extra ;; at end
		
		for i in range(len(allseqs)):
			for j in range(i + 1, len(allseqs)):
				n1 = allseqs[i].name
				n2 = allseqs[j].name
				key = n1 + ";;" + n2
				if not key in seen_keys:
					relstr += n1 + "\t" + n2 + "\t" + "Paralogs" + ";;"
			
		
		f = open(self.relations_filenames[0], 'w')
		f.write(relstr)
		f.close()
		
		
		
		
		
class OrthoMCLClusterPredictor(HyPPO_ClusterPredictor):
	
	cached_clusters = []
	
	def get_clusters_filename(self, gene_family, no_iter, speciestree_filename = None):
		return join(gene_family.workdir, os.path.basename(gene_family.seq_file) + ".orthomcl_" + str(no_iter) + ".clusters")
	
	def infer_clusters(self, hyppo_scores, gene_family, no_iter, speciestree_filename = None):
		
		workdir_seqs = join(gene_family.workdir, "in")
		workdir_out = join(gene_family.workdir, "out")
		
		
		if len(self.cached_clusters) == 0:
			gfiles = []
			for g in self.all_gene_families:
				gfiles.append(g.seq_file)
			seqs_by_species = Sequence.combineSequenceFilesBySpecies(files = gfiles, species_separator = self.other_args["species_separator"], species_index = int(self.other_args["species_index"]))
			
			
			if not os.path.exists(workdir_seqs):
				os.mkdir(workdir_seqs)
			for key in seqs_by_species:
				outfile = join(workdir_seqs, key + ".fasta")
				Sequence.outputSequencesToFasta(seqs_by_species[key], outfile, "", False, False)
		
			cmd = "/u/lafonman/src/orthomcl-pipeline/bin/orthomcl-pipeline -i " + workdir_seqs + " -o " + workdir_out + " -m /u/lafonman/src/orthomcl-pipeline/orthomcl.conf --nocompliant --yes"
			print("EXEC " + cmd)
			os.system(cmd)
		
			seen_genes = set()
			clusters = []
			clfile = join(workdir_out, "groups/groups.txt")
			f = open(clfile, 'r')
			for line in f:
				line = line.replace("\n", "")
				if line != "":
					gz = line.split(":")[1].split()
					cluster = set()
					for g in gz:
						gname = g.split("|")[1]
						cluster.add(gname)
						seen_genes.add(gname)
					clusters.append(cluster)
			
			f.close()
			
			for s in gene_family.sequences:
				name = s.name
				if not name in seen_genes:
					cl = set()
					cl.add(name)
					clusters.append(cl)
					
			self.cached_clusters = clusters
		else:
			print("USING CACHED CLUSTERS")
			clusters = self.cached_clusters
		
		
		#restrict clusters to current family
		f_set = set()
		for s in gene_family.sequences:
			f_set.add(s.name)
		f_clusters = []
		for cl in clusters:
			inter = f_set.intersection(cl)
			if len(inter) > 0:
				f_clusters.append(inter)
		
			
		
		cl_outfile = self.get_clusters_filename(gene_family, no_iter, speciestree_filename, other_args)
		write_clusters(cl_outfile, f_clusters)
		
		
		
		
		