from sequences import Sequence, Distances
import sys
import random
import hyppo_classes
import os
#import numpy
import filecmp
import math
from os import listdir
from os.path import isfile, join
from shutil import copyfile
from timeit import default_timer as timer

#eg python3 hyppo.py --infiles=./testdir/alignment_01_1.0000000000_TRUE.phy,./testdir/alignment_02_1.0000000000_TRUE.phy,./testdir/alignment_03_1.0000000000_TRUE.phy,./testdir/alignment_04_1.0000000000_TRUE.phy,./testdir/alignment_05_1.0000000000_TRUE.phy,./testdir/alignment_06_1.0000000000_TRUE.phy,./testdir/alignment_07_1.0000000000_TRUE.phy,./testdir/alignment_08_1.0000000000_TRUE.phy,./testdir/alignment_09_1.0000000000_TRUE.phy,./testdir/alignment_10_1.0000000000_TRUE.phy

#eg --modes=ScoresPctID,MaxScoreClusterPredictor,BottomUpSpeciesTreeMaker,BottomUpClusterPredictor,GreedyInterClusterPredictor

#modes = "hyppo_classes.ScoresPctID,hyppo_classes.MaxScoreClusterPredictor,hyppo_classes.BottomUpSpeciesTreeMaker,hyppo_classes.BottomUpClusterPredictor,hyppo_classes.GreedyInterClusterPredictor"
modes = "hyppo_classes.ScoresPctID,hyppo_classes.MaxScoreClusterPredictor,hyppo_classes.BottomUpSpeciesTreeMaker,None,hyppo_classes.GreedyInterClusterPredictor"



'''
THE PIPELINE

hyppo proceed in four main steps (with step 3 split into 2 substeps)
1) compute a score between every pair of genes.  
2) compute the set of orthology clusters from these scores
3.a) compute a species tree from these clusters (skipped if a species tree is provided)
3.b) re-compute orthology clusters given the constructed species tree (optional)
4) compute orthologies between genes from distinct clusters

By default, step 3.b is not performed, as the added value of the species tree is not fully understood at the moment. 
Nonetheless, if this step is performed (through a command line argument), it is possible to loop 
through steps 3.a and 3.b until convergence, as the new clusters found at step 3.b can be used to obtain a 
better species tree.  A maximum number of iterations can be specified.

Each of the five above tasks is performed by a module in hyppo.  
The user can therefore specify how each of these tasks should be performed, using the command line arguments.
To do this, the user needs to provide the name of the module that should execute the task.



HOW TO ADD A MODULE

hyppo works with five main superclasses, each with its own responsibility in the pipeline.
HyPPO_Scores: computes pairwise scores between genes
HyPPO_ClusterPredictor: computes orthology clusters
HyPPO_SpeciesTreeMaker computes a species tree from clusters
HyPPO_InterClusterPredictor: computes orthologs between clusters

Each of these classes is abstract, and the pipeline uses instances for classes that 
inherit from these.  
For instance by default, hyppo uses an instance of hyppo_classes.ScoresPctID as the actual HyPPO_Scores that is used, 
which simply assumes the the given sequences are aligned and computes the % id between genes. 

If the used wishes to provide his/her own implementation of scores computation, one must:
1) create a class that inherits HyPPO_Scores and implements all the necessary functions (get_scores_filename and load_scores.
   Suppose that this class is called my_classes.MyClass
2) specify that my_classes.MyClass should be used to compute scores.  
   For this, call the proram with the argument --scores_mode=my_classes.MyClass

'''




infiles = ""
nb_iter = 1
workdir = ""
outdir = ""
fullmode = ""
speciestree_filename = ""
mode_string = ""
species_separator = "__"
species_index = "0"
timefile = ""

compute_pctid = False

arg_scores_mode = ""
arg_cluster_init_mode = ""
arg_cluster_sp_mode = ""

other_args = {}
other_args_str = ""

rename_genes = False

for arg in sys.argv[1:]:
  if arg.startswith("--modes="):
    modes = arg.split("=")[1]
  if arg.startswith("--mode_string="):
    mode_string = arg.split("=")[1]
  if arg.startswith("--infiles="):
    infiles = arg.split("=")[1]
  if arg.startswith("--nbiter="):
    nb_iter = int(arg.split("=")[1])
  if arg.startswith("--workdir="):
    workdir = arg.split("=")[1]
  if arg.startswith("--outdir="):
    outdir = arg.split("=")[1]
  if arg.startswith("--fullmode="):
    fullmode = arg.split("=")[1]
  if arg.startswith("--speciestree="):
    speciestree_filename = arg.split("=")[1]
  if arg.startswith("--scores_mode="):
    arg_scores_mode = arg.split("=")[1]
  if arg.startswith("--cluster_init_mode="):
    arg_cluster_init_mode = arg.split("=")[1]
  if arg.startswith("--cluster_sp_mode="):
    arg_cluster_sp_mode = arg.split("=")[1]
  if arg.startswith("--rename_genes"): 
    rename_genes = True
  if arg.startswith("--compute_pctid"):
    compute_pctid = True
  if arg.startswith("--other_args="):
    other_args_str = arg.split("=", 1)[1]
  if arg.startswith("--species_separator="):
    species_separator = arg.split("=")[1]
  if arg.startswith("--species_index="):
    species_index = arg.split("=")[1]
  if arg.startswith("--timefile="):
    timefile = arg.split("=")[1]
  if arg in ["--use_cache", "--use_clustal", "--use_speciestree_cache"]:
    other_args[arg.replace("--", "")]  = "1"

	
if timefile != "":
	tfile = open(timefile, 'a')
	tfile.write("INPUT=" + (" ".join(sys.argv[:])) + "\n")
	tfile.close()
	tstart = timer()
	


other_args["species_separator"] = species_separator
other_args["species_index"] = species_index
other_args["mode_string"] = mode_string

if compute_pctid:
	other_args["compute_pctid"] = "1"

if other_args_str != "":
	argz = other_args_str.split(";;")
	for a in argz:
		az = a.split("=")
		other_args[az[0]] = az[1]
  
  

  
if infiles == "":
	print("--infiles must be specified.")
	sys.exit()

if mode_string == "":
	mode_string = fullmode
	
if workdir == "":
	workdir = os.path.dirname(infiles.split(",")[0])
	print("--workdir not set.  Will use " + workdir)
else:
	if not os.path.exists(workdir):
		os.makedirs(name=workdir, exist_ok=True)

		
if outdir == "":
	outdir = workdir
	print("--outdir not set.  Will output to " + workdir)
else:
	if not os.path.exists(outdir):
		os.makedirs(name=outdir, exist_ok=True)
		
def get_class( kls ):
	parts = kls.split('.')
	module = ".".join(parts[:-1])
	m = __import__( module )
	for comp in parts[1:]:
		m = getattr(m, comp)            
	return m
 		
mz = modes.split(",")
if arg_scores_mode != "":
	mz[0] = arg_scores_mode
if arg_cluster_init_mode != "":
	mz[1] = arg_cluster_init_mode
if arg_cluster_sp_mode != "":
	mz[3] = arg_cluster_sp_mode

#here's an annoying problem: when having multiple input files, some (distinct) genes will have the same name
#(this happens with simphy)
#so here we copy the alginment files, but rename the genes with its file index
newinfiles = ""
files_list = infiles.split(",")

newfile_to_old_file = {}

if rename_genes:
	print("Copying alignment files, ensuring gene name uniqueness...")
	for i in range(len(files_list)):
		f = files_list[i]
		filename, ext = os.path.splitext(f)
		newfile = join(workdir, os.path.basename(f).replace(ext, "_" + str(i) + ".fa"))
		sequences = Sequence.readSequences(f)
		Sequence.outputSequencesToFasta(sequences, newfile, str(i), True)
		
		if newinfiles != "":
			newinfiles += ","
		newinfiles += newfile
		
		newfile_to_old_file[newfile] = f
	
	infiles = newinfiles
else:
	for i in range(len(files_list)):
		f = files_list[i]
		newfile_to_old_file[f] = f
	



#############################################################################################
# If a fullmode is set, we just call the fullmode class, which, by definition of full, takes care of every step
#############################################################################################
if fullmode != "":
	hyppo_predictor = get_class(fullmode)()
	
	hyppo_predictor.other_args = other_args
	hyppo_predictor.predict_orthologs(infiles.split(","), workdir, speciestree_filename)
	
	outfile = join(outdir, mode_string + ".clusters")
	print("Copying " + hyppo_predictor.clusters_filenames[0] + " to " + outfile)
	copyfile(hyppo_predictor.clusters_filenames[0], outfile)
	
	outfile = join(outdir, mode_string + ".relations")
	print("Copying " + hyppo_predictor.relations_filenames[0] + " to " + outfile)
	copyfile(hyppo_predictor.relations_filenames[0], outfile)
	
	
else:
	print("The following algorithms will be used:")
	print("Scores: " + mz[0])
	print("Clusters init: " + mz[1])
	
	
	if speciestree_filename == "":
		print("Species tree: " + mz[2])
	else:
		print("Species tree: " + speciestree_filename)
	
	print("Clusters from species tree: " + mz[3])
	print("Inter clusters: " + mz[4] + "\n")
	
	#initialize the objects for every step
	hyppo_scores = get_class(mz[0])()
	hyppo_scores.other_args = other_args
	
	if mz[1] != "None":
		hyppo_clusters_init = get_class(mz[1])()
		hyppo_clusters_init.other_args = other_args
	else:
		hyppo_clusters_init = None
		
	hyppo_speciestree = get_class(mz[2])()
	hyppo_speciestree.other_args = other_args
	
	if mz[3] != "None":
		hyppo_clusters = get_class(mz[3])()
		hyppo_clusters.other_args = other_args
	else:
		hyppo_clusters = None
		
	hyppo_intercluster = get_class(mz[4])()
	hyppo_intercluster.other_args = other_args
	
	if speciestree_filename != "":
		print("Species tree is set.  Will jump directly to Clusters inference using this species tree")
	
	#############################################################################################
	# Step 0: just read the scores from the provided alignment files
	#############################################################################################
	print("Step 0: Reading scores")

	gene_families = []
	for f in infiles.split(","):
		if not isfile(f):
			print("ERROR: file " + f + " does not exist.  Program will exit")
			sys.exit()
		g = hyppo_classes.GeneFamily()
		g.workdir = workdir
		
		if "use_clustal" in other_args:
			clfile = f + ".clustal.fa"
			cmd = "clustalo -i " + f + " --dealign -o " + clfile
			print("EXEC " + cmd)
			if not os.path.isfile(clfile):
				os.system(cmd)
			else:
				print("Actually, file exists and we'll use it.")
			
			g.load_from_alignment_file(clfile)
			
			
			newfile_to_old_file[clfile] = newfile_to_old_file[f]
		else:
			g.load_from_alignment_file(f)
		hyppo_scores.load_scores(g)
		gene_families.append(g)

	if hyppo_clusters_init != None:
		hyppo_clusters_init.all_gene_families = gene_families
	
	if hyppo_clusters != None:
		hyppo_clusters.all_gene_families = gene_families

	#############################################################################################
	# Step 1: sptree is set: if we already have the species tree, we use the cluster-with-species inference if it is set, if not, we use the init-cluster-inferer
	#############################################################################################
	if speciestree_filename != "":
		
		if hyppo_clusters != None:
			inferer = hyppo_clusters
		else:
			inferer = hyppo_clusters_init
		
		print("Step 1: Inferring clusters using given species tree.")
		for g in gene_families:
			print("Inferring clusters for " + g.seq_file)
			basefile = newfile_to_old_file[g.seq_file]
			inferer.infer_clusters(hyppo_scores, g, 1, speciestree_filename)
			clusters_file = inferer.get_clusters_filename(g, 1, speciestree_filename)
			outfile_clusters = join(outdir, os.path.basename(basefile) + "." + mode_string + ".clusters")
			copyfile(clusters_file, outfile_clusters)
			g.clusters_filename = outfile_clusters

		sptreefile = speciestree_filename
	else:
		sptreefile = ""
		
		#############################################################################################
		# Step 1: Infer initial set of clusters
		#############################################################################################
		
		print("Step 1: inferring initial set of clusters")
		
		for g in gene_families:
			print("Inferring clusters for " + g.seq_file)
			hyppo_clusters_init.infer_clusters(hyppo_scores, g, 0, None)
			
			basefile = newfile_to_old_file[g.seq_file]
			init_clusters_file = hyppo_clusters_init.get_clusters_filename(g, 0, "")
			g.init_outfile_clusters_filename = join(outdir, os.path.basename(basefile) + "." + mode_string + ".clusters")
			g.clusters_filename = init_clusters_file
		
		#############################################################################################
		# Step 2, 3: Species tree - clusters loop
		#############################################################################################
		
		for i in range(0, nb_iter):
			print("Step 2: Inferring species tree (iter=" + str(i + 1) + "/" + str(nb_iter) + ")")
			sptreefile = hyppo_speciestree.infer_species_tree(hyppo_scores, gene_families, i + 1)
			
			#if the inferred species tree hasn't changed, no point in going on - we have "converged" and can stop here
			sametree = False
			if i > 0 and filecmp.cmp(prevsptreefile, sptreefile):
				sametree = True
				
			print("Step 3: Inferring clusters using species tree (iter=" + str(i + 1) + "/" + str(nb_iter) + ")")
			
			for g in gene_families:
				if hyppo_clusters != None:
					print("Reinferring clusters for " + g.seq_file + "(iter=" + str(i + 1) + "/" + str(nb_iter) + ")")
					hyppo_clusters.infer_clusters(hyppo_scores, g, i + 1, sptreefile)
					clusters_file = hyppo_clusters.get_clusters_filename(g, 1, speciestree_filename)
					g.clusters_filename = clusters_file
					
					#if we inferred the final set of clusters, 
					if i == nb_iter - 1 or sametree:
						basefile = newfile_to_old_file[g.seq_file]
						outfile_clusters = join(outdir, os.path.basename(basefile) + "." + mode_string + ".clusters")
						
						print("Copying " + clusters_file + " to " + outfile_clusters)
						copyfile(clusters_file, outfile_clusters)
				else:
					if i == nb_iter - 1 or sametree:
						print("Using initial clusters in step 3.  outfile=" + g.init_outfile_clusters_filename)
						copyfile(g.clusters_filename, g.init_outfile_clusters_filename)

			prevsptreefile = sptreefile
			if sametree:
				print("Species tree hasn't changed - we can stop here.")
				break
				
				
			
	#at this point, g.clusters should be set
	
	#############################################################################################
	# Infer inter cluster relations
	#############################################################################################
	print("Step 4: Inferring inter-cluster relations")
	for g in gene_families:
		print("Inferring relations for " + g.seq_file)
		hyppo_intercluster.infer_intercluster_relations(hyppo_scores, g, sptreefile)
		#copy final relations to outdir
		relsfile = hyppo_intercluster.get_relations_filename(g, sptreefile)
		basefile = newfile_to_old_file[g.seq_file]
		outfile = join(outdir, os.path.basename(basefile) + "." + mode_string + ".relations")
		print("Copying " + relsfile + " to " + outfile)
		copyfile(relsfile, outfile)
	
		
if timefile != "":
	ttime = timer() - tstart
	tfile = open(timefile, 'a')
	tfile.write("TIME=" + str(ttime) + "\n")
	tfile.close()
	
		
		