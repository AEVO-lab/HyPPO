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


#these are the default modes
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
Please refer to the argument documentation and to the hyppo_classes.py file for more information on 
the available modules.

Consult the hyppo_classes.py file for more information on creating new modules.

'''




infiles = ""
nb_iter = 1
workdir = ""
outdir = ""
fullmode = ""
speciestree_filename = ""
mode_string = "default"
species_separator = "__"
species_index = "0"
timefile = ""

compute_pctid = False

arg_scores_mode = ""
arg_cluster_init_mode = ""
arg_cluster_sp_mode = ""
arg_sptree_mode = ""
arg_intercluster_mode = ""

other_args = {}
other_args_str = ""

help_arg = False

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
  if arg.startswith("--sptree_mode="):
    arg_sptree_mode = arg.split("=")[1]
  if arg.startswith("--intercluster_mode="):
    arg_intercluster_mode = arg.split("=")[1]
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
  if arg in ["--help", "-h"]:
    help_arg = True

	
help_params = []
help_params.append(["--infiles=", "Set of input sequence files (in fasta or phylip format), separated by a comma.  e.g. --infiles=seq1.fa,seq2.fa"])
help_params.append(["--workdir=", "Directory for temporary files.  HyPPO will not delete nor cleanup.  Useful for caching long tasks such as alignements."])
help_params.append(["--mode_string=", "Recommended to be set by user.  Suffix of the temporary and out files.  The generated files will have this string preceding the extension.  Helps avoid overwriting files when running different modes.  Default:empty."])
help_params.append(["--outdir=", "Output directory."])
help_params.append(["--scores_mode=", "Set the hyppo class to use to compute pairwise scores.  Default: hyppo_classes.ScoresPctID.  Recommended: hyppo_classes.MAFFTPctID if MAFFT is installed."])
help_params.append(["--cluster_init_mode=", "Set the hyppo class to compute the initial set of clusters (without the species tree).  Default: hyppo_classes.MaxScoreClusterPredictor."])
help_params.append(["--sptree_mode=", "Set the hyppo class to compute the species tree from the clusters.  Default: hyppo_classes.BottomUpSpeciesTreeMaker."])
help_params.append(["--cluster_sp_mode=", "Set the hyppo class to compute the set of clusters using the species tree.  Default: None"])
help_params.append(["--intercluster_mode=", "Set the hyppo class to compute the inter cluster relations.  Default: GreedyInterClusterPredictor"])
help_params.append(["--fullmode=", "Use a special class that implements all four steps in its own manner.  Not set by default.  If set, overrides all step-specific class specified.  e.g. --fullmode=hyppo_classes.OMAPredictor."])
help_params.append(["--speciestree=", "Filename of the species tree, in newick format, if 'true' species tree is known.  Default: not set."])
help_params.append(["--species_separator=", "The character in the gene names that is used to separate the gene name from the species name.  Default:__ (double-underscore)"])
help_params.append(["--species_index=", "The position of the species name in the gene name, after separation by the separator.  Default:0"])
help_params.append(["--use_cache", "If set, the argument 'use_cache' is passed to the hyppo classes.  These classes may or may not use it.  Most classes that compute alignements use it to avoid re-computing the alignment."])

help_params.append(["--nbiter=", "maximum number of iterations of the species tree - cluster loop to perform.  Only used if cluster_sp_mode is set.  Default: 10"])
help_params.append(["--timefile=", "Filename of a file in which the total time taken is output.  Default: not set"])
help_params.append(["--rename_genes", "When this flag is set, every gene in the work dir is renamed uniquely.  We append the index of the gene family to the name of every gene.  Useful if input sequence files have gene names in common."])
help_params.append(["--other_args=", "Other arguments that are not used directly, but are sent to each hyppo class.  Refer to the hyppo classes for specific usage.  The format is --other_args=param1=value1;;param2=value2;;param3=value3"])



helpstr = "HyPPO (Hybrid Prediction of Paralogs and Orthologs \n\n"
helpstr += "--infiles is the only required argument, and consists in a list of sequence filenames.  Currently, fasta and phylip are supported.  The extension of the sequence filenames must be in {.fa, .fasta, .fst, .phy, .phylip}\n\n"
helpstr += "The user can specify which module to use for each step of the pipeline.  Please refer to the argument list below, and to the hyppo classes documentation for more informaiton.\n\n"
helpstr += "The gene names in the sequence files are expected to contain the name of the species that contains them.  There must be a string separating the gene name from the species name.  For instance, a gene name could be HUMAN__POP1.  The default separator is '__' (double underscore), and the position of the species in this string is 0 by default.  This can be changed using the species_separator and species_index arguments.  For instance, if your gene string is POP1_HUMAN, set --species_separator='_' and --species_index=1.  What really matters is that hyppo finds the species name at the given index - the rest of the name does not matter.  For instance, POP1_HUMAN_SOME_RANDOM|ANNOTATION with the same separator/index will work.\n\nThe full list of arguments is below.  Arguments ending with a = need to be assigned a value, the others are flags that do not need a value.\n\n"

for hp in help_params:
	helpstr += hp[0].ljust(25) + hp[1] + "\n"

if help_arg:
	print(helpstr)
	sys.exit()
	
	
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
if arg_sptree_mode != "":
	mz[2] = arg_sptree_mode
if arg_cluster_sp_mode != "":
	mz[3] = arg_cluster_sp_mode
if arg_intercluster_mode != "":
	mz[4] = arg_intercluster_mode

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
	
		
		