# HyPPO
Hybrid Prediction of Paralogs and Orthologs

HyPPO is a program that takes as input as set of gene sequences (DNA or AA) and predicts the orthology clusters (primary orthologs) and the inter-cluster orthologs (secondary orthologs).

To use the default settings, run

\> python3 hyppo.py --infiles=myfile1.fa,myfile2.fa,myfile3.fa

where the files, separated by a comma, are aligned fasta of phylip sequence files.  
The gene names must have the format SPECIES__GENE to identify their species (the format can be changed, see below).

It is highly recommended to specify a --workdir=[workdir] directory, as HyPPO can create many temporary files.
An --outdir=[outdir] directory can also be specified. 
The modules to be used by HyPPO can also be specified, see below for more information.


# Prerequisites, compiling and running

The program is divided in modules, each of which can be implemented independently.  
The core of HyPPO requires python 3 to run (it should in principle work with python 2, but this was not tested).
The default modules provided are implemented in C++, and must be compiled.

These modules are all part of the OCR program.  OCR stands for Orthology Cluster Recovery - though the program does more than that now.
OCR was developed using QtCreator, which generates the makefile provided in the directory. 
Note that there are no dependencies with the Qt framework.
Normally, entering the ./OCR directory and running make should do the job.

\> cd ./OCR
\> make

should create the OCR binary. 
** This binary should then be accessible through the PATH environment variable.

If make does not work and you have Qt installed, you can try to regenreate the makefile, as follows:

\> cd ./OCR
\> qmake OCR.pro config+=RELEASE
\> make

If this still fails, contact me at mlafond2@uottawa.ca and I will create the standard 
configure - make - make install commands accessible.

# The HyPPO pipeline

HyPPO proceed in four main steps (with step 3 split into 2 substeps)

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


# Output format

If k files are given as input, HyPPO outputs 2k files, two per input file.
Say the --infiles contains myfile.fa.  Then the output dir (cwd by default) will contain 

myfile.fa.default.clusters

myfile.fa.default.relations

The "default" suffix can be changed using the --mode_string argument. 
The .clusters file contains a partition of the genes present in the myfile.fa file. 
These are the orthology clusters.
The file format is as follows:

\<GROUP\> \
gene1 \
gene2 \
... \
\<\/GROUP\> \
\<GROUP\> \
gene3 \
gene4 \
... \
\<\/GROUP\> \
... 

where each \<GROUP\> corresponds to an orthology cluster.

The .relations file has two lines, one for the gene tree that displays the relations, and one that has the relations.
The format is 

TREE=(the tree in newick format) \
RELATIONS=g1[tab]g2[tab]reltype;;g1[tab]g3[tab]reltype;;...

That is, the relations are all on one big line, separated by ;;

Each relation is represented as a tab-separated (the \\t character) string with 3 elements g1,g2,reltype.
g1 and g2 are the two genes, reltype is one of "Orthologs" or "Paralogs".

# Command line arguments

This is the help message output when running python3 hyppo.py --help.

<pre>

--infiles is the only required argument, and consists in a list of sequence filenames.  Currently, fasta and phylip are supported.  The 
extension of the sequence filenames must be in {.fa, .fasta, .fst, .phy, .phylip}

The user can specify which module to use for each step of the pipeline.  Please refer to the argument list below, and to the hyppo 
classes documentation for more informaiton.

The gene names in the sequence files are expected to contain the name of the species that contains them.  There must be a string 
separating the gene name from the species name.  For instance, a gene name could be HUMAN__POP1.  The default separator is '__' (double 
underscore), and the position of the species in this string is 0 by default.  This can be changed using the species_separator and 
species_index arguments.  For instance, if your gene string is POP1_HUMAN, set --species_separator='_' and --species_index=1.  What 
really matters is that hyppo finds the species name at the given index - the rest of the name does not matter.  For instance, 
POP1_HUMAN_SOME_RANDOM|ANNOTATION with the same separator/index will work.

--infiles=               Set of input sequence files (in fasta or phylip format), separated by a comma.  e.g. --infiles=seq1.fa,seq2.fa
--workdir=               Directory for temporary files.  HyPPO will not delete nor cleanup.  Useful for caching long tasks such as alignements.
--mode_string=           Recommended to be set by user.  Suffix of the temporary and out files.  The generated files will have this string preceding the extension.  Helps avoid overwriting files when running different modes.  Default:empty.
--outdir=                Output directory.
--scores_mode=           Set the hyppo class to use to compute pairwise scores.  Default: hyppo_classes.ScoresPctID.  Recommended: hyppo_classes.MAFFTPctID if MAFFT is installed.
--cluster_init_mode=     Set the hyppo class to compute the initial set of clusters (without the species tree).  Default: hyppo_classes.MaxScoreClusterPredictor.
--sptree_mode=           Set the hyppo class to compute the species tree from the clusters.  Default: hyppo_classes.BottomUpSpeciesTreeMaker.
--cluster_sp_mode=       Set the hyppo class to compute the set of clusters using the species tree.  Default: None
--intercluster_mode=     Set the hyppo class to compute the inter cluster relations.  Default: GreedyInterClusterPredictor
--fullmode=              Use a special class that implements all four steps in its own manner.  Not set by default.  If set, overrides all step-specific class specified.  e.g. --fullmode=hyppo_classes.OMAPredictor.
--speciestree=           Filename of the species tree, in newick format, if 'true' species tree is known.  Default: not set.
--species_separator=     The character in the gene names that is used to separate the gene name from the species name.  Default:__ (double-underscore)
--species_index=         The position of the species name in the gene name, after separation by the separator.  Default:0
--use_cache              If set, the argument 'use_cache' is passed to the hyppo classes.  These classes may or may not use it.  Most classes that compute alignements use it to avoid re-computing the alignment.
--nbiter=                maximum number of iterations of the species tree - cluster loop to perform.  Only used if cluster_sp_mode is set.  Default: 10
--timefile=              Filename of a file in which the total time taken is output.  Default: not set
--rename_genes           When this flag is set, every gene in the work dir is renamed uniquely.  We append the index of the gene family to the name of every gene.  Useful if input sequence files have gene names in common.
--other_args=            Other arguments that are not used directly, but are sent to each hyppo class.  Refer to the hyppo classes for specific usage.  The format is --other_args=param1=value1;;param2=value2;;param3=value3

</pre>
