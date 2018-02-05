# HyPPO
Hybrid Prediction of Paralogs and Orthologs

HyPPO is a program that takes as input as set of gene sequences (DNA or AA) and predicts the orthology clusters (primary orthologs) and the inter-cluster orthologs (secondary orthologs).

To use the default settings, run

> python3 hyppo.py --infiles=myfile1.fa,myfile2.fa,myfile3.fa

where the files, separated by a comma, are aligned fasta of phylip sequence files.  
The gene names must have the format SPECIES__GENE to identify their species (the format can be changed, see below).

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

