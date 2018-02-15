This directory contains the scripts that were used to :

- generate the trees using simphy (exec_simphy.py)
- multiply the duplication branch lengths (exec_simphy.py, using the sgutils program)
- generate the sequences for each tree (exec_simphy.py)
- run the experiments (hyppo_simphy_util.py)
- get some statistics from the results (hyppo_simphy_util.py and util.py)

The file data_commands.txt contains the commands that were executed to 

The swisstree_experiments.tar.gz contains all the files used and generated for our real data experiments. 
The simphy_files directory contains the config files for simphy.  
If you need the files for the simulations, Manuel Lafond will happily provide them
(there were too many files for github).
Please email me at  mlafond2@uottawa.ca

These scripts are not intended for general usage, and henceforth are not very user friendly. 
It contains hard-coded paths and has essentially no documentation. 
The user shold in principle be able to replicate the results (although 
different trees will be obtained due to randomness) after changing the hard-coded paths in the scripts.
Do not forget to compile sgutils.


If you want to replicate the results, Manuel Lafond will be more than happy to offer his full collaboration.
Once again, please email me at  mlafond2@uottawa.ca

