
import sys
import random
import os
from os import listdir
from os.path import isfile, join


conffile = ""
simphydir = ""
simphy_outdir = ""
nbreplicates = 0
mode = "all"
amino_acids = False

#python3 exec_simphy.py --simphydir=./simphy --conf=DL000001.conf --replicates=10 --simphyoutdir=DL000001
#python3 exec_simphy.py --simphydir=./simphy --conf=DLPROT05.conf --replicates=10 --simphyoutdir=DLPROT05 --use_aa

for arg in sys.argv[1:]:
  if arg.startswith("--conf="):
    pz = arg.split("=")
    conffile = pz[1]
  if arg.startswith("--simphydir="):
    pz = arg.split("=")
    simphydir = pz[1]
  if arg.startswith("--simphyoutdir="):
    pz = arg.split("=")
    simphy_outdir = pz[1]
  if arg.startswith("--replicates="):
    pz = arg.split("=")
    nbreplicates = int(pz[1])
  if arg.startswith("--mode="):
    pz = arg.split("=")
    mode = pz[1]
  if arg.startswith("--use_aa"):
    amino_acids = True
 
#execute simphy
if (mode == "all" or mode == "simphy"):
	cmd = "$HOME/bin/simphy/simphy_lnx64 -I " + join(simphydir, conffile) + " -o " + join(simphydir, simphy_outdir)
	print("EXEC " + cmd)
	os.system(cmd)

	#use indelible wrapper to generate control file, for usage with indelible 

	ctrl_file = join(simphydir, "indelible_control.txt")
	if amino_acids:
		ctrl_file = join(simphydir, "indelible_control_aa.txt")
	cmd = join(simphydir, "INDELIble_control_gen.pl") + " " + join(simphydir, simphy_outdir) + " " + ctrl_file + " 22 1"
	print("EXEC " + cmd)
	os.system(cmd)

cwd = os.getcwd()

for i in range(nbreplicates):
	str_i = str(i + 1)
	if len(str_i) == 1:
		str_i = "0" + str_i
		
	workdir = join(join(simphydir, simphy_outdir), str_i)
		
	#makes the mult directory along with multiplied edge lengths, and the summary.txt file, and the control file for indelible
	if (mode == "all" or mode == "sgutils"):
		cmd = "./sgutils/sgutils -d " + workdir
		print("EXEC " + cmd)
		os.system(cmd)
	
	multdir = join(workdir, "mult")
	
	#goto dir and execute indelible
	print("chdir: " + multdir)
	os.chdir(multdir)
	
	if (mode == "all" or mode == "indelible"):
		cmd = "indelible"
		print("EXEC " + cmd)
		os.system(cmd)
	
	#then calculate distances
	print("chdir: " + cwd)
	os.chdir(cwd)
	
	if (mode == "all" or mode == "distances"):
		cmd = "python3 distances_maker.py --indir=" + multdir + " --outdir=" + multdir
		print("EXEC " + cmd)
		os.system(cmd)
		
		cmd = "python3 distances_maker.py --mode=normalize --indir=" + multdir + " --outdir=" + multdir
		print("EXEC " + cmd)
		os.system(cmd)
		
		cmd = "python3 distances_maker.py --mode=bbh --indir=" + multdir + " --outdir=" + multdir
		print("EXEC " + cmd)
		os.system(cmd)
	

	if (mode == "all" or mode == "nni"):
		spfile = join(multdir, "speciestree_nobranch.nw")
		
		for n in [1,3,5,10]:
			spout = join(multdir, "speciestree_nobranch.nni" + str(n) + ".nw")
			cmd = "./sgutils/sgutils -m random_nni -n " + str(n) + " -f " + spfile + " -o " + spout
			print("EXEC " + cmd)
			os.system(cmd)	

