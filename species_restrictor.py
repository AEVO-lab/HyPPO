from sequences import Sequence, Distances
import sys
import random
import os
from os import listdir
from os.path import isfile, join

overwriteExisting = False
indir = ""
outdir = ""
speciesSetFile = "species_set.txt"

for arg in sys.argv[1:]:
	if arg.startswith("--indir="):
		pz = arg.split("=")
		indir = pz[1]
	if arg.startswith("--outdir="):
		pz = arg.split("=")
		outdir = pz[1]
	if arg.startswith("--overwrite"):
		overwriteExisting = True 
	if arg.startswith("--species"):
		pz = arg.split("=")
		speciesSetFile = pz[1]
	

species = set()
	
#read the species we restrict to
with open(speciesSetFile) as f:
	for line in f:
		line = line.replace("\n", "").replace("\r", "")
		if len(line) > 0:
			species.add(line)

print("Will restrict to species")
print(species)

#then restrict every distance file to the species
list = listdir(indir)

for f in list:
	distFile = join(indir, f)
	if isfile(distFile) and distFile.endswith(".needle"):
		outfilename = join(outdir, f).replace(".needle", ".rneedle")
		
		if (not isfile(outfilename) or overwriteExisting): 
			print("Restricting from " + distFile + " to " + outfilename)
			outfile = open(outfilename, 'w')
			with open(distFile) as fdist:
				for line in fdist:
					line = line.replace("\n", "").replace("\r", "")
					
					#format of a line is GENE1__SPECIES1 GENE2__SPECIES2 DIST (DIST2)
					if len(line) > 0:
						pz = line.split(" ")
						if len(pz) >= 2:
							p1 = pz[0].split("__")
							p2 = pz[1].split("__")
							if p1[1] in species and p2[1] in species:
								outfile.write(line + "\n")
					
			outfile.close()
			
		else:
			print("Skipped " + outfile + " (already exists)")
			
		