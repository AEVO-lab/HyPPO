import numpy
import sys
import os
from os.path import isfile, join

fname = "" #"/u/lafonman/Projects/simgraphs/simphy/DL0000005/01/mult/alignment_10_50.0000000000_TRUE.bbh"

outfilename = ""
species_separator = "__"
species_index = 0


for arg in sys.argv[1:]:
	if arg.startswith("--infile="):
		pz = arg.split("=")
		fname = pz[1]
	
	if arg.startswith("--outfile="):
		pz = arg.split("=")
		outfilename = pz[1]
	if arg.startswith("--nooverwrite"):
		overwriteExisting = False
		
	if arg.startswith("--species_separator="):
		species_separator = arg.split("=")[1]
	if arg.startswith("--species_index="):
		species_index = int(arg.split("=")[1])
	
		
if fname == "":
	print("--infile must be specified")
	sys.exit()
		

		
class GreedyClusterer:
	
	
	original_matrix = []
	original_labels = []
	matrix = []

	labels = []
	
	species_indices = {}
	label_indices = {}
	
	def __init__(self, matrix, labels):
		self.matrix = matrix
		self.labels = labels
		self.original_labels = labels[:]
		
		self.original_matrix = [row[:] for row in matrix]
		
		for i in range(len(self.labels)):
			self.label_indices[self.labels[i]] = i
			
			sp = self.get_species(self.labels[i])
			if sp not in self.species_indices:
				self.species_indices[sp] = []
			self.species_indices[sp].append(i)
			
		
		
		
	def get_species(self, label):
		#TODO: remove evil globals
		global species_separator
		global species_index
	
		pz = label.split(species_separator)
		return pz[species_index]
		
	def get_fav_in_species(self, i, sp):
		indices = self.species_indices[sp]
		
		max = -1
		jmax = 0
		for j in indices:
			if self.original_matrix[i][j] > max:
				max = self.original_matrix[i][j]
				jmax = j

		return jmax
				
	def find_clusters(self):
		
		done = False

		###################################################
		#1st pass: merge all direct favorite clusters with no common species
		while not done:
			
			imax = 0
			jmax = 0
			max = -1
			for i in range(len(self.matrix)):
				
				favj = self.find_favorite_cluster(i, True)
				
				if favj >= 0 and not self.have_common_species(i, favj) and self.matrix[i][favj] > max:
					imax = i
					jmax = favj
					max = self.matrix[i][favj]

			if max < 0:
				done = True
			else:	
				#print("Merging " + labels[rmax] + " " + labels[cmax])
				self.merge_elements(imax, jmax)
			
		###################################################
		#2nd pass: find more clusters to join
		done = False
		
		while not done:
			
			imax = 0
			jmax = 0
			max = -1
			for i in range(len(self.matrix)):
				
				favj = self.find_favorite_cluster(i, False)
				
				if favj >= 0 and self.matrix[i][favj] > max and self.have_enough_favorites(i, favj):
					imax = i
					jmax = favj
					max = self.matrix[i][favj]

			if max < 0:
				done = True
			else:	
				print("Phase 2: Merging " + str(imax) + " " + str(jmax))
				self.merge_elements(imax, jmax)
			

		clusters = []
		
		for lbl in self.labels:
			clusters.append(lbl.split(";;"))
		return clusters
	
	
	def have_enough_favorites(self, i, j):
		nbfav = 0
		nbrels = 0
		l1 = self.labels[i].split(";;")
		l2 = self.labels[j].split(";;")
		for k1 in range(len(l1)):
			for k2 in range(len(l2)):
				lbl1 = l1[k1]
				lbl2 = l2[k2]
				i = self.label_indices[lbl1]
				j = self.label_indices[lbl2]
				sp1 = self.get_species(lbl1)
				sp2 = self.get_species(lbl2)
				
				fav_of_i = self.get_fav_in_species(i, sp2)
				if (fav_of_i == j):
					nbfav += 1
					
				fav_of_j = self.get_fav_in_species(j, sp1)
				if (fav_of_j == i):
					nbfav += 1
				
				nbrels+=2
				
		pct = float(nbfav)/float(nbrels)
		
		return (pct > 0.7)
	
	
	def find_favorite_cluster(self, i, allow_same_species):
		
		maxj = -1
		max = -1
		for j in range(len(self.matrix)):
			if i != j and self.matrix[i][j] > max:
				maxj = j
				max = self.matrix[i][j]
		return maxj
	
			
	def have_common_species(self, i, j):
		l1 = self.labels[i].split(";;")
		s1 = set()
		
		for l in l1:
			sx = self.get_species(l)
			s1.add(sx)
			
		l2 = self.labels[j].split(";;")
		for l in l2:
			sx = self.get_species(l)
			if sx in s1:
				return True
				
		return False
					
	def merge_elements(self, index1, index2):
		i = min(index1, index2)
		j = max(index1, index2)
		
		newlabel = self.labels[i] + ";;" + self.labels[j]
		
		newrow = [0]*len(self.matrix)
		
		for k in range(len(self.matrix)):
			#if self.have_common_species(i, k) or self.have_common_species(j, k):
			#	newrow[k] = -1
			#else:
			newrow[k] = max(self.matrix[i][k], self.matrix[j][k])
			
			self.matrix[k].append(newrow[k])
		
		newrow.append(0)
		
		self.matrix.append(newrow)
		self.labels.append(newlabel)
		
		#remove i-th and j-th column.  Note that i < j
		for r in range(len(self.matrix)):
			del self.matrix[r][j]
			del self.matrix[r][i]
		
		del self.matrix[j]
		del self.matrix[i]
		
		del self.labels[j]
		del self.labels[i]
		
		
		
		

matrix = []

def read_matrix(fname):
	l_matrix = []
	l_labels = []

	with open(fname) as f:

		nblines = 0
		for line in f:
			
			line = line.replace("\n", "")
			if nblines == 0:
				l_labels = line.split(';')[1:-1]
				ncols = len(l_labels)
			else:
				
				pz = line.split(';')[1:-1]	#lines end with an extra ;
				arr = []
				
				for p in range(len(pz)):
					arr.append( float(pz[p]) )
				
				l_matrix.append( arr )
			nblines += 1
			
		#if debug:
		#	print(str(len(l_labels)) + " labels")
		#	sys.exit()

	return [l_matrix, l_labels]
		
p = read_matrix(fname)
matrix = p[0]
labels = p[1]

		
cluster_guy = GreedyClusterer(matrix, labels)

clusters = cluster_guy.find_clusters()

trout = ""	#trout because it's more fun than strout
for c in clusters:
	trout += "<GROUP>\n" + "\n".join(c) + "\n</GROUP>\n"

if outfilename == "":
	print(trout)
else:
	fout = open(outfilename, 'w')
	fout.write(trout)
	fout.close()
		