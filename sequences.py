import sys
import os


class Sequence:
	'''
	class Sequence
	
	Stores a DNA/AA sequence. 
	has many static util functionalities.
	'''
	
	name = ""		#name of the sequence
	seq = ""		#the sequence, without gaps
	alignedSeq = "" #the sequence, exactly as given or as read from an alignement file (may include gaps) 

	def __init__(self, _name, _seq):
		self.name = _name 
		self.alignedSeq = _seq
		self.seq = _seq.replace("-", "").replace("*", "")
	
	@staticmethod
	def saveUnalignedSequencesToFastaFile(seqoutfile, sequences):
		'''
		Unaligns sequences and outputs them in fasta format
		sequences is a list of Sequence objects
		'''
		with open( seqoutfile, 'w') as file:
			for i in range(len(sequences)):
				file.write(">" + sequences[i].name + "\n" + sequences[i].seq + "\n")
	
	
	@staticmethod
	def readSequences(seqfile):
		'''
		Reads given alignment file and returns list of Sequence objects.
		seqfile must end with {".fa", ".fasta", ".fst", ".phy", ."phylip"}
		'''
		filename, ext = os.path.splitext(seqfile)
		
		if ext == ".fa" or ext == ".fasta" or ext == ".fst":
			return Sequence.readSequencesFromFastaFile(seqfile)
		elif ext == ".phy" or ext == ".phylip":
			return Sequence.readSequencesFromPhylipFile(seqfile)
	
		return []
	
	@staticmethod
	def readSequencesFromFastaFile(seqfile):
		'''
		Reads fasta file, returns Sequence
		'''
		sequences = []
		curseq = ""
		curseqname = ""
		
		with open(seqfile) as f:
			for line in f:
				line = line.replace("\n", "").replace("\r", "")
				if len(line) > 0:
					if line[0] == ">":
						if curseqname != "":
							#curseq = curseq.replace("-", "")
							s = Sequence( curseqname, curseq )
							sequences.append(s)
						curseqname = line.lstrip("> ")
						curseq = ""
					else:
						curseq = curseq + line
		if curseqname != "":
			#curseq = curseq.replace("-", "")
			s = Sequence( curseqname, curseq )
			sequences.append(s)

		return sequences
		
		
	@staticmethod
	def readSequencesFromPhylipFile(seqfile):
		'''
		Reads phylip file, returns Sequence
		'''
		sequences = []
		curseq = ""
		curseqname = ""
		
		nblines = 0
		with open(seqfile) as f:
			for line in f:
				line = line.replace("\n", "").replace("\r", "")
				
				if len(line) > 0:
					if nblines == 0:
						nblines += 1 #we ignore the first line
					else:	
						pz = [x for x in line.split(' ') if x] #magically splits and removes empty substrings
						
						if len(pz) >= 2:
							s = Sequence( pz[0], pz[1] )
							sequences.append(s)
						nblines += 1
		

		return sequences
		
	@staticmethod
	def combineSequenceFilesBySpecies(files, species_separator = "__", species_index = 0):
		'''
		Reads all sequences in all given files, returns a dictionary where 
		the key is the species name, the value is a list of Sequence objects for 
		the genes in this species.  
		Gene names must contain species name, as species by the separator and index.
		'''
		seqs_by_species = {}
		
		for f in files:
		
			seqs = Sequence.readSequences(f)
			
			for s in seqs:
				pz = s.name.split(species_separator)
				sp = pz[species_index]
				
				if sp not in seqs_by_species:
					seqs_by_species[sp] = [s]
				else:
					seqs_by_species[sp].append(s)
					
				
		
		return seqs_by_species
		
	@staticmethod
	def outputSequencesToFasta(sequences, filename, name_suffix = "", aligned = False, convertToAA = False, name_prefix = ""):
		'''
		outputs list of Sequence objects to fasta format.
		sequences: list of Sequence objects
		filename: where to output
		name_suffix: if needed, append a suffix to each sequence name
		aligned: do you want your sequences aligned or not (if they were not aligned to begin with, doesn't matter)?
		convertToAA: if True, will convert each codon to its corresponding AA before outputting.
		name_prefix: add a prefix to every name
		'''
		fout = open(filename, 'w')
		for s in sequences:
			seq = s.seq
			
			if convertToAA:
				seq = Sequence.getAminoAcidSequence(seq)
			elif aligned:
				seq = s.alignedSeq
			fout.write(">" + name_prefix + s.name + name_suffix + "\n" + seq + "\n")
		fout.close()
	
	
	
	@staticmethod
	def getAminoAcidSequence(dnaseq):
		'''
		Converts a DNA sequence to an AA sequence by translating each codon.
		'''
		codes = "tttf ttcf ttal ttgl cttl ctcl ctal ctgl atti atci atai atgm gttv gtcv gtav gtgv tcts tccs tcas tcgs cctp cccp ccap ccgp actt acct acat acgt gcta gcca gcaa gcga taty tacy taax tagx cath cach caaq cagq aatn aacn aaak aagk gatd gacd gaae gage tgtc tgcc tgax tggw cgtr cgcr cgar cggr agts agcs agar aggr ggtg ggcg ggag gggg ".upper()
		
		code_dict = {}
		
		c = 0
		while c + 5 <= len(codes):
			code_dict[codes[c:c+3]] = codes[c+3:c+4]
			c += 5
		
		
		aaseq = ""
		c = 0
		while c + 3 < len(dnaseq):
			trip = dnaseq[c:c+3]
			aaseq += code_dict[trip]
			
			c += 3
			
		return aaseq
		
	
	
class Distances:
	'''
	Mostly static util functions for sequence distances
	'''
	sequences = []
	

		
	@staticmethod	
	def getIdentityPct(s1, s2):
		'''
		Returns nbequal base pairs divided by max lengths (not including gaps)
		'''
		nbch = 0
		nbeq = 0
		
		#TODO: what if unequal lengths
		for i in range( max( len(s1), len(s2) ) ):
			if s1[i] != "-" or s2[i] != "-":
				nbch += 1
				if s1[i] == s2[i]:
					nbeq += 1
					
		#len1 = len(s1.lstrip("-").rstrip("-"))
		#len2 = len(s2.lstrip("-").rstrip("-"))
		#minlen = min(len1, len2)
		
		return float(nbeq)/float(nbch)

		
	@staticmethod
	def computeIdentityPct(s1, s2, convertToAA = False):
		'''
		Uses needleman-wunsch to compute maxscore between 2 sequences
		'''
		if convertToAA:
			s1 = Sequence.getAminoAcidSequence(s1)
			s2 = Sequence.getAminoAcidSequence(s2)
		
		scores = [[0 for col in range(len(s2))] for row in range(len(s1))]
		
		gapsInserted = 0
		
		for i in range(len(s1)):
			for j in range(len(s2)):
				if i == 0 or j == 0:
					scores[i][j] = 0
				else:
					ssame = scores[i - 1][j - 1]
					if s1[i] == s2[j]:
						ssame = scores[i - 1][j - 1] + 1
						
					
					
					scores[i][j] = max( scores[i - 1][j], scores[i][j - 1], ssame)
		return float(scores[len(s1) - 1][len(s2) - 1])/float(min(len(s1), len(s2)))
	
	
	@staticmethod
	def getPairwisePctID(sequences, verbose, run_nw_algorithm = False, print_progress = False):
		'''
		Computes pairwise identity pct, returns a 2D array, indices corresponding to index in sequences
		'''

		scores = [[0 for col in range(len(sequences))] for row in range(len(sequences))]
		

		for i in range(len(sequences)):
			for j in range(len(sequences)):
				if i > j:
					scores[i][j] = scores[j][i]	#avoid double-computing scores
				else:
					if run_nw_algorithm:
						scores[i][j] = Distances.computeIdentityPct( sequences[i].seq, sequences[j].seq, True )
						if print_progress:
							print("Computed " + str(i) + "," + str(j) + "  (nbseqs=" + str(len(sequences)) + ")")
					else:
						scores[i][j] = Distances.getIdentityPct( sequences[i].alignedSeq, sequences[j].alignedSeq )
						
			if verbose:
				print("Computed alignment " + str(i) + " of " + str(len(sequences)))	
		return scores
	
	@staticmethod
	def getScoresFromEdgeList(sequences, filename):
		'''
		Loads a file in the edge list format, returns 2D array
		'''
		scores = [[0 for col in range(len(sequences))] for row in range(len(sequences))]
		
		seq_index = {}
		for i in range(len(sequences)):
			seq_index[sequences[i].name] = i
		
		f = open(filename)
		for line in f:
			line = line.replace("\n", "")
			if line != "":
				pz = line.split(";")
				i1 = seq_index[pz[0]]
				i2 = seq_index[pz[1]]
				scores[i1][i2] = float(pz[2])
				scores[i2][i1] = float(pz[2])
		f.close()
		
		return scores
	
	@staticmethod
	def getMatrixString(sequences, distances):
		'''
		Converts 2D array of scores in a matrix string, ready to be output
		'''
		matrixstr = ";"

		for i in range(len(sequences)):
			matrixstr += sequences[i].name.replace("\n", "") + ";"
		matrixstr += "\n"

		for i in range(len(sequences)):
			matrixstr += sequences[i].name + ";"
			for j in range(len(sequences)):
				
				#matrixstr += str(pct) + ";"
				matrixstr += str(distances[i][j]) + ";"
			matrixstr += "\n"
			
		return matrixstr		
		
		
	

	@staticmethod
	def readMatrixFile(filename, header_mode = "csv", sep = ";"):
		'''
		Reads a file containing scores in a matrix format, returns a dict with two keys.
		dict["labels"] has the list of all names found in the file
		dict["matrix"] has a 2D array, indices correspond to those from the labels
		'''
		matrix = []
		labels = []
		with open(filename) as f:

			nblines = 0
			for line in f:
				line = line.replace("\n", "")
				if nblines == 0:
					#if header_mode == "csv":
					#	labels = line.split(sep)
					#	if labels[-1] == "":	#lines end with an extra ;
					#		labels = labels[1:-1]
					#elif header_mode == "count":
					#	pass
					pass
				else:
					
					pzx = line.split(sep)	
					pz = []
					for p in pzx:
						if p != "":
							pz.append(p)
					
					arr = []
					
					labels.append(pz[0])
					
					
					for p in range(1, len(pz)):
						arr.append( float(pz[p]) )
					
					matrix.append( arr )
				nblines += 1
		
		return {"labels" : labels, "matrix" : matrix}
		
