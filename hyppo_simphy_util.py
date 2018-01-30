from sequences import Sequence, Distances
import sys
import random
import hyppo_classes
import os
import mwmatching
import numpy
import math
from os import listdir
from os.path import isfile, join

'''
Run default
python3 hyppo_simphy_util.py --simphydir=./simphy/simphy/DL0000005/ --workdir=./work0000005_s_m1/ --mult=1 --replicates=10 --directories=1:10
python3 hyppo_simphy_util.py --simphydir=./simphy/simphy/DL0000005/ --workdir=./work0000005_s_m1/ --mult=1 --replicates=10 --directories=1:10 --mode=results

python3 hyppo_simphy_util.py --simphydir=./simphy/simphy/DL0000005/ --workdir=./work0000005_s_m1/ --mult=1 --replicates=10 --directories=1:10 --mode_string=default_init --cluster_sp_mode=None 
'''

'''
Run OMA 
python3 hyppo_simphy_util.py --simphydir=./simphy/simphy/DL0000005/ --mult=1 --replicates=10 --workdir=./work0000005_s_m1/ --directories=1:10 --fullmode=hyppo_classes.OMAPredictor --use_species_tree --split_into_families
python3 hyppo_simphy_util.py --simphydir=./simphy/simphy/DLPROT05/ --mult=1 --replicates=10 --workdir=./workPROT05_s_m1/ --directories=1:10 --fullmode=hyppo_classes.OMAPredictor --use_species_tree --split_into_families
python3 hyppo_simphy_util.py --simphydir=./simphy/simphy/DLPROT05/ --mult=1 --replicates=10 --workdir=./work/workPROT05_s_m1/ --directories=1:10 --fullmode=hyppo_classes.OMAPredictor --use_species_tree --split_into_families --other_args=seqtype=AA

python3 hyppo_simphy_util.py --simphydir=./simphy/simphy/DL0000005/ --mult=2 --replicates=10 --workdir=./work/work0000005_s_m2/ --directories=1:10 --fullmode=hyppo_classes.OMAPredictor --use_species_tree --split_into_families
python3 hyppo_simphy_util.py --simphydir=./simphy/simphy/DL000001/ --mult=2 --replicates=10 --workdir=./work/work000001_s_m2/ --directories=1:10 --fullmode=hyppo_classes.OMAPredictor --use_species_tree --split_into_families
python3 hyppo_simphy_util.py --simphydir=./simphy/simphy/DLFAST/ --mult=2 --replicates=10 --workdir=./work/workFAST_s_m2/ --directories=1:10 --fullmode=hyppo_classes.OMAPredictor --use_species_tree --split_into_families &
python3 hyppo_simphy_util.py --simphydir=./simphy/simphy/DL0000005/ --mult=2 --replicates=10 --workdir=./work/work0000005_s_m2/ --directories=1:10 --fullmode=hyppo_classes.OMAPredictor --use_species_tree --split_into_families &
'''


# OMA results
'''
python3 hyppo_simphy_util.py --simphydir=./simphy/simphy/DL0000005/ --mult=1 --replicates=10 --workdir=./work/work0000005_s_m1/ --directories=1:10 --fullmode=hyppo_classes.OMAPredictor --use_species_tree --split_into_families --mode=results

'''
# python3 hyppo_simphy_util.py --simphydir=./simphy/simphy/DLFAST/ --mult=2 --replicates=10 --workdir=./workFAST/ --directories=1:10 --fullmode=hyppo_classes.OMAPredictor --use_speciestree --mode=results
# default results
# python3 hyppo_simphy_util.py --simphydir=./simphy/simphy/DLFAST/ --mult=2 --replicates=10 --workdir=./workFAST/ --directories=1:10 --mode=results


# Run orthoMCL
# python3 hyppo_simphy_util.py --simphydir=./simphy/simphy/DL000001/ --mult=1 --replicates=10 --workdir=./work/workmcl/ --directories=1:10 --cluster_init_mode=hyppo_classes.OrthoMCLClusterPredictor --use_species_tree


mult = "2"
simphydir = ""
replicates = 10
workdir = ""
outdir = ""
nbiter = 10
fullmode = ""
use_species_tree = False
directories = ""
mode = "exec"
split_into_families = False
use_true_clusters = False
scores_mode = ""
cluster_init_mode = ""
cluster_sp_mode = ""
other_arg_str = ""
mode_string_custom = ""

args_to_pass_directly = []

for arg in sys.argv[1:]:
  if arg.startswith("--simphydir="):
    simphydir = arg.split("=")[1]
  if arg.startswith("--mode="):
    mode = arg.split("=")[1]
  if arg.startswith("--mult="):
    mult = arg.split("=")[1]
  if arg.startswith("--nbiter="):
    nbiter = int(arg.split("=")[1])
  if arg.startswith("--replicates="):
    replicates = int(arg.split("=")[1])
  if arg.startswith("--workdir="):
    workdir = arg.split("=")[1]
  if arg.startswith("--outdir="):
    outdir = arg.split("=")[1]
  if arg.startswith("--fullmode="):
    fullmode = arg.split("=")[1]
  if arg.startswith("--use_species_tree"):
    use_species_tree = True
  if arg.startswith("--use_true_clusters"):
    use_true_clusters = True
  if arg in  ["--compute_pctid", "--use_cache", "--use_clustal", "--use_speciestree_cache"]:
    args_to_pass_directly.append(arg)
  if arg.startswith("--timefile="):
    args_to_pass_directly.append(arg)
  if arg.startswith("--split_into_families"):
    split_into_families = True
  if arg.startswith("--directories="):
    directories = arg.split("=")[1]
  if arg.startswith("--scores_mode="):
    scores_mode = arg.split("=")[1]
  if arg.startswith("--cluster_init_mode="):
    cluster_init_mode = arg.split("=")[1]
  if arg.startswith("--cluster_sp_mode="):
    cluster_sp_mode = arg.split("=")[1]
  if arg.startswith("--other_args="):
    other_arg_str = arg.split("=", 1)[1]
  if arg.startswith("--mode_string="):
    mode_string_custom = arg.split("=")[1]
	
mult = mult + ".0000000000"	




if simphydir == "":
	print("--simphydir must be set.")
	sys.exit()

if workdir == "":
	print("--workdir not set.  Will use simphydir")
	workdir = simphydir
else:
	if not os.path.exists(workdir):
		os.makedirs(workdir)
 
 
if outdir == "":
	outdir = workdir


	
modestr = "default"
if use_true_clusters:
	modestr = "true_clusters"
if fullmode != "":
	modestr = fullmode
if mode_string_custom != "":
	modestr = mode_string_custom
	
	
	
class SummaryInfo:
	
	nbLoci = 0
	species_tree_str = 0
	locus_info = {}
	
	cur_locus = ""
	speciestree_file = ""
	
	def __init__(self):
		pass
		
	def set_species_tree(self, species_tree_str):
		self.species_tree = species_tree_str
		
	def set_cur_locus(self, name):
		self.cur_locus = name
		self.locus_info[name] = {}
		
	def set_cur_locus_relations(self, relations_str):
		self.set_locus_info(self.cur_locus, "RELATIONS", relations_str)
		
	def set_cur_locus_clusters(self, clusters_str):
		self.set_locus_info(self.cur_locus, "CLUSTERS", clusters_str)
		
	def set_locus_info(self, locus_name, info_key, info_value):
		self.locus_info[locus_name][info_key]  = info_value
		
	def parse_summary_file(self, summary_filename):	
	
		indir = os.path.dirname(summary_filename)
		with open(summary_filename) as summary_file:
					
			for line in summary_file:
				line = line.replace("\n", "").replace("\r", "")
				
				if line != "":
					pz = line.split("=")
					key = pz[0]
					value = pz[1]
				
				#-------------------------------------------------------------------------------------------------------
				if key == "SPECIESTREE":
					self.set_species_tree(value)
					with open(join(indir, "speciestree.nw"), "w") as sp_outfile:
						sp_outfile.write(value)
					self.speciestree_file = join(indir, "speciestree.nw")
				#-------------------------------------------------------------------------------------------------------
				#Here a locus is actually a gene tree.  
				elif key == "LOCUS":	
					self.set_cur_locus(value)
				#-------------------------------------------------------------------------------------------------------
				else:
					
					if key == "CLUSTERS":
						#make a clusters file called true_clusters_[LOCUS].clusters
						self.set_cur_locus_clusters(value)
						clusters = value.split(";;")
					
						with open( join(indir, "true_clusters_" + self.cur_locus + ".clusters"), "w" ) as outclusters:
							for c in clusters:
								outclusters.write("<GROUP>\n")
								outclusters.write(c.replace(" ", "\n") + "\n")						
								outclusters.write("</GROUP>\n")
					elif key == "RELATIONS":
						self.set_cur_locus_relations(value)

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
	
	
def get_clusters_from_string(strcl, name_suffix = ""):
	clusters = []
	
	pz = strcl.split(";;")
	for p in pz:
		if p != "":
			cz = p.split()
			cluster = set()
			for c in cz:
				cluster.add(c + name_suffix)
			clusters.append(cluster)
	
	return clusters
				

def compare_clusters(cl1, cl2):
	
	n1 = len(cl1)
	edges = []
	edge_weights = {}
	
	
	for c in range(len(cl1)):
		
		for d in range(len(cl2)):
			
			inter = len(cl1[c].intersection(cl2[d]))
			if inter > 0:
				edges.append( (c, (d + 1)*10*n1, inter) )
				edge_weights[ (c, (d + 1)*10*n1) ] = inter
			
				
	match = mwmatching.maxWeightMatching(edges)
	
	
	w = 0
	for c in range(len(cl1)):
		mate = match[c]
		
		if mate >= 0:
			w += edge_weights[ (c, mate) ]
			
	return w
		
  
def get_relations_from_string(relations_string):
	rels = {}
	vals = relations_string.split(";;")
	
	for v in vals:
		pz = v.split("\t")
		
		if len(pz) == 3:
			rels[pz[0] + ";;" + pz[1]] = pz[2]
			rels[pz[1] + ";;" + pz[0]] = pz[2]
	
	return rels
  
  
def get_relations_from_file(filename, waitForRelationsString = True):
	rels = {}
	done = False
	with open(filename) as f:
		for line in f:
			if done:
				break
			line = line.replace("\n", "").replace("\r", "")
			
			if waitForRelationsString:
				if line.startswith("RELATIONS="):
					pz = line.split("=")
					rels = get_relations_from_string(pz[1])
					done = True
			else:
				rels = get_relations_from_string(line)
				done = True
	
	return rels
  


def get_relations_comparison(true_rels, my_rels, replicate_index):
	
	nbTotal = 0
	nbCommon = 0
	nbDiff = 0
	
	nbOKOrtho = 0
	nbBadOrtho = 0
	nbOrtho = 0
	nbUnknownOrtho = 0
	nbOKPara = 0
	nbBadPara = 0
	nbUnknownPara = 0
	nbPara = 0
	
	
	for key in true_rels:
		r1 = true_rels[key]
		
		if key in my_rels:
			r2 = my_rels[key]
		else:
			px = key.split(";;")
			key2 = px[0] + str(replicate_index) + ";;" + px[1] + str(replicate_index)
			
			
			
			if not key2 in my_rels:
				print("WARNING: " + key2 + " NOT IN my_rels")
				r2 = "Unknown"
			else:
				r2 = my_rels[key2]
		
		if (r1 == r2):
			nbCommon += 1
		else:
			nbDiff += 1
		
		nbTotal += 1
		
		if r1 == "Orthologs":
			nbOrtho += 1
			if (r1 == r2):
				nbOKOrtho += 1
			elif r2 == "Paralogs":
				nbBadPara += 1
			elif r2 == "Unknown":
				nbUnknownOrtho += 1
		elif r1 == "Paralogs":
			nbPara += 1
			if (r1 == r2):
				nbOKPara += 1
			elif r2 == "Orthologs":
				nbBadOrtho += 1
			elif r2 == "Unknown":
				nbUnknownPara += 1
	
	vals = {}
	vals["nbtotal"] = nbTotal
	vals["nbcommon"] = nbCommon
	vals["nbdiff"] = nbDiff
	vals["nbokortho"] = nbOKOrtho
	vals["nbbadortho"] = nbBadOrtho
	vals["nbunknownortho"] = nbUnknownOrtho
	vals["nbortho"] = nbOrtho
	vals["nbokpara"] = nbOKPara
	vals["nbbadpara"] = nbBadPara
	vals["nbunknownpara"] = nbUnknownPara
	vals["nbpara"] = nbPara
	
	vals["pctcommon"] = float(nbCommon)/float(nbTotal)
	
	if (nbOKOrtho + nbBadOrtho == 0):
		vals["precision_ortho"] = 1
	else:
		vals["precision_ortho"] = float(nbOKOrtho) / float(nbOKOrtho + nbBadOrtho)
		
	if (nbOrtho == 0):
		vals["recall_ortho"] = 1
	else:
		vals["recall_ortho"] = float(nbOKOrtho) / float(nbOrtho)
		
	if (nbOKPara + nbBadPara == 0):
		vals["precision_para"] = 1
	else:
		vals["precision_para"] = float(nbOKPara) / float(nbOKPara + nbBadPara)
		
	if (nbPara == 0):
		vals["recall_para"] = 1
	else:
		vals["recall_para"] = float(nbOKPara) / float(nbPara)
		
	if (nbOrtho + nbPara == 0):
		vals["accuracy"] = 1
	else:
		vals["accuracy"] = float(nbOKPara + nbOKOrtho) / float(nbOrtho + nbPara)
  
	return vals

def add_relations(r1, r2, gene_suffix = ""):
	for key in r2:
		if gene_suffix != "":
			pz = key.split(";;")
			newkey = pz[0] + gene_suffix + ";;" + pz[1] + gene_suffix
		else:
			newkey = key
		r1[newkey] = r2[key]
	return r1
	
def add_clusters(c1, c2, gene_suffix = ""):
	for c in c2:
		if gene_suffix != "":
			cp = set()
			for g in c:
				cp.add(g + gene_suffix)
			c1.append(cp)
		else:
			c1.append(c)
	return c1

	
def get_results_header_line(namelen):
	return "NAME".ljust(namelen + 5) + "PREC_ORTHO".ljust(15) + "RECALL_ORTHO".ljust(15) + "PREC_PARA".ljust(15) + "RECAL_PARA".ljust(15) + "GOOD/NBRELS".ljust(15) + " ACC".ljust(15) + "CLSCORE/NBGENES".ljust(15) + "CLSCORE".ljust(15)
	
#def get_results_line(name, namelen, nbcommonrels, nbrels, clusters_score, nbgenes):
	#pctcommon = float(nbcommonrels)/float(nbrels)
	#str_pctr = "{0:.5f}".format(pctcommon)
	#pctcl = float(clusters_score)/float(nbgenes)
	#str_pctc = "{0:.5f}".format(pctcl)
	#line = name.ljust(namelen + 5) + (str(int(nbcommonrels/2)) + "/" + str(int(nbrels/2))).ljust(15) + str_pctr.ljust(15)
	#line += (str(clusters_score) + "/" + str(nbgenes)).ljust(15) + str_pctc
	
	
	
	#return line
	
def get_results_line(name, namelen, precision_ortho, recall_ortho, precision_para, recall_para, nbcommonrels, nbrels, accuracy, clusters_score, nbgenes):
	
	pctcl = float(clusters_score)/float(nbgenes)
	str_pctc = "{0:.5f}".format(pctcl)
	
	line = name.ljust(namelen + 5) + "{0:.5f}".format(precision_ortho).ljust(15) + "{0:.5f}".format(recall_ortho).ljust(15)
	line += "{0:.5f}".format(precision_para).ljust(15) + "{0:.5f}".format(recall_para).ljust(15)
	line += (str(int(nbcommonrels/2)) + "/" + str(int(nbrels/2))).ljust(15) + "{0:.5f}".format(accuracy).ljust(15)
	line += (str(clusters_score) + "/" + str(nbgenes)).ljust(15) + str_pctc
	return line

def process_directory(simphydir, fullmode, workdir, replicates, mult, nbiter, outdir, use_species_tree, results_outfile = ""):
	
	global scores_mode
	global cluster_init_mode
	global cluster_sp_mode
	global other_arg_str
	global modestr
	global args_to_pass_directly
	
	filestr = ""
	files = []
	maxfilename_length = 0	#used for display later on

	for i in range(replicates):
		if filestr != "":
			filestr += ","
		
		stri = str(i + 1)
		if i + 1 < 10:
			stri = "0" + stri
		f = join(simphydir, "mult/alignment_" + stri + "_" + mult + "_TRUE.phy")
		
		if len(os.path.basename(f)) > maxfilename_length:
			maxfilename_length = len(os.path.basename(f))
		
		files.append(f)
		filestr += f

	summary_info = SummaryInfo()
	summary_info.parse_summary_file(join(simphydir, "mult/summary.txt"))

	inferred_relations = {}
	inferred_clusters = []
	if fullmode == "":
		strarg = ""
		
		cmd = "python3 hyppo.py --infiles=" + filestr + " --workdir=" + workdir + " --nbiter=" + str(nbiter) + " --outdir=" + outdir + " --rename_genes"
		if use_species_tree:
			cmd += " --speciestree=" + summary_info.speciestree_file
		
		if scores_mode != "":
			cmd += " --scores_mode=" + scores_mode
		
		if cluster_init_mode != "":
			cmd += " --cluster_init_mode=" + cluster_init_mode
			
		if cluster_sp_mode != "":
			cmd += " --cluster_sp_mode=" + cluster_sp_mode
		for a in args_to_pass_directly:
			cmd += " " + a
		if other_arg_str != "":
			cmd += " --other_args=\"" + other_arg_str + "\""
		
		
		#if we use true clusters, we pass a special true clusters string to hyppo
		if use_true_clusters:
			for i in range(replicates):
				stri = str(i + 1)
				if i + 1 < 10:
					stri = "0" + stri
				if strarg != "":
					strarg += ","
				strarg += files[i] + "=" + join(simphydir, "mult", "true_clusters_" + stri + ".clusters") + ","
				strarg += files[i] + "_suffix=" + str(i)
		
		
		cmd += " --mode_string=" + modestr
		
		if strarg != "":
			cmd += " --args=\"" + strarg + "\""
		
		print("EXEC " + cmd)
		os.system(cmd)
		
		for i in range(replicates):
			f = files[i]
			relsfile = join(outdir, os.path.basename(f) + "." + modestr + ".relations")
			tmp_rel = get_relations_from_file(relsfile, True)
			clustersfile = join(outdir, os.path.basename(f) + "." + modestr + ".clusters")
			tmp_cl = get_clusters_from_file(clustersfile)
			inferred_relations = add_relations(inferred_relations, tmp_rel)
			inferred_clusters = add_clusters(inferred_clusters, tmp_cl)
			
	else:
		if not split_into_families:
			cmd = "python3 hyppo.py --fullmode=" + fullmode + " --infiles=" + filestr + " --nbiter=" + str(nbiter) + " --workdir=" + workdir + " --outdir=" + outdir + " --rename_genes" 
			if use_species_tree:
				cmd += " --speciestree=" + summary_info.speciestree_file
			if other_arg_str != "":
				cmd += " --other_args=\"" + other_arg_str + "\""
			
			for a in args_to_pass_directly:
				cmd += " " + a
			
			
			print("EXEC " + cmd)
			os.system(cmd)
			
			relsfile = join(outdir, fullmode + ".relations")
			inferred_relations = get_relations_from_file(relsfile, False)
			clustersfile = join(outdir, fullmode + ".clusters")
			inferred_clusters = get_clusters_from_file(clustersfile)
		else:
			for i in range(replicates):
				f = files[i]
				modestr_i = modestr + str(i + 1)
								
				workdir_i = join(workdir, str(i + 1))
				outdir_i = workdir_i
				cmd = "python3 hyppo.py --fullmode=" + fullmode + " --mode_string=" + modestr_i + " --nbiter=" + str(nbiter) + " --infiles=" + f + " --workdir=" + workdir_i + " --outdir=" + outdir_i
				
				if use_species_tree:
					cmd += " --speciestree=" + summary_info.speciestree_file
				if other_arg_str != "":
					cmd += " --other_args=\"" + other_arg_str + "\""
				for a in args_to_pass_directly:
					cmd += " " + a
				
				print("EXEC " + cmd)
				os.system(cmd)
		
				relsfile = join(outdir_i, modestr_i + ".relations")
				tmp_rel = get_relations_from_file(relsfile, False)
				inferred_relations = add_relations(inferred_relations, tmp_rel, str(i))
				
				clustersfile = join(outdir_i, modestr_i + ".clusters")
				tmp_cl = get_clusters_from_file(clustersfile)
				inferred_clusters = add_clusters(inferred_clusters, tmp_cl, str(i))
				
	
	print("\n\nRESULTS \n")
	
	outstr = ""
	outstr += get_results_header_line(maxfilename_length) + "\n"

	totals = {}
	totals["nbcommon"] = 0
	totals["nbtotal"] = 0
	totals["nbgenes"] = 0
	totals["cluster_scores"] = 0

	#in fullmode, we get ALL the relations in one file
	#if fullmode != "":
	#	relsfile = join(workdir, fullmode + ".relations")
	#	inferred_relations = get_relations_from_file(relsfile, False)
	#	clustersfile = join(workdir, fullmode + ".clusters")
	#	inferred_clusters = get_clusters_from_file(clustersfile)

	
	
	for i in range(replicates):
		f = files[i]
		
		stri = str(i + 1)
		if i + 1 < 10:
			stri = "0" + stri
		
		#in not fullmode, we get the relations per family
		#if fullmode == "":
		#	relsfile = join(outdir, os.path.basename(f) + ".default.relations")
		#	inferred_relations = get_relations_from_file(relsfile, True)
		#	clustersfile = join(outdir, os.path.basename(f) + ".default.clusters")
		#	inferred_clusters = get_clusters_from_file(clustersfile)
		
		true_relations = get_relations_from_string(summary_info.locus_info[stri]["RELATIONS"])
		
		vals = get_relations_comparison(true_relations, inferred_relations, i)
		
		#totals["nbcommon"] = totals["nbcommon"] + vals["nbcommon"]
		#totals["nbtotal"] = totals["nbtotal"] + vals["nbtotal"]
		
		
		
		#handle clusters
		true_clusters = get_clusters_from_string(summary_info.locus_info[stri]["CLUSTERS"], str(i))
		score = compare_clusters(true_clusters, inferred_clusters)
		
		nbgenes = 0
		for c in true_clusters:
			nbgenes += len(c)
		pctcl = float(score)/float(nbgenes)
		
		for v in vals:
			if v not in totals:
				totals[v] = 0
			totals[v] = totals[v] + vals[v]
		totals["cluster_scores"]  = totals["cluster_scores"] + score
		totals["nbgenes"] = totals["nbgenes"] + nbgenes
		
		
		#outstr += get_results_line(os.path.basename(f), maxfilename_length, vals["nbcommon"], vals["nbtotal"], score, nbgenes) + "\n"
		outstr += get_results_line(os.path.basename(f), maxfilename_length, vals["precision_ortho"], vals["recall_ortho"], vals["precision_para"], vals["recall_para"], vals["nbcommon"], vals["nbtotal"], vals["accuracy"], score, nbgenes) + "\n"
	
	#hopfully nothing will divide by 0
	totals["precision_ortho"] = float(totals["nbokortho"]) / float(totals["nbokortho"] + totals["nbbadortho"])
	totals["recall_ortho"] = float(totals["nbokortho"]) / float(totals["nbortho"])
	totals["precision_para"] = float(totals["nbokpara"]) / float(totals["nbokpara"] + totals["nbbadpara"])
	totals["recall_para"] = float(totals["nbokpara"]) / float(totals["nbpara"])
	totals["accuracy"] = float(totals["nbokpara"] + totals["nbokortho"]) / float(totals["nbortho"] + totals["nbpara"])
	
	
	#pcttotal = float(totals["nbcommon"])/float(totals["nbtotal"])
	
	#outstr += get_results_line("TOTAL", maxfilename_length, totals["nbcommon"], totals["nbtotal"], totals["cluster_scores"], totals["nbgenes"]) + "\n"
	outstr += get_results_line("TOTAL", maxfilename_length, totals["precision_ortho"], totals["recall_ortho"], totals["precision_para"], totals["recall_para"], totals["nbcommon"], totals["nbtotal"], totals["accuracy"], totals["cluster_scores"], totals["nbgenes"]) + "\n"
	
	
	print(outstr)
	
	if results_outfile != "":
		fout = open(results_outfile, 'w')
		fout.write(outstr)
		fout.close()
		
'''
if mode == "results":
	pctrel_per_family = 0
	pctcl_per_family = 0
	nb_families = 0
	pctrel_per_file = 0
	pctcl_per_file = 0
	nb_files = 0
	nbrels_correct = 0
	clscore_sum = 0
	nbgenes = 0
	nbrels = 0
	
	files = os.listdir(workdir)
	for file in files:
		if file.endswith("." + modestr + ".results"):
			
			f = open(join(workdir, file), 'r')
			for line in f:
				line = line.replace("\n", "")
				if line != "" and not line.startswith("TOTAL"):
					pz = line.split()
					
					pctrel_per_family += float(pz[2])
					pctcl_per_family += float(pz[4])
					
					nb_families += 1
				elif line.startswith("TOTAL"):
					pz = line.split()
					pctrel_per_file += float(pz[2])
					pctcl_per_file += float(pz[4])
					nb_files += 1
					
					qz = pz[1].split("/")
					nbrels_correct += int(qz[0])
					nbrels += int(qz[1])
					
					rz = pz[3].split("/")
					clscore_sum += int(rz[0])
					nbgenes += int(rz[1])
			
			f.close()

	print("AVG CLUSTER SCORE PER FAMILY\t" + str(float(pctcl_per_family)/float(nb_families)))
	print("AVG CLUSTER SCORE PER DATASET\t" + str(float(pctcl_per_file)/float(nb_files)))
	print("TOTAL CLUSTER SCORE\t\t" + str(clscore_sum) + "/" + str(nbgenes) + "\t" + str(float(clscore_sum) / float(nbgenes)))
		
	print("AVG RELATIONS PER FAMILY\t" + str(float(pctrel_per_family)/float(nb_families)))
	print("AVG RELATIONS PER DATASET\t" + str(float(pctrel_per_file)/float(nb_files)))
	print("TOTAL RELATIONS\t\t\t" + str(nbrels_correct) + "/" + str(nbrels) + "\t" + str(float(nbrels_correct) / float(nbrels)))
'''

if mode == "results":
	precision_ortho_per_family = 0
	precision_para_per_family = 0
	recall_ortho_per_family = 0
	recall_para_per_family = 0
	clscore_per_family = 0
	acc_per_family = 0
	
	precision_ortho_per_file = 0
	precision_para_per_file = 0
	recall_ortho_per_file = 0
	recall_para_per_file = 0
	clscore_per_file = 0
	acc_per_file = 0

	nb_families = 0
	clscore_sum = 0
	nbgenes = 0
	nbrels = 0
	nbrels_correct = 0
	nb_files = 0
	
	
	files = os.listdir(workdir)
	for file in files:
		if file.endswith("." + modestr + ".results"):
			
			f = open(join(workdir, file), 'r')
			for line in f:
				line = line.replace("\n", "")
				if line != "" and not line.startswith("TOTAL") and not line.startswith("NAME"):
					pz = line.split()
					
					precision_ortho_per_family += float(pz[1])
					precision_para_per_family += float(pz[3])
					recall_ortho_per_family += float(pz[2])
					recall_para_per_family += float(pz[4])
					acc_per_family += float(pz[6])
					clscore_per_family += float(pz[8])
					
					nb_families += 1
				elif line.startswith("TOTAL"):
					pz = line.split()
					precision_ortho_per_file += float(pz[1])
					precision_para_per_file += float(pz[3])
					recall_ortho_per_file += float(pz[2])
					recall_para_per_file += float(pz[4])
					acc_per_file += float(pz[6])
					clscore_per_file += float(pz[8])
					
					
					qz = pz[5].split("/")
					nbrels_correct += int(qz[0])
					nbrels += int(qz[1])
					
					rz = pz[7].split("/")
					clscore_sum += int(rz[0])
					nbgenes += int(rz[1])
					
					nb_files += 1
			
			f.close()

	
	print("AVG CLUSTER SCORE PER FAMILY\t" + str(float(clscore_per_family)/float(nb_families)))
	print("AVG CLUSTER SCORE PER DATASET\t" + str(float(clscore_per_file)/float(nb_files)))
	print("TOTAL CLUSTER SCORE\t\t" + str(clscore_sum) + "/" + str(nbgenes) + "\t" + str(float(clscore_sum) / float(nbgenes)))
		
	print("AVG PREC ORTHO,PREC PARA PER FAMILY\t" + str(float(precision_ortho_per_family)/float(nb_families)) + "   " + str(float(precision_para_per_family)/float(nb_families)))
	print("AVG RECALL ORTHO,RECALL PARA PER FAMILY\t" + str(float(recall_ortho_per_family)/float(nb_families)) + "   " + str(float(recall_para_per_family)/float(nb_families)))
	print("AVG ACCURACY PER FAMILY\t" + str(float(acc_per_family)/float(nb_families)))
	
	print("AVG PREC ORTHO,PREC PARA PER DATASET\t" + str(float(precision_ortho_per_file)/float(nb_files)) + "   " + str(float(precision_para_per_file)/float(nb_files)))
	print("AVG RECALL ORTHO,RECALL PARA PER DATASET\t" + str(float(recall_ortho_per_file)/float(nb_files)) + "   " + str(float(recall_para_per_file)/float(nb_files)))
	print("AVG ACCURACY PER DATASET\t" + str(float(acc_per_file)/float(nb_files)))
	
	
	print("TOTAL RELATIONS\t\t\t" + str(nbrels_correct) + "/" + str(nbrels) + "\t" + str(float(nbrels_correct) / float(nbrels)))


else:
	if directories != "":
		dz = directories.split(":")
		start = int(dz[0])
		end = int(dz[1])
		
		
		
		i = start
		while i <= end:
			strdir = str(i)
			if i < 10:
				strdir = "0" + strdir
		
			l_simphydir = join(simphydir, strdir)
			l_workdir = join(workdir, strdir)
			l_outdir = join(outdir, strdir)
			resout = join(workdir, strdir + "." + modestr + ".results")
			process_directory(l_simphydir, fullmode, l_workdir, replicates, mult, nbiter, l_outdir, use_species_tree, resout)
		
			i += 1
			
	else:
		process_directory(simphydir, fullmode, workdir, replicates, mult, nbiter, outdir, use_species_tree)