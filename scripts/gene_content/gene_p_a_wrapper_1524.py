#!/usr/bin/py
import ete2
import pandas
import pickle
import cPickle
import os
import numpy as np
from Bio import Phylo
import matplotlib
matplotlib.use('Agg') #this will force to plot figures
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.insert(0, '/ebio/abt6_projects9/Pseudomonas_diversity/code/run_scripts/tree_comparisons')

from gene_p_a_matrix import *

def calc_abundance(is_patho, not_patho, gene1):
	p1 = sum(roary[is_patho].ix[gene1])
	percent_path = float(p1)/len(is_patho)
	p2 = sum(roary[not_patho].ix[gene1])
	percent_not = float(p2)/len(not_patho)
	return percent_path, percent_not




#the goal of this script is to assign gene content and trait values to a phylogeny

#Tree of strains
t=Phylo.read("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/phylogeny/RAxML_bestTree.1524_only_ML_all_sites_gap5_5_2018", "newick")#("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/phylogeny/RAxML_1821_bipartitions.boot100_divided_reordered.nwk1524.nwk", "newick")
t.ladderize()

#mapping of gene cluster to order
path='/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/'
sorted_genelist = load_sorted_clusters(path)

#using michael's effectors:
eff=[line.strip().split(',') for line in open("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/effectors/ROARY_effectors/hop-effectors/hop-effector-all-molecules-reference_40p-ident_e1e-5_60ident_60len/raw.csv").readlines()]
eff_pd=pandas.DataFrame([line[1:] for line in eff[1:]], columns=eff[0][1:], index=[line[0].replace('plate','p') for line in eff[1:]]).astype(float)>0
eff_fin=eff_pd.astype(int)
idx = eff_fin.sum(axis=0).sort_values(ascending=False).index
roary_sorted = eff_fin[idx]
roary_sorted = roary_sorted.ix[[x.name for x in t.get_terminals()]].transpose()
roary_sorted = roary_sorted.astype(float)
plot_pres_abs(roary_sorted,t, "/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/phylogeny/effector_gene_matrix_mgiolai.pdf")

#presence absence pattern for entire genome
my_genes=pickle.load(open("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/geneCluster/dt_genePresence.cpk"))
for rec in my_genes:
    temp=[gene for gene in my_genes[rec]]
    my_genes[rec]=temp
 
# Load roary
#each column in a different genotype, each row a different gene
# Sort the matrix by the sum of strains presence

roary=pandas.DataFrame.from_dict(my_genes)
roary=roary.astype(float)
idx = roary.sum(axis=1).sort_values(ascending=False).index
roary_sorted = roary.ix[idx]

#need to rename columns in roary_sorted to match phylogeny
temp=roary_sorted.columns
temp=[rec.replace("plate", "p").strip(".annotation") for rec in temp]
roary_sorted.columns=temp
roary_sorted = roary_sorted[[x.name for x in t.get_terminals()]]
#roary_sorted = roary_sorted.astype(float)

#create figure 
plot_pres_abs(roary_sorted,t, "/ebio/abt6_projects9/Pseudomonas_diversity/all_gene_matrix.pdf",red=True)
matrix_pdf="/ebio/abt6_projects9/Pseudomonas_diversity/poop.pdf"
reduced_mat=roary_sorted.T[roary_sorted.T.columns[0:30000]]
plt.imshow(reduced_mat, cmap=plt.cm.Blues,
        vmin=0, vmax=1,
        aspect='auto',
        interpolation='none')

plt.savefig(matrix_pdf, dpi=600)

 #60000 cuts off
 #20000 fine
 #30000 fine
 #35000 cuts off at 32000
 #40000 cuts off at 25000
 #50000 cuts off at 15000


#now determine which orthology groups are special to OTU5
classification={}
file="/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/ForTalia_Strains_OTUs_272717.txt"
with open(file) as f:
    for line in f:
        (key, val)=line.split()
        classification[key]=val

patho_dict={}
for rec in classification:
    if(classification[rec]=='1'):
        patho_dict[rec]="r"
    else:
        patho_dict[rec]="black"
is_patho=[rec for rec in patho_dict.keys() if patho_dict[rec]=='r']
not_patho=[rec for rec in patho_dict.keys() if patho_dict[rec]!='r' and rec != "Strain"]


#now look at abundance of genes in these different groups
path_roary=roary[is_patho]
not_roary = roary[not_patho]
sum_not = not_roary.sum(1)
sum_is = path_roary.sum(1)

absent_out = [rec for rec in sum_not.index if sum_not[rec]/float(len(not_patho))<0.1]
present_OTU5 = [rec for rec in sum_is.index if sum_is[rec]/float(len(is_patho))>0.9]
#the following are the genes conserved 
high_in_low_out = [rec for rec in present_OTU5 if rec in absent_out]

#now I want to determine the identity of those clusters:
#mapping = pickle.load(open('/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/geneID_to_geneSeqID.cpk'))
annot = pickle.load(open('/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/geneID_to_description.cpk'))

ids = [key for key in mapping if mapping[key] in high_in_low_out]
final_annot = {int(key.split("_")[1]):value for key,value in annot.iteritems()}
abs_annotations = [final_annot[rec]['annotation'] for rec in high_in_low_out]
hyp = [rec for rec in abs_annotations if 'hypothetical_protein' in rec]



