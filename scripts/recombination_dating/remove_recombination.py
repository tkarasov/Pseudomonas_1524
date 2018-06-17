#!/usr/bin/py
#this file takes the output from clonalframe, removes the recomination events in the fasta files

#commandline for clonalframe


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
import os
import sys
import ete2
from ete3 import Tree
from Bio import Phylo
sys.path.insert(0, '/ebio/abt6_projects9/Pseudomonas_diversity/code/run_scripts/tree_comparisons')
sys.path.insert(0, '/ebio/abt6_projects9/Pseudomonas_diversity/code/run_scripts/demography_comparisons')
from strain_distribution import *
from gene_p_a_matrix import *
from collapse_tree import *

def mylayout(node):
    node.img_style["size"] = 0

    
def children_dict(clonal_tree):
    clonal_dict={}
    for node in clonal_tree.traverse():
        if node.is_leaf():
            clonal_dict[node.name]=node.name
        else:
            all_descendants=node.get_descendants()
            clonal_dict[node.name]=[rec.name for rec in all_descendants if "NODE" not in rec.name]
    return clonal_dict

def remove_regions(recom_events, strain, mutable):
    #my events
    strain_events=recom_events[strain]
    strain_seq=mutable[strain]
    for rec in strain_events:
        start=int(rec[0])-1
        end=int(rec[1])-1
        strain_seq[start:end]="-"*(end-start)
    mutable[strain]=strain_seq
    return mutable



#run clonalframe then next part extracts sequences from fasta which have not undergone recobmination. importation_status.txt has first and last events
#cd /ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe
import_file=[line.strip().split() for line in open("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe/clonalframe_107.out.importation_status.txt").readlines()]
all_concat=list(SeqIO.parse(open("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe/all_concat.fasta"), "fasta"))
mutable={rec.name:rec.seq.tomutable() for rec in all_concat}

clonal_tree=Tree("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe/clonalframe_107.out.labelled_tree.newick", format=1)
clonal_dict=children_dict(clonal_tree)
recom_events={}
for rec in clonal_tree.get_leaves():
    recom_events[rec.name]=[]
    
for rec in import_file:
    if 'Node' not in rec[0]:
        affected=clonal_dict[rec[0]]
        if type(affected)!=str:
            for sample in affected:`` 
                recom_events[sample].append([rec[1], rec[2]])
        else:
            recom_events[affected].append([rec[1], rec[2]])


for strain in recom_events:
    mutable=remove_regions(recom_events, strain, mutable)

#now write to file:
for rec in mutable:
    mutable[rec].name=rec

for rec in all_concat:
    rec.seq=Seq(str(mutable[rec.name]))

SeqIO.write(all_concat, "concat_107_no_recom.fasta", "fasta")

my_xml=[]
for rec in all_concat:
    temp=str(rec.seq)
    my_seq=Seq(temp, IUPAC.ambiguous_dna)
    rec.seq=my_seq

SeqIO.write(all_concat, "concat_107_no_recom.xml", "nexus")

#run /ebio/abt6_projects9/Pseudomonas_diversity/code/run_scripts/recombination_dating/raxml_no_recombine_107.sh


raxml_tree=Tree("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe/RAxML_bestTree.concat_107_no_recom", format=1)
R=raxml_tree.get_midpoint_outgroup()
raxml_tree.set_outgroup(R)
ts = ete2.TreeStyle()
ts.mode = 'r'
ts.show_leaf_name = False
ts.layout_fn = mylayout
raxml_tree.render("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe/raxml_107.pdf", tree_style=ts)#,tree_style=circular_style)
samp=[]
for rec in raxml_tree.get_leaves():
    samp.append(raxml_tree.get_distance(R, rec))
    
    
