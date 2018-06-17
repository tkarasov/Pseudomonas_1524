#usr/bin/py
from io import StringIO
from Bio.Phylo.NewickIO import Parser
from ete3 import *
#reorder
#output from raxml bootstrapping was in bootstraps where ete3 seems to only accept percentage. Changed this in ape, printed out then read in again
raxml=Tree('/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/phylogeny/RAxML_bipartitions.boot100_divided.nwk', format=2)

panx=Tree("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/all16s_verified_25_7_2017_1524.newick", format=1)
keep=[rec for rec in raxml.get_leaf_names() if rec.replace('plate', 'p').split('.annotation')[0] in panx.get_leaf_names()]
raxml.prune(keep)
for node in raxml:
    if node.is_leaf():
        temp=node.name.replace("plate",'p').split(".annotation")[0]
        node.name=temp


 descendants={}
 for node in raxml.traverse():
     if node.is_leaf()==False:
        desc=[leaf.name for leaf in node.get_leaves()]
        descendants[desc]=node.name
