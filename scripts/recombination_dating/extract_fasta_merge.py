
from Bio import SeqIO
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
import os
import subprocess


#this script prepares the core genes for a group of strains for clonalframeML. At present it is set up to do this for ten strains


def search_strains(filename, strains):
    rec=list(SeqIO.parse(open(str(core_gene_folder)+"/"+str(filename)), "fasta"))
    keep=[fasta for fasta in rec if fasta.id.split(".annotation")[0] in strains]
    print filename
    if len(keep)==len(strains):
        for fasta in keep:
            rec_dict[fasta.id.split(".annotation")[0]].append(fasta)
            
    return rec_dict



tree_old=Tree("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/all16s_verified_25_7_2017_1524.newick", format=1)

#collapsed=Tree("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/all16s_verified_25_7_2017_1524.newick0.0032231227493.nwk", format=1)
collapsed=Tree("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/phylogeny/RAxML_1821_bipartitions.boot100_divided_reordered.nwk0.0001_1524_only_OTU5.nwk", format=1)

classification={}
file="/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/ForTalia_Strains_OTUs_272717.txt"
with open(file) as f:
    for line in f:
        (key, val)=line.split()
        classification[key]=val

my_otu=[key for key in classification.keys() if classification[key]=='1']
#strains=[rec.name.split("_")[0] for rec in collapsed.get_leaves() if rec.name.split("_")[0] in my_otu]
#strains=[rec.replace('p','plate') for rec in strains]

#strains of interest:
#strains=("plate3.A5", "plate4.C9", "plate4.B2", "plate4.G8", "plate5.E6", "plate3.B12", "plate5.B1", "plate3.A9", "plate5.E1",  "plate5.A11")
strains=("plate25.C11", "plate11.F10", "plate2.H5", "plate26.B10", "plate6.E4", "plate11.E1", "plate8.A7", "plate7.E11", "plate11.C12", "plate3.F12", "plate1.D11", "plate8.B3", "plate9.B7", "plate7.D8", "plate12.B6", "plate6.C2", "plate13.B6", "plate4.A6", "plate3.A3", "plate4.B3", "plate20.G1", "plate9.A7", "plate24.B1", "plate5.E12", "plate12.D9", "plate5.D12", "plate3.A6", "plate20.D1", "plate3.D10", "plate20.C1", "plate20.E1", "plate5.D10", "plate22.A3", "plate5.A4", "plate7.D9", "plate24.H1", "plate2.C2", "plate3.C11", "plate8.C3", "plate3.B6", "plate13.H5", "plate1.A3", "plate12.G6", "plate12.G7", "plate4.G8", "plate3.F1", "plate2.H4", "plate26.G1", "plate23.B2", "plate23.A2", "plate5.B8", "plate13.H3", "plate26.B8", "plate6.F8", "plate6.B2", "plate22.D2", "plate26.G2", "plate22.E2", "plate26.D2", "plate5.F8", "plate3.C8", "plate22.G3", "plate1.H4", "plate11.H6", "plate11.H9", "plate12.B10", "plate24.H7", "plate8.G7", "plate26.D3", "plate11.G10", "plate9.B2", "plate9.D7", "plate27.F4", "plate3.F6", "plate23.C2", "plate26.F10", "plate11.A2", "plate1.A5", "plate23.A3", "plate8.H7", "plate3.B9", "plate13.F3", "plate26.H8", "plate24.G2", "plate2.D1", "plate9.C3", "plate8.E5", "plate26.B7", "plate4.H9", "plate22.A8", "plate13.C11", "plate6.H1", "plate9.G12", "plate13.D10", "plate22.C1", "plate25.B12", "plate7.F9", "plate5.H2", "plate11.H8", "plate21.C10", "plate13.H11", "plate24.G8", "plate23.B3", "plate6.A11", "plate23.A8", "plate7.E7", "plate23.D3")

#pull out core genes for concatenation
core_gene_folder="/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/coregenes_alignments"

rec_dict={}
for strain in strains:
    rec_dict[strain]=[]

#go through core
for filename in os.listdir(core_gene_folder):
    if filename.endswith("_na_aln.fa"):
        search_strains(filename, strains)
        
#now concatenate sequences
concat_dict={}
for strain in strains:
    concat_dict[strain]=''
    for fasta in rec_dict[strain]:
        concat_dict[strain]=concat_dict[strain]+str(fasta.seq)

#write out concatenated fasta            
handle=open("all_concat.fasta", "w")
for rec in concat_dict:
    handle.write(">"+rec+"\n"+concat_dict[rec]+"\n")

handle.close()
#SeqIO.write(list(SeqIO.parse(open("all_concat.fasta"), "fasta")), "all_concat.phy", "phylip")

#now write out tree
for node in tree_old.traverse():
    if node.is_leaf():
        temp=node.name.replace('p','plate')
        node.name=temp
        

tree_old.prune(strains, preserve_branch_length=T)
write.tree(tree_old, "all_concat.newick", formatrue=1)



for leaf in collapsed.get_leaves():
    temp=leaf.name.replace('p', 'plate').split("_")[0]
    leaf.name=temp
collapsed.write(outfile="concat_107.newick", format=1)

test=subprocess(Popen(["/ebio/abt6_projects9/Pseudomonas_diversity/Programs/bin/ClonalFrameML",  "/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe/concat_107.newick", "/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe/all_concat.fasta"]), stdout=subprocess.PIPE))
output = test.communicate()[0]
