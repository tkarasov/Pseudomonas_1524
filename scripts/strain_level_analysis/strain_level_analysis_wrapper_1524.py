#!/usr/bin/py
import ete2
import pandas
import pickle
import cPickle
import os
import numpy as np
from Bio import Phylo
import matplotlib
import csv
matplotlib.use('Agg') #this will force to plot figures
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.insert(0, '/ebio/abt6_projects9/Pseudomonas_diversity/code/run_scripts/tree_comparisons')
sys.path.insert(0, '/ebio/abt6_projects9/Pseudomonas_diversity/code/run_scripts/demography_comparisons')
from strain_distribution import *
from gene_p_a_matrix import *
from collapse_tree import *
#this is a wrapper script that will generate many of the figures on strain distribution in the paper

def mylayout(node):
    node.img_style["size"] = 0

def circle_layout(node):
    node.img_style["size"] = 0
    #global node_counted
    if node.name in node_counted:
        colors=['Black']
        F= faces.CircleFace(100*node_counted[node.name]/float(sum(node_counted.values())),'Black')
        #F.border.width = None
        F.hz_align=1
        F.opacity = 0.5
        faces.add_face_to_node(F,node, 0, position="aligned")

def circle_layout_seasons(node):
    node.img_style["size"] = 0
    global node_counted
    global tree
    eyach1=('11_12_2015', 'Eyach')
    eyach2=('23_3_2016', 'Eyach')
    det1=('15_12_2015', 'Det-2')
    det2=('31_3_2016', 'Det-2')
    eyach1_dict=node_count_specific_time_location(tree, strain_dict, eyach1[0], eyach1[1])[0]
    eyach2_dict=node_count_specific_time_location(tree, strain_dict, eyach2[0], eyach2[1])[0]
    det1_dict=node_count_specific_time_location(tree, strain_dict, det1[0],det1[1])[0]
    det2_dict=node_count_specific_time_location(tree, strain_dict, det2[0], det2[1])[0]
    if node.name in node_counted:
        colors=['Black']
        try:
            Fe1= faces.CircleFace(100*eyach1_dict[node.name]/sum(eyach1_dict.values()),'Black')
            Fe1.hz_align=1
            Fe1.margin_left = 10
            faces.add_face_to_node(Fe1,node, 0, position='aligned')
        except KeyError:
            Fe1= faces.CircleFace(0,'Black')
            Fe1.hz_align=1
            Fe1.margin_left = 10
            faces.add_face_to_node(Fe1,node, 0,  position='aligned')
        try:
            Fe2= faces.CircleFace(100*eyach2_dict[node.name]/sum(eyach2_dict.values()),'Black')
            Fe2.hz_align=1
            Fe2.margin_left = 10
            faces.add_face_to_node(Fe2,node, 1,  position='aligned')
        except KeyError:
            Fe2= faces.CircleFace(0,'Black')
            Fe2.hz_align=1
            Fe2.margin_left = 10
            faces.add_face_to_node(Fe2,node, 1,  position='aligned')

            
################
#collapse tree
intree="RAxML_bestTree.1524_only_ML_all_sites_gap5_5_2018" #/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/phylogeny/RAxML_bestTree.1524_only_ML_gap5_5_2018" #RAxML_bestTree.1524_only_ML_gap5"
#panx tree
#intree ="/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/" #RAxML_bestTree.1524_only_ML_gap5"

alignment_all_sites="/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/phylogeny/all_concat_1524.fasta"
SNP_alignment="/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/phylogeny/all_concat_1524_snp.fasta"



dist=0.0001
percent_SNP=1#0.40879766490179903 #len(list(SeqIO.parse(open(SNP_alignment), "fasta"))[0].seq)*1.0/len(list(SeqIO.parse(open(alignment_all_sites), "fasta"))[0].seq)
adj_dist=0.001/percent_SNP
output=intree+str(dist)+"_1524.nwk"
t=Tree(intree, format=1)
t.set_outgroup(t.get_midpoint_outgroup()) 
for node in t.traverse():
    if node.is_leaf:
        temp=node.name.split(".annotation")[0].replace('plate','p')
        node.name=temp

ts = ete2.TreeStyle()
ts.mode = 'r'
ts.show_leaf_name = False
ts.layout_fn = mylayout
t.render(output+"original.pdf", tree_style=ts)#,tree_style=circular_style)
t.write(format=1, outfile=intree+"1524.nwk")
collapse(t, adj_dist) 

#Now just write to file
t.write(format=1, outfile=output)
ts = ete2.TreeStyle()
ts.mode = 'r'
ts.show_leaf_name = False
ts.layout_fn = mylayout
t.render(output+".pdf", tree_style=ts)

#read in collapsed tree 
tree=ete2.Tree("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/phylogeny/RAxML_bestTree.1524_only_ML_all_sites_gap5_5_20180.0001_1524.nwk", format=1)

'''midpoint root the tree'''
tree.set_outgroup(tree.get_midpoint_outgroup()) 
ts = ete2.TreeStyle()
ts.mode = 'r'
ts.show_leaf_name = False
ts.layout_fn = mylayout
ts.branch_vertical_margin = 4
#ts.force_topology = True
ts.show_scale = True
#tree_old.convert_to_ultrametric()
#tree_old.render("basic_phylogeny.png", w=183, units="mm",tree_style=ts)

ts = ete2.TreeStyle()
ts.mode = 'r'
ts.show_leaf_name = False
ts.layout_fn = mylayout
ts.branch_vertical_margin = 4
#ts.force_topology = True
ts.show_scale = True
#tree.convert_to_ultrametric()
tree.render("basic_phylogeny_collapsed.pdf", w=600, dpi=300, units="mm",tree_style=ts)


################
'''Now for incorporating meta data'''
#Original meta number of plant

# read in strain meta file
f= open("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/sequenced_strain_meta2_amend_30_10_2017.csv", 'rU')
reader = csv.reader(f, delimiter=',', dialect=csv.excel_tab)


#for now I need to change the id of the strain
strain_dict={}
for line in reader:
    temp=line[0].replace('plate_', 'p')+"."+line[1]
    strain_dict[temp]=line[1:]
    
#create figure with circle for the abundance of a strain overall
node_counted=node_count(tree, strain_dict)[0]
node_perc=node_count(tree, strain_dict)[1]

#global node_counted #format also requires node_counted
nstyle = ete2.NodeStyle()
nstyle["shape"] = "sphere"
nstyle["size"] = 1
nstyle["fgcolor"] = "black"
ts = ete2.TreeStyle()
ts.mode = 'r'
ts.show_leaf_name = False
ts.layout_fn = circle_layout
ts.branch_vertical_margin = 4
#ts.force_topology = True
ts.show_scale = False

#ts.layout_fn = phyparts_pie_layout(mynode, node_perc)#phyparts_pie_layout #LAYOUT PROBLEM HERE DUE TO NON_GLOBAL VARIABLE
for n in tree.traverse():
    n.set_style(nstyle)

#tree.convert_to_ultrametric()
tree.render("strain_pie_12_2017_poo_chart.pdf", w=600, dpi=300, units="mm", tree_style=ts) 

#create figure with circle for the abundance of a strain in each season and in each site. First thing is to make a node counted for each time point/location
coll_options=set([(line[1], line[2]) for line in strain_dict.values() if "NA" not in line[1] and "NA" not in line[2] and 'Collection' not in line[1]])
node_counted=node_count(tree, strain_dict)[0]
node_perc=node_count(tree, strain_dict)[1]
#global tree   #the formatting of figure requires tree, which is not global until now
ts.show_leaf_name = False
ts.layout_fn = circle_layout_seasons
ts.branch_vertical_margin = 4
#ts.force_topology = True
ts.show_scale =True
#tree.convert_to_ultrametric()
tree.render("strain_pie_12_2017_poo_poo_chart.pdf", w=183, units="mm",tree_style=ts) 

#generate bar chart of strain composition per LEAF!!!
node_info=build_node_info(tree, strain_dict)[0]
node_leaf=build_node_info(tree, strain_dict)[1]
counts=build_plant(node_info) #counts will count the number of each strain per plant

#node classification: is a node "pathogenic" or not
#Ill use p26.F7 and p2.B6 to classify pathogenic vs. not pathogenic
node1=[node.name for node in tree.get_leaves() if 'p26.F7' in node][0]
node2=[node.name for node in tree.get_leaves() if 'p2.B6' in node][0]
ancestor_leaves=[node.name for node in tree_old.get_common_ancestor(node1, node2).get_leaves()]



################
#classification of OTU5

#Read in Juliana Classification
#tree.unroot()
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
not_patho=[rec for rec in patho_dict.keys() if patho_dict[rec]!='r']
patho_node={}
for rec in tree.get_leaves():
    patho_node[rec.name]=patho_dict[rec.name.split("_")[0]]

#now output list and tree of only this classification
class_5=[rec for rec in classification if classification[rec]=='1']
tree.prune(rec.name for rec in tree.get_leaves() if rec.name.split("_")[0] in class_5)
R=tree.get_midpoint_outgroup()
tree.set_outgroup(R)
tree.write(format=1, outfile='/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/phylogeny/RAxML_bestTree.1524_only_ML_gap5_only_OTU5.nwk')
samp=[]
for rec in tree.get_leaves():
    samp.append(rec.name.split("_")[0])
outfile=open('/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/phylogeny/otu5_strains.txt','w')
strain_list=(', ').join([rec.replace('p', 'plate') for rec in samp])
print>> outfile,strain_list
outfile.close()


#visualize_bar(counts, tree)
#tree=ete2.Tree("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/phylogeny/RAxML_1821_bipartitions.boot100_divided_reordered.nwk0.0001_1524.nwk", format=1)
#visualize_bar_amended([line[1] for line in counts], tree)

################
#I want to look at the composition and probabilities of similarity
#leaf_counts=build_leaf(node_leaf)

tree=t


#first let's ask within a leaf what is the probability of two strains being the same.
#first I need a list of all of the strains
sort_counts=sorted([line[1] for line in counts], key=lambda lin:sum(lin.values()))
sort_count_df=DataFrame(columns=tree.get_leaf_names())
i=0
for line in sort_counts:
    sort_count_df.loc[i]=[0]*len(sort_count_df.columns)
    for key in line.keys():
        sort_count_df[key][i]=line[key]
    i=i+1
    


sort_count_df.to_csv("plant_strain_count.csv", index=True,sep=',')

#write pathogen_strains to file for R
sort_count_df[[rec for rec in sort_count_df.keys() if rec.split('_')[0] in is_patho]].to_csv("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/phylogeny/pathogen_nodes.csv", sep=',')


#can get probability in plant from sort_count_df, as well as probability within plant and number of trials

wit_plant= (within_plant(sort_count_df)) #When in a plant, its frequncy
acr_plant=(across_plant(sort_count_df)) #percentage of plants with given strain
sort_wit=[[rec, wit_plant[rec]] for rec in [line[0] for line in  acr_plant] if rec in wit_plant]
corr_wit_btwn=[[rec[0], rec[1], numpy.mean(wit_plant[rec[0]])] for rec in [line[0:2] for line in  acr_plant] if rec[0] in wit_plant]
corr_wit_btwn2=[[[rec,val,acr[1]] for rec in wit_plant for val in wit_plant[rec] for acr in acr_plant if acr[0]==rec]][0]
plt.scatter(x=[line[2] for line in corr_wit_btwn2],y=[line[1] for line in corr_wit_btwn2])
new=[]
for rec in sort_wit:
    for line in rec[1]:
        new.append([rec[0], line])


#plot the ranked abundances colored according to their pathogenic identity
acr_plant_f=pandas.DataFrame(acr_plant)
acr_plant_f['path']=pandas.Series([patho_node[line[0]] for line in acr_plant], index=acr_plant_f.index)
acr_plant_f.columns=names=["strain", "freq", "color"]
plt.clf()
plt.figure()
ax=sns.barplot(data=acr_plant_f, x="strain", y="freq",palette=acr_plant_f['color'])
widthbars=[.75]*len(ax.patches)
sns.set_style("whitegrid")
ax.set(xlabel='Bacterial Strains in order of abundance', ylabel='Fraction of plants in which strain found', xticklabels=[])
plt.savefig("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/Ta1524/phylogeny/percent_plant_in.pdf")
'''
win=[]
spring=[]

for leaf in tree.get_leaves():
    strains=[leaf.name.strip("--").split("_")]
    winter=[rec for rec in strains if '2015' in strain_dict[rec[0]][1]]
    print winter
    spr=[rec for rec in strains if '2016' in strain_dict[rec[0]][1]]
    if len(winter)>0:
        win.append(len(winter[0]))
    if len(spr)>0:
        spring.append(len(spr[0]))    

'''