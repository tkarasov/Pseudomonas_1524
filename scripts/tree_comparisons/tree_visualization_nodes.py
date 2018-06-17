#!/usr/bin/py
#the goal of this script is to assign gene content and trait values to a phylogeny
import ete2
import pandas as pd
import pickle
import os
import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
import seaborn as sns

def mylayout(node):
    # If node is a leaf, split the name and paste it back together to remove the underscore
    if node.is_leaf():
        # species name
        temp = node.name.split('_')
        sp = temp[0] + ' ' + temp[1]
        temp2 = ete.faces.TextFace(sp, fgcolor = colorGuide[node.name], fsize = 18, fstyle = 'italic')
        ete.faces.add_face_to_node(temp2, node, column=0)
    # make a circle for SLA, weighted by SLA values
    sla = ete.CircleFace(radius = slaGuide[node.name]*15, color = colorGuide[node.name], style = 'circle')
    sla.margin_left = 10
    sla.hz_align = 1
    ete.faces.add_face_to_node(sla, node, column = 0, position = 'aligned')
    # same with toughness
    toughness = ete.CircleFace(radius = toughGuide[node.name]*15, color = colorGuide[node.name], style = 'circle')
    toughness.margin_left = 40
    toughness.hz_align = 1
    ete.faces.add_face_to_node(toughness, node, column = 1, position = 'aligned')

def load_sorted_clusters(path):
    '''
    load gene clusters and sort 1st by abundance and then by clusterID
    '''
    geneClusterPath='%s%s'%(path,'protein_faa/diamond_matches/')
    diamond_geneCluster_dt=load_pickle(geneClusterPath+'orthamcl-allclusters_final.cpk')
    from operator import itemgetter
    # sort by decreasing abundance (-v[0], minus to achieve decreasing)
    # followed by increasing clusterID GC_00001
    return sorted(diamond_geneCluster_dt.iteritems(),
                   key=lambda (k,v): (-itemgetter(0)(v),k), reverse=False)

def load_pickle(filename):
    f = open(filename,"rb")
    p = cPickle.load(f)
    f.close()
    return(p)


 # load geneID_to_description_dict
 '''PTL261 /ebio/abt6_projects7/pangenome/data/Talia/PTL261'''
    path='/ebio/abt6_projects7/pangenome/data/Talia/PTL261'
    geneID_to_description_dict=pickle.load(open(path+'/geneID_to_description.cpk')) #keys are the strain and gene id. values are gene annotations

     # if disable_RNA_clustering==0:
        # load RNAID_to_description_file
       # geneID_to_description_dict.update(pickle.load(open(path+'\RNAID_to_description.cpk'))
    output_path='%s%s'%(path,'/geneCluster/')
    visualzition_path='%s%s'%(path,'vis/')
    gene_diversity_Dt=load_pickle(path+'/gene_diversity.cpk') #ID is GC cluster number

    ## sorted clusters
    sorted_genelist= load_sorted_clusters(path+"/") #this is a list with GC and all the genes that map to that cluster

    ## prepare geneId_Dt_to_locusTag

    #geneId_Dt_to_locusTag=defaultdict(list)
    #geneId_Dt_to_locusTag={v:k for k,v in locusTag_to_geneId_Dt.items()}

    ## load gain/loss event count dictionary
    dt_geneEvents= load_pickle(path+'/dt_geneEvents.cpk') #keys are a number probably related to GC number somehow

    #write_file_lst_json.write('['); begin=0
    ## sorted_genelist: [(clusterID, [ count_strains,[memb1,...],count_genes]),...]
    for gid, (clusterID, gene) in enumerate(sorted_genelist):
        strain_count, gene_list, gene_count = gene
        if begin==0:
            begin=1
        else:
            write_file_lst_json.write(',\n')


my_genes=pickle.load(open("dt_genePresence.cpk"))












    


cd /Users/tkarasov/Dropbox/germany_pathogen_collections/data_files_rmarkdown#/ebio/abt6_projects9/Pseudomonas_diversity/data/genome_processing/run_0034_Hiseq/PTL266/geneCluster
my_genes=pickle.load(open("dt_genePresence.cpk"))
for rec in my_genes:
    my_genes[rec]=[int(i) for i in my_genes[rec]]
#t= Phylo.read("tree_result.newick", 'newick')

#I need to build a data file that's compatible with roary and make into spreadsheat
build=pandas.DataFrame.from_dict(my_genes, orient="index").transpose()
build.to_csv(header=True, path_or_buf="./gene_content.txt", sep="\t")
old_t = Phylo.read("tree_result_no_p5_G10_p4_G7_collapse_0.001.newick", 'newick')
#this phylogeny needs to have one node removed
t = old_t.common_ancestor([rec for rec in old_t.get_terminals() if rec.name!='p4_H7'])


#this next part of code works very well
# Max distance to create better plots
mdist = max([t.distance(t.root, x) for x in t.get_terminals()])

# Load roary
#each column in a different genotype, each row a different gene
# Sort the matrix by the sum of strains presence
roary=pd.read_table("./gene_content.txt")
idx = roary.sum(axis=1).sort_values(ascending=False).index
roary_sorted = roary.ix[idx]
roary_sorted = roary_sorted[[x.name for x in t.get_terminals()]]

# Plot presence/absence matrix against the tree
with sns.axes_style('whitegrid'):
        fig = plt.figure(figsize=(17, 10))

        ax1=plt.subplot2grid((1,40), (0, 10), colspan=30)
        a=ax1.matshow(roary_sorted.T, cmap=plt.cm.Blues,
                   vmin=0, vmax=1,
                   aspect='auto',
                   interpolation='none',
                    )
        ax1.set_yticks([])
        ax1.set_xticks([])
        ax1.axis('off')

        ax = fig.add_subplot(1,2,1)
        ax=plt.subplot2grid((1,40), (0, 0), colspan=10, axisbg='white')

        fig.subplots_adjust(wspace=0, hspace=0)

        ax1.set_title('Presence/absence matrix\n(%d gene clusters)'%roary_sorted.shape[0])

        Phylo.draw(t, axes=ax, 
                       show_confidence=False,
                       label_func=lambda x: None,
                       xticks=([],), yticks=([],),
                       ylabel=('',), xlabel=('',),
                       xlim=(-mdist*0.1,mdist+mdist*0.1),
                       axis=('off',),
                       title=('Tree\n(%d genotypes)'%roary_sorted.shape[1],),
                       do_show=False,
                      )
        plt.savefig('pangenome_0.001_cluster_matrix.pdf', dpi=300)
        plt.clf()

    # Plot the pangenome pie chart
    plt.figure(figsize=(10, 10))

    core     = roary[(roary.sum(axis=1) >= roary.shape[1]*0.99) & (roary.sum(axis=1) <= roary.shape[1]     )].shape[0]
    softcore = roary[(roary.sum(axis=1) >= roary.shape[1]*0.95) & (roary.sum(axis=1) <  roary.shape[1]*0.99)].shape[0]
    shell    = roary[(roary.sum(axis=1) >= roary.shape[1]*0.15) & (roary.sum(axis=1) <  roary.shape[1]*0.95)].shape[0]
    cloud    = roary[roary.sum(axis=1)  < roary.shape[1]*0.15].shape[0]

    total = roary.shape[0]
    









    
    
# load data
traits = pd.read_csv('/Users/Nate/Documents/FIU/Research/Invasion_TraitPhylo/Data/plantTraits.csv')
SERCphylo = ete.Tree('/Users/Nate/Documents/FIU/Research/SERC_Phylo/SERC_Nov1-2013.newick.tre')
 
#### TRAIT CLEANUP ####
# put an underscore in trait species
traits['species'] = traits['species'].map(lambda x: x.replace(' ', '_'))
# pull out the relevant traits and only keep complete cases
traits = traits[['species', 'introduced', 'woody', 'SLA', 'seedMass', 'toughness']]
traits = traits.dropna()
 
# next, prune down the traits data
traitsPrune = traits[traits['species'].isin(SERCphylo.get_leaf_names())]
 
# prune the phylogeny so only species with traits are kept
SERCphylo.prune(traitsPrune['species'], preserve_branch_length = True)
 
# basic phylogenetic plot
SERCphylo.show()

# guide for color
cols = [['black', 'red'][x] for x in traitsPrune['introduced']]
colorGuide = dict(zip(traitsPrune['species'], cols))
# weights (scaled to 1)
slaGuide = dict(zip(traitsPrune['species'], traitsPrune['SLA']/traitsPrune['SLA'].max()))
toughGuide = dict(zip(traitsPrune['species'], traitsPrune['toughness']/traitsPrune['toughness'].max()))
seedGuide = dict(zip(traitsPrune['species'], traitsPrune['seedMass']/traitsPrune['seedMass'].max()))
# set the base style of the phylogeny with thick lines
for n in SERCphylo.traverse():
style = ete.NodeStyle()
style['hz_line_width'] = 2
style['vt_line_width'] = 2
style['size'] = 0
n.set_style(style)
