#the goal of this script is to assign gene content and trait values to a phylogeny
import ete2
import pandas
import pickle
import cPickle
import os
import numpy as np
from Bio import Phylo
#import matplotlib
#matplotlib.use('Agg') #this will force to plot figures
import matplotlib.pyplot as plt
import seaborn as sns


def load_sorted_clusters(path): #Taken from Wei
    '''
    load gene clusters and sort 1st by abundance and then by clusterID
    '''
    geneClusterPath='%s%s'%(path,'protein_faa/diamond_matches/')
    geneCluster_dt=load_pickle(geneClusterPath+'allclusters_postprocessed.cpk')
    from operator import itemgetter
    # sort by decreasing abundance (-v[0], minus to achieve decreasing)
    # followed by increasing strain count
    return sorted(geneCluster_dt.iteritems(),
                key=lambda (k,v): (-itemgetter(0)(v),itemgetter(2)(v)), reverse=False)
    #return sorted(geneCluster_dt.iteritems(),
    #            key=lambda (k,v): (-itemgetter(0)(v),itemgetter(2)(v)), reverse=False)

def load_pickle(filename):
    f = open(filename,"rb")
    p = cPickle.load(f)
    f.close()
    return(p)


# Plot presence/absence matrix against the tree. This requires the roary_sorted file and the phylogenetic tree
def plot_pres_abs(roary_sorted,t, matrix_pdf, red=False):
    #this next part of code works very well, except does not do well with entire pan genome
    if red==True:
        reduced_mat=roary_sorted.T[roary_sorted.T.columns[0:30000]]
    else:
        reduced_mat=roary_sorted.T
    # Max distance to create better plots
    mdist = max([t.distance(t.root, x) for x in t.get_terminals()])
    with sns.axes_style('whitegrid'):
        fig = plt.figure(figsize=(17, 10))
        ax1=plt.subplot2grid((1,40), (0, 10), colspan=30)
        a=ax1.matshow(reduced_mat, cmap=plt.cm.Blues,
        vmin=0, vmax=1,
        aspect='auto',
        interpolation='none')
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
        plt.savefig(matrix_pdf, dpi=600)
        plt.clf()

def subset_my_genes(sorted_genelist, effectors, my_genes):
    #this function is meant to subset the gene clusters to be visualized. It modifies my_genes to include only the desired subset
    subsetted={}
    GC=[line[0] for line in sorted_genelist]
    #first I need to find the indeces in the sorted_genelist for the desired effectors
    inds=[GC.index(line[0].strip(".faa")) for line in effectors]
    for rec in my_genes:
       new_subset=[my_genes[rec][i] for i in inds]
       subsetted[rec]=new_subset
    return subsetted


