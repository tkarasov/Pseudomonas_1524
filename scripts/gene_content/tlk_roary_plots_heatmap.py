#!/usr/bin/env python
# Modified version of roary script that cannot be disseminated for non-academic usage:
#original copyright (C) <2015> EMBL-European Bioinformatics Institute
#cmd  python tlk_roary_plots.py RAxML_bestTree.1524_only_ML_all_sites_gap5_5_2018 raw_transpose.csv --format pdf
#stored on cluster in /ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/effector_identification/toxin_effectors

#cmd  python tlk_roary_plots.py RAxML_bestTree.1524_only_ML_all_sites_gap5_5_2018 DC3000_print.csv --format pdf #this gives the DC3000 heatmap
#python tlk_roary_plots_heatmap.py RAxML_bestTree.1524_only_ML_all_sites_gap5_5_2018 viridiflava_print.csv --format pdf #this gives the DC3000 heatmap



def get_options():
    import argparse

    # create the top-level parser
    description = "Create plots from roary outputs"
    parser = argparse.ArgumentParser(description = description,
                                     prog = 'roary_plots.py')

    parser.add_argument('tree', action='store',
                        help='Newick Tree file', default='accessory_binary_genes.fa.newick')
    parser.add_argument('spreadsheet', action='store',
                        help='Roary gene presence/absence spreadsheet', default='gene_presence_absence.csv')

    parser.add_argument('--labels', action='store_true',
                        default=False,
                        help='Add node labels to the tree (up to 10 chars)')
    parser.add_argument('--format',
                        choices=('png',
                                 'tiff',
                                 'pdf',
                                 'svg'),
                        default='png',
                        help='Output format [Default: pdf]')
    parser.add_argument('-N', '--skipped-columns', action='store',
                        type=int,
                        default=1,
                        help='First N columns of Roary\'s output to exclude [Default: 14]')
    
    parser.add_argument('--version', action='version',
                         version='%(prog)s '+__version__)

    return parser.parse_args()


def colorbar_index(ncolors, cmap):
    cmap = cmap_discretize(cmap, ncolors)
    mappable=a #a=roary presence absence map
    mappable = cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(-0.5, ncolors+0.5)
    colorbar = plt.colorbar(mappable, fraction=0.9, pad=0.9, aspect=8, orientation='vertical') #pad=distance, aspect=heigth
    colorbar.set_ticks(np.linspace(0, ncolors, ncolors))
    #colorbar.set_ticklabels(range(ncolors))
    colorbar.set_label('% bit-score reference', size=11)
    colorbar.set_ticklabels(['0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'])
    

def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet. 
        N: number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    if type(cmap) == str:
        cmap = plt.get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki])
                       for i in xrange(N+1) ]
    # Return colormap object.
    return mcolors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)




if __name__ == "__main__":
    options = get_options()
    import matplotlib
    matplotlib.use("pdf")
    import matplotlib.pyplot as plt
    import seaborn as sns
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors

    sns.set_style('white')

    import os
    import pandas as pd
    import numpy as np
    from Bio import Phylo

    t = Phylo.read("RAxML_bestTree.1524_only_ML_all_sites_gap5_5_2018", 'newick')
    t.ladderize()

    # Max distance to create better plots
    mdist = max([t.distance(t.root, x) for x in t.get_terminals()])

    # Load roary
    roary = pd.read_table(options.spreadsheet,
                         sep=',',
                         low_memory=False)

    print "YES"
    # Set index (group name)
    roary.set_index('gene_ID', inplace=True)
    # Drop the other info columns
    roary.drop(list(roary.columns[:1-1]), axis=1, inplace=True)

    # Transform it in a presence/absence matrix (1/0)
    roary.replace('.{2,100}', 1, regex=True, inplace=True)
    roary.replace(np.nan, 0, regex=True, inplace=True)

    # Sort the matrix by the sum of strains presence ***DONT DO THIS FOR HRP_HRC
    #idx = roary.sum(axis=1).sort_values(ascending=False).index
    #roary_sorted = roary.ix[idx]

    #Do not sort the matrix
    #roary_sorted=roary

    #sort the matrix
    idx = roary.sum(axis=1).sort_values(ascending=False).index
    roary_sorted = roary#roary.ix[idx]

    # Pangenome frequency plot
    plt.figure(figsize=(7, 5))

    plt.hist(roary.sum(axis=1), roary.shape[1],
             histtype="stepfilled", alpha=.7)

    plt.xlabel('No. of genomes')
    plt.ylabel('No. of genes')

    sns.despine(left=True,
                bottom=True)
    plt.savefig('pangenome_frequency.%s'%"pdf", dpi=300)
    plt.clf()

    # Sort the matrix according to tip labels in the tree
    roary_sorted = roary_sorted[[x.name.replace('p', 'p') for x in t.get_terminals()]]

    # Plot presence/absence matrix against the tree
    with sns.axes_style('whitegrid'):
        fig = plt.figure(figsize=(17, 10))

        set_cmap='Blues'#'hot_r' #CHANGE THIS TO CHANGE COLORS
        #set_cmap='viridis'
        ax1=plt.subplot2grid((1,40), (0, 10), colspan=28)
        a=ax1.matshow(roary_sorted.T, cmap=plt.get_cmap(set_cmap), #color has to be the same as in colorbar below
                   vmin=0, vmax=100, #edited vmax from 50 to 100
                   aspect='auto',
                   interpolation='none',
                    )

        ax1.set_yticks([])
        ax1.set_xticks([])
        ax1.set_xticklabels(['']+roary_sorted.index)
        ax1.axis('off')

        #Colorbar
        ax3 = plt.subplot2grid((1,40), (0, 39), colspan=1)
        ax3.set_visible(False)
        colorbar_index(ncolors=11, cmap=plt.get_cmap(set_cmap))

        
        ax = fig.add_subplot(1,2,1)
        ax=plt.subplot2grid((1,40), (0, 0), colspan=10, axisbg='white')

        fig.subplots_adjust(wspace=0, hspace=0)
        ax1.set_title('Roary matrix\n(%d gene clusters)'%roary.shape[0])        
        
        

        if options.labels:
            fsize = 12 - 0.1*roary.shape[1]
            if fsize < 7:
                fsize = 7
            with plt.rc_context({'font.size': fsize}):
                Phylo.draw(t, axes=ax, 
                           show_confidence=False,
                           label_func=lambda x: str(x)[:10],
                           xticks=([],), yticks=([],),
                           ylabel=('',), xlabel=('',),
                           xlim=(-mdist*0.1,mdist+mdist*0.45-mdist*roary.shape[1]*0.001),
                           axis=('off',),
                           title=('Tree\n(%d strains)'%roary.shape[1],), 
                           do_show=False,
                          )
        else:
            Phylo.draw(t, axes=ax, 
                       show_confidence=False,
                       label_func=lambda x: None,
                       xticks=([],), yticks=([],),
                       ylabel=('',), xlabel=('',),
                       xlim=(-mdist*0.1,mdist+mdist*0.1),
                       axis=('off',),
                       title=('Tree\n(%d strains)'%roary.shape[1],),
                       do_show=False,
                      )




        plt.savefig('pangenome_matrix.%s'%options.format, dpi=300)
        plt.clf()

    # Plot the pangenome pie chart
    plt.figure(figsize=(10, 10))

    core     = roary[(roary.sum(axis=1) >= roary.shape[1]*0.99) & (roary.sum(axis=1) <= roary.shape[1]     )].shape[0]
    softcore = roary[(roary.sum(axis=1) >= roary.shape[1]*0.95) & (roary.sum(axis=1) <  roary.shape[1]*0.99)].shape[0]
    shell    = roary[(roary.sum(axis=1) >= roary.shape[1]*0.15) & (roary.sum(axis=1) <  roary.shape[1]*0.95)].shape[0]
    cloud    = roary[roary.sum(axis=1)  < roary.shape[1]*0.15].shape[0]

    total = roary.shape[0]
    
    def my_autopct(pct):
        val=int(round(pct*total/100.0))
        return '{v:d}'.format(v=val)

    a=plt.pie([core, softcore, shell, cloud],
          labels=['core\n(%d <= strains <= %d)'%(roary.shape[1]*.99,roary.shape[1]),
                  'soft-core\n(%d <= strains < %d)'%(roary.shape[1]*.95,roary.shape[1]*.99),
                  'shell\n(%d <= strains < %d)'%(roary.shape[1]*.15,roary.shape[1]*.95),
                  'cloud\n(strains < %d)'%(roary.shape[1]*.15)],
          explode=[0.1, 0.05, 0.02, 0], radius=0.9,
          colors=[(0, 0, 1, float(x)/total) for x in (core, softcore, shell, cloud)],
          autopct=my_autopct)
    plt.savefig('pangenome_pie.%s'%options.format, dpi=300)
    plt.clf()
