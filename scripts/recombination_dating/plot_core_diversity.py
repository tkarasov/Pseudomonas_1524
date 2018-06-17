from Bio import SeqIO
import numpy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import random
from matplotlib.pyplot import cm
import pickle

'''the goal of this script is to take an alignment file with and without recombination events and plot diversity along the alignment file. Run this on alignment files that have come out after clonalframeml corrections have been made'''

def make_binary(alignment):
    pd_alignment=pd.DataFrame(alignment)
    pd_SNP=pd.DataFrame(columns=pd_alignment.columns, index=["SNP"])
    pd_SNP.ix["SNP"]=pd_alignment.apply(is_SNP, axis=0)
    return pd_SNP

def is_SNP(pd_column):
    ref=set([rec for rec in pd_column if rec!="-"])
    is_other=[rec for rec in pd_column if rec!=ref and rec !="-"]
    #print is_other
    if len(ref)>1:
        return 1
    elif len(ref)==0:
        return 'nan'
    else:
        return 0


def slidingWindow(SNP_matrix,sequence,winSize=1000,step=500 ):
    """Returns a generator that will iterate through
    the defined chunks of input sequence.  Input sequence
    must be iterable."""
 
    # Verify the inputs
    try: it = iter(sequence)
    except TypeError:
        raise Exception("**ERROR** sequence must be iterable.")
    if not ((type(winSize) == type(0)) and (type(step) == type(0))):
        raise Exception("**ERROR** type(winSize) and type(step) must be int.")
    if step > winSize:
        raise Exception("**ERROR** step must not be larger than winSize.")
    if winSize > len(sequence):
        print len(sequence)
        print sequence
        raise Exception("**ERROR** winSize must not be larger than sequence length.")
 
    # Pre-compute number of chunks to emit
    numOfChunks = ((len(sequence)-winSize)/step)+1
 
    # Do the work
    for i in range(0,numOfChunks*step,step):
        position= numpy.mean(sequence[i:i+winSize])
        SNP_yes=sum([rec for rec in list(SNP_matrix[[z for z in range(i, i+winSize)]].ix['SNP']) if rec!='nan'])
        actual_sites=(len([rec for rec in list(SNP_matrix[[z for z in range(i, i+winSize)]].ix['SNP']) if rec!='nan'])*1.0)
        if actual_sites==0:
            SNP_yes=0
            actual_sites=1
        yield [position, SNP_yes/actual_sites]

def generate_windows(alignment):
    pd_SNP=make_binary(alignment)
    windows=slidingWindow( pd_SNP,sequence=pd_SNP.columns,winSize=1000,step=100)
    all_seg=[]
    for window in windows:
        all_seg.append(list(window))
    return all_seg


def single_genome( alignment_full, number=5):
    samples=random.sample(range(len(alignment_full)), number)
    ref=samples[0]
    rest=samples[1:]
    fig=plt.figure()
    fig = plt.figure(figsize=(18, 12))
    fig.patch.set_facecolor('w')
    fig.add_axes=True
    color=iter(cm.rainbow(numpy.linspace(0,1,number)))
    i=1
    for other in rest:
        c=next(color)
        align=[alignment_full[ref], alignment_full[other]]
        seg=generate_windows(align)
        x1=[rec[0] for rec in seg]
        y1=[rec[1] for rec in seg]
        mean_div=numpy.mean(y1)
        sub1=plt.subplot(2, 2, i)
        plt.plot(x1, y1, color=c)
        plt.axhline(y=mean_div, color='r', linestyle='-')
        plt.xlabel("Position in Core Genome (bp)")
        print ref
        print other
        plt.ylabel('Pairwise nucleotide diversity', fontsize=16)
        plt.xlim(0, max(x1))
        i=i+1
    fig.patch.set_facecolor('w')
    fig.add_axes=True
    fig.tight_layout()
    plt.savefig("pairwise_divergence_statistics.pdf", facecolor=fig.get_facecolor(), edgecolor='black', transparent=True)
    print "Reference is:" + str(ref)
    print "Rest are:" + str(rest)


def compare_all(alignment_full):
    ref=64
    keep_scores={}
    for other in range(len(alignment_full)):
        if other!=ref:
            align=[alignment_full[ref], alignment_full[other]]
            seg=generate_windows(align)
            keep_scores[alignment_full[ref].name]=seg
    return keep_scores
            


    
def plot_segregating_sites():
    fig=plt.figure()
    x1=[rec[0] for rec in all_seg]
    x2=[rec[0] for rec in all_seg_recom]
    y1=[rec[1] for rec in all_seg]
    y2=[rec[1] for rec in all_seg_recom]
    plt.plot(x1, y1, color='k')
    plt.plot(x2, y2, color='g')
    plt.xlim(0, max(x1))
    plt.xlabel('Position in Core Genome (bp)', fontsize=18)
    plt.ylabel('Fraction Sites Segregating (Sn)', fontsize=16)
    fig.savefig("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe/segregating_sites_107.pdf")

'''first look at full data'''    

alignment_full=list(SeqIO.parse(open("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe/all_concat.fasta"), "fasta"))
alignment_norecom=list(SeqIO.parse(open("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe/concat_107_no_recom.fasta"), "fasta"))
all_seg=generate_windows(alignment_full)
all_seg_recom=generate_windows(alignment_norecom)

#plot segretgating sites plot
plot_segregating_sites()



'''now do pairwise comparisons'''
single_genome( alignment_full, number=5)
#ones used in the paper are: 
#ref 59 'plate11.C12'
#plot 1: 36 'plate9.C3'
#plot 2: 66 'plate6.F8'
#plot 3: 17 'plate7.F9'
#plot 4: 101 'plate1.A5'




#distance on tree vs gene sharing vs 

def compare_gene_content(strain1, strain2):
    num_dif=len([i for i in xrange(len(strain1)) if strain1[i] != strain2[i]])
    return num_dif



intree="/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe/RAxML_bestTree.concat_107_no_recom"
tree=Tree(intree, format=1)
distance=pd.DataFrame(columns=[rec for rec in tree.get_leaf_names() if 'plate' in rec], index=[rec for rec in tree.get_leaf_names() if 'plate' in rec])
pres_abs=pickle.load(open("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/geneCluster/dt_genePresence.cpk"))
pres_abs1={k: pres_abs[str(k+'.annotation')] for k in distance.keys()}
pres_abs_diff=pd.DataFrame(columns=pres_abs1.keys(), index=pres_abs1.keys())
i=0
for col in distance.columns:
    print i
    i=i+1
    for row in distance.index:
        if col!=row:
            distance[col][row]=tree.get_distance(col, row)
            pres_abs_diff[col][row]=compare_gene_content(pres_abs1[col], pres_abs1[row])
            print  pres_abs_diff[col][row]
        else:
            distance[col][row]=0
            pres_abs_diff[col][row]=0

pres_abs_diff.to_pickle('./pres_abs.cpk')
distance.to_pickle('./distance.cpk')

fig=plt.figure()
strain='plate11.F10'
x1=numpy.array(distance[strain])
y1=[ pres_abs_diff[strain].ix[rec] for rec in distance[strain].index]
plt.scatter(x1, y1, color='k')
plt.xlim(0, max(x1))
plt.ylim(0, max((y1)))
plt.xlabel('Distance on Tree', fontsize=18)
plt.ylabel('Number of differing orthogroups', fontsize=16)
fig.savefig("/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe/distance_vs_orthogroup.pdf")

