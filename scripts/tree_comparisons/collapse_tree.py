#!/usr/bin/py
from ete2 import Tree
import sys
from ete2 import Tree, TreeStyle
from Bio import SeqIO
import os
import numpy as np

circular_style = TreeStyle()
circular_style.mode = "c" # draw tree in circular mode
circular_style.scale = 20
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = False
ts.show_branch_support = False

#script stole from https://www.biostars.org/p/97409/

def mean(array):
    return sum(array)/float(len(array))


def cache_distances(tree):
    ''' precalculate distances of all nodes to the root''' 
    node2rootdist = {tree:0}
    for node in tree.iter_descendants('preorder'):
        node2rootdist[node] = node.dist + node2rootdist[node.up]
    return node2rootdist


def collapse(tree, min_dist):
    # cache the tip content of each node to reduce the number of times the tree is traversed
    #Because the collapsing is happening from the inner nodes, we could iterate over the collapsing several times, but we will not.
    num_tips=len(tree.get_leaves())
    node2tips = tree.get_cached_content()
    root_distance = cache_distances(tree)
    #re-run the collapsing until the number of nodes settles (until two iterations previously the same value)
#    while num_tips_old>num_tips_iteration_new:
    num_tips_old=num_tips
    num_tips=len(tree.get_leaves())
    collapse_for_loop(node2tips, root_distance, tree, min_dist)
    num_tips_iteration_new=len(tree.get_leaves())
    print "Was "+ str(num_tips)+ ", is now "+str(num_tips_iteration_new)+ " tips"

             
def collapse_for_loop(node2tips, root_distance, tree, min_dist):
    for node in tree.get_descendants('preorder'):
        if not node.is_leaf():
            max_distance_to_tips = max([root_distance[tip]-root_distance[node]
                                         for tip in node2tips[node]])
            if max_distance_to_tips < min_dist:
                    #rename first tip, average distance of all other tips, remove rest of tips
                    new_tip=  '_'.join([tip.name for tip in node2tips[node]])
                    new_tip_dist=np.mean([root_distance[tip]-root_distance[node]
                                         for tip in node2tips[node]])
                    node.add_features(collapsed=True)
                    #support_vec=np.mean([node.support for support in node.get_children()])
                    for ch in node.get_children():
                        ch.detach()

                    node.add_child(name=new_tip, dist=new_tip_dist)
                    #node.is_leaf()=True
                    #node.img_style['draw_descendants'] = Fals

    return tree


def mean(array):
    return sum(array)/float(len(array))


def cache_distances(tree):
    ''' precalculate distances of all nodes to the root''' 
    node2rootdist = {tree:0}
    for node in tree.iter_descendants('preorder'):
        node2rootdist[node] = node.dist + node2rootdist[node.up]
    return node2rootdist


def collapse_rename(tree, min_dist):
    #does the same thing as collapse but renames nodes with species identification
    # cache the tip content of each node to reduce the number of times the tree is traversed
    node2tips = tree.get_cached_content()
    root_distance = cache_distances(tree)

    for node in tree.get_descendants('preorder'):
        if not node.is_leaf():
            avg_distance_to_tips = mean([root_distance[tip]-root_distance[node]
                                         for tip in node2tips[node]])

            if avg_distance_to_tips < min_dist:

                node.name += ' COLLAPSED avg_d:%g {%s}' %(avg_distance_to_tips,
                                                 ','.join([str(tip.name + '_'+ species_dict[tip.name]) for tip in node2tips[node]]))
                node.add_features(collapsed=True)

                node.img_style['draw_descendants'] = False


def rename_singular(tree, min_dist):
    #does the same thing as collapse but renames nodes with single identification
    # cache the tip content of each node to reduce the number of times the tree is traversed
    node2tips = tree.get_cached_content()
    root_distance = cache_distances(tree)

    for node in tree.get_descendants('preorder'):
        if not node.is_leaf():
            avg_distance_to_tips = mean([root_distance[tip]-root_distance[node]
                                         for tip in node2tips[node]])

            if avg_distance_to_tips < min_dist:

                node.name = [(tip.name) for tip in node2tips[node]][0]
                node.add_features(collapsed=True)

                node.img_style['draw_descendants'] = False


def calc_core_gene_alignment_length(core_gene_folder, SNP_alignment):
    #calculates the length of SNP alignment
    i=0
    fasta_sequences = SeqIO.parse(open(SNP_alignment),'fasta')
    for fasta in fasta_sequences:
        i=i+1
        if i>1:
            pass
        else:
            SNP_length=len(fasta.seq)
             
    #calculates of all core gene alignment
    tot=[]
    for input_file in os.listdir(core_gene_folder):
            i=0
            if input_file.endswith("fa"):
                fasta_sequences = SeqIO.parse(open(core_gene_folder+input_file),'fasta')
                for fasta in fasta_sequences:
                    i=i+1
                    if i>1:
                            pass
                    else:
                            gene_length=len(fasta.seq)
                            tot.append(gene_length)
    return float(SNP_length)/sum(tot)
                    
    
def pull_single_node(tree):
#when I need to simplify the tree to pull 16s sequences for just one sequence per node. Note that I exclude the genomes not sequenced here.
    tree=t #t is the collapsed tree
    is_mine=[]
   
    for node in tree.get_leaves():
        strains=[node.name.split( ) for node in tree.get_leaves()]
    keep=[]
    count_node=[]
    for rec in strains:
        is_mine=[]
        if 'NODE' in rec[0]:
            plates= rec[3].split("-")
            is_mine=[thing for thing in plates if 'plate' in thing]
            #print rec[0]
            if len(is_mine)==0:
                print "This doesn't work"+ str(is_mine)
                #pass
                #val=rec[3].split("-")[0].strip("{")
                #keep.append(val)
            else:
                print is_mine[0]
                val=is_mine[0].strip("{")
                keep.append([val, len(is_mine)])
        else:
            if "plate" in rec[0]:
                is_mine=["here"]
            keep.append([rec[0], len(is_mine)])

    keep_mine=[rec for rec in keep if rec[1]!=0]

    new_tree=tree.copy()
    for node in new_tree.get_leaves():
        is_mine=[]
        if 'NODE' in node.name:
            plates= node.name.split("-")
            is_mine=[thing for thing in plates if 'plate' in thing]
            #print rec[0]
            if len(is_mine)==0:
                print "This doesn't work"+ str(is_mine)
                #pass
                #val=rec[3].split("-")[0].strip("{")
                #keep.append(val)
            else:
                print is_mine[0]
                val=is_mine[0].split("{")[1].strip("{")
                node.name=val
        else:
            if "plate" in node.name:
                is_mine=["here"]
            node.name=(node.name)

            
    my_node=new_tree.get_common_ancestor(["plate23.A5.annotation","NZ_AP014637","plate2.C11.annotation"])
    Psyringae_leaves=[rec.name for rec in my_node.get_leaves() if 'plate' in rec.name]
    P_other_leaves=[rec.name for rec in new_tree.get_leaves() if rec.name not in Psyringae_leaves]
    for rec in keep_mine:
        if rec[0] in Psyringae_leaves:
            rec.append("P_syringae")
        if rec[0] in P_other_leaves:
            rec.append("P_other")
        rec[0]=rec[0].strip(".annotation")
    outfile=open("strain_per_node_count.txt", "w")
    for item in keep_mine:
        print>> outfile, '\t'.join([str(rec) for rec in item])
    outfile.close()

def output_strain_assignment(t):
    is_mine=[]
    orig_tree=Tree(intree, format=1)
    for node in orig_tree.get_leaves():
        strains=[node.name.split( ) for node in tree.get_leaves()]
    keep=[]
    count_node=[]
    for rec in strains:
        is_mine=[]
        if 'NODE' in rec[0]:
            plates= rec[3].split("-")
            is_mine=[thing.strip("{").strip("}") for thing in plates if 'plate' in thing]
            #print rec[0]
            if len(is_mine)==0:
                print "This doesn't work"+ str(is_mine)
                #pass
                #val=rec[3].split("-")[0].strip("{")
                #keep.append(val)
            else:
                print is_mine[0]
                val=is_mine[0].strip("{")
                keep.append([val, is_mine])
        else:
            if "plate" in rec[0]:
                is_mine=["here"]
            keep.append([rec[0], len(is_mine)])

    keep_mine=[rec for rec in keep if rec[1]!=0]

    #Now assign create table in which every isolate is first column, strain identification is second, pathogenicity identification is third column
    my_node=orig_tree.get_common_ancestor(["plate23.A5.annotation","NZ_AP014637","plate2.C11.annotation"])
    Psyringae_leaves=[rec.name for rec in my_node.get_leaves() if 'plate' in rec.name]
    P_other_leaves=[rec.name for rec in orig_tree.get_leaves() if rec.name not in Psyringae_leaves]
    for row in keep_mine:
        if type(row[1])==int:
            row[1]=row[0]

    strain_assignment=[]
    for row in keep_mine:
        if type(row[1])==str:
            strain=row[0]
            isolate=row[1]
            strain_assignment.append([isolate, strain, complexN])
        else:
            for rec in row[1]:
                strain=row[0]
                isolate=rec
                if rec in Psyringae_leaves:
                    complexN="P_syringae"
                if rec in P_other_leaves:
                    complexN="P_other"
                strain_assignment.append([isolate, strain, complexN])

    outfile=open("strain_assignment_node_species.txt", "w")
    for item in strain_assignment:
        print>> outfile, '\t'.join([str(rec) for rec in item])
    outfile.close()
    
#print t.get_ascii(attributes=["dist", 'name'])

'''if __name__ == "__main__":
    intree="/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/all16s_verified_25_7_2017_1524.newick"
    core_gene_folder="/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/coregenes_alignments/"
    SNP_alignment="/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/geneCluster/SNP_whole_matrix.aln"

    
#sys.argv[1]#"tree_result_no_p5_G10_p4_G7.newick"#sys.argv[1] #file path and name
    dist=0.0001#float(sys.argv[2]) #desired distance for collapsing
   # output=intree+str(adj_dist)+".nwk"#"tree_result_no_p5_G10_p4_G7_collapse_0.001.newick" #sys.argv[3] #output file name
    percent_SNP=calc_core_gene_alignment_length(core_gene_folder, SNP_alignment)
    adj_dist=0.0001/percent_SNP
    output=intree+str(adj_dist)+".nwk"#"tree_result_no_p5_G10_p4_G7_collapse_0.001.newick" #sys.argv[3] #output file name
    print adj_dist


    t=Tree(intree, format=1)
 #   t.render(output+"original.pdf", tree_style=ts)#,tree_style=circular_style)
    collapse(t, adj_dist) 
    #print t.get_ascii(attributes=['name'])
    for n in t.search_nodes(collapsed=True):
        for ch in n.get_children():
             ch.detach()
             
    #Now just write to file
    t.write(format=1, outfile=output)
    t.render(output+".pdf", tree_style=ts)
    
#rename nodes --not needed unless putting species on
#first read in species identification tree
#spec="/Users/tkarasov/Dropbox/germany_pathogen_collections/data_files_rmarkdown/run_0034_species_assignment.txt"
#species_ident=[line.strip().split(',') for line in open(spec).readlines()]
#species_dict={}
#for line in species_ident:
#    species_dict[line[0]]=line[1]


#collapse_rename(t, min_dist)
#t.show()

'''
