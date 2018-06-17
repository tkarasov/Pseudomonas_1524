#!/usr/bin/py

#the goal of this script is to take in a user defined tree, and a number of nodes to be kept, and to write out a new tree with those nodes only,
#cd /Users/tkarasov/Dropbox/germany_pathogen_collections/data_files_rmarkdown

#**********To be acceptable to Megan, need to add a root manually to the map file and to the tre file.

from ete2 import Tree
tre="strain_tree.nwk"#sys.argv[1]
nod="or_strains_reduced.txt"#sys.argv[2]

#Load tree from Newick file
old_tree=Tree(tre, format=1)

#Now take nodes to subset
keep_nodes=[line.strip().split()[0] for line in open(nod).readlines()]
old_tree.prune(keep_nodes)

#I think it needs to be re-rooted now...
#old_tree.set_outgroup("p3.G10")

#Now traverse the tree renaming
mapping=[]
synonym=[]
i=1
for node in old_tree.traverse():
    if node.is_leaf():
        curr=node.name
        mapping.append([i, curr])
        synonym.append([curr, i])
        node.name=i
        print "Leaf  "+ str(i)
    if not node.is_leaf():
        node.name=i
        curr_min=node.get_leaves()
        curr='_'.join([str(rec.name) for rec in curr_min])
        mapping.append([i, curr])
        print "Node "+ str(i)
    i=i+1


#The mapping list of lists is the megan mapping file, and the tree is the tre file. The third file 
#mini.map
with open("mini.map", 'w') as file:
    file.writelines('\t'.join([str(val) for val in rec])+'\n' for rec in mapping)

#mini.tre
old_tree.write(outfile="mini.tre", format=8)

#synonyms2min.map
with open("synonyms2min.map", 'w') as file:
    file.writelines('\t'.join([str(val) for val in i])+'\n' for i in synonym)
