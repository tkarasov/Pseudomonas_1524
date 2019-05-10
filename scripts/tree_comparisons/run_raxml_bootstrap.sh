#!/bin/sh
#
#  Reserve 20 CPUs for this job
#$ -pe parallel 20
#  
#$ -l h_vmem=32G
#
#  The name shown in the qstat output and in the output file(s). The
#  default is to use the script name.
#$ -N $letter$number
#
#  The path used for the standard output stream of the job
#$ -o $HOME/tmp/stdout_of_job
#
# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
#$ -j y
#
#  Use /bin/bash to execute this script
#$ -S /bin/bash
#
#  Run job from current working directory
#$ -cwd
#
#  Send email when the job begins, ends, aborts, or is suspended
#$ -m beas

cd /ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/phylogeny
alignment=/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/geneCluster/SNP_whole_matrix.aln

#only bootstrap+ML tree estimate, 100 bootstraps. This is very slow 10-20minutes per bootstrap
/usr/bin/raxmlHPC-PTHREADS-SSE3  -T 20 -p 12345 -x 12345 -# 100  -s $alignment -n branches -c 25 -m GTRGAMMA -p 344312987 -f a




#­f b draw bipartition information on a tree provided with ­t (typically the bestknownML tree) based on multiple trees (e.g., from a bootstrap) in a file  specified by ­z
# /usr/bin/raxmlHPC-PTHREADS-SSE3  -T 20 -p 12345 -x 12345 -# 100  -s $alignment -n branches -c 25 -m GTRGAMMA -p 344312987 -f a

bootstrap=/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/phylogeny/RAxML_bootstrap.branches
original=/ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/strain_tree.nwk


#now we can analyze the quality of the partition from the boostraps above
/usr/bin/raxmlHPC -f b -t $original -z $boostrap -m GTRGAMMA -n test
#/usr/bin/raxmlHPC-PTHREADS-SSE3  -m GTRGAMMA  -n test  -z /ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/phylogeny/RAxML_bootstrap.branches  -f b -t /ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/strain_tree.nwk
/usr/bin/raxmlHPC-PTHREADS-SSE3  -m GTRGAMMA  -n raxml_tree_bootstrapped.tree  -z /ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/phylogeny/RAxML_bootstrap.branches  -f b -t /ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/phylogeny/RAxML_bestTree.branches
