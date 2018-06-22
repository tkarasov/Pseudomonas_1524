#!/bin/sh
#
#  Reserve 20 CPUs for this job
#$ -pe parallel 20
#  Request 640G of RAM
#$ -l h_vmem=32G
#  
#$ -o $HOME/tmp/stdout_of_job
#
#$ -j y
#
#  Use /bin/bash to execute this script
#$ -S /bin/bash
#
#  Run job from current working directory
#$ -cwd

cd /ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe
/usr/bin/raxmlHPC-PTHREADS-SSE3  -T 20 -p 654329  -s concat_107_no_recom.fasta  -m GTRGAMMA -n concat_107_no_recom
