#!/bin/sh
#
#  Reserve 20 CPUs for this job
#$ -pe parallel 20
#  Request 640G of RAM
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



cd /ebio/abt6_projects9/Pseudomonas_diversity/data/post_assembly_analysis/pan_genome/data/vis/clonalframe
/usr/bin/raxmlHPC-PTHREADS-SSE3  -T 20 -p 654329  -s concat_107_no_recom.fasta  -m GTRGAMMA -n concat_107_no_recom
