#!/bin/bash   
#SBATCH -p normal       # partition (queue)
#SBATCH -N 1            # number of nodes
#SBATCH -n 1            # number of tasks
#SBATCH --mem 10000       # memory pool for all cores
#SBATCH -t 0-2:00          # time (D-HH:MM)
#SBATCH -o logs/slurm.%N.%j.out     # STDOUT
#SBATCH -e logs/slurm.%N.%j.err     # STDERR
#SBATCH --job-name enrichGO
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=Perla.TroncosoRey@earlham.ac.uk # send-to address

source R-3.2.2_GSEA
Rscript enrichmentGO.R MultiPEN-Rankings_lambda0.1.txt
