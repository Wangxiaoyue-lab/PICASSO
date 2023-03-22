#!/bin/bash
#SBATCH -J test
#SBATCH -p normal
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --mem 200g

~/miniconda3/envs/r1/bin/Rscript /public/home/luoliheng/SINGLE/Sample_pipline/qc.R
