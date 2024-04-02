#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --job-name=nnSVG
#SBATCH --output=nnSVG_FFPE_CA.out
#SBATCH --mail-user=du.m@wehi.edu.au
#SBATCH --mail-type=END,FAIL

module load R/4.2.0
module load ImageMagick
Rscript nnSVG.R
