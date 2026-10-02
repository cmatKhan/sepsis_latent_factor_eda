#!/bin/bash
#SBATCH --job-name=wgcna_diagnostics
#SBATCH --output=logs/wgcna_diagnostics-%j.out
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --container=/ref/mblab/containers/chasem/r-wgcna_1.74--aa9d6f2e8e081486.sif

set -euo pipefail

# --recompute-dendro is resumable (fits that already have a dendro are
# skipped -- see compute_wgcna_diagnostics.R's header), so if this job
# hits the time limit partway through, just resubmit the same script.
Rscript R/compute_wgcna_diagnostics.R --datasets datasets.txt --db results/stability.sqlite --recompute-dendro
