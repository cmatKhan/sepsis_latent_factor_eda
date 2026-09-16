# Image for the ingest slurm pipeline (R/create_ingest_slurm_bundle.R and
# the job families it stages: ingest_core, fgsea_grid, gprofiler_grid,
# projectr_within_grid/projectr_cross_grid -- see that script's header and
# config/ingest_slurm_config.yml). One combined image so a single
# `libPaths`/container entry can cover every stage instead of juggling the
# separate bioconductor-*/r-arrow_r-gprofiler2/r-projectr images referenced
# there today.
#
# Build (from the repo root, with BuildKit's extra --build-context so the
# local projectR checkout doesn't need to live inside this repo's own
# build context):
#
#   DOCKER_BUILDKIT=1 docker build \
#     --build-context projectr=/home/chase/projects/projectR \
#     -t chasem/ingest-slurm:latest -f Dockerfile .
#
# Convert to a Singularity/Apptainer .sif for the cluster (see
# R/README.md's "Packaging a bundle for the cluster" section):
#
#   apptainer build ingest-slurm.sif docker-daemon://chasem/ingest-slurm:latest

FROM rocker/tidyverse:4.6.0

# ---- system dependencies -----------------------------------------------
# libcurl/openssl/xml2: arrow + gprofiler2 (network calls to g:Profiler) +
# various Bioconductor package installs. libglpk: `clue`'s solve_LSAP.
RUN apt-get update && apt-get install -y --no-install-recommends \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev \
      libglpk-dev \
      zlib1g-dev \
      && rm -rf /var/lib/apt/lists/*

# ---- CRAN packages ------------------------------------------------------
# dplyr/tidyr already present via rocker/tidyverse.
RUN install2.r --error --skipinstalled \
      here \
      optparse \
      yaml \
      DBI \
      RSQLite \
      rslurm \
      arrow \
      matrixStats \
      clue \
      mclust \
      gprofiler2 \
      msigdbr \
      remotes \
      BiocManager

# ---- Bioconductor packages ------------------------------------------------
# CoGAPS transitively pulls in fgsea/DESeq2 (see
# config/ingest_slurm_config.yml's comment on the ingest_core/fgsea image),
# so listing it here covers both job families' needs.
RUN R -e 'BiocManager::install(c("fgsea", "BiocParallel", "DESeq2", "CoGAPS"), update = FALSE, ask = FALSE)'

# ---- projectR, from the user's local fork --------------------------------
# The released Bioconductor projectR is NOT sufficient here -- this project
# depends on unreleased changes in the local checkout at
# /home/chase/projects/projectR (see R/README.md's container notes and
# config/ingest_slurm_config.yml's `projectr:` entry). Requires the
# `projectr` build context above; install_local() also picks up projectR's
# own DESCRIPTION-declared dependencies.
COPY --from=projectr . /tmp/projectR
RUN R -e 'remotes::install_local("/tmp/projectR", dependencies = TRUE, upgrade = "never")' \
    && rm -rf /tmp/projectR

# ---- sanity check ---------------------------------------------------------
# Fails the build immediately if any package doesn't import cleanly, rather
# than surfacing as a mysterious "missing" package deep into a slurm job
# (see config/ingest_slurm_config.yml's note about a leaked host
# R_LIBS_USER previously causing this same failure mode).
RUN R -e 'library(here); library(optparse); library(yaml); library(DBI); library(RSQLite); \
          library(rslurm); library(arrow); library(matrixStats); library(clue); library(mclust); \
          library(gprofiler2); library(msigdbr); library(fgsea); library(BiocParallel); \
          library(DESeq2); library(CoGAPS); library(projectR); \
          message("all packages imported cleanly")'
