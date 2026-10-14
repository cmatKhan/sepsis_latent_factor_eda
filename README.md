# Sepsis latent-factor stability

[![Docs](https://img.shields.io/badge/docs-gh--pages-blue)](https://cmatkhan.github.io/sepsis_latent_factor_eda/)

This provides a [targets](https://books.ropensci.org/targets/) workflow 
configured with [crew.cluster](https://wlandau.github.io/crew.cluster/) to 
conduct variancePartition analysis, along with a number of latent factor
decomposition methods, to examine repeat measure Sepsis datasets. It may 
be extended with both new data, and new analysis methods. Results are 
ingested into a database and visualized with a Shiny app.

## Installation

For more in depth setup instructions, especially for running the targets
workflow on a SLURM cluster, see the
[setup docs](https://cmatkhan.github.io/sepsis_latent_factor_eda/setup.html).  

The following are general setup instructions, especially suitable for running 
the shiny app from previously generated data or local development.

This is a R project, not a package, with an environment managed by `renv`. It 
is expected that you'll clone this repo and then work inside the project 
as opposed to installing it as a package and calling it into some R session.

1. clone the repo to your local

```
git clone https://github.com/cmatKhan/sepsis_latent_factor_eda.git
```

2. Next, open the R project. There are multiple ways of doing this, either 
from in Rstudio or by navigating to the .Rproj file in the project 
directory and clicking it.

3. Install the environment. **You need R 4.6.1** (the version recorded in
`renv.lock`) and a toolchain that can build R packages from source (a
C/C++/Fortran compiler, `make`, and the usual system libraries such as libcurl,
openssl and libxml2). Then start R **in the project directory** (opening the
`.Rproj` does this). The first start downloads and installs renv itself and
prints a message that the project library is empty; that is expected. In the 
console, run:

```r
options(renv.install.timeout = 6 * 3600)   # many packages build from source; the default is 1 hour
renv::restore()                            # answer y when asked to proceed
```

`renv::restore()` installs every package at the version in `renv.lock` into the
project's own library (`renv/library/`), leaving your other R libraries alone.
Expect it to take a long time on a fresh machine (hours, since roughly 300
packages build from source) and minutes where renv's cache already holds the
builds. If it stops, run it again: finished builds are cached and it continues
from there. When it finishes, check the result with `renv::status()`, which
should report that the project is in a consistent state.
Two packages need extra steps (CoGAPS is rebuilt with OpenMP; projectR comes
from a fork branch that `restore()` fetches from GitHub). See
[docs/setup.qmd](docs/setup.qmd). The Prerequisites and Set up sections below
are the SLURM cluster (HTCF) route and are not needed to browse results or to
work on a laptop.

## Launching the app with public data

To launch the app, you need the results of the targets pipeline. A version 
of this is available [here](https://huggingface.co/datasets/cmatkhan/sepsis_latent_factors_db_data).
Download that data and either mount it to a directory called `results` in your 
project, or extract it. If you're using linux or mac, you probably already 
have [squash](https://github.com/plougher/squashfs-tools).
If not, you'll need to install it.

```
# first, make a directory called results in the project
mkdir results
```
Then, mount the results.squash file to it:

```
# if you have sudo privileges
sudo mount -t squashfs -o loop results_no_rose.squash results 

# if you do not have sudo privileges, install squashfuse and then do
squashfuse results_no_rose.squash results 
```

When finished, stop the app and unmount

```
# if you have sudo privileges
sudo umount results

# if you used squashfuse
fusermount -u results
```

then, either open `app/app.R` and click "Run App" in rstudio, or from the 
console enter `shiny::runApp("app")`
