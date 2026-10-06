# Two renv start-up costs that every crew worker (and tar_make()'s callr
# process) pays on each R start, measured on HTCF:
#  - synchronized.check: scans every R file in the project to compare the
#    library with renv.lock -- ~1 min on /scratch (BeeGFS). Run
#    renv::status() by hand instead.
#  - sandbox: re-creates a sandboxed copy of the system library under
#    ~/.cache (Ceph) -- 39 s per worker when 20 start at once.
options(
  renv.config.synchronized.check = FALSE,
  renv.config.sandbox.enabled = FALSE
)
source("renv/activate.R")
