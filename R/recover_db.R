# Pure-R salvage for a `database disk image is malformed` stability DB --
# no sqlite3 CLI required (its `.recover` command is more thorough if it's
# available on your cluster -- see R/README.md's "if the DB is corrupted"
# section -- but this covers the common case where corruption is confined
# to a handful of tables/pages, not the master schema itself).
#
# Strategy: bootstrap a brand-new, empty DB with a clean schema
# (open_stability_db() on a path that doesn't exist yet), then for every
# table the CORRUPTED db still knows about, try a plain `SELECT *` and copy
# whatever succeeds into the new DB. A table that errors is reported and
# skipped (its rows are lost) rather than aborting the whole recovery.
#
# Usage:
#   Rscript R/recover_db.R <corrupted_db_path> <new_db_path>
# `stability_artifacts/` is NOT touched/copied here -- artifacts live
# next to the DB file already; once you're satisfied with new_db_path,
# move it into place yourself (see printed instructions at the end) so
# resolve_artifact()'s "artifacts live next to the DB" assumption holds.

library(here)
source(here("R/lib/ingest/db.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript R/recover_db.R <corrupted_db_path> <new_db_path>")
old_path <- args[[1]]
new_path <- args[[2]]
if (file.exists(new_path)) stop(new_path, " already exists -- refusing to overwrite. Pick a new path.")

message("Opening corrupted DB (read side): ", old_path)
old_con <- DBI::dbConnect(RSQLite::SQLite(), old_path)

message("Integrity check (informational -- recovery proceeds regardless):")
ic <- tryCatch(DBI::dbGetQuery(old_con, "PRAGMA integrity_check"),
               error = function(e) data.frame(integrity_check = paste("check itself failed:", conditionMessage(e))))
print(ic)

tables <- tryCatch(DBI::dbListTables(old_con), error = function(e) {
  stop("Could not even list tables in the corrupted DB (", conditionMessage(e), ") -- ",
       "this needs the sqlite3 CLI's `.recover` (or a lower-level tool); pure-R recovery can't help here.")
})
message("\nFound ", length(tables), " table(s): ", paste(tables, collapse = ", "))

message("\nBootstrapping a fresh, empty DB with a clean schema at: ", new_path)
new_con <- open_stability_db(new_path)

recovered <- character(0); failed <- character(0)
for (tbl in tables) {
  res <- tryCatch({
    d <- DBI::dbGetQuery(old_con, sprintf("SELECT * FROM %s", tbl))
    if (nrow(d) > 0) DBI::dbWriteTable(new_con, tbl, d, append = TRUE, row.names = FALSE)
    nrow(d)
  }, error = function(e) e)
  if (inherits(res, "error")) {
    message("  [", tbl, "] FAILED to recover: ", conditionMessage(res))
    failed <- c(failed, tbl)
  } else {
    message("  [", tbl, "] recovered ", res, " row(s)")
    recovered <- c(recovered, tbl)
  }
}

DBI::dbDisconnect(old_con)
DBI::dbDisconnect(new_con)

message("\n===== summary =====")
message("Recovered cleanly: ", paste(recovered, collapse = ", "))
if (length(failed) > 0) {
  message("FAILED (rows lost, table left empty in new DB): ", paste(failed, collapse = ", "))
  message("If any of these are `fits`/`factors`/`ingests`, you'll likely want to purge and re-ingest")
  message("whichever dataset(s) they belonged to -- see R/purge_dataset.R -- rather than trust a")
  message("partially-recovered fits table.")
}
message("\nNext: verify integrity, then swap the new DB into place, e.g.:")
message("  mv ", old_path, " ", old_path, ".corrupt-backup")
message("  mv ", new_path, " ", old_path)
