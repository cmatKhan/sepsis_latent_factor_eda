# Connection to the stability DB (R/db/schema.sql + views.sql). The DB is a
# derived output: open_stability_db() creates a new file from the schema and
# refuses one built from a different schema version -- delete it and rebuild
# with tar_make() (docs: Database). Artifacts live under
# stability_artifacts/<dataset_id>/ next to the DB, referenced by paths
# relative to the DB's directory.

#' Bump when schema.sql or views.sql changes.
SCHEMA_VERSION <- 1L

#' Split a schema file into statements
#'
#' Statements are separated by lines holding only `;;`, since a statement may
#' itself contain semicolons.
#'
#' @param path Path to a `.sql` file.
#' @return Character vector of non-empty statements.
read_sql_statements <- function(path) {
  sql <- paste(readLines(path), collapse = "\n")
  parts <- trimws(strsplit(sql, "\n;;\n", fixed = TRUE)[[1]])
  parts[nzchar(gsub("--[^\n]*", "", parts) |> trimws())]
}

#' Open the stability DB
#'
#' Creates the file from `schema.sql` and `views.sql` if it's new and stamps
#' `PRAGMA user_version`; otherwise refuses a DB whose version isn't
#' `SCHEMA_VERSION`. Sets `foreign_keys = ON`, WAL journaling,
#' `synchronous = NORMAL` and a 30 s busy timeout.
#'
#' @param db_path Path to the SQLite file.
#' @param schema_dir Directory holding `schema.sql` and `views.sql`.
#' @return A DBI connection.
open_stability_db <- function(db_path, schema_dir = "R/db") {
  dir.create(dirname(db_path), recursive = TRUE, showWarnings = FALSE)
  is_new <- !file.exists(db_path) || file.size(db_path) == 0
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  DBI::dbExecute(con, "PRAGMA foreign_keys = ON")
  # WAL + synchronous = NORMAL: no full fsync on every commit (costly on
  # network filesystems). If db_path's filesystem can't do WAL's shared
  # memory locking, SQLite fails loudly here.
  DBI::dbExecute(con, "PRAGMA journal_mode = WAL")
  DBI::dbExecute(con, "PRAGMA synchronous = NORMAL")
  DBI::dbExecute(con, "PRAGMA busy_timeout = 30000")

  if (is_new) {
    DBI::dbExecute(con, "BEGIN")
    for (f in c("schema.sql", "views.sql")) {
      for (stmt in read_sql_statements(file.path(schema_dir, f))) DBI::dbExecute(con, stmt)
    }
    DBI::dbExecute(con, sprintf("PRAGMA user_version = %d", SCHEMA_VERSION))
    DBI::dbExecute(con, "COMMIT")
  } else {
    v <- DBI::dbGetQuery(con, "PRAGMA user_version")$user_version
    if (!identical(as.integer(v), SCHEMA_VERSION)) {
      DBI::dbDisconnect(con)
      stop("Stability DB ", db_path, " has schema version ", v, ", but this code writes version ",
           SCHEMA_VERSION, ". Delete it (and its stability_artifacts/ directory) and rebuild with ",
           "targets::tar_make(); there are no migrations.", call. = FALSE)
    }
  }
  con
}

#' A dataset's artifact directory
#'
#' @param db_path Path to the SQLite file.
#' @param dataset_id Dataset id.
#' @return `<db dir>/stability_artifacts/<dataset_id>`.
artifacts_dir <- function(db_path, dataset_id) {
  file.path(dirname(db_path), "stability_artifacts", dataset_id)
}

#' An artifact's path as the DB stores it
#'
#' @param dataset_id Dataset id.
#' @param fname File name.
#' @return `stability_artifacts/<dataset_id>/<fname>`, relative to the DB's directory.
artifact_rel <- function(dataset_id, fname) file.path("stability_artifacts", dataset_id, fname)

#' Resolve stored artifact paths
#'
#' @param path Paths as stored in the DB (relative to the DB's directory);
#'   absolute paths pass through.
#' @param db_path Path to the SQLite file.
#' @return Absolute paths; `NA` for `NA` or empty input.
resolve_artifact <- function(path, db_path) {
  out <- rep(NA_character_, length(path))
  ok <- !is.na(path) & nzchar(path)
  is_abs <- ok & startsWith(path, "/")
  out[is_abs] <- path[is_abs]
  rel <- ok & !is_abs
  out[rel] <- file.path(normalizePath(dirname(db_path)), path[rel])
  out
}
