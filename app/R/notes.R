# Missing data. Anything the app expects from the database or the artifact files that
# isn't there -- a method that hasn't been run yet, a table from an older schema, an
# artifact file that wasn't copied -- is logged once to the console and shown as a notice.
# It never stops the app.

.notes_seen <- new.env(parent = emptyenv())

app_log <- function(level, ...) {
  message(sprintf("[stability-app %s] %s: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, paste0(...)))
}

#' Log that something expected is missing (once per message and detail).
note_missing <- function(what, detail = NULL) {
  key <- paste(what, detail)
  if (is.null(.notes_seen[[key]])) {
    .notes_seen[[key]] <- TRUE
    app_log("missing", what, if (!is.null(detail)) paste0(" (", detail, ")"))
  }
  invisible(NULL)
}

#' An inline notice for a card whose data is missing.
missing_notice <- function(what, detail = NULL) {
  div(class = "alert alert-warning py-2 my-2", role = "alert",
      strong("Not available: "), what, if (!is.null(detail)) span(" ", detail))
}

#' A query that returns NULL (and logs) instead of failing, e.g. when the table is absent.
safe_query <- function(con, sql, params = NULL, what = sql) {
  tryCatch(DBI::dbGetQuery(con, sql, params = params),
           error = function(e) { note_missing(what, conditionMessage(e)); NULL })
}

#' Inside a render function: stop with a message (and log) when `x` is NULL or has no rows.
need_data <- function(x, what, hint = NULL) {
  if (is.null(x) || (is.data.frame(x) && nrow(x) == 0) || length(x) == 0) {
    note_missing(what, hint)
    shiny::validate(shiny::need(FALSE, paste0("Not available: ", what, if (!is.null(hint)) paste0(" (", hint, ")"))))
  }
  invisible(x)
}

#' Open the stability DB read-only. A database in WAL mode on a read-only file system
#' (e.g. a mounted snapshot of the cluster results) can't be opened normally, because
#' SQLite wants to create its -shm file; it is then opened as an immutable snapshot, which
#' means the file must not change while the app runs.
connect_read_only <- function(path) {
  path <- normalizePath(path)
  con <- tryCatch(DBI::dbConnect(RSQLite::SQLite(), path, flags = RSQLite::SQLITE_RO, synchronous = NULL),
                  error = function(e) NULL)
  works <- !is.null(con) && tryCatch({ DBI::dbGetQuery(con, "SELECT 1 FROM sqlite_master LIMIT 1"); TRUE },
                                     error = function(e) FALSE)
  if (works) return(con)
  if (!is.null(con)) try(DBI::dbDisconnect(con), silent = TRUE)
  app_log("note", "could not open ", path, " read-only in place; opening it as an immutable snapshot ",
          "(restart the app after replacing the file)")
  DBI::dbConnect(RSQLite::SQLite(), paste0("file:", path, "?immutable=1"),
                 flags = RSQLite::SQLITE_RO, synchronous = NULL)
}
