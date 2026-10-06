# Builds docs/reference.qmd -- the "Function reference" page -- from the
# roxygen comments in the project's R code. Run by Quarto before every
# render (`pre-render` in docs/_quarto.yml), from the docs/ directory.
#
# This is a project, not a package, so there is no roxygenise() / man/;
# roxygen2::parse_file() reads the blocks without evaluating any code.
# Functions with no roxygen block are listed as undocumented.

# Quarto runs this from docs/, where the project's .Rprofile (and so renv)
# isn't loaded: activate the project library first.
if (!requireNamespace("roxygen2", quietly = TRUE)) {
  owd <- setwd("..")
  source("renv/activate.R")
  setwd(owd)
}

root <- normalizePath("..")
sections <- c(
  "R/lib" = "Configuration, data paths, matrices and method discovery",
  "R/lib/ingest" = "Ingestion computations",
  "R/db" = "Database connection",
  "R/targets" = "Target factories, controllers, ingestion and DB writers",
  "R/ingest_jobs" = "Enrichment and projectR jobs",
  "R/methods" = "Methods",
  "R/de" = "Differential expression"
)

md_escape <- function(x) gsub("\n", " ", trimws(x))

tag_value <- function(block, tag) {
  v <- roxygen2::block_get_tag_value(block, tag)
  if (is.null(v)) NA_character_ else paste(v, collapse = "\n\n")
}

# Name of the object a block documents: `name <- ...` / `name = ...`.
block_name <- function(block) {
  call <- block$call
  if (is.call(call) && length(call) >= 2 && as.character(call[[1]]) %in% c("<-", "=")) {
    return(deparse(call[[2]]))
  }
  NA_character_
}

# Top-level function definitions in a file, by name.
defined_functions <- function(file) {
  exprs <- parse(file, keep.source = FALSE)
  names <- vapply(exprs, function(e) {
    if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") && is.call(e[[3]]) &&
        identical(e[[3]][[1]], as.name("function"))) deparse(e[[2]]) else NA_character_
  }, character(1))
  names[!is.na(names)]
}

document_block <- function(block, name) {
  out <- c(sprintf("#### `%s()` {.unnumbered}", name), "")
  title <- tag_value(block, "title")
  if (!is.na(title)) out <- c(out, paste0("**", md_escape(title), "**"), "")
  desc <- tag_value(block, "description")
  if (!is.na(desc)) out <- c(out, desc, "")
  params <- roxygen2::block_get_tags(block, "param")
  if (length(params)) {
    out <- c(out, "| Argument | Description |", "|---|---|",
             vapply(params, function(p) sprintf("| `%s` | %s |", p$val$name,
                                                gsub("\\|", "\\\\|", md_escape(p$val$description))),
                    character(1)), "")
  }
  ret <- tag_value(block, "return")
  if (!is.na(ret)) out <- c(out, paste0("**Returns:** ", md_escape(ret)), "")
  out
}

lines <- c("---", 'title: "Function reference"', "---", "",
           "Generated from the roxygen comments in the code by `docs/build_reference.R`;",
           "edit the comments, not this page.", "")
n_documented <- 0L
undocumented <- character(0)
for (dir in names(sections)) {
  files <- sort(list.files(file.path(root, dir), pattern = "\\.R$", full.names = TRUE))
  if (length(files) == 0) next
  lines <- c(lines, sprintf("## %s (`%s/`)", sections[[dir]], dir), "")
  for (f in files) {
    rel <- sub(paste0(root, "/"), "", f, fixed = TRUE)
    blocks <- roxygen2::parse_file(f, env = NULL)
    names <- vapply(blocks, block_name, character(1))
    fns <- defined_functions(f)
    keep <- !is.na(names) & names %in% fns
    if (!any(keep) && length(fns) == 0) next
    lines <- c(lines, sprintf("### `%s`", rel), "")
    for (i in which(keep)) {
      lines <- c(lines, document_block(blocks[[i]], names[[i]]))
      n_documented <- n_documented + 1L
    }
    missing <- setdiff(fns, names[keep])
    if (length(missing)) {
      undocumented <- c(undocumented, paste0("`", rel, "::", missing, "()`"))
      lines <- c(lines, paste0("*Undocumented:* ", paste0("`", missing, "()`", collapse = ", ")), "")
    }
  }
}
writeLines(lines, file.path(root, "docs", "reference.qmd"))
message(sprintf("reference.qmd: %d documented functions, %d undocumented",
                n_documented, length(undocumented)))
if (length(undocumented)) message("undocumented: ", paste(undocumented, collapse = ", "))
