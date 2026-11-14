# Remove macOS / Word / editor metadata junk from the CHARM repo.
# Run from the CHARM folder:   source("clean_junk_files.R")
# It first LISTS what it would delete and asks you to type "yes" before doing anything.
# What counts as junk:  ._*  (macOS AppleDouble files)   .DS_Store   __pycache__/   ~$*  (Word lock files)   *.bak
# Never touches: .git/, .Rproj.user/, or any other file.

local({
root <- getwd()
stopifnot(file.exists("app.R"), dir.exists(".git"))           # safety: must be the CHARM repo folder

all_paths <- list.files(root, recursive = TRUE, all.files = TRUE, include.dirs = TRUE, full.names = FALSE)
all_paths <- all_paths[!grepl("(^|/)\\.git(/|$)|(^|/)\\.Rproj\\.user(/|$)", all_paths)]
base <- basename(all_paths)
is_dir <- dir.exists(file.path(root, all_paths))

junk_file <- !is_dir & (grepl("^\\._", base) | base == ".DS_Store" | grepl("^~\\$", base) | grepl("\\.bak$", base))
junk_dir  <-  is_dir & base == "__pycache__"
files <- all_paths[junk_file]
dirs  <- all_paths[junk_dir]

cat(sprintf("Found %d junk files and %d __pycache__ folders.\n", length(files), length(dirs)))
if (length(files) + length(dirs) == 0) { cat("Nothing to clean.\n"); return(invisible()) }
print(head(c(dirs, files), 60))
if (length(files) + length(dirs) > 60) cat("... (", length(files) + length(dirs) - 60, "more )\n")

ans <- if (isTRUE(get0("CLEAN_CONFIRM"))) "yes" else readline("Delete these? Type 'yes' to continue: ")
if (!identical(tolower(trimws(ans)), "yes")) { cat("Cancelled, nothing was deleted.\n"); return(invisible()) }

# 1. stop git from tracking any junk file that was committed earlier (keeps history, removes from the next commit)
tracked <- tryCatch(system2("git", c("ls-files"), stdout = TRUE), error = function(e) character(0))
tracked_junk <- tracked[grepl("(^|/)(\\._[^/]*|\\.DS_Store|~\\$[^/]*|[^/]*\\.bak)$", tracked) | grepl("(^|/)__pycache__/", tracked)]
if (length(tracked_junk) > 0) {
  cat(sprintf("%d junk files were tracked by git: removing them from the index.\n", length(tracked_junk)))
  for (chunk in split(tracked_junk, ceiling(seq_along(tracked_junk) / 100)))
    system2("git", c("rm", "--cached", "-q", "--ignore-unmatch", "--", shQuote(chunk)))
}

# 2. delete from disk
unlink(file.path(root, files))
unlink(file.path(root, dirs), recursive = TRUE)
cat(sprintf("Deleted %d files and %d folders.\n", length(files), length(dirs)))

# 3. make sure they cannot come back into git
gi <- if (file.exists(".gitignore")) readLines(".gitignore", warn = FALSE) else character(0)
want <- c("._*", ".DS_Store", "__pycache__/", "~$*", "*.bak")
add <- setdiff(want, trimws(gi))
if (length(add) > 0) {
  writeLines(c(gi, "", "# macOS / Word / backup junk", add), ".gitignore")
  cat("Added to .gitignore:", paste(add, collapse = "  "), "\n")
}
cat("Done. Review with 'git status' before committing.\n")

})
