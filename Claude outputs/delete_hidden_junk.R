# Delete macOS/FUSE junk files from the CHARM repository
#
# Targets only:
#   ._*            AppleDouble metadata files (created when a Mac/FUSE mount
#                  writes to a non-Mac drive)
#   .fuse_hidden*  leftovers from files deleted while still open on a FUSE mount
#
# Everything else that starts with a dot (.git, .gitignore, .Rproj.user,
# .Rhistory, .RData, .github, .Rbuildignore, .nojekyll, .dockerignore) is kept.
#
# Run once with dry_run <- TRUE to see the list, then set it to FALSE.

root    <- "C:/Users/alexandre.afonso/Documents/GitHub/CHARM"
dry_run <- TRUE

all_files <- list.files(
  root,
  recursive    = TRUE,
  all.files    = TRUE,   # include dot files
  full.names   = TRUE,
  include.dirs = FALSE,
  no..         = TRUE
)

junk <- all_files[grepl("^(\\._|\\.fuse_hidden)", basename(all_files))]

cat(sprintf("Found %d junk files (%.1f MB)\n",
            length(junk), sum(file.size(junk), na.rm = TRUE) / 1024^2))
print(table(ifelse(startsWith(basename(junk), "._"), "AppleDouble (._*)", ".fuse_hidden*")))

if (dry_run) {
  print(utils::head(sub(paste0("^", root, "/"), "", junk), 50))
  message("Dry run only - set dry_run <- FALSE to delete.")
} else {
  ok <- suppressWarnings(file.remove(junk))
  cat(sprintf("Deleted %d of %d files\n", sum(ok), length(junk)))
  if (any(!ok)) {
    message("Could not delete (probably still open/locked):")
    print(junk[!ok])
  }
}
