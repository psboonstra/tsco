files <- c(
  "DESCRIPTION",
  "NAMESPACE",
  list.files("R", pattern = "\\.R$", full.names = TRUE),
  list.files("tests/testthat", pattern = "\\.R$", full.names = TRUE)
)

# Exclude scratch/check files if still present
files <- files[!grepl("check_", files)]

cat(
  paste0(
    "\n\n", strrep("=", 80), "\n",
    files,
    "\n", strrep("=", 80), "\n",
    vapply(files, function(f) paste(readLines(f, warn = FALSE), collapse = "\n"), character(1)),
    collapse = "\n"
  ),
  file = "dev/package-review.txt"
)
