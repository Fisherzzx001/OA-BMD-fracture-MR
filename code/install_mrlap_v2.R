options(repos = "https://cloud.r-project.org")
pkgdir <- "C:/Users/Administrator/WorkBuddy/2026-09-08-22-17-40/OA_B_PROJECT/code/ref/pkg"

cat("=== install GenomicSEM CRAN deps ===\n")
deps <- c("e1071", "gdata", "lavaan", "doParallel", "splitstackshape",
          "R.utils", "mgsub", "simsalapar")
install.packages(deps)
for (p in deps) cat(p, ":", as.character(requireNamespace(p, quietly = TRUE)), "\n")

cat("\n=== install GenomicSEM from local tarball (deps=FALSE) ===\n")
r1 <- tryCatch({
  install.packages(file.path(pkgdir, "GenomicSEM.tar.gz"), repos = NULL,
                   type = "source", dependencies = FALSE)
  "OK"
}, error = function(e) paste("ERR:", conditionMessage(e)))
cat("GenomicSEM install:", as.character(r1), "\n")

cat("\n=== install MRlap from local tarball (deps=FALSE) ===\n")
r2 <- tryCatch({
  install.packages(file.path(pkgdir, "MRlap.tar.gz"), repos = NULL,
                   type = "source", dependencies = FALSE)
  "OK"
}, error = function(e) paste("ERR:", conditionMessage(e)))
cat("MRlap install:", as.character(r2), "\n")

cat("\n=== final verify ===\n")
for (p in c("MRlap", "GenomicSEM", "OpenMx", "TwoSampleMR", "ieugwasr",
            "lavaan", "e1071")) {
  cat(p, ":", as.character(requireNamespace(p, quietly = TRUE)), "\n")
}
cat("ALL_DONE\n")
