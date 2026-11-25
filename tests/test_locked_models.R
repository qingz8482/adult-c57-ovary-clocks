#!/usr/bin/env Rscript
# Smoke test of the locked integrated clocks on synthetic input (no biological data):
# 1. predict_frozen.R runs on a genes x samples ln1p(CPM) CSV;
# 2. its predictions equal a manual replay from the exported text files
#    (raw-scale intercept + coefficients on ln1p(CPM); PCA: gene weights = loadings x PC coefficients / scale);
# 3. missing genes are imputed with the training mean (zero contribution);
# 4. the exported coefficients match Supplementary Tables 3a/3b.
suppressPackageStartupMessages(library(data.table))
fa <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo <- normalizePath(file.path(dirname(fa), ".."))
rel <- file.path(repo, "locked_models", "hvg_exclusion_vbpca_20261012")
rscript <- file.path(R.home("bin"), "Rscript")
tmp <- tempfile("locked_smoke_"); dir.create(tmp)
set.seed(20261125)

manual <- function(dir, x) {
  # Exported gene-level coefficients are in weeks per ln1p(CPM) unit and intercept.txt is the matching
  # raw-scale intercept, so prediction = intercept + sum(coefficient * x); unmeasured genes take the centre.
  co <- fread(file.path(dir, "coefficients.csv"))
  xi <- x[co$gene, , drop = FALSE]
  miss <- is.na(xi); xi[miss] <- matrix(co$center, nrow(xi), ncol(xi))[miss]
  b0 <- as.numeric(readLines(file.path(dir, "intercept.txt")))
  if (file.exists(file.path(dir, "pca_loadings.csv"))) {  # gene weights = loadings x PC coefficients / scale
    L <- as.matrix(fread(file.path(dir, "pca_loadings.csv"))[, -1])
    pc <- fread(file.path(dir, "pc_coefficients.csv"))$coefficient
    stopifnot(max(abs(as.numeric(L %*% pc) / co$scale - co$coefficient)) < 1e-10)
  }
  b0 + as.numeric(t(xi) %*% co$coefficient)
}

for (m in c("HVG-ElasticNet-exclude-single-age", "VersionB-PCA-glmnet")) {
  dir <- file.path(rel, m)
  sc <- fread(file.path(dir, "scaler.csv"))
  stopifnot(nrow(sc) == 3000L, identical(sc$gene, readLines(file.path(dir, "genes.txt"))))
  n <- 6L
  x <- matrix(pmax(0, sc$center + sc$scale * rnorm(nrow(sc) * n)), nrow(sc))
  rownames(x) <- sc$gene; colnames(x) <- paste0("synthetic_", seq_len(n))
  x[sample(nrow(x), 30), 1] <- NA  # unmeasured genes in one profile
  inp <- file.path(tmp, paste0(m, "_input.csv")); out <- file.path(tmp, paste0(m, "_pred.csv"))
  fwrite(data.table(gene = rownames(x), x), inp)
  status <- system2(rscript, c(shQuote(file.path(rel, "scripts", "predict_frozen.R")),
                               shQuote(file.path(dir, "model.rds")), shQuote(inp), shQuote(out)),
                    env = c("OPENBLAS_NUM_THREADS=1", "OMP_NUM_THREADS=1"))
  stopifnot(status == 0L, file.exists(out))
  pred <- fread(out)
  err <- max(abs(pred$pred - manual(dir, x)))
  stopifnot(nrow(pred) == n, err < 1e-8, abs(pred$missing_gene_fraction[1] - 0.01) < 1e-12)
  cat(sprintf("%s: predict_frozen.R matches text-export replay (max |diff| = %.2e weeks)\n", m, err))
}

# Exported coefficients versus the published Supplementary Tables 3a/3b.
en <- fread(file.path(rel, "HVG-ElasticNet-exclude-single-age", "coefficients.csv"))[coefficient != 0]
st3a <- fread(file.path(repo, "metadata", "manuscript_tables", "ST3a_EN_nonzero_coefficients.tsv"), skip = 1)
stopifnot(nrow(en) == 79L, nrow(st3a) == 79L, setequal(en$gene, st3a$gene))
stopifnot(max(abs(en$coefficient - as.numeric(st3a[[2]][match(en$gene, st3a$gene)]))) < 1e-9)
pc <- fread(file.path(rel, "VersionB-PCA-glmnet", "pc_coefficients.csv"))
stopifnot(identical(which(pc$coefficient != 0), c(1L, 2L, 5L, 8L)))
cat("Exported coefficients match Supplementary Table 3a (79 EN genes); PCA/EN active PCs = PC1, PC2, PC5, PC8\n")
unlink(tmp, recursive = TRUE)
cat("locked-model smoke test passed\n")
