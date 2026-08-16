suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(patchwork)
  library(MASS)
  library(glmnet)
  library(matrixStats)
  library(Matrix)
})

options(stringsAsFactors = FALSE, scipen = 999)
SEED <- 20260719L
set.seed(SEED)

script_arg <- commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
audit_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
module_dir <- normalizePath(file.path(audit_dir, ".."), mustWork = TRUE)
if (basename(audit_dir) != "3_gtex_age_signal_adequacy") stop("Unexpected GTEx adequacy directory: ", audit_dir)

paudit <- function(...) file.path(audit_dir, ...)
pmod <- function(...) file.path(module_dir, ...)
dir_create <- function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE)
invisible(lapply(paudit(c("results", "figures", "figures/source_data", "reports", "derived_data", "logs", "external_data")), dir_create))

write_tsv <- function(x, path) {
  dir_create(dirname(path))
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

zsafe <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

safe_num <- function(x) suppressWarnings(as.numeric(x))
safe_name <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("(^_+|_+$)", "", x)
}

fmt <- function(x, digits = 3) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")
}

p_string <- function(x) {
  ifelse(!is.finite(x), "NA", ifelse(x < 0.001, formatC(x, digits = 2, format = "e"), formatC(x, digits = 3, format = "f")))
}

palette_age <- c(
  "20-29" = "#4C78A8", "30-39" = "#72A0C1", "40-49" = "#9ECAC3",
  "50-59" = "#E3B65F", "60-69" = "#E07A5F", "70-79" = "#A8434C"
)
palette_contract <- c(
  neutral_dark = "#272727", neutral_mid = "#767676", neutral_light = "#D8D8D8",
  signal_blue = "#3182BD", signal_teal = "#33B5A5", accent_red = "#D24B40",
  accent_orange = "#E28E2C", elasticnet = "#2878B5", pca_glmnet = "#D65F5F"
)

theme_nature_contract <- function(base_size = 7, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.title = element_text(size = base_size + 0.7, face = "bold"),
      panel.grid = element_blank()
    )
}
theme_set(theme_nature_contract())

save_pub_r <- function(plot, stem, width_mm = 183, height_mm = 120, dpi = 600) {
  dir_create(dirname(stem))
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  grDevices::svg(paste0(stem, ".svg"), width = w, height = h, pointsize = 7, family = "sans")
  print(plot); grDevices::dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h, pointsize = 7, family = "sans")
  print(plot); grDevices::dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = w, height = h, units = "in", res = dpi)
  print(plot); grDevices::dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in", res = dpi, compression = "lzw")
  print(plot); grDevices::dev.off()
  invisible(stem)
}

save_base_bundle <- function(stem, draw_fun, width_mm = 183, height_mm = 120, dpi = 600) {
  dir_create(dirname(stem))
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  grDevices::svg(paste0(stem, ".svg"), width = w, height = h, pointsize = 7, family = "sans")
  draw_fun(); grDevices::dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h, pointsize = 7, family = "sans")
  draw_fun(); grDevices::dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = w, height = h, units = "in", res = dpi)
  draw_fun(); grDevices::dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in", res = dpi, compression = "lzw")
  draw_fun(); grDevices::dev.off()
  invisible(stem)
}

read_gmt <- function(path, min_size = 10L, max_size = 500L) {
  ln <- readLines(path, warn = FALSE)
  spl <- strsplit(ln, "\t", fixed = TRUE)
  nm <- vapply(spl, `[`, character(1), 1L)
  sets <- lapply(spl, function(z) unique(z[-c(1L, 2L)]))
  names(sets) <- make.unique(nm)
  sets[lengths(sets) >= min_size & lengths(sets) <= max_size]
}

empirical_summary <- function(observed, null, tail = c("greater", "less", "two.sided")) {
  tail <- match.arg(tail)
  p <- switch(tail,
    greater = (1 + sum(null >= observed, na.rm = TRUE)) / (length(null) + 1),
    less = (1 + sum(null <= observed, na.rm = TRUE)) / (length(null) + 1),
    two.sided = (1 + sum(abs(null - median(null, na.rm = TRUE)) >= abs(observed - median(null, na.rm = TRUE)), na.rm = TRUE)) / (length(null) + 1)
  )
  data.table(
    observed = observed, null_median = median(null, na.rm = TRUE),
    null_q025 = unname(quantile(null, 0.025, na.rm = TRUE)),
    null_q975 = unname(quantile(null, 0.975, na.rm = TRUE)),
    empirical_p = p, monte_carlo_se = sqrt(p * (1 - p) / (length(null) + 1))
  )
}

max_pairwise_smd <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(factor(g[ok]))
  lev <- levels(g)
  if (length(lev) < 2L) return(NA_real_)
  vals <- c()
  for (i in seq_len(length(lev) - 1L)) for (j in (i + 1L):length(lev)) {
    a <- x[g == lev[i]]; b <- x[g == lev[j]]
    sp <- sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / max(1, length(a) + length(b) - 2))
    vals <- c(vals, if (is.finite(sp) && sp > 0) abs(mean(a) - mean(b)) / sp else NA_real_)
  }
  vals <- vals[is.finite(vals)]
  if (!length(vals)) NA_real_ else max(vals)
}

cramers_v <- function(tab) {
  if (nrow(tab) < 2L || ncol(tab) < 2L) return(NA_real_)
  cs <- suppressWarnings(chisq.test(tab, correct = FALSE))
  sqrt(unname(cs$statistic) / (sum(tab) * max(1, min(dim(tab)) - 1)))
}

association_numeric <- function(x, age_factor, age_order) {
  ok <- is.finite(x) & !is.na(age_factor)
  if (sum(ok) < 10L || sd(x[ok]) == 0) return(list(p_omnibus = NA_real_, p_ordered = NA_real_, rho = NA_real_))
  a <- tryCatch(anova(lm(x[ok] ~ droplevels(factor(age_factor[ok]))))[1, "Pr(>F)"], error = function(e) NA_real_)
  s <- suppressWarnings(cor.test(x[ok], age_order[ok], method = "spearman", exact = FALSE))
  list(p_omnibus = as.numeric(a), p_ordered = s$p.value, rho = unname(s$estimate))
}

association_categorical <- function(x, age_factor) {
  ok <- !is.na(x) & x != "" & !is.na(age_factor)
  tab <- table(factor(age_factor[ok]), factor(x[ok]))
  if (nrow(tab) < 2L || ncol(tab) < 2L) return(list(p = NA_real_, v = NA_real_))
  tst <- suppressWarnings(chisq.test(tab, simulate.p.value = any(tab < 5), B = 5000))
  list(p = tst$p.value, v = cramers_v(tab))
}

age_levels <- c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")

load_core_inputs <- function() {
  gtex <- readRDS(pmod("derived_data/gtex/gtex_ovary_prepared.rds"))
  sa <- fread(pmod("derived_data/gtex/raw/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"))
  sp <- fread(pmod("derived_data/gtex/raw/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt"))
  meta <- copy(gtex$metadata)
  meta[, age_bracket := factor(as.character(age_bracket), levels = age_levels, ordered = TRUE)]
  meta[, `:=`(age_order = as.integer(age_bracket) - 1L, age_decade = as.integer(age_bracket) - 1L)]
  sa <- sa[match(meta$sample_id, SAMPID)]
  sp <- sp[match(meta$donor_id, SUBJID)]
  proj <- fread(pmod("derived_data/stageB/gtex_frozen_zero_shot_scores.tsv"))
  arts <- readRDS(pmod("derived_data/mouse_models/authoritative_artifacts.rds"))
  stagea <- fread(pmod("derived_data/stageA/mouse_human_gene_effect_master.tsv"))
  list(gtex = gtex, sa = sa, sp = sp, meta = meta, proj = proj, arts = arts, stagea = stagea)
}
