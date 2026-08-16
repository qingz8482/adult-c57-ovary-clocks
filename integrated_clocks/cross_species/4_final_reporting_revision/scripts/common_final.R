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
SEED <- 20260720L
set.seed(SEED)

script_arg <- commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
caller_dir <- dirname(script_path)
final_dir <- if (basename(caller_dir) == "scripts") {
  normalizePath(file.path(caller_dir, ".."), mustWork = TRUE)
} else {
  normalizePath(caller_dir, mustWork = TRUE)
}
module_dir <- normalizePath(file.path(final_dir, ".."), mustWork = TRUE)
base_dir <- normalizePath(file.path(module_dir, "3_gtex_age_signal_adequacy"), mustWork = TRUE)
if (basename(final_dir) != "4_final_reporting_revision") stop("Unexpected final reporting directory: ", final_dir)
if (basename(base_dir) != "3_gtex_age_signal_adequacy") stop("Unexpected GTEx adequacy directory: ", base_dir)

pfinal <- function(...) file.path(final_dir, ...)
pbase <- function(...) file.path(base_dir, ...)
pmod <- function(...) file.path(module_dir, ...)
dir_create <- function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE)
invisible(lapply(pfinal(c("results", "figures", "figures/source_data", "reports", "derived_data", "logs")), dir_create))

write_tsv <- function(x, path) {
  dir_create(dirname(path))
  fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}

zsafe <- function(x) {
  x <- as.numeric(x); s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

fmt <- function(x, d = 3) ifelse(is.finite(x), formatC(x, digits = d, format = "f"), "NA")
p_string <- function(x) ifelse(!is.finite(x), "NA", ifelse(x < .001, formatC(x, digits = 2, format = "e"), formatC(x, digits = 3, format = "f")))
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

age_levels <- c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")
palette_age <- c("20-29"="#4C78A8", "30-39"="#72A0C1", "40-49"="#9ECAC3", "50-59"="#E3B65F", "60-69"="#E07A5F", "70-79"="#A8434C")
palette_contract <- c(neutral_dark="#272727", neutral_mid="#767676", neutral_light="#D8D8D8", signal_blue="#3182BD", signal_teal="#33B5A5", accent_red="#D24B40", accent_orange="#E28E2C")

theme_final <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(axis.line = element_line(linewidth=.35), axis.ticks = element_line(linewidth=.35),
      strip.background = element_blank(), strip.text = element_text(face="bold"),
      plot.title = element_text(face="bold", size=base_size+.7), panel.grid = element_blank())
}
theme_set(theme_final())

save_pub_r <- function(plot, stem, width_mm=183, height_mm=120, dpi=600) {
  dir_create(dirname(stem)); w <- width_mm/25.4; h <- height_mm/25.4
  grDevices::svg(paste0(stem,".svg"), width=w, height=h, pointsize=7, family="sans"); print(plot); dev.off()
  grDevices::cairo_pdf(paste0(stem,".pdf"), width=w, height=h, pointsize=7, family="sans"); print(plot); dev.off()
  ragg::agg_png(paste0(stem,".png"), width=w, height=h, units="in", res=dpi); print(plot); dev.off()
  ragg::agg_tiff(paste0(stem,".tiff"), width=w, height=h, units="in", res=dpi, compression="lzw"); print(plot); dev.off()
  invisible(stem)
}

read_gmt <- function(path, min_size=10L, max_size=500L) {
  spl <- strsplit(readLines(path, warn=FALSE), "\t", fixed=TRUE)
  nm <- vapply(spl, `[`, character(1), 1L)
  sets <- lapply(spl, function(z) unique(z[-c(1L,2L)])); names(sets) <- make.unique(nm)
  sets[lengths(sets)>=min_size & lengths(sets)<=max_size]
}

load_inputs <- function() {
  gtex <- readRDS(pmod("derived_data/gtex/gtex_ovary_prepared.rds"))
  qc <- readRDS(pbase("derived_data/gtex_age_signal_qc_objects.rds"))
  multi <- readRDS(pbase("derived_data/gtex_age_multivariate_objects.rds"))
  de <- readRDS(pbase("derived_data/gtex_age_de_objects.rds"))
  audit <- copy(qc$audit)
  audit[, age_bracket := factor(as.character(age_bracket), levels=age_levels, ordered=TRUE)]
  list(gtex=gtex, qc=qc, multi=multi, de=de, audit=audit)
}

weighted_median <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; w <- w[ok]; o <- order(x); x <- x[o]; w <- w[o]
  x[which(cumsum(w) >= sum(w)/2)[1]]
}

emp_p_two <- function(obs, null) {
  null <- null[is.finite(null)]
  (1 + sum(abs(null) >= abs(obs))) / (length(null) + 1)
}

material_passport <- function(type, status="ANALYZED") c(
  "## Material Passport", "",
  paste0("- Artifact type: ", type),
  paste0("- Verification status: ", status),
  "- Data scope: GTEx V11 ovary, 193 donors, frozen upstream adequacy results",
  "- Frozen boundary: no clock, Stage A/B, primary-model, or decision changes",
  paste0("- Revision seed: ", SEED),
  "- Plotting backend: R exclusively", ""
)
