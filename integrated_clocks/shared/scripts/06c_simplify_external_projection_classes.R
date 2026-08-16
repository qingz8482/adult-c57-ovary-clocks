#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- script_arg[grepl("^--file=", script_arg)]
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else file.path(getwd(), "scripts")
if (!dir.exists(script_dir)) {
  script_dir <- file.path(getwd(), "adult_c57_ovary_clocks", "bulk-pseudobulk_clock", "scripts")
}

project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
results_dir <- Sys.getenv("RESULTS_ROOT", unset = file.path(project_dir, "results"))
results_dir <- normalizePath(results_dir, mustWork = TRUE)
if (grepl("results_groupmapped_v2", results_dir, fixed = TRUE) || grepl("_v2", basename(results_dir), fixed = TRUE)) {
  stop("Refusing to simplify v1 external projection classes against a v2 results root: ", results_dir, call. = FALSE)
}

projection_dir <- file.path(results_dir, "external_projection")
out_dir <- file.path(projection_dir, "simplified_classes")
plot_dir <- file.path(out_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

allowed_projection_classes <- c(
  "reference_control", "genetic_model", "diet", "microbiome",
  "intervention", "stress", "developmental_OOD", "other"
)
allowed_projection_roles <- c("reference_control", "external_treatment", "OOD_supplement")

read_required <- function(filename) {
  path <- file.path(projection_dir, filename)
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
  fread(path, sep = "\t", na.strings = c("", "NA", "NaN"))
}

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- NA_character_
  x
}

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}

ensure_cols <- function(dt, cols) {
  for (cc in cols) if (!cc %in% names(dt)) dt[, (cc) := NA]
  dt
}

clean_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nchar(x) == 0, "plot", x)
}

theme_projection <- function() {
  theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      strip.background = element_rect(fill = "grey92", color = "grey70"),
      legend.position = "right"
    )
}

save_plot <- function(p, stem, width = 12, height = 7) {
  ggsave(file.path(plot_dir, paste0(stem, ".pdf")), p, width = width, height = height, limitsize = FALSE)
  ggsave(file.path(plot_dir, paste0(stem, ".png")), p, width = width, height = height, dpi = 300, limitsize = FALSE)
}

label_from_text <- function(x, fallback = "other") {
  x <- safe_chr(x)
  x[is.na(x)] <- fallback
  x
}

first_nonempty <- function(...) {
  pieces <- list(...)
  out <- rep(NA_character_, length(pieces[[1]]))
  for (x in pieces) {
    x <- safe_chr(x)
    take <- is.na(out) & !is.na(x)
    out[take] <- x[take]
  }
  out
}

add_simplified_classes <- function(dt) {
  dt <- copy(dt)
  needed <- c(
    "sample_key", "dataset_folder", "projection_class", "projection_role",
    "projection_reason", "matched_reference_control_flag", "matched_reference_control_source",
    "used_in_deployment_training_flag", "primary_external_projection_include",
    "primary_effect_interpretation_include", "reference_control_group_id",
    "reference_control_group_label", "gse279446_cohort", "treatment", "genotype",
    "sample_context", "sample_title", "sample_id", "accession", "matched_control_method",
    "relative_delta_age_week", "projection_group", "comparison_id", "condition",
    "model", "injection_age", "post_injection_days", "age_group"
  )
  dt <- ensure_cols(dt, needed)
  for (cc in c("matched_reference_control_flag", "primary_external_projection_include", "primary_effect_interpretation_include")) {
    dt[, (cc) := as_bool(get(cc))]
  }

  original_class <- safe_chr(dt$projection_class)
  original_role <- safe_chr(dt$projection_role)
  original_control_source <- safe_chr(dt$matched_reference_control_source)
  original_projection_group <- safe_chr(dt$projection_group)
  original_comparison_id <- safe_chr(dt$comparison_id)
  text <- tolower(paste(
    safe_chr(dt$dataset_folder), safe_chr(dt$sample_key), safe_chr(dt$projection_class),
    safe_chr(dt$projection_role), safe_chr(dt$projection_reason), safe_chr(dt$reference_control_group_id),
    safe_chr(dt$reference_control_group_label), safe_chr(dt$treatment), safe_chr(dt$genotype),
    safe_chr(dt$sample_context), safe_chr(dt$sample_title), safe_chr(dt$sample_id), safe_chr(dt$accession),
    sep = " "
  ))
  dt[, recode_text := text]
  control_text <- tolower(paste(
    original_class, original_role, original_control_source,
    safe_chr(dt$treatment), safe_chr(dt$genotype), safe_chr(dt$sample_context),
    safe_chr(dt$sample_title), safe_chr(dt$sample_id), safe_chr(dt$accession),
    sep = " "
  ))

  dt[, projection_class_original := original_class]
  dt[, projection_class := "other"]
  dt[, projection_role := "external_treatment"]
  dt[, projection_group := first_nonempty(original_projection_group, treatment, sample_title, projection_class_original)]
  dt[, comparison_id := first_nonempty(original_comparison_id, reference_control_group_id, matched_control_method)]
  dt[, genotype_class := fifelse(grepl("\\bwt\\b|control", text), "WT", NA_character_)]
  dt[, cohort := first_nonempty(gse279446_cohort, bio_age_bin)]
  dt[, control_group_label := first_nonempty(reference_control_group_label, reference_control_group_id)]
  dt[, plot_include_primary := primary_external_projection_include %in% TRUE]
  dt[, plot_include_ood_supplement := FALSE]

  is_control <- dt$matched_reference_control_flag %in% TRUE |
    grepl("reference_control|internal_control|control_meta|same_study.*control|wt.*control|control_diet|young.*donor.*microbiota|wt chow fertile", control_text)
  is_ood <- dt$primary_external_projection_include %in% FALSE |
    grepl("ood|developmental|prepubertal|younger_than_7|outside current age-bin|age_bin", text)
  is_stress <- grepl("stress|cums", text)
  is_microbiome <- grepl("fmt|microbiota|donor", text)
  is_diet <- grepl("hfhs|diet|high-fat|high fat|high-sucrose|high sucrose|vegetable-oil-free|vof|metabolic", text)
  is_genetic <- grepl("bnc1|zeb1_200h|zeb1_200m|rnf111|abcg5|abcg8|ko\\b|cko\\b|transgenic|heterozygous|homozygous|mutant", text)
  is_intervention <- grepl("intervention|treatment|drug|nmn|vehicle|water control|plus", text)

  dt[is_control, `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    projection_group = "reference control"
  )]
  dt[!is_control & is_genetic, projection_class := "genetic_model"]
  dt[!is_control & !is_genetic & is_diet, projection_class := "diet"]
  dt[!is_control & !is_genetic & !is_diet & is_microbiome, projection_class := "microbiome"]
  dt[!is_control & !is_genetic & !is_diet & !is_microbiome & is_stress, projection_class := "stress"]
  dt[!is_control & !is_genetic & !is_diet & !is_microbiome & !is_stress & is_intervention, projection_class := "intervention"]
  dt[projection_class != "reference_control", projection_role := "external_treatment"]
  dt[is_ood, `:=`(
    projection_class = "developmental_OOD",
    projection_role = "OOD_supplement",
    projection_group = "age-bin OOD",
    plot_include_primary = FALSE,
    plot_include_ood_supplement = TRUE
  )]

  dt[grepl("heterozygous|zeb1_200h", text), genotype_class := "heterozygous"]
  dt[grepl("homozygous|zeb1_200m", text), genotype_class := "homozygous"]
  dt[grepl("transgenic|bnc1", text), genotype_class := "transgenic"]
  dt[grepl("\\bko\\b|cko\\b|abcg5|abcg8|rnf111", text), genotype_class := "KO"]
  dt[projection_class == "reference_control" & is.na(genotype_class), genotype_class := "WT"]

  ds <- dt$dataset_folder

  dt[ds == "GSE194196" & grepl("WT-[147]|WT ovary|\\bWT\\b", text), `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    projection_group = "WT control",
    control_group_label = "WT-1/WT-4/WT-7",
    comparison_id = "GSE194196_Bnc1_vs_WT",
    genotype_class = "WT"
  )]
  dt[ds == "GSE194196" & grepl("bnc1|chun", text), `:=`(
    projection_class = "genetic_model",
    projection_role = "external_treatment",
    projection_group = "Bnc1tr/tr",
    comparison_id = "GSE194196_Bnc1_vs_WT",
    genotype_class = "transgenic"
  )]

  dt[ds == "GSE252791" & grepl("\\bWT\\b|WT ovary|control", text), `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    projection_group = "WT control",
    control_group_label = "WT-1/WT-4/WT-7",
    comparison_id = "GSE252791_Rnf111_vs_WT",
    genotype_class = "WT"
  )]
  dt[ds == "GSE252791" & grepl("rnf111|heterozygous|mutant|ko", text) & !grepl("\\bWT\\b", text), `:=`(
    projection_class = "genetic_model",
    projection_role = "external_treatment",
    projection_group = "Rnf111 heterozygous KO",
    comparison_id = "GSE252791_Rnf111_vs_WT",
    genotype_class = "heterozygous"
  )]

  dt[ds == "GSE279446", cohort := fifelse(grepl("young", recode_text), "young", "adult")]
  dt[ds == "GSE279446" & grepl("zeb1_wt|wt_ovary|wt ovary", text), `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    genotype_class = "WT",
    projection_group = paste(cohort, "WT"),
    comparison_id = paste0("GSE279446_", cohort, "_Zeb1_vs_WT")
  )]
  dt[ds == "GSE279446" & grepl("zeb1_200h", text), `:=`(
    projection_class = "genetic_model",
    projection_role = "external_treatment",
    genotype_class = "heterozygous",
    projection_group = paste(cohort, "Zeb1_200H"),
    comparison_id = paste0("GSE279446_", cohort, "_Zeb1_vs_WT")
  )]
  dt[ds == "GSE279446" & grepl("zeb1_200m", text), `:=`(
    projection_class = "genetic_model",
    projection_role = "external_treatment",
    genotype_class = "homozygous",
    projection_group = paste(cohort, "Zeb1_200M"),
    comparison_id = paste0("GSE279446_", cohort, "_Zeb1_vs_WT")
  )]

  dt[ds == "GSE279573" & grepl("wt chow fertile", text), `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    genotype_class = "WT",
    projection_group = "WT chow fertile",
    control_group_label = "WT chow fertile",
    comparison_id = "GSE279573_KO_or_rescue_vs_WT"
  )]
  dt[ds == "GSE279573" & grepl("ko chow infertile", text), `:=`(
    projection_class = "genetic_model",
    projection_role = "external_treatment",
    genotype_class = "KO",
    projection_group = "KO chow infertile",
    comparison_id = "GSE279573_KO_or_rescue_vs_WT"
  )]
  dt[ds == "GSE279573" & grepl("vegetable-oil-free|vof", text), `:=`(
    projection_class = "diet",
    projection_role = "external_treatment",
    genotype_class = "KO",
    projection_group = "KO vegetable-oil-free rescue",
    comparison_id = "GSE279573_KO_or_rescue_vs_WT"
  )]

  dt[ds == "GSE240547" & grepl("control_diet|control diet", text), `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    projection_group = "control diet",
    control_group_label = "control diet",
    comparison_id = "GSE240547_HFHS_vs_control_diet"
  )]
  dt[ds == "GSE240547" & grepl("hfhs|high-fat|high fat|high-sucrose|high sucrose", text), `:=`(
    projection_class = "diet",
    projection_role = "external_treatment",
    projection_group = "HFHS diet",
    comparison_id = "GSE240547_HFHS_vs_control_diet"
  )]

  dt[ds == "PRJNA1076509" & grepl("young female donor|fmt-yf|young_female", text), `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    projection_group = "young female donor FMT",
    control_group_label = "young female donor FMT",
    comparison_id = "PRJNA1076509_estropausal_vs_young_FMT"
  )]
  dt[ds == "PRJNA1076509" & grepl("estropausal|fmt-of", text), `:=`(
    projection_class = "microbiome",
    projection_role = "external_treatment",
    projection_group = "estropausal female donor FMT",
    comparison_id = "PRJNA1076509_estropausal_vs_young_FMT"
  )]

  prjna863443 <- ds == "PRJNA863443"
  dt[prjna863443 & grepl("natural ageing reference|natural_aging_reference", text), `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    comparison_id = "PRJNA863443_natural_aging_reference",
    genotype_class = "WT"
  )]
  dt[prjna863443 & grepl("^ctl_| ctl_|condition ctl|\\bctl\\b|vehicle", text), `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    genotype_class = "WT"
  )]
  dt[prjna863443 & grepl("^vcd_| vcd_|condition vcd|\\bvcd\\b", text), `:=`(
    projection_class = "intervention",
    projection_role = "external_treatment",
    genotype_class = "WT"
  )]
  dt[prjna863443 & grepl("foxl2\\+/\\+|foxl2_wt", text), `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    genotype_class = "WT"
  )]
  dt[prjna863443 & grepl("foxl2\\+/-|foxl2_het", text), `:=`(
    projection_class = "genetic_model",
    projection_role = "external_treatment",
    genotype_class = "heterozygous"
  )]
  restore_prjna_group <- prjna863443 & !is.na(original_projection_group)
  restore_prjna_comparison <- prjna863443 & !is.na(original_comparison_id)
  dt[restore_prjna_group, projection_group := original_projection_group[restore_prjna_group]]
  dt[restore_prjna_comparison, comparison_id := original_comparison_id[restore_prjna_comparison]]
  dt[prjna863443 & comparison_id == "PRJNA863443_natural_aging_reference", control_group_label := projection_group]
  dt[prjna863443 & projection_role == "reference_control" & is.na(control_group_label), control_group_label := projection_group]

  dt[ds == "GSE228975", `:=`(
    projection_class = "developmental_OOD",
    projection_role = "OOD_supplement",
    projection_group = "age-bin OOD",
    comparison_id = NA_character_,
    plot_include_primary = FALSE,
    plot_include_ood_supplement = TRUE
  )]

  dt[ds == "CRA003645" & projection_role == "reference_control", projection_group := "postpubertal WT reference"]
  dt[ds == "CRA003645" & projection_class == "developmental_OOD", projection_group := fifelse(grepl("3w|3 w", recode_text), "3-week prepubertal ovary", "0-day developmental ovary")]
  dt[ds == "CRA003645" & is.na(comparison_id), comparison_id := "CRA003645_developmental_vs_postpubertal_WT"]

  dt[ds == "CRA017492" & grepl("young 2 months control", text), `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    projection_group = "young control",
    control_group_label = "young 2 months control",
    comparison_id = "CRA017492_NMN_vs_aged_water_with_young_reference"
  )]
  dt[ds == "CRA017492" & grepl("water control", text), `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    projection_group = "aged water control",
    control_group_label = "aged 12 months water control",
    comparison_id = "CRA017492_NMN_vs_aged_water_with_young_reference"
  )]
  dt[ds == "CRA017492" & grepl("nmn|plus", text), `:=`(
    projection_class = "intervention",
    projection_role = "external_treatment",
    projection_group = "aged NMN",
    control_group_label = "aged 12 months water control",
    comparison_id = "CRA017492_NMN_vs_aged_water_with_young_reference"
  )]

  dt[ds == "GSE273724" & projection_role == "reference_control", `:=`(
    projection_group = "reference control",
    comparison_id = "GSE273724_intervention_vs_control"
  )]
  dt[ds == "GSE273724" & projection_role != "reference_control", `:=`(
    projection_class = "intervention",
    projection_group = label_from_text(treatment, "intervention"),
    comparison_id = "GSE273724_intervention_vs_control"
  )]

  dt[ds == "GSE329024" & projection_role == "reference_control", `:=`(
    projection_group = fcase(
      grepl("2m[-_ ]*8w", recode_text), "2M-8W control",
      grepl("2m[-_ ]*12w", recode_text), "2M-12W control",
      grepl("6m[-_ ]*8w", recode_text), "6M-8W control",
      grepl("6m[-_ ]*12w", recode_text), "6M-12W control",
      default = NA_character_
    ),
    comparison_id = fcase(
      grepl("2m[-_ ]*8w", recode_text), "GSE329024_CUMS_2M_8W_vs_control",
      grepl("2m[-_ ]*12w", recode_text), "GSE329024_CUMS_2M_12W_vs_control",
      grepl("6m[-_ ]*8w", recode_text), "GSE329024_CUMS_6M_8W_vs_control",
      grepl("6m[-_ ]*12w", recode_text), "GSE329024_CUMS_6M_12W_vs_control",
      default = NA_character_
    )
  )]
  dt[ds == "GSE329024" & projection_role != "reference_control", `:=`(
    projection_class = "stress",
    projection_group = fcase(
      grepl("2m[-_ ]*8w", recode_text), "2M-8W CUMS stress",
      grepl("2m[-_ ]*12w", recode_text), "2M-12W CUMS stress",
      grepl("6m[-_ ]*8w", recode_text), "6M-8W CUMS stress",
      grepl("6m[-_ ]*12w", recode_text), "6M-12W CUMS stress",
      default = NA_character_
    ),
    comparison_id = fcase(
      grepl("2m[-_ ]*8w", recode_text), "GSE329024_CUMS_2M_8W_vs_control",
      grepl("2m[-_ ]*12w", recode_text), "GSE329024_CUMS_2M_12W_vs_control",
      grepl("6m[-_ ]*8w", recode_text), "GSE329024_CUMS_6M_8W_vs_control",
      grepl("6m[-_ ]*12w", recode_text), "GSE329024_CUMS_6M_12W_vs_control",
      default = NA_character_
    )
  )]

  dt[ds == "GSE235021" & grepl("polished rice", recode_text), `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    projection_group = fifelse(age_week < 30, "6-month polished rice", "9-month polished rice"),
    control_group_label = fifelse(
      age_week < 30,
      "GSE235021 6-month polished-rice controls",
      "GSE235021 9-month polished-rice controls"
    ),
    comparison_id = fifelse(
      age_week < 30,
      "GSE235021_rice_diets_vs_polished_6m",
      "GSE235021_rice_diets_vs_polished_9m"
    ),
    genotype_class = "WT"
  )]
  dt[ds == "GSE235021" & !grepl("polished rice", recode_text), `:=`(
    projection_class = "diet",
    projection_role = "external_treatment",
    projection_group = label_from_text(treatment, "rice diet"),
    comparison_id = fifelse(
      age_week < 30,
      "GSE235021_rice_diets_vs_polished_6m",
      "GSE235021_rice_diets_vs_polished_9m"
    ),
    genotype_class = "WT"
  )]

  dt[ds == "GSE154890", `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    projection_group = paste0(round(age_week, 1), "-week natural-age reference"),
    comparison_id = "GSE154890_natural_age_external_placement",
    genotype_class = "WT"
  )]

  dt[ds == "GSE249151" & grepl("sesame oil control|vehicle control", recode_text), `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    projection_group = "sesame-oil vehicle control",
    control_group_label = "GSE249151 sesame-oil vehicle controls",
    comparison_id = "GSE249151_VCD_vs_vehicle",
    genotype_class = "WT"
  )]
  dt[ds == "GSE249151" & grepl("\\bVCD\\b", recode_text, ignore.case = TRUE), `:=`(
    projection_class = "intervention",
    projection_role = "external_treatment",
    projection_group = label_from_text(treatment, "VCD intervention"),
    comparison_id = "GSE249151_VCD_vs_vehicle",
    genotype_class = "WT"
  )]

  dt[ds == "GSE106339" & grepl("goat|knockout", safe_chr(genotype), ignore.case = TRUE), `:=`(
    projection_class = "genetic_model",
    projection_role = "external_treatment",
    projection_group = "3-week GOAT KO",
    comparison_id = "GSE106339_prepubertal_GOAT_KO_vs_WT",
    genotype_class = "KO",
    plot_include_primary = TRUE,
    plot_include_ood_supplement = FALSE
  )]
  dt[ds == "GSE106339" & toupper(trimws(safe_chr(genotype))) == "WT", `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    projection_group = "3-week WT",
    control_group_label = "GSE106339 3-week WT controls",
    comparison_id = "GSE106339_prepubertal_GOAT_KO_vs_WT",
    genotype_class = "WT",
    plot_include_primary = TRUE,
    plot_include_ood_supplement = FALSE
  )]

  dt[ds == "GSE128240" & grepl("(^|[^a-z])control([^a-z]|$)", recode_text) &
       !grepl("cyclophosphamide|ctx", recode_text), `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    projection_group = "same-study control",
    control_group_label = "GSE128240 same-study controls",
    comparison_id = "GSE128240_CTX_vs_control",
    genotype_class = "WT"
  )]
  dt[ds == "GSE128240" & grepl("cyclophosphamide|ctx", recode_text), `:=`(
    projection_class = "intervention",
    projection_role = "external_treatment",
    projection_group = "cyclophosphamide",
    comparison_id = "GSE128240_CTX_vs_control",
    genotype_class = "WT"
  )]

  dt[ds == "GSE98497" & grepl("(^|[^a-z0-9])wt([^a-z0-9]|$)", recode_text) &
       !grepl("cko|knockout", recode_text), `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    projection_group = "WT control",
    control_group_label = "GSE98497 WT controls",
    comparison_id = "GSE98497_Mtor_cKO_vs_WT",
    genotype_class = "WT"
  )]
  dt[ds == "GSE98497" & grepl("cko|knockout", recode_text), `:=`(
    projection_class = "genetic_model",
    projection_role = "external_treatment",
    projection_group = "Mtor oocyte-conditional KO",
    comparison_id = "GSE98497_Mtor_cKO_vs_WT",
    genotype_class = "KO"
  )]

  dt[ds == "PRJNA787028" & grepl("2 mice per cage|2pc", recode_text), `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    projection_group = "2 mice per cage",
    control_group_label = "PRJNA787028 two-mice-per-cage controls",
    comparison_id = "PRJNA787028_high_vs_low_housing_density",
    genotype_class = "WT"
  )]
  dt[ds == "PRJNA787028" & grepl("8 mice per cage|8pc", recode_text), `:=`(
    projection_class = "stress",
    projection_role = "external_treatment",
    projection_group = "8 mice per cage",
    comparison_id = "PRJNA787028_high_vs_low_housing_density",
    genotype_class = "WT"
  )]

  prjna1002222 <- ds == "PRJNA1002222"
  prjna1002222_wt <- prjna1002222 & toupper(trimws(safe_chr(dt$genotype))) == "WT"
  prjna1002222_ko <- prjna1002222 & (
    grepl("cd38", safe_chr(dt$genotype), ignore.case = TRUE) |
      grepl("__Cd38KO_", safe_chr(dt$sample_key), fixed = TRUE)
  )
  prjna1002222_age <- round(dt$age_week * 12 / 52.1775)
  dt[prjna1002222_wt & prjna1002222_age == 2, `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    projection_group = "WT 2M young reference",
    control_group_label = "WT 2M young reference",
    comparison_id = "PRJNA1002222_WT_natural_aging_reference_2M",
    genotype_class = "WT"
  )]
  for (age_month in c(8, 12)) {
    age_label <- paste0(age_month, "M")
    dt[prjna1002222_wt & prjna1002222_age == age_month, `:=`(
      projection_class = "reference_control",
      projection_role = "reference_control",
      projection_group = paste0("WT_", age_label),
      control_group_label = paste0("WT_", age_label),
      comparison_id = paste0("PRJNA1002222_Cd38KO_vs_WT_", age_label),
      genotype_class = "WT"
    )]
    dt[prjna1002222_ko & prjna1002222_age == age_month, `:=`(
      projection_class = "genetic_model",
      projection_role = "external_treatment",
      projection_group = paste0("Cd38KO_", age_label),
      control_group_label = paste0("WT_", age_label),
      comparison_id = paste0("PRJNA1002222_Cd38KO_vs_WT_", age_label),
      genotype_class = "KO"
    )]
  }

  dt[projection_class == "reference_control", projection_role := "reference_control"]
  dt[projection_class == "developmental_OOD", `:=`(
    projection_role = "OOD_supplement",
    plot_include_primary = FALSE,
    plot_include_ood_supplement = TRUE
  )]
  dt[projection_class != "developmental_OOD" & is.na(plot_include_primary), plot_include_primary := TRUE]
  dt[projection_class != "developmental_OOD" & is.na(plot_include_ood_supplement), plot_include_ood_supplement := FALSE]
  dt[is.na(projection_group) | projection_group == "", projection_group := projection_class]
  dt[is.na(cohort) | cohort == "", cohort := "not_applicable"]
  dt[is.na(control_group_label) & projection_role == "reference_control", control_group_label := projection_group]
  dt[is.na(genotype_class) | genotype_class == "", genotype_class := "not_applicable"]

  dt[, recode_text := NULL]
  dt[]
}

validate_simplified <- function(tables) {
  checks <- list()
  add_check <- function(table, check, passed, detail = "") {
    checks[[length(checks) + 1L]] <<- data.table(
      table = table,
      check = check,
      passed = isTRUE(passed),
      detail = as.character(detail)
    )
  }

  for (nm in names(tables)) {
    dt <- tables[[nm]]
    add_check(nm, "projection_class_original_present", "projection_class_original" %in% names(dt) && all(!is.na(dt$projection_class_original)), paste0("missing=", sum(is.na(dt$projection_class_original))))
    add_check(nm, "projection_class_present", "projection_class" %in% names(dt) && all(!is.na(dt$projection_class)), paste0("missing=", sum(is.na(dt$projection_class))))
    add_check(nm, "projection_role_present", "projection_role" %in% names(dt) && all(!is.na(dt$projection_role)), paste0("missing=", sum(is.na(dt$projection_role))))
    add_check(nm, "projection_class_allowed", all(unique(dt$projection_class) %in% allowed_projection_classes), paste(setdiff(unique(dt$projection_class), allowed_projection_classes), collapse = ";"))
    add_check(nm, "projection_role_allowed", all(unique(dt$projection_role) %in% allowed_projection_roles), paste(setdiff(unique(dt$projection_role), allowed_projection_roles), collapse = ";"))
    add_check(nm, "no_more_than_8_projection_classes", uniqueN(dt$projection_class) <= 8, paste0("n=", uniqueN(dt$projection_class)))
    add_check(nm, "no_more_than_3_projection_roles", uniqueN(dt$projection_role) <= 3, paste0("n=", uniqueN(dt$projection_role)))
    primary_missing_group <- dt[(plot_include_primary %in% TRUE) & (is.na(projection_group) | projection_group == "")]
    add_check(nm, "primary_rows_have_projection_group", nrow(primary_missing_group) == 0, paste0("rows=", nrow(primary_missing_group)))
    if ("relative_delta_age_week" %in% names(dt)) {
      rel_missing_comp <- dt[is.finite(suppressWarnings(as.numeric(relative_delta_age_week))) & (is.na(comparison_id) | comparison_id == "")]
      add_check(nm, "finite_relative_delta_rows_have_comparison_id", nrow(rel_missing_comp) == 0, paste0("rows=", nrow(rel_missing_comp)))
    }
    zeb1_wt_bad <- dt[dataset_folder == "GSE279446" & grepl("zeb1_wt|wt_ovary|wt ovary", tolower(paste(sample_key, sample_title, treatment, genotype, sep = " "))) & genotype_class == "heterozygous"]
    add_check(nm, "GSE279446_Zeb1_WT_not_heterozygous", nrow(zeb1_wt_bad) == 0, paste0("rows=", nrow(zeb1_wt_bad)))
    gse279446_ids <- unique(dt[dataset_folder == "GSE279446" & cohort %in% c("adult", "young"), .(cohort, comparison_id)])
    adult_ids <- unique(gse279446_ids[cohort == "adult"]$comparison_id)
    young_ids <- unique(gse279446_ids[cohort == "young"]$comparison_id)
    add_check(nm, "GSE279446_adult_young_comparison_ids_separate", length(intersect(adult_ids, young_ids)) == 0, paste(intersect(adult_ids, young_ids), collapse = ";"))
    gse228_bad <- dt[dataset_folder == "GSE228975" & plot_include_primary %in% TRUE]
    add_check(nm, "GSE228975_plot_include_primary_false", nrow(gse228_bad) == 0, paste0("rows=", nrow(gse228_bad)))
  }
  rbindlist(checks, use.names = TRUE)
}

write_summary <- function(pred_all, validation, warnings) {
  class_counts <- pred_all[, .(n_samples = uniqueN(sample_key)), by = projection_class][order(projection_class)]
  role_counts <- pred_all[, .(n_samples = uniqueN(sample_key)), by = projection_role][order(projection_role)]
  lines <- c(
    "# Simplified external projection classes",
    "",
    "## Why simplify projection_class",
    "",
    "The original projection_class mixed broad biology, dataset-specific control names, genotype details, and OOD labels. The simplified files preserve the old values in projection_class_original and reserve projection_class for broad categories.",
    "",
    "## Final projection_class labels",
    "",
    paste0("- ", allowed_projection_classes),
    "",
    "## Final projection_role labels",
    "",
    paste0("- ", allowed_projection_roles),
    "",
    "## How to use projection_group and comparison_id",
    "",
    "- projection_group is the readable within-dataset plotting label.",
    "- comparison_id identifies the matched-control comparison used for relative delta-age plots.",
    "",
    "## Dataset-specific recoding summary",
    "",
    "- GSE194196: Bnc1tr/tr is genetic_model; WT-1/WT-4/WT-7 are reference_control.",
    "- GSE252791: Rnf111 heterozygous KO is genetic_model; WT controls are reference_control.",
    "- GSE279446: Zeb1_WT is WT reference_control; Zeb1_200H is heterozygous; Zeb1_200M is homozygous; adult and young cohorts keep separate comparison_id values.",
    "- GSE279573: WT chow fertile is reference_control; KO chow infertile is genetic_model; KO vegetable-oil-free rescue is diet.",
    "- GSE240547: control diet is reference_control; HFHS is diet.",
    "- PRJNA1076509: young female donor FMT is reference_control; estropausal donor FMT is microbiome.",
    "- PRJNA863443: natural ageing, CTL, and Foxl2+/+ are reference_control; VCD is intervention; Foxl2+/- is genetic_model.",
    "- GSE106339 GOAT-KO samples are genetic_model and WT samples are reference_control and enter the 288-sample analysis-domain set; their prepubertal contrast remains descriptive.",
    "- GSE228975 and other developmental rows are developmental_OOD/OOD_supplement.",
    "",
    "## Replacement plots",
    "",
    "- Use 01_projection_cohort_flow_simplified instead of detailed projection-class flow plots.",
    "- Use 02_raw_delta_by_projection_class_simplified for broad-category raw delta-age.",
    "- Use 03_raw_delta_by_dataset_group for dataset-specific absolute placement.",
    "- Use 04_relative_delta_by_comparison, 05_genetic_model_relative_delta, and 06_diet_microbiome_relative_delta for matched-control effect estimates.",
    "",
    "## Caution",
    "",
    "Raw delta-age plots show absolute placement in the clock; relative delta-age plots show matched-control effect estimates.",
    "",
    "## Sample counts",
    "",
    "Projection classes:",
    paste0("- ", class_counts$projection_class, ": ", class_counts$n_samples, " unique samples"),
    "",
    "Projection roles:",
    paste0("- ", role_counts$projection_role, ": ", role_counts$n_samples, " unique samples"),
    "",
    "## Validation",
    "",
    paste0("- Checks: ", nrow(validation)),
    paste0("- Failures: ", sum(validation$passed != TRUE)),
    paste0("- Warnings: ", length(warnings))
  )
  writeLines(lines, file.path(out_dir, "external_projection_simplified_summary.md"))
}

pred_all <- read_required("external_projection_predictions_all.tsv")
pred_primary <- read_required("external_projection_predictions.tsv")
rel_all <- read_required("external_projection_relative_to_controls_all.tsv")
rel_primary <- read_required("external_projection_relative_to_controls.tsv")
availability <- read_required("external_projection_dataset_control_availability.tsv")
control_audit <- read_required("external_projection_control_matching_audit.tsv")

original_projection_class_levels <- sort(unique(pred_all$projection_class))

pred_all_s <- add_simplified_classes(pred_all)
pred_primary_s <- pred_all_s[plot_include_primary %in% TRUE]
rel_all_s <- add_simplified_classes(rel_all)
rel_primary_s <- rel_all_s[plot_include_primary %in% TRUE]

fwrite(pred_all_s, file.path(out_dir, "external_projection_predictions_all_simplified.tsv"), sep = "\t")
fwrite(pred_primary_s, file.path(out_dir, "external_projection_predictions_simplified.tsv"), sep = "\t")
fwrite(rel_all_s, file.path(out_dir, "external_projection_relative_to_controls_all_simplified.tsv"), sep = "\t")
fwrite(rel_primary_s, file.path(out_dir, "external_projection_relative_to_controls_simplified.tsv"), sep = "\t")

dictionary <- unique(rbindlist(list(
  pred_all_s[, .(dataset_folder, projection_class_original, projection_class, projection_role, projection_group, comparison_id, genotype_class, cohort, control_group_label)],
  rel_all_s[, .(dataset_folder, projection_class_original, projection_class, projection_role, projection_group, comparison_id, genotype_class, cohort, control_group_label)]
), use.names = TRUE, fill = TRUE))
setorder(dictionary, dataset_folder, projection_class, projection_role, projection_group)
fwrite(dictionary, file.path(out_dir, "external_projection_class_recode_dictionary.tsv"), sep = "\t")

group_summary <- pred_all_s[, .(
  n_samples = uniqueN(sample_key),
  n_model_rows = .N,
  mean_age_week = mean(age_week, na.rm = TRUE),
  mean_predicted_age_week = mean(predicted_age_week, na.rm = TRUE),
  mean_delta_age_week = mean(delta_age_week, na.rm = TRUE),
  median_delta_age_week = median(delta_age_week, na.rm = TRUE)
), by = .(dataset_folder, projection_class, projection_role, projection_group, comparison_id, model_family, selector)]

rel_summary <- rel_all_s[, .(
  n_relative_rows = .N,
  n_relative_samples = uniqueN(sample_key),
  mean_relative_delta_age_week = mean(relative_delta_age_week, na.rm = TRUE),
  median_relative_delta_age_week = median(relative_delta_age_week, na.rm = TRUE)
), by = .(dataset_folder, projection_class, projection_role, projection_group, comparison_id, model_family, selector)]

group_summary <- merge(
  group_summary,
  rel_summary,
  by = c("dataset_folder", "projection_class", "projection_role", "projection_group", "comparison_id", "model_family", "selector"),
  all.x = TRUE
)
setorder(group_summary, dataset_folder, comparison_id, projection_group, model_family)
fwrite(group_summary, file.path(out_dir, "external_projection_simplified_group_summary.tsv"), sep = "\t")

validation <- validate_simplified(list(
  predictions_all = pred_all_s,
  predictions = pred_primary_s,
  relative_to_controls_all = rel_all_s,
  relative_to_controls = rel_primary_s
))
fwrite(validation, file.path(out_dir, "external_projection_simplified_validation_audit.tsv"), sep = "\t")

sample_group_counts <- unique(pred_all_s[, .(sample_key, dataset_folder, projection_group, projection_class, projection_role)])
small_groups <- sample_group_counts[, .(n_samples = uniqueN(sample_key)), by = .(dataset_folder, projection_group, projection_class, projection_role)][n_samples < 2]
warnings <- character()
if (nrow(small_groups) > 0) {
  warnings <- paste0(
    "Projection group with fewer than 2 samples: ",
    small_groups$dataset_folder, " / ", small_groups$projection_group,
    " (", small_groups$projection_class, ", ", small_groups$projection_role, "), n=", small_groups$n_samples
  )
}
if (length(warnings) == 0) warnings <- "No simplified-class warnings."
writeLines(warnings, file.path(out_dir, "external_projection_simplified_warnings.txt"))

write_summary(pred_all_s, validation, if (identical(warnings, "No simplified-class warnings.")) character() else warnings)

flow_dt <- unique(pred_all_s[, .(sample_key, dataset_folder, projection_class, projection_role, plot_include_primary, plot_include_ood_supplement)])
flow_dt[, panel := fifelse(plot_include_primary %in% TRUE, "Primary", "OOD/supplement")]
p <- ggplot(flow_dt, aes(projection_class, fill = projection_role)) +
  geom_bar() +
  facet_wrap(~panel, scales = "free_y") +
  labs(x = "Simplified projection class", y = "Unique samples", fill = "Projection role") +
  theme_projection()
save_plot(p, "01_projection_cohort_flow_simplified", width = 12, height = 7)

raw_primary <- pred_all_s[plot_include_primary %in% TRUE]
p <- ggplot(raw_primary, aes(projection_class, delta_age_week, color = projection_role, shape = projection_role)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_boxplot(aes(group = projection_class), outlier.shape = NA, alpha = 0.12, color = "grey45") +
  geom_point(position = position_jitter(width = 0.15, height = 0), alpha = 0.8, size = 2) +
  facet_wrap(~model_family) +
  labs(x = "Simplified projection class", y = "Raw delta age (weeks)", color = "Projection role", shape = "Projection role") +
  theme_projection()
save_plot(p, "02_raw_delta_by_projection_class_simplified", width = 12, height = 7)

for (ds in sort(unique(raw_primary$dataset_folder))) {
  dt <- raw_primary[dataset_folder == ds]
  if (nrow(dt) == 0) next
  p <- ggplot(dt, aes(projection_group, delta_age_week, color = projection_class, shape = projection_role)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_boxplot(aes(group = projection_group), outlier.shape = NA, alpha = 0.10, color = "grey55") +
    geom_point(position = position_jitter(width = 0.15, height = 0), alpha = 0.85, size = 2) +
    facet_wrap(~model_family, ncol = 1) +
    labs(x = "Projection group", y = "Raw delta age (weeks)", color = "Projection class", shape = "Projection role") +
    theme_projection()
  save_plot(p, paste0("03_raw_delta_by_dataset_group_", clean_filename(ds)), width = max(8, min(15, 1.4 * uniqueN(dt$projection_group) + 5)), height = 8)
}

rel_plot <- rel_all_s[is.finite(suppressWarnings(as.numeric(relative_delta_age_week))) & plot_include_primary %in% TRUE]
if (nrow(rel_plot) > 0) {
  p <- ggplot(rel_plot, aes(projection_group, relative_delta_age_week, color = projection_class, shape = projection_role)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_boxplot(aes(group = projection_group), outlier.shape = NA, alpha = 0.10, color = "grey55") +
    geom_point(position = position_jitter(width = 0.15, height = 0), alpha = 0.85, size = 1.9) +
    facet_grid(model_family ~ comparison_id, scales = "free_x", space = "free_x") +
    labs(x = "Projection group", y = "Matched-control relative delta age (weeks)", color = "Projection class", shape = "Projection role") +
    theme_projection() +
    theme(strip.text.x = element_text(size = 9), axis.text.x = element_text(size = 9, angle = 35, hjust = 1))
  save_plot(p, "04_relative_delta_by_comparison", width = 18, height = 9)
}

genetic_rel <- rel_plot[comparison_id %in% c(
  "GSE194196_Bnc1_vs_WT",
  "GSE252791_Rnf111_vs_WT",
  "GSE279446_adult_Zeb1_vs_WT",
  "GSE279446_young_Zeb1_vs_WT",
  "GSE279573_KO_or_rescue_vs_WT",
  "PRJNA863443_Foxl2_het_vs_WT_young",
  "PRJNA863443_Foxl2_het_vs_WT_mid_age"
) & projection_class %in% c("reference_control", "genetic_model", "diet")]
if (nrow(genetic_rel) > 0) {
  p <- ggplot(genetic_rel, aes(projection_group, relative_delta_age_week, color = projection_class, shape = projection_role)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_boxplot(aes(group = projection_group), outlier.shape = NA, alpha = 0.10, color = "grey55") +
    geom_point(position = position_jitter(width = 0.15, height = 0), alpha = 0.85, size = 2) +
    facet_grid(model_family ~ comparison_id, scales = "free_x", space = "free_x") +
    labs(x = "Projection group", y = "Matched-control relative delta age (weeks)", color = "Projection class", shape = "Projection role") +
    theme_projection() +
    theme(strip.text.x = element_text(size = 9), axis.text.x = element_text(size = 9, angle = 35, hjust = 1))
  save_plot(p, "05_genetic_model_relative_delta", width = 16, height = 8)
}

diet_micro_rel <- rel_plot[comparison_id %in% c(
  "GSE240547_HFHS_vs_control_diet",
  "GSE279573_KO_or_rescue_vs_WT",
  "PRJNA1076509_estropausal_vs_young_FMT",
  "PRJNA863443_VCD_vs_CTL_3m_30d",
  "PRJNA863443_VCD_vs_CTL_3m_90d",
  "PRJNA863443_VCD_vs_CTL_10m_30d",
  "PRJNA863443_VCD_vs_CTL_10m_90d"
) & projection_class %in% c("reference_control", "diet", "microbiome", "genetic_model", "intervention")]
if (nrow(diet_micro_rel) > 0) {
  p <- ggplot(diet_micro_rel, aes(projection_group, relative_delta_age_week, color = projection_class, shape = projection_role)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_boxplot(aes(group = projection_group), outlier.shape = NA, alpha = 0.10, color = "grey55") +
    geom_point(position = position_jitter(width = 0.15, height = 0), alpha = 0.85, size = 2) +
    facet_grid(model_family ~ comparison_id, scales = "free_x", space = "free_x") +
    labs(x = "Projection group", y = "Matched-control relative delta age (weeks)", color = "Projection class", shape = "Projection role") +
    theme_projection() +
    theme(strip.text.x = element_text(size = 9), axis.text.x = element_text(size = 9, angle = 35, hjust = 1))
  save_plot(p, "06_diet_microbiome_relative_delta", width = 14, height = 8)
}

heat <- unique(pred_all_s[, .(sample_key, projection_group, projection_class, projection_role)])
heat <- heat[, .(n_samples = uniqueN(sample_key)), by = .(projection_group, projection_class, projection_role)]
heat[, class_role := paste(projection_class, projection_role, sep = " / ")]
p <- ggplot(heat, aes(class_role, projection_group, fill = n_samples)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_samples), size = 3) +
  scale_fill_gradient(low = "#EEF5F3", high = "#287C6F") +
  labs(x = "Projection class / role", y = "Projection group", fill = "Unique samples") +
  theme_projection() +
  theme(axis.text.x = element_text(size = 9, angle = 35, hjust = 1))
save_plot(p, "07_category_group_audit_heatmap", width = 12, height = max(7, 0.32 * uniqueN(heat$projection_group) + 3))

annotation_cols <- c(
  "sample_key", "model_family", "selector", "dataset_folder", "projection_class_original",
  "projection_class", "projection_role", "projection_group", "comparison_id",
  "genotype_class", "cohort", "control_group_label", "plot_include_primary",
  "plot_include_ood_supplement"
)
annotation_lookup <- unique(pred_all_s[, ..annotation_cols])

availability_summary <- unique(pred_all_s[, .(
  projection_class = paste(sort(unique(projection_class)), collapse = ";"),
  projection_role = paste(sort(unique(projection_role)), collapse = ";"),
  projection_group = paste(sort(unique(projection_group)), collapse = ";"),
  comparison_id = paste(sort(unique(na.omit(comparison_id))), collapse = ";"),
  plot_include_primary = any(plot_include_primary %in% TRUE),
  plot_include_ood_supplement = any(plot_include_ood_supplement %in% TRUE)
), by = dataset_folder])
availability_s <- merge(availability, availability_summary, by = "dataset_folder", all.x = TRUE)

control_join_cols <- intersect(c("sample_key", "model_family"), names(control_audit))
control_anno_cols <- setdiff(annotation_cols, c("selector", setdiff(c("sample_key", "model_family"), control_join_cols)))
control_audit_s <- merge(
  control_audit,
  unique(annotation_lookup[, ..control_anno_cols]),
  by = control_join_cols,
  all.x = TRUE
)
fwrite(availability_s, file.path(out_dir, "external_projection_dataset_control_availability_simplified.tsv"), sep = "\t")
fwrite(control_audit_s, file.path(out_dir, "external_projection_control_matching_audit_simplified.tsv"), sep = "\t")

pca_scores_path <- file.path(projection_dir, "pca_projection_scores.tsv")
if (file.exists(pca_scores_path)) {
  pca_scores <- fread(pca_scores_path, sep = "\t", na.strings = c("", "NA", "NaN"))
  pca_scores_s <- merge(pca_scores, annotation_lookup, by = intersect(c("sample_key", "model_family", "selector"), names(pca_scores)), all.x = TRUE)
  fwrite(pca_scores_s, file.path(out_dir, "pca_projection_scores_simplified.tsv"), sep = "\t")
}

pca_contrib_path <- file.path(projection_dir, "pca_projection_pc_contributions.tsv")
if (file.exists(pca_contrib_path)) {
  pca_contrib <- fread(pca_contrib_path, sep = "\t", na.strings = c("", "NA", "NaN"))
  pca_contrib_s <- merge(pca_contrib, annotation_lookup, by = intersect(c("sample_key", "model_family", "selector"), names(pca_contrib)), all.x = TRUE)
  fwrite(pca_contrib_s, file.path(out_dir, "pca_projection_pc_contributions_simplified.tsv"), sep = "\t")
}

failures <- validation[passed != TRUE]
cat("Original projection_class levels:", length(original_projection_class_levels), "\n")
cat("Simplified projection_class levels:", uniqueN(pred_all_s$projection_class), "\n")
cat("Projection_role levels:", uniqueN(pred_all_s$projection_role), "\n")
cat("Output directory:", out_dir, "\n")
if (nrow(failures) == 0) {
  cat("Validation: PASS\n")
} else {
  cat("Validation: FAIL\n")
  print(failures)
  stop("Simplified external projection validation failed.", call. = FALSE)
}
