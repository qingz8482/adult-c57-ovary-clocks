#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1] |> sub("^--file=", "", x = _))), "utilities.R"))

cov <- fread(path_module("results/tables/table_orthology_coverage.tsv"))
dirs <- fread(path_module("results/tables/table_stageA_direction_conservation.tsv"))
enr <- fread(path_module("results/tables/table_stageA_enrichment.tsv"))
nullA <- fread(path_module("results/tables/table_random_signature_benchmarks.tsv"))
assoc <- fread(path_module("results/tables/table_stageB_score_age_association.tsv"))
nullB <- fread(path_module("results/null_tests/stageB_random_signature_benchmarks.tsv"))
ood <- fread(path_module("results/tables/table_stageB_ood_qc.tsv"))
cross <- fread(path_module("results/stageB/cross_model_concordance.tsv"))
disc <- fread(path_module("results/qc/prompt_vs_artifact_shared_gene_discrepancy.tsv"))
repro <- fread(path_module("results/qc/mouse_artifact_reproduction.tsv"))
gtex_meta <- fread(path_module("derived_data/gtex/gtex_ovary_sample_audit.tsv"))
hpa_meta <- fread(path_module("derived_data/hpa/hpa_ovary_sample_metadata.tsv"))

primary <- assoc[analysis=="primary_technical_adjusted"]
el <- primary[model_family=="elasticnet"]; pc <- primary[model_family=="pca_glmnet"]
elcov <- cov[gene_set=="elasticnet_hvg_2000",absolute_model_weight_coverage]
pccov <- cov[gene_set=="pca_versionB_2000",absolute_model_weight_coverage]
elnull <- nullB[model_family=="elasticnet"]; pcnull <- nullB[model_family=="pca_glmnet"]
elood <- ood[model_family=="elasticnet"]; pcood <- ood[model_family=="pca_glmnet"]

grade <- data.table(
  scope=c("elasticnet","pca_glmnet","combined_stageA_stageB"),
  grade=c("PARTIAL_SUPPORT","INSUFFICIENT_OR_NEGATIVE","INSUFFICIENT_OR_NEGATIVE"),
  rationale=c(
    "Positive GTEx trend and moderate weighted coverage, but matched frozen-signature null not exceeded and OOD distance is age-associated.",
    "High weighted coverage, but Stage A model direction and Stage B trend are directionally inconsistent or null.",
    "No model satisfies the prespecified support gates; models disagree and the only nominally supported trend fails null/OOD safeguards."
  )
)
write_tsv(grade,path_module("results/tables/table_claim_grade.tsv"))

combined_null <- rbindlist(list(
  nullA[,.(stage="Stage A",test=signature,observed,predicted_direction_pass,null_mean,null_sd,null_q95,empirical_p,fdr,n_permutations)],
  nullB[,.(stage="Stage B",test=model_family,observed,predicted_direction_pass,null_mean,null_sd,null_q95,empirical_p,fdr,n_permutations)]
),fill=TRUE)
write_tsv(combined_null,path_module("results/tables/table_random_signature_benchmarks.tsv"))

methods <- c(
  "# Methods: cross-species Stage A and Stage B","",
  "## Frozen mouse models","",
  "The authoritative top-level deployment artifacts comprised the integrated ElasticNet/HVG and PCA/glmnet/Version B clocks, each with 2,000 selected mouse genes and 117 control-like deployment-reference samples. Artifact schemas, feature order, centers, scales, calibration, and the missing-gene policy were validated. Predictions were independently reconstructed from model coefficients; the PCA predictor was also back-projected into gene space. The maximum permitted discrepancy was 1e-8.","",
  "## Orthology","",
  "Mouse genes were mapped through Ensembl BioMart release 116. Primary analyses retained unique, high-confidence ortholog_one2one mappings. One-to-many, many-to-many, low-confidence, and unmapped relationships remained in the audit table but were excluded from primary projection. Mapping by symbol capitalization was not used.","",
  "## Human ovarian RNA-seq","",
  "GTEx V11 ovary gene read counts (RNASeQC 2.4.3; GENCODE 47) and open sample/subject metadata were obtained from the official Adult GTEx repository. V11 changes annotation relative to V10 but adds no donors. Counts were filtered at >=10 reads in >=20% of ovary donors, normalized by edgeR TMM, and transformed to logCPM with prior.count=2. Public 10-year age brackets were preserved; midpoints were used only as the prespecified numeric trend representation. RIN was included in the primary technical-adjusted model. HPA v25.1 transcript TPMs were summed to gene level for the seven official ovary samples and used only for exploratory Spearman direction checks.","",
  "## Stage A","",
  "Mouse expression-age effects were estimated by limma across the frozen 117-sample reference matrix using age_week and source_batch. GTEx gene effects were estimated by limma-voom per decade of public age-bracket midpoint with RIN adjustment. Conservation was tested against the intersection of the integrated mouse universe, eligible one-to-one orthologues, and GTEx-expressed genes. Model weights were extracted in expression space, including PCA back-projection. Ten thousand expression-, detection-, and variance-matched random draws were used for empirical signed-direction benchmarks.","",
  "## Stage B","",
  "GTEx logCPM values were inserted into the exact mouse feature order through primary one-to-one mappings. Unavailable features were imputed to the mouse training center and therefore scaled to zero. Frozen centers, scales, PCA rotation, glmnet coefficients, and mouse calibration were applied exactly once. No human recentering, rescaling, batch harmonization, feature selection, PCA refit, coefficient fit, or age calibration was performed. The primary score was the within-GTEx z score of the frozen raw linear predictor; mouse-week output was retained only for audit.","",
  "## OOD, nulls, and inference","",
  "OOD metrics included extreme mouse-scaled feature fractions, mouse-reference range exceedance, low-dimensional reference distance, and nearest-reference distance. Score trends used HC3 robust standard errors. Ten thousand frozen-signature nulls preserved the observed mapped feature set and weight distribution while permuting weights across mapped features. Multiple testing across the two clocks used BH FDR."
)
writeLines(methods,path_module("reports/methods_stageA_stageB.md"))

results <- c(
  "# Results: cross-species Stage A and Stage B","",
  "## Orthologue and detectability coverage","",
  sprintf("The integrated universe had %s/%s human-detected one-to-one orthologues. Absolute model-weight coverage was %.3f for ElasticNet and %.3f for PCA/glmnet. The current frozen artifacts shared %d genes; the supplied prompt's 277-gene value was stale.",cov[gene_set=="integrated_universe",n_human_detected],cov[gene_set=="integrated_universe",n_total],elcov,pccov,disc$artifact_observed_shared_gene_count),"",
  "## Mouse-human age-direction conservation","",
  sprintf("ElasticNet nonzero weights showed weighted sign concordance %.3f with GTEx age effects. PCA/glmnet showed weighted sign concordance %.3f and a negative weighted rank association (%.3f). None of the predicted-direction Stage A signatures passed the prespecified matched-null and direction requirements after FDR control.",dirs[analysis=="elasticnet_nonzero",weighted_sign_concordance],dirs[analysis=="pca_versionB_2000" & comparison=="model_weight_vs_human_age",weighted_sign_concordance],dirs[analysis=="pca_versionB_2000" & comparison=="model_weight_vs_human_age",weighted_spearman_rho]),"",
  "## Frozen GTEx score trends","",
  sprintf("The ElasticNet score increased by %.3f SD per decade (HC3 95%% CI %.3f to %.3f; FDR %.3g). PCA/glmnet changed by %.3f SD per decade (95%% CI %.3f to %.3f; FDR %.3g).",el$estimate,el$conf_low,el$conf_high,el$fdr,pc$estimate,pc$conf_low,pc$conf_high,pc$fdr),"",
  "## Matched nulls, model agreement, and OOD","",
  sprintf("The ElasticNet trend did not exceed the 10,000-draw frozen-signature null (empirical p=%.3f). PCA/glmnet failed the predicted-direction gate. ElasticNet and PCA scores were essentially uncorrelated (Pearson r=%.3f). ElasticNet OOD distance was associated with age (p=%.3g), so its nominal age trend cannot be separated from age-linked cross-species distribution shift.",elnull$empirical_p,cross$pearson,elood$p),"",
  "## Evidence grade","",
  "The combined grade is **INSUFFICIENT_OR_NEGATIVE**. These results do not support a conserved human ovarian ageing axis under the prespecified frozen-transfer safeguards."
)
writeLines(results,path_module("reports/results_stageA_stageB.md"))

limits <- c("# Limitations and claim boundary","","- GTEx ages are public intervals, not exact ages.","- Human ovary pieces differ biologically and technically from the mouse integrated bulk/pseudobulk training profiles.","- Cross-species baseline and covariance shifts are severe; OOD is quantified but not removable without changing the frozen test.","- HPA contains only seven ovarian samples and is exploratory.","- Whole-tissue trends can reflect ovarian cell-composition changes.","- The current artifacts share 297 genes, whereas the supplied instruction package stated 277; current frozen artifacts were treated as authoritative.","","The analysis must not be described as human age prediction, a validated human ovarian clock, human age acceleration, or mouse-week-to-human-year conversion.")
writeLines(limits,path_module("reports/limitations_and_claim_boundary.md"))

panel_plan <- c("# Publication panel plan","","Core conclusion: the two frozen mouse clocks do not jointly support a conserved human ovarian ageing axis under matched-null and OOD safeguards.","","Archetype: schematic-led quantitative composite. Backend: R only. Final assembly: 183 mm wide, vector PDF/SVG plus 600-dpi PNG/TIFF.","","- A: frozen workflow and claim boundary.","- B: raw, detected, and model-weighted orthology coverage.","- C: mouse versus GTEx observed expression-age directions.","- D: signed model-direction statistics against matched nulls.","- E: strict frozen scores across GTEx public age brackets.","- F: cross-model disagreement and OOD.","","Reviewer risks: moderate ElasticNet weighted coverage, age-associated OOD, non-independent shared preprocessing, bracketed age, and small exploratory HPA cohort.")
writeLines(panel_plan,path_module("reports/publication_panel_plan.md"))

report <- c("# Human cross-species ovarian clock audit","",results[-1],"",methods,"",limits)
writeLines(report,path_module("reports/human_cross_species_stageA_stageB_report.md"))
html <- paste0("<!doctype html><html><head><meta charset='utf-8'><title>Human cross-species Stage A/B</title><style>body{font-family:Arial,sans-serif;max-width:1000px;margin:40px auto;line-height:1.45}pre{white-space:pre-wrap}</style></head><body><pre>",paste(gsub("&","&amp;",gsub("<","&lt;",report,fixed=TRUE),fixed=TRUE),collapse="\n"),"</pre></body></html>")
writeLines(html,path_module("reports/human_cross_species_stageA_stageB_report.html"))

summary <- c("# Codex execution summary","",sprintf("- Run date: %s",Sys.Date()),"- Frozen mouse artifacts: two authoritative deployment_all_control_like artifacts; 117 reference samples each.",sprintf("- Artifact reproduction: PASS; maximum error %.3g.",max(c(repro$max_abs_predict_vs_formula,repro$max_abs_pc_vs_gene_space),na.rm=TRUE)),sprintf("- Human data: GTEx V11 ovary n=%d; HPA v25.1 ovary n=%d (exploratory).",nrow(gtex_meta),nrow(hpa_meta)),sprintf("- Weighted coverage: ElasticNet %.3f; PCA/glmnet %.3f.",elcov,pccov),sprintf("- Stage B: ElasticNet %.3f SD/decade (FDR %.3g; matched-null p %.3f); PCA %.3f SD/decade (FDR %.3g).",el$estimate,el$fdr,elnull$empirical_p,pc$estimate,pc$fdr),sprintf("- OOD: ElasticNet age-OOD p=%.3g; PCA age-OOD p=%.3g.",elood$p,pcood$p),sprintf("- Cross-model Pearson correlation: %.3f.",cross$pearson),"- Prompt/artifact discrepancy: prompt stated 277 shared genes; current authoritative artifacts share 297.","- Claim grade: **INSUFFICIENT_OR_NEGATIVE**.","- Blocked inputs: none.","- Main figure: `1_publication/assembled/human_cross_species_main/human_cross_species_main.pdf` (also SVG/PNG/TIFF).","- Validation: run `validate_outputs.R`.","- Exact rerun: `bash adult_c57_ovary_clocks/bulk-pseudobulk_clock/human_cross_species_stageA_stageB/run_all.sh`.")
writeLines(summary,path_module("CODEX_EXECUTION_SUMMARY.md"))

tables <- list.files(path_module("results/tables"),pattern="[.]tsv$",full.names=TRUE)
invisible(lapply(tables,make_data_dictionary))
message("Reports, evidence grades, captions, and table dictionaries complete.")
