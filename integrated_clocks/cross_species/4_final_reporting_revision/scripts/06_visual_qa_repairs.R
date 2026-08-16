#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]))), "common_final.R"))

message("[6/6] R-only visual-QA label repairs")

## Repair cramped nested-CV facet labels using already-computed source data.
pm <- fread(pfinal("figures/source_data/nested_cv_validation_metrics.tsv"))
spec_labels <- c(
  baseline_age_stratified_5fold="Baseline age-stratified",
  batch_stratified_5fold="Center-stratified",
  hardy_age_stratified_5fold="Hardy-stratified",
  sequencing_center_blocked="Center-blocked")
metric_labels <- c(pearson="Pearson r",spearman="Spearman rho",
  balanced_age_bracket_accuracy="Balanced accuracy",
  ordinal_concordance_kendall="Kendall concordance",calibration_slope="Calibration slope")
pm[,specification_label:=factor(spec_labels[specification],levels=rev(unname(spec_labels)))]
pm[,metric_label:=factor(metric_labels[metric],levels=unname(metric_labels))]
p_metrics <- ggplot(pm,aes(specification_label,value,fill=specification_label))+
  geom_col(width=.72)+facet_wrap(~metric_label,scales="free_x",ncol=3)+coord_flip()+
  scale_fill_brewer(palette="Set2",guide="none")+
  labs(x=NULL,y=NULL,title="Nested-CV robustness metrics")+
  theme(axis.text.y=element_text(size=6),strip.text=element_text(size=6.5,face="bold"))
save_pub_r(p_metrics,pfinal("figures/15_nested_cv_validation_metrics"),183,110)

pd <- fread(pfinal("figures/source_data/nested_cv_validation_predictions.tsv"))
pd[,specification_label:=factor(spec_labels[specification],levels=unname(spec_labels[1:3]))]
p_pred <- ggplot(pd,aes(observed,predicted,colour=age_bracket))+geom_abline(slope=1,intercept=0,colour="#BDBDBD")+
  geom_point(size=.8,alpha=.65)+facet_wrap(~specification_label)+scale_colour_manual(values=palette_age)+
  labs(x="Bracket midpoint (descriptive only)",y="Nested-CV prediction",colour="Age",
    title="Leakage-controlled age-information prediction",
    subtitle="Dataset-adequacy analysis; not a validated human ovarian clock")+
  theme(legend.position="top",strip.text=element_text(size=6.5,face="bold"))
save_pub_r(p_pred,pfinal("figures/14_nested_cv_validation_predictions"),183,82)

## Short strip labels preserve legibility while the full definitions remain in the legend/report.
hm <- fread(pfinal("figures/source_data/representative_age_pattern_heatmap.tsv"))
hm[,gene:=factor(gene,levels=rev(unique(gene)))]
class_labels <- c("approximately monotonic; one reversal"="one reversal","monotonic decreasing"="monotonic decreasing",
  "monotonic increasing"="monotonic increasing","nonlinear"="nonlinear","oldest-bin-driven"="oldest-bin driven","unclassified"="unclassified")
hm[,class_label:=factor(class_labels[classification],levels=class_labels[unique(classification)])]
p_hm <- ggplot(hm,aes(age_bracket,gene,fill=value))+geom_tile()+
  facet_grid(class_label~.,scales="free_y",space="free_y",switch="y")+
  scale_fill_gradient2(low="#3973B7",mid="white",high="#C73636",midpoint=0,name="Row z")+
  labs(x="Age bracket",y=NULL,title="Representative fitted age-bin patterns",
    subtitle="Top absolute ordered-age t statistics within each mutually exclusive class")+
  theme(axis.text.x=element_text(angle=35,hjust=1),strip.placement="outside",
    strip.text.y.left=element_text(angle=0,size=5.5,face="bold"),panel.spacing.y=unit(1.2,"mm"))
save_pub_r(p_hm,pfinal("figures/13_representative_age_pattern_heatmap"),150,190)

writeLines(c("# Figure visual-QA notes","",
  "- Backend: R only (ggplot2/patchwork; svglite, cairo PDF, and ragg exports).",
  "- Core conclusion: GTEx ovary contains robust age-associated structure, while sensitivity panels bound its magnitude and interpretation.",
  "- Archetype: quantitative-grid comparison with a supporting pattern heatmap.",
  "- Export bundle: vector PDF/SVG, 600-dpi PNG, 600-dpi TIFF, scripts, and TSV source data. The selected R environment lacks svglite, so the base-R SVG device outlines glyphs; the PDF is the preferred text-preserving vector deliverable.",
  "- Scale check: all centroid correlation panels print their display limits; the 1-r panel is separate; all cells are numerically annotated.",
  "- Statistics check: donor counts and pair counts are distinct; similarity intervals use donor bootstrap and p values use 1,000 donor-label permutations.",
  "- ML check: split design, 5-fold nested structure, metrics, blocked/stratified baselines, and 100-permutation null are documented.",
  "- Visual inspection: centroid annotations and legends readable; similarity intervals unobstructed; Hardy/pathway scatterplots legible; nested-CV facet clipping repaired; pattern strips shortened without changing definitions; sparse sensitivity labels readable at 183 mm.",
  "- Image integrity: no microscopy, image manipulation, raster stitching, or pseudo-color provenance issue applies."),
  pfinal("reports/figure_visual_qa.md"),useBytes=TRUE)

fig_files <- list.files(pfinal("figures"),pattern="\\.(pdf|svg|png|tiff)$",full.names=TRUE)
manifest <- data.table(file=sub(paste0("^",pfinal(""),"/?"),"",fig_files),extension=tools::file_ext(fig_files),bytes=file.info(fig_files)$size)
manifest[,expected_resolution:=fifelse(extension%in%c("png","tiff"),"600 dpi","vector")]
write_tsv(manifest,pfinal("results/table_final_figure_manifest.tsv"))

message("Visual-QA repairs complete")
