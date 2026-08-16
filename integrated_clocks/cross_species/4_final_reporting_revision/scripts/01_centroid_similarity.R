#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]))), "common_final.R"))

message("[1/5] Centroid heatmaps and donor-level similarity inference")
x <- load_inputs(); audit <- x$audit
cent <- fread(pbase("results/table_age_bin_centroid_correlations_distances.tsv"))
boot_cent <- fread(pbase("results/table_age_bin_centroid_bootstrap_draws.tsv"))

metrics <- c("raw_pearson", "raw_spearman", "residualized_pearson", "residualized_spearman", "residualized_one_minus_pearson")
metric_figure_stems <- c(raw_pearson="02_centroid_heatmap_raw_pearson",raw_spearman="03_centroid_heatmap_raw_spearman",
  residualized_pearson="04_centroid_heatmap_residualized_pearson",residualized_spearman="05_centroid_heatmap_residualized_spearman",
  residualized_one_minus_pearson="06_centroid_heatmap_residualized_one_minus_pearson")
metric_titles <- c(raw_pearson="Raw Pearson centroid correlation", raw_spearman="Raw Spearman centroid correlation",
  residualized_pearson="Residualized Pearson centroid correlation", residualized_spearman="Residualized Spearman centroid correlation",
  residualized_one_minus_pearson="1 - residualized Pearson correlation")

make_matrix_source <- function(metric) {
  metric_value <- metric
  d <- cent[metric == metric_value, .(age_bin_a, age_bin_b, observed, bootstrap_ci_low, bootstrap_ci_high, bootstrap_replicates, n_a, n_b)]
  revd <- copy(d); setnames(revd, c("age_bin_a","age_bin_b"), c("age_bin_b","age_bin_a"))
  diagval <- if (grepl("one_minus", metric)) 0 else 1
  dg <- data.table(age_bin_a=age_levels, age_bin_b=age_levels, observed=diagval,
    bootstrap_ci_low=diagval, bootstrap_ci_high=diagval, bootstrap_replicates=NA_integer_,
    n_a=audit[, .N, by=age_bracket]$N[match(age_levels, audit[, .N, by=age_bracket]$age_bracket)],
    n_b=audit[, .N, by=age_bracket]$N[match(age_levels, audit[, .N, by=age_bracket]$age_bracket)])
  out <- rbindlist(list(d, revd, dg), fill=TRUE)
  out[, `:=`(metric=metric, age_bin_a=factor(age_bin_a, levels=age_levels), age_bin_b=factor(age_bin_b, levels=rev(age_levels)))]
  out
}

heat_src <- rbindlist(lapply(metrics, make_matrix_source))
write_tsv(heat_src, pfinal("figures/source_data/centroid_correlation_heatmaps.tsv"))

plots <- list()
for (m in metrics) {
  d <- heat_src[metric == m]
  if (grepl("one_minus", m)) {
    lim <- c(0, ceiling(max(d$observed, na.rm=TRUE)*1000)/1000)
    br <- pretty(lim, n=4)
    p <- ggplot(d, aes(age_bin_a, age_bin_b, fill=observed)) + geom_tile(colour="white", linewidth=.45) +
      geom_text(aes(label=sprintf("%.4f", observed)), size=2.15) +
      scale_fill_gradient(low="white", high=unname(palette_contract["accent_red"]), limits=lim, breaks=br, name="1 - r") +
      labs(x=NULL, y=NULL, title=metric_titles[m], subtitle=sprintf("Displayed scale: %.3f to %.3f", lim[1], lim[2]))
  } else {
    lo <- floor(min(d$observed, na.rm=TRUE)*1000)/1000; lim <- c(lo, 1)
    br <- unique(c(lo, pretty(lim, n=3), 1)); br <- br[br>=lo & br<=1]
    p <- ggplot(d, aes(age_bin_a, age_bin_b, fill=observed)) + geom_tile(colour="white", linewidth=.45) +
      geom_text(aes(label=sprintf("%.4f", observed)), size=2.15) +
      scale_fill_gradient(low="#EFF3F8", high=unname(palette_contract["signal_blue"]), limits=lim, breaks=br, name="Correlation") +
      labs(x=NULL, y=NULL, title=metric_titles[m], subtitle=sprintf("Displayed scale: %.3f to 1.000", lo))
  }
  p <- p + coord_fixed() + theme(axis.text.x=element_text(angle=40,hjust=1), plot.subtitle=element_text(size=5.8), legend.position="right")
  plots[[m]] <- p
  save_pub_r(p, pfinal("figures", metric_figure_stems[[m]]), 100, 82)
}
panel <- (plots$raw_pearson | plots$raw_spearman) /
  (plots$residualized_pearson | plots$residualized_spearman) /
  (plots$residualized_one_minus_pearson | plot_spacer()) + plot_annotation(tag_levels="a") &
  theme(plot.tag=element_text(face="bold",size=8))
save_pub_r(panel, pfinal("figures/01_centroid_correlation_heatmaps_complete"), 183, 215)

## Sample-level dependence-aware inference -------------------------------------
sim <- fread(pbase("results/table_sample_level_age_bin_similarity.tsv"))
am <- audit[, .(sample_id, age_bracket, sequencing_center)]
sim <- merge(sim, am, by.x="sample_a", by.y="sample_id"); setnames(sim, c("age_bracket","sequencing_center"), c("age_a","center_a"))
sim <- merge(sim, am, by.x="sample_b", by.y="sample_id"); setnames(sim, c("age_bracket","sequencing_center"), c("age_b","center_b"))
sim[, age_bin_separation := abs(as.integer(age_a)-as.integer(age_b))]
sim[, pair_class := factor(fcase(age_bin_separation==0,"within bracket", age_bin_separation==1,"adjacent brackets",
  age_bin_separation %in% 2:3,"intermediate (2-3)", age_bin_separation>=4,"distant (4+)"),
  levels=c("within bracket","adjacent brackets","intermediate (2-3)","distant (4+)"))]
metric_cols <- c("raw_pearson","raw_spearman","residualized_pearson","residualized_spearman")

donor_counts <- audit[, .(donors=.N), by=age_bracket][order(age_bracket)]
pair_counts <- sim[, .(pairs=.N), by=pair_class][order(pair_class)]
write_tsv(donor_counts, pfinal("results/table_similarity_donor_counts.tsv"))
write_tsv(pair_counts, pfinal("results/table_similarity_pair_counts.tsv"))

obs_group <- rbindlist(lapply(metric_cols, function(m) sim[, .(metric=m, mean=mean(get(m)), median=median(get(m))), by=pair_class]))
classes <- levels(sim$pair_class)

B <- 1000L; set.seed(SEED)
boot_rows <- vector("list", B)
ids <- audit$sample_id
age_by_id <- setNames(as.character(audit$age_bracket), ids)
for (b in seq_len(B)) {
  w <- setNames(integer(length(ids)), ids)
  for (a in age_levels) {
    ii <- ids[as.character(audit$age_bracket)==a]
    tb <- table(sample(ii, length(ii), replace=TRUE)); w[names(tb)] <- as.integer(tb)
  }
  pw <- w[sim$sample_a] * w[sim$sample_b]
  rr <- list()
  for (m in metric_cols) for (cl in classes) {
    ii <- sim$pair_class==cl & pw>0
    rr[[length(rr)+1L]] <- data.table(iteration=b, metric=m, pair_class=cl,
      mean=sum(sim[[m]][ii]*pw[ii])/sum(pw[ii]), median=weighted_median(sim[[m]][ii],pw[ii]))
  }
  boot_rows[[b]] <- rbindlist(rr)
  if (b %% 200L == 0L) message("  donor bootstrap ",b,"/",B)
}
boot <- rbindlist(boot_rows); boot[, pair_class:=factor(pair_class,levels=classes)]
write_tsv(boot, pfinal("results/table_similarity_donor_bootstrap_draws.tsv"))

## Permute donor age labels within sequencing center, then reclassify fixed similarities.
perm_rows <- vector("list", B); set.seed(SEED+1L)
for (b in seq_len(B)) {
  lab <- setNames(as.integer(audit$age_bracket), ids)
  for (s in unique(audit$sequencing_center)) {
    ii <- ids[audit$sequencing_center==s]; lab[ii] <- sample(lab[ii])
  }
  sep <- abs(lab[sim$sample_a]-lab[sim$sample_b])
  cls <- factor(fcase(sep==0,"within bracket",sep==1,"adjacent brackets",sep%in%2:3,"intermediate (2-3)",sep>=4,"distant (4+)"),levels=classes)
  rr <- list()
  for (m in metric_cols) for (cl in classes) {
    z <- sim[[m]][cls==cl]
    rr[[length(rr)+1L]] <- data.table(iteration=b,metric=m,pair_class=cl,mean=mean(z),median=median(z))
  }
  perm_rows[[b]] <- rbindlist(rr)
  if (b %% 200L == 0L) message("  donor-label permutation ",b,"/",B)
}
perm <- rbindlist(perm_rows); perm[,pair_class:=factor(pair_class,levels=classes)]
write_tsv(perm, pfinal("results/table_similarity_donor_label_permutation_draws.tsv"))

## Group summaries and contrasts against within-bracket similarity.
group_summary <- merge(obs_group, boot[, .(
  mean_ci_low=quantile(mean,.025), mean_ci_high=quantile(mean,.975),
  median_ci_low=quantile(median,.025), median_ci_high=quantile(median,.975)
),by=.(metric,pair_class)],by=c("metric","pair_class"))
group_summary <- merge(group_summary,pair_counts,by="pair_class")
write_tsv(group_summary,pfinal("results/table_sample_similarity_dependence_aware_summary.tsv"))

contrast_rows <- list()
for (m in metric_cols) for (cl in classes[-1]) for (stat in c("mean","median")) {
  obs_w <- obs_group[metric==m & pair_class==classes[1],get(stat)]
  obs_c <- obs_group[metric==m & pair_class==cl,get(stat)]
  bd <- boot[metric==m & pair_class==cl,get(stat)] - boot[metric==m & pair_class==classes[1],get(stat)]
  pd <- perm[metric==m & pair_class==cl,get(stat)] - perm[metric==m & pair_class==classes[1],get(stat)]
  obs_diff <- obs_c-obs_w
  contrast_rows[[length(contrast_rows)+1L]] <- data.table(metric=m,comparison=paste(cl,"minus within bracket"),statistic=stat,
    observed_difference=obs_diff,bootstrap_ci_low=quantile(bd,.025),bootstrap_ci_high=quantile(bd,.975),
    empirical_permutation_p=emp_p_two(obs_diff,pd),n_permutations=B,bootstrap_replicates=B,
    inference_unit="donor; age labels permuted by donor within sequencing center; donor bootstrap stratified by age bracket")
}
contrasts <- rbindlist(contrast_rows)
write_tsv(contrasts,pfinal("results/table_sample_similarity_dependence_aware_contrasts.tsv"))

plotdat <- melt(group_summary,id.vars=c("metric","pair_class","pairs"),measure.vars=c("mean","median"),variable.name="statistic",value.name="estimate")
lo <- melt(group_summary,id.vars=c("metric","pair_class"),measure.vars=c("mean_ci_low","median_ci_low"),variable.name="tmp",value.name="ci_low")
lo[,statistic:=sub("_ci_low$","",tmp)]; lo[,tmp:=NULL]
hi <- melt(group_summary,id.vars=c("metric","pair_class"),measure.vars=c("mean_ci_high","median_ci_high"),variable.name="tmp",value.name="ci_high")
hi[,statistic:=sub("_ci_high$","",tmp)]; hi[,tmp:=NULL]
plotdat <- merge(plotdat,lo,by=c("metric","pair_class","statistic")); plotdat <- merge(plotdat,hi,by=c("metric","pair_class","statistic"))
plotdat[,expression:=fifelse(grepl("^raw",metric),"Raw expression","Technically residualized")]
plotdat[,correlation:=fifelse(grepl("pearson",metric),"Pearson","Spearman")]
p_sim <- ggplot(plotdat,aes(pair_class,estimate,colour=statistic,group=statistic))+
  geom_line(linewidth=.45,position=position_dodge(width=.22))+geom_point(size=1.5,position=position_dodge(width=.22))+
  geom_errorbar(aes(ymin=ci_low,ymax=ci_high),width=.08,linewidth=.35,position=position_dodge(width=.22))+
  facet_grid(correlation~expression)+scale_colour_manual(values=c(mean=unname(palette_contract["signal_blue"]),median=unname(palette_contract["accent_red"])))+
  labs(x=NULL,y="Similarity (donor-bootstrap 95% CI)",colour=NULL,title="Age-bin sample similarity",
    subtitle="Inference uses donor bootstrap and donor-label permutations; pair counts are descriptive")+
  theme(axis.text.x=element_text(angle=30,hjust=1),legend.position="top")
save_pub_r(p_sim,pfinal("figures/07_sample_similarity_dependence_aware"),183,125)
write_tsv(plotdat,pfinal("figures/source_data/sample_similarity_dependence_aware.tsv"))

message("Centroid and similarity revision complete")
