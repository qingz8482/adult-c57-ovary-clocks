#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]))), "common_final.R"))

message("[2/5] Hardy stability, effect sizes, clock features, and age-pattern classes")
x <- load_inputs(); gtex <- x$gtex; audit <- x$audit
ordered <- fread(pbase("results/table_ordered_age_trend_full.tsv"))
hardy0 <- fread(pbase("results/table_ordered_age_trend_donor_death_sensitivity_full.tsv"))
omni <- fread(pbase("results/table_age_omnibus_nonlinear_full.tsv"))
path0 <- fread(pbase("results/table_reactome_preranked_pathway_enrichment.tsv"))
stagea <- fread(pmod("derived_data/stageA/mouse_human_gene_effect_master.tsv"))

primary <- ordered[analysis_set=="primary_all_193"]
qcpass <- ordered[analysis_set=="sensitivity_predefined_technical_qc_pass"]
gene_sd_hardy <- rowSds(gtex$logcpm[,is.finite(audit$donor_death_hardy_scale),drop=FALSE])
hardy <- hardy0[, .(human_gene_id_stable,human_gene_symbol,
  coefficient=logFC, standardized_beta=logFC/gene_sd_hardy[match(human_gene_id_stable,rownames(gtex$logcpm))],
  signed_t=signed_quasi_likelihood_t,p=PValue,fdr=FDR,expression=logCPM,detection=human_detection_fraction)]
prim <- primary[,.(human_gene_id_stable,human_gene_symbol,coefficient=ordered_age_coefficient_log2fc_per_decade,
  standardized_beta=standardized_coefficient_per_decade,signed_t=signed_quasi_likelihood_t,p=nominal_p,fdr=bh_fdr,
  expression=expression_level_logcpm,detection=detection_prevalence)]
qcp <- qcpass[,.(human_gene_id_stable,human_gene_symbol,coefficient=ordered_age_coefficient_log2fc_per_decade,
  standardized_beta=standardized_coefficient_per_decade,signed_t=signed_quasi_likelihood_t,p=nominal_p,fdr=bh_fdr,
  expression=expression_level_logcpm,detection=detection_prevalence)]
specs <- list(primary=prim,qc_pass=qcp,hardy_adjusted=hardy)

## Gene-level coefficient stability --------------------------------------------
stability <- list(); overlap <- list()
for (nm in c("qc_pass","hardy_adjusted")) {
  d <- merge(prim,specs[[nm]],by="human_gene_id_stable",suffixes=c("_primary","_comparison"))
  psig <- d$fdr_primary<.05
  stability[[nm]] <- data.table(comparison=paste("primary vs",nm),n_genes=nrow(d),
    pearson_coefficient=cor(d$coefficient_primary,d$coefficient_comparison),
    spearman_coefficient=cor(d$coefficient_primary,d$coefficient_comparison,method="spearman"),
    sign_concordance=mean(sign(d$coefficient_primary)==sign(d$coefficient_comparison)),
    primary_fdr_genes=sum(psig),proportion_primary_fdr_retaining_direction=mean(sign(d$coefficient_primary[psig])==sign(d$coefficient_comparison[psig])),
    proportion_primary_fdr_retaining_fdr_lt_0_05=mean(d$fdr_comparison[psig]<.05))
  for (th in c(.10,.25,.50)) {
    a <- d$fdr_primary<.05 & abs(d$coefficient_primary)>th
    b <- d$fdr_comparison<.05 & abs(d$coefficient_comparison)>th
    overlap[[paste(nm,th)]] <- data.table(comparison=paste("primary vs",nm),abs_log2fc_threshold=th,
      primary_n=sum(a),comparison_n=sum(b),overlap_n=sum(a&b),union_n=sum(a|b),jaccard=sum(a&b)/max(1,sum(a|b)),
      overlap_direction_concordance=if(sum(a&b))mean(sign(d$coefficient_primary[a&b])==sign(d$coefficient_comparison[a&b]))else NA_real_)
  }
}
stability <- rbindlist(stability); overlap <- rbindlist(overlap)
write_tsv(stability,pfinal("results/table_hardy_gene_coefficient_stability.tsv"))
write_tsv(overlap,pfinal("results/table_hardy_joint_fdr_effect_overlap.tsv"))

## Effect-size threshold counts, explicitly separated by direction -------------
effect_rows <- list()
for (nm in names(specs)) {
  d <- specs[[nm]]
  defs <- list(
    fdr_lt_0_05=d$fdr<.05,
    fdr_and_abs_log2fc_gt_0_10=d$fdr<.05&abs(d$coefficient)>.10,
    fdr_and_abs_log2fc_gt_0_25=d$fdr<.05&abs(d$coefficient)>.25,
    fdr_and_abs_log2fc_gt_0_50=d$fdr<.05&abs(d$coefficient)>.50,
    fdr_and_abs_standardized_beta_gt_0_10=d$fdr<.05&abs(d$standardized_beta)>.10,
    fdr_and_abs_standardized_beta_gt_0_25=d$fdr<.05&abs(d$standardized_beta)>.25)
  for (dn in names(defs)) {
    z <- defs[[dn]]
    effect_rows[[length(effect_rows)+1L]] <- data.table(specification=nm,criterion=dn,total=sum(z,na.rm=TRUE),
      positive=sum(z&d$coefficient>0,na.rm=TRUE),negative=sum(z&d$coefficient<0,na.rm=TRUE),total_tested=nrow(d))
  }
}
effect_counts <- rbindlist(effect_rows)
write_tsv(effect_counts,pfinal("results/table_effect_size_thresholds_by_specification.tsv"))

## Hardy-adjusted Reactome enrichment and pathway stability --------------------
symbols <- gtex$age_effects$human_gene_symbol[match(rownames(gtex$logcpm),gtex$age_effects$human_gene_id_stable)]
gmt <- read_gmt(pbase("external_data/ReactomePathways.gmt"),15L,500L)
indices <- lapply(gmt,function(gs)which(symbols%in%gs)); indices <- indices[lengths(indices)>=10L&lengths(indices)<=500L]
st_h <- hardy$signed_t[match(rownames(gtex$logcpm),hardy$human_gene_id_stable)]; st_h[!is.finite(st_h)]<-0
cam <- as.data.table(cameraPR(st_h,indices,use.ranks=TRUE,sort=FALSE),keep.rownames="pathway")
cam[,`:=`(specification="hardy_adjusted",normalized_enrichment_score=qnorm(pmax(PValue/2,1e-300),lower.tail=FALSE)*fifelse(Direction=="Up",1,-1))]
cam[,leading_edge_genes:=vapply(pathway,function(pw){ii<-indices[[pw]];zz<-st_h[ii];o<-if(Direction[pathway==pw][1]=="Up")order(zz,decreasing=TRUE)else order(zz);paste(head(unique(symbols[ii][o]),20),collapse=";")},character(1))]
setnames(cam,c("NGenes","Direction","PValue","FDR"),c("n_genes","direction","p","fdr"))
path_sel <- rbindlist(list(path0[specification%in%c("primary_edgeR_technical_adjusted","qc_pass_edgeR_technical_adjusted")],cam),fill=TRUE)
path_sel[,specification:=fcase(specification=="primary_edgeR_technical_adjusted","primary",specification=="qc_pass_edgeR_technical_adjusted","qc_pass",default=specification)]
write_tsv(path_sel,pfinal("results/table_reactome_hardy_adjusted_and_reference.tsv"))

path_stab_rows <- list()
for (nm in c("qc_pass","hardy_adjusted")) {
  a <- path_sel[specification=="primary"]
  b <- path_sel[specification==nm]
  d <- merge(a,b,by="pathway",suffixes=c("_primary","_comparison"))
  led_j <- mapply(function(x,y){xx<-strsplit(x,";",fixed=TRUE)[[1]];yy<-strsplit(y,";",fixed=TRUE)[[1]];length(intersect(xx,yy))/max(1,length(union(xx,yy)))},d$leading_edge_genes_primary,d$leading_edge_genes_comparison)
  path_stab_rows[[nm]] <- data.table(comparison=paste("primary vs",nm),n_pathways=nrow(d),
    pearson_NES=cor(d$normalized_enrichment_score_primary,d$normalized_enrichment_score_comparison),
    spearman_NES=cor(d$normalized_enrichment_score_primary,d$normalized_enrichment_score_comparison,method="spearman"),
    direction_concordance=mean(d$direction_primary==d$direction_comparison),
    primary_fdr_pathways=sum(d$fdr_primary<.05),proportion_primary_fdr_retaining_fdr=mean(d$fdr_comparison[d$fdr_primary<.05]<.05),
    median_leading_edge_jaccard=median(led_j),mean_leading_edge_jaccard=mean(led_j))
}
path_stab <- rbindlist(path_stab_rows)
write_tsv(path_stab,pfinal("results/table_hardy_pathway_stability.tsv"))

## Feature-set comparison with expression/detection/variance matched nulls ------
univ <- data.table(human_gene_id_stable=rownames(gtex$logcpm),symbol=symbols,
  mean_expression=rowMeans(gtex$logcpm),detection=rowMeans(gtex$counts_filtered>=10),variance=rowVars(gtex$logcpm),
  age_t=prim$signed_t[match(rownames(gtex$logcpm),prim$human_gene_id_stable)],
  age_fdr=prim$fdr[match(rownames(gtex$logcpm),prim$human_gene_id_stable)],
  age_beta=prim$coefficient[match(rownames(gtex$logcpm),prim$human_gene_id_stable)])
univ[,`:=`(expr_bin=pmax(1L,pmin(10L,ceiling(frank(mean_expression)/.N*10))),
  detect_bin=pmax(1L,pmin(5L,ceiling(frank(detection)/.N*5))),var_bin=pmax(1L,pmin(10L,ceiling(frank(variance)/.N*10))))]
univ[,stratum:=paste(expr_bin,detect_bin,var_bin,sep="_")]
sets <- list(all_testable_genes=univ$human_gene_id_stable,
  elasticnet_hvg_2000=unique(stagea[elasticnet_hvg_2000==TRUE&!is.na(human_gene_id_stable),human_gene_id_stable]),
  elasticnet_nonzero=unique(stagea[elasticnet_nonzero==TRUE&!is.na(human_gene_id_stable),human_gene_id_stable]),
  pca_versionB_2000=unique(stagea[pca_versionB_2000==TRUE&!is.na(human_gene_id_stable),human_gene_id_stable]),
  deployment_shared_artifact_297=unique(stagea[deployment_shared_artifact==TRUE&!is.na(human_gene_id_stable),human_gene_id_stable]))
sets <- lapply(sets,intersect,y=univ$human_gene_id_stable)
pools <- split(seq_len(nrow(univ)),univ$stratum); set.seed(SEED+2L)
fsum <- list(); fnull <- list()
metric_fun <- function(d)c(median_abs_t=median(abs(d$age_t)),mean_abs_t=mean(abs(d$age_t)),positive_fraction=mean(d$age_t>0),
  sign_imbalance=abs(mean(d$age_t>0)-.5),fdr_significant_proportion=mean(d$age_fdr<.05),
  effect_qualified_proportion=mean(d$age_fdr<.05&abs(d$age_beta)>.25))
for (nm in names(sets)) {
  d <- univ[human_gene_id_stable%in%sets[[nm]]]; obs <- metric_fun(d)
  if(nm=="all_testable_genes") {nullm<-matrix(NA_real_,0,length(obs));ep<-NA_real_;nullmed<-NA_real_} else {
    nullm <- matrix(NA_real_,1000L,length(obs),dimnames=list(NULL,names(obs)))
    for(b in 1:1000){sel<-vapply(d$stratum,function(s){pool<-pools[[s]];sample(pool,1L)},integer(1));nullm[b,]<-metric_fun(univ[sel])}
    ep <- (1+sum(nullm[,"mean_abs_t"]>=obs["mean_abs_t"]))/(nrow(nullm)+1);nullmed<-median(nullm[,"mean_abs_t"])
    fnull[[nm]]<-data.table(feature_set=nm,iteration=1:1000,as.data.table(nullm))
  }
  fsum[[nm]] <- as.data.table(c(
    list(feature_set=nm,n_genes=nrow(d)), as.list(obs),
    list(matched_null_median_mean_abs_t=nullmed,
      empirical_p_more_age_informative=ep,
      classification=ifelse(nm=="all_testable_genes","reference_universe",
        ifelse(ep<.05&obs["mean_abs_t"]>nullmed,"more_age_informative",
          ifelse(obs["mean_abs_t"]<nullmed,"less_age_informative","equivalent_or_marginal"))))
  ))
}
fsum<-rbindlist(fsum,fill=TRUE);fnull<-rbindlist(fnull,fill=TRUE)
write_tsv(fsum,pfinal("results/table_clock_feature_age_signal_matched_variance.tsv"))
write_tsv(fnull,pfinal("results/table_clock_feature_matched_null_draws.tsv"))

## Mutually exclusive age-pattern classes, with and without the sparse oldest bin.
parse_means <- function(x)as.numeric(strsplit(x,";",fixed=TRUE)[[1]])
classify <- function(v,nonlinear=FALSE,tol=.02) {
  d<-diff(v);net<-tail(v,1)-v[1];sgn<-if(net>=0)1 else -1
  opp<-sum(d*sgn < -tol); same<-sum(d*sgn > tol)
  oldest<-abs(tail(d,1)) > 2*(median(abs(head(d,-1)))+tol) && abs(tail(d,1))>.10
  if(oldest)return("oldest-bin-driven")
  if(all(d>=-tol)&&net>tol)return("monotonic increasing")
  if(all(d<=tol)&&net< -tol)return("monotonic decreasing")
  if(opp==1&&same>=length(d)-2)return("approximately monotonic; one reversal")
  if(nonlinear||opp>=2)return("nonlinear")
  "unclassified"
}
pat <- merge(omni,prim[,.(human_gene_id_stable,primary_t=signed_t,primary_fdr=fdr)],by="human_gene_id_stable")
pat <- pat[primary_fdr<.05]
means <- t(vapply(pat$age_bin_means_logcpm,parse_means,numeric(6)))
pat[,classification_all_bins:=vapply(seq_len(.N),function(i)classify(means[i,],nonlinear_effect[i]),character(1))]
pat[,classification_exclude_70_79:=vapply(seq_len(.N),function(i)classify(means[i,1:5],nonlinear_effect[i]),character(1))]
pat[,classification_changed_after_excluding_70_79:=classification_all_bins!=classification_exclude_70_79]
write_tsv(pat[,.(human_gene_id_stable,human_gene_symbol,primary_t,primary_fdr,nonlinear_fdr,age_bin_means_logcpm,
  classification_all_bins,classification_exclude_70_79,classification_changed_after_excluding_70_79)],
  pfinal("results/table_age_pattern_classification_full.tsv"))
psummary <- rbindlist(list(
  pat[,.(N=.N),by=.(classification=classification_all_bins)][,analysis:="all_six_bins"],
  pat[,.(N=.N),by=.(classification=classification_exclude_70_79)][,analysis:="exclude_70_79"]
), use.names=TRUE)
psummary[,proportion:=N/sum(N),by=analysis]
write_tsv(psummary,pfinal("results/table_age_pattern_classification_summary.tsv"))

sel <- pat[order(-abs(primary_t)),head(.SD,8),by=classification_all_bins]
hm <- rbindlist(lapply(seq_len(nrow(sel)),function(i)data.table(human_gene_id_stable=sel$human_gene_id_stable[i],gene=sel$human_gene_symbol[i],
  classification=sel$classification_all_bins[i],age_bracket=age_levels,value=as.numeric(scale(means[match(sel$human_gene_id_stable[i],pat$human_gene_id_stable),])))))
hm[,gene:=factor(gene,levels=rev(unique(gene)))]
p_hm <- ggplot(hm,aes(age_bracket,gene,fill=value))+geom_tile()+facet_grid(classification~.,scales="free_y",space="free_y")+
  scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",limits=c(-2.2,2.2),name="Row z")+
  labs(x="Age bracket",y=NULL,title="Representative fitted age-bin patterns",subtitle="Top absolute ordered-age t statistics within each mutually exclusive class")+
  theme(axis.text.x=element_text(angle=35,hjust=1),axis.text.y=element_text(size=5.2))
save_pub_r(p_hm,pfinal("figures/13_representative_age_pattern_heatmap"),150,190)
write_tsv(hm,pfinal("figures/source_data/representative_age_pattern_heatmap.tsv"))

## Stability figures ------------------------------------------------------------
scatter <- rbindlist(lapply(c("qc_pass","hardy_adjusted"),function(nm){d<-merge(prim,specs[[nm]],by="human_gene_id_stable",suffixes=c("_primary","_comparison"));d[,comparison:=nm];d}))
p_coef <- ggplot(scatter,aes(coefficient_primary,coefficient_comparison,colour=fdr_primary<.05))+geom_abline(slope=1,intercept=0,colour="#BDBDBD")+
  geom_point(size=.35,alpha=.35)+facet_wrap(~comparison)+scale_colour_manual(values=c(`FALSE`="#BDBDBD",`TRUE`=unname(palette_contract["accent_red"])))+
  labs(x="Primary log2FC per decade",y="Sensitivity log2FC per decade",colour="Primary FDR<.05",title="Gene-level age-effect stability")
save_pub_r(p_coef,pfinal("figures/08_hardy_gene_coefficient_stability"),145,75)

pathplot <- rbindlist(lapply(c("qc_pass","hardy_adjusted"),function(nm){d<-merge(path_sel[specification=="primary"],path_sel[specification==nm],by="pathway",suffixes=c("_primary","_comparison"));d[,comparison:=nm];d}))
p_path <- ggplot(pathplot,aes(normalized_enrichment_score_primary,normalized_enrichment_score_comparison,colour=fdr_primary<.05))+geom_abline(slope=1,intercept=0,colour="#BDBDBD")+
  geom_point(size=.55,alpha=.55)+facet_wrap(~comparison)+scale_colour_manual(values=c(`FALSE`="#BDBDBD",`TRUE`=unname(palette_contract["signal_blue"])))+
  labs(x="Primary Reactome NES",y="Sensitivity Reactome NES",colour="Primary FDR<.05",title="Pathway-level age-signal stability")
save_pub_r(p_path,pfinal("figures/09_hardy_pathway_stability"),145,75)

write_tsv(scatter,pfinal("figures/source_data/hardy_gene_coefficient_stability.tsv"))
write_tsv(pathplot,pfinal("figures/source_data/hardy_pathway_stability.tsv"))
message("Hardy/feature/pattern revision complete")
