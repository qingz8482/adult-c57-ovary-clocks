#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a<-commandArgs(FALSE);here<-dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE)));out<-dirname(here);dev<-file.path(dirname(out),'model_development')
new<-readRDS(file.path(out,'inputs/model_gene_expression.rds'));old<-readRDS(file.path(dev,'shared_inputs/prepared_expression.rds'))
keys<-intersect(colnames(new$expression),colnames(old$expression));stopifnot(length(keys)>0)
xx<-new$expression[,keys,drop=FALSE];yy<-old$expression[rownames(xx),keys,drop=FALSE]
stopifnot(identical(is.na(xx),is.na(yy)));e<-max(abs(xx-yy),na.rm=TRUE);stopifnot(e<1e-10)
fwrite(data.table(check='expression agrees with development preparation for overlapping profiles',n_profiles=length(keys),max_absolute_error=e,pass=TRUE),file.path(out,'provenance/training_overlap_expression_replay.tsv'),sep='\t')
# Recompute every saved test from the plotted samples and compare all numeric statistics.
d<-fread(file.path(out,'tables/external_projection_relative_to_controls.tsv'));s<-fread(file.path(out,'tables/external_projection_contrast_statistics.tsv'))
source_contrasts<-fread(file.path(out,'inputs/prespecified_contrasts.tsv'))
ns<-s[model=='HVG-ElasticNet',.(n_control=sum(n_control),n_experimental=sum(n_experimental)),by=contrast_id]
ns<-merge(ns,source_contrasts,by='contrast_id');stopifnot(nrow(ns)==40,all(ns$n_control==ns$n_reference_or_control_samples),all(ns$n_experimental==ns$n_experimental_or_target_samples))
fwrite(ns[,.(contrast_id,n_control,n_experimental,n_reference_or_control_samples,n_experimental_or_target_samples)],file.path(out,'provenance/contrast_membership_count_replay.tsv'),sep='\t')
errs<-lapply(seq_len(nrow(s)),function(i){r<-s[i];v<-d[analysis_contrast_id==r$analysis_contrast_id&model==r$model];tt<-t.test(v[role=='Target',delta_age_week],v[role=='Reference',delta_age_week]);data.table(analysis_contrast_id=r$analysis_contrast_id,model=r$model,p_error=abs(tt$p.value-r$welch_p),t_error=abs(unname(tt$statistic)-r$t_statistic),ci_error=max(abs(tt$conf.int-c(r$ci95_low,r$ci95_high))))})
errs<-rbindlist(errs);stopifnot(max(errs$p_error)<1e-12,max(errs$t_error)<1e-10,max(errs$ci_error)<1e-9)
fwrite(errs,file.path(out,'provenance/statistics_replay_checks.tsv'),sep='\t')
pc<-fread(file.path(out,'tables/pca_projection_pc_contributions.tsv'));p<-fread(file.path(out,'external_projection_predictions.tsv'))[model=='VersionB-PCA/glmnet']
reconstructed<-pc$model_intercept+rowSums(as.matrix(pc[,grep('^PC[0-9]+$',names(pc)),with=FALSE]));er<-max(abs(reconstructed-p$predicted_age_week[match(pc$sample_key,p$sample_key)]));stopifnot(er<1e-8)
fwrite(data.table(check='PCA contributions sum to prediction',max_absolute_error=er,pass=TRUE),file.path(out,'provenance/pca_contribution_replay.tsv'),sep='\t')
message('PASS: overlapping training inputs, all statistics, and PCA contribution reconstruction.')
