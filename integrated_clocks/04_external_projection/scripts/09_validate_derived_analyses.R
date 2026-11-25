#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
a<-commandArgs(FALSE);here<-dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE)));out<-dirname(here)
rd<-function(f)fread(file.path(out,f))
p<-rd('external_projection_predictions_all.tsv');tt<-rd('tables/all_reference_contrast_statistics.tsv');mem<-rd('tables/all_reference_contrast_membership.tsv')
replay<-rbindlist(lapply(seq_len(nrow(tt)),function(i){
 z<-tt[i];m<-mem[analysis_contrast_id==z$analysis_contrast_id];d<-p[model==z$model]
 x<-d[sample_key%in%m[!contrast_role%in%c('control','placement_reference'),sample_key],delta_age_week]
 y<-d[sample_key%in%m[contrast_role%in%c('control','placement_reference'),sample_key],delta_age_week]
 nx<-length(x);ny<-length(y);vx<-var(x)/nx;vy<-var(y)/ny
 statistic<-(mean(x)-mean(y))/sqrt(vx+vy);df<-(vx+vy)^2/(vx^2/(nx-1)+vy^2/(ny-1));pv<-2*pt(-abs(statistic),df)
 ci<-mean(x)-mean(y)+c(-1,1)*qt(.975,df)*sqrt(vx+vy)
 data.table(analysis_contrast_id=z$analysis_contrast_id,model=z$model,t_error=abs(statistic-z$t_statistic),df_error=abs(df-z$df),p_error=abs(pv-z$welch_p),ci_error=max(abs(ci-c(z$ci95_low,z$ci95_high))),n_match=nx==z$n_experimental&&ny==z$n_control)
}))
stopifnot(nrow(replay)==90,all(replay$n_match),max(replay$t_error)<1e-9,max(replay$df_error)<1e-9,max(replay$p_error)<1e-9,max(replay$ci_error)<1e-9)
ref<-rd('tables/frozen_deployment_reference_predictions.tsv')
distance<-p[,{mo<-.BY$model;refd<-ref[model==mo,projection_distance];.(max_percentile_error=max(abs(projection_distance_percentile-vapply(projection_distance,function(z)mean(refd<=z),numeric(1)))))},by=model]
stopifnot(max(distance$max_percentile_error)<1e-12)
r<-rd('external_projection_relative_to_controls.tsv');r[,expected_relative:=delta_age_week-mean(delta_age_week[role=='Reference']),by=.(analysis_contrast_id,model)]
stopifnot(max(abs(r$expected_relative-r$relative_delta_age_week))<1e-9)
original<-rd('tables/external_projection_relative_to_controls.tsv');current<-r[analysis_family=='prespecified']
cmp<-merge(original[,.(sample_key,analysis_contrast_id,model,original=relative_delta_age_week)],current[,.(sample_key,analysis_contrast_id,model,current=relative_delta_age_week)],by=c('sample_key','analysis_contrast_id','model'))
stopifnot(nrow(cmp)==nrow(original),max(abs(cmp$original-cmp$current))<1e-9)
fwrite(replay,file.path(out,'provenance/all_90_welch_formula_replay.tsv'),sep='\t')
fwrite(distance,file.path(out,'provenance/ood_percentile_replay.tsv'),sep='\t')
fwrite(data.table(check=c('90 Welch statistics and confidence intervals independently replayed','training-reference empirical distance percentiles replayed','relative means centered within contrast/model','prespecified relative deltas identical'),pass=TRUE),file.path(out,'provenance/derived_analysis_validation.tsv'),sep='\t')
message('PASS: all 90 Welch tests, OOD percentiles and relative deltas independently verified.')
