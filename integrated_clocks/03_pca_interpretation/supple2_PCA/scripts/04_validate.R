#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(glmnet)})
setDTthreads(1);glmnet.control(fdev=0,devmax=1)
arg<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),fixed=TRUE)
out<-normalizePath(file.path(dirname(arg),'..'));dd<-file.path(out,'derived_data');dev<-normalizePath(file.path(out,'../../..'))
read<-function(n)fread(file.path(dd,paste0(n,'.tsv')))
checks<-list();check<-function(n,ok,error=NA_real_){checks[[length(checks)+1]]<<-data.table(check=n,passed=isTRUE(ok),error=error);if(!isTRUE(ok))stop(n)}
m<-readRDS(file.path(dev,'locked_models/count_tuning_20261011/VersionB-PCA-glmnet/model.rds'))
inp<-readRDS(file.path(dev,'shared_inputs/prepared_expression.rds'))
meta<-as.data.table(inp$metadata)[match(m$training_keys,sample_key)]
z<-t(inp$expression[m$representation$genes,m$training_keys,drop=FALSE]);z<-sweep(sweep(z,2,m$representation$center,'-'),2,m$representation$scale,'/')
scores<-z%*%m$representation$rotation
fit_audit<-read('panel_c_outer_fit_audit');inner<-read('panel_c_inner_fold_membership');tune<-read('panel_c_inner_tuning_audit')
lopo<-read('panel_c_lopo_sample_metrics');full<-read('panel_c_full_conditional_oof_predictions')
check('70_conditional_outer_refits',nrow(fit_audit)==70&&!anyDuplicated(fit_audit[,.(variant,outer_test_study)]))
check('910_inner_study_splits',nrow(inner)==910)
check('outer_test_studies_excluded_from_all_inner_roles',all(mapply(function(o,t,s)!o%in%c(t,strsplit(s,';',fixed=TRUE)[[1]]),inner$outer_test_study,inner$inner_test_study,inner$train_studies)))
check('inner_test_studies_excluded_from_inner_training',all(mapply(function(t,s)!t%in%strsplit(s,';',fixed=TRUE)[[1]],inner$inner_test_study,inner$train_studies)))
check('tuning_complete_for_13_inner_studies',nrow(tune)==70*6*100&&all(tune$n_studies==13))
replay<-rbindlist(lapply(seq_len(nrow(fit_audit)),function(i){
 a<-fit_audit[i];tr<-which(meta$study!=a$outer_test_study);te<-which(meta$study==a$outer_test_study)
 keep<-strsplit(a$keep_pcs,';',fixed=TRUE)[[1]]
 ww<-1/as.numeric(table(meta$study[tr])[meta$study[tr]]);ww<-ww/mean(ww)
 path<-10^seq(m$config$lambda_log10_max,m$config$lambda_log10_min,length.out=m$config$lambda_values)
 fit<-glmnet(scores[tr,keep,drop=FALSE],meta$age_week[tr],weights=ww,alpha=a$alpha,lambda=path,
  standardize=FALSE,intercept=TRUE,control=list(thresh=m$config$glmnet_thresh,maxit=m$config$glmnet_maxit))
 pp<-drop(predict(fit,scores[te,keep,drop=FALSE],s=a$lambda))
 saved<-if(a$variant=='full_active_PC_refit')full$prediction[match(meta$sample_key[te],full$sample_key)] else {
  q<-lopo[omitted_pc==sub('omit_','',a$variant)];q$lopo_prediction[match(meta$sample_key[te],q$sample_key)]}
 # Independently recalculate the candidate-specific 1SE selections and alpha choice.
 tt<-tune[variant==a$variant & outer_test_study==a$outer_test_study]
 winners<-rbindlist(lapply(unique(tt$alpha),function(aa){
  t<-tt[alpha==aa];minrow<-t[which.min(mean_mae)];allowed<-t[mean_mae<=minrow$mean_mae+minrow$se];best<-allowed[which.max(lambda)]
  data.table(alpha=aa,lambda=best$lambda,mae=best$mean_mae)
 }));setorder(winners,mae,-alpha);best<-winners[1]
 data.table(variant=a$variant,study=a$outer_test_study,max_prediction_error=max(abs(pp-saved)),
  lambda_error=abs(best$lambda-a$lambda),alpha_error=abs(best$alpha-a$alpha),inner_mae_error=abs(best$mae-a$inner_MAE))
}))
fwrite(replay,file.path(dd,'conditional_oof_independent_replay.tsv'),sep='\t')
check('70_outer_prediction_replays',max(replay$max_prediction_error)<1e-8,max(replay$max_prediction_error))
check('70_alpha_lambda_selections_recomputed',max(replay$lambda_error)<1e-8&&max(replay$alpha_error)==0&&max(replay$inner_mae_error)<1e-10)
check('paired_error_differences',max(abs(lopo$paired_delta_abs_error-(abs(lopo$lopo_prediction-lopo$age_week)-abs(lopo$full_prediction-lopo$age_week))))<1e-10)
summary<-read('panel_c_lopo_delta_mae_bootstrap_ci')
err<-lopo[,.(mae=mean(paired_delta_abs_error)),by=omitted_pc]
check('sample_delta_MAE_recomputed',max(abs(summary$delta_mae-err$mae[match(summary$omitted_pc,err$omitted_pc)]))<1e-10)
con<-read('panel_d_unaggregated_pc_contributions');bins<-read('panel_d_sample_bins');agg<-read('panel_d_age_binned_pc_contributions')
check('544_exact_PC_contributions',nrow(con)==136*4&&!anyDuplicated(con[,.(sample_key,pc)]))
check('20_nonoverlapping_age_bins',nrow(bins)==136&&!anyDuplicated(bins$sample_key)&&uniqueN(bins$age_bin)==20&&all(diff(bins$age_week)>=0))
expected<-merge(con,bins[,.(sample_key,age_bin)],by='sample_key')[,.(expected=median(contribution)),by=.(age_bin,pc)]
both<-merge(agg,expected,by=c('age_bin','pc'))
check('age_bin_medians_recomputed',max(abs(both$contribution-both$expected))<1e-10)
slots<-read('gsea_panel_slots')[panel%in%c('f','f2','f3')];selected<-slots[status=='selected']
check('GSEA_panels_have_12_documented_slots',nrow(slots)==12)
check('10_distinct_selected_pathways',nrow(selected)==10&&!anyDuplicated(selected$gene_set_name))
steroid<-read('panel_f3_all_tested_steroid_terms')
check('absent_significant_steroid_terms_not_substituted',all(steroid$fdr_within_pc>=.05)&&slots[panel=='f3'&category=='steroid / cholesterol',all(status=='no_terms_pass_FDR_0.05')])
for(pn in c('f','f2','f3')) {
 s<-read(paste0('panel_',pn,'_selected_terms'));c<-read(paste0('panel_',pn,'_enrichment_curves'));q<-read(paste0('panel_',pn,'_peak_rank_validation'))
 check(paste0(pn,'_significance_threshold'),all(s$fdr_within_pc<.05))
 check(paste0(pn,'_exact_enrichment_reconstruction'),all(q$peak_exact_match)&&max(q$es_absolute_discrepancy)<1e-12)
 check(paste0(pn,'_running_sum_closes_at_zero'),max(abs(c[,tail(running_es,1),by=gene_set_name]$V1))<1e-10)
 a<-read(paste0('panel_',pn,'_selection_audit'))
 check(paste0(pn,'_selection_redundancy_threshold'),all(a[selected==TRUE,max_jaccard_at_evaluation]<.5))
}
qc<-rbindlist(checks);fwrite(qc,file.path(dd,'independent_validation.tsv'),sep='\t')
message(nrow(qc),' independent checks passed.')
