#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
arg<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),fixed=TRUE)
out<-normalizePath(file.path(dirname(arg),'..')); development<-normalizePath(file.path(out,'../../..'))
source(file.path(development,'scripts/core.R'))
dd<-file.path(out,'derived_data');write_tsv<-function(x,n)fwrite(x,file.path(dd,paste0(n,'.tsv')),sep='\t')
m<-readRDS(file.path(development,'locked_models/count_tuning_20261011/VersionB-PCA-glmnet/model.rds'))
inp<-readRDS(file.path(development,'shared_inputs/prepared_expression.rds'))
meta<-as.data.table(inp$metadata)[match(m$training_keys,sample_key)]
x<-inp$expression[m$representation$genes,m$training_keys,drop=FALSE]
stopifnot(identical(meta$sample_key,m$training_keys),m$recipe=='VersionB-PCA',m$calibration=='none',all(is.finite(x)))
rep<-m$representation;cfg<-m$config;hp<-m$hyperparameters
z<-sweep(sweep(t(x),2,rep$center,'-'),2,rep$scale,'/')
scores<-z%*%rep$rotation
candidate_pcs<-colnames(scores);pcs<-candidate_pcs[m$pc_coefficients!=0]
beta<-setNames(m$pc_coefficients,candidate_pcs)[pcs]
stopifnot(identical(pcs,c('PC1','PC2','PC5','PC8')))
y<-meta$age_week;studies<-sort(unique(meta$study));deployment_lambda<-hp$lambda

# a: the actual grouped, fold-local feature/PCA tuning path for the winning
# N/PC/alpha candidate. Nonzero counts refer to the all-training deployment fit.
cv<-fread(file.path(development,'VersionB-PCA/glmnet/tuning_audit/final_r00_f0/count_tuning_lambda_summary.csv'))
cv_dt<-cv[N==hp$N & pc==hp$pc & alpha==hp$alpha]
stopifnot(nrow(cv_dt)==100,all(cv_dt$n_studies==length(studies)))
setorder(cv_dt,-lambda)
ix<-vapply(cv_dt$lambda,function(v)which.min(abs(m$glmnet$lambda-v)),integer(1))
stopifnot(max(abs(m$glmnet$lambda[ix]-cv_dt$lambda))<1e-8)
cv_dt[,`:=`(log_lambda=log(lambda),se_mae=se,lower=mean_mae-se,upper=mean_mae+se,n_nonzero=as.integer(m$glmnet$df[ix]))]
winner<-cv_dt[which.min(abs(lambda-deployment_lambda))]
stopifnot(abs(winner$mean_mae-hp$inner_MAE)<1e-10)
threshold<-min(cv_dt$mean_mae)+cv_dt$se[which.min(cv_dt$mean_mae)]
stopifnot(abs(max(cv_dt[mean_mae<=threshold,lambda])-deployment_lambda)<1e-8)
write_tsv(cv_dt,'panel_a_cv_error_path')

# b: use the saved glmnet path, not a newly fitted path with different weights.
path_beta<-as.matrix(m$glmnet$beta)
path_dt<-as.data.table(as.table(path_beta));setnames(path_dt,c('pc','lambda_index','coefficient'))
path_dt[,lambda_index:=as.integer(sub('^s','',lambda_index))+1L]
path_dt[,`:=`(lambda=m$glmnet$lambda[lambda_index],log_lambda=log(m$glmnet$lambda[lambda_index]),selected=pc%in%pcs)]
write_tsv(path_dt,'panel_b_coefficient_paths')

# c: conditional leave-one-active-PC-out diagnostic. The deployment gene
# selection/scaling/rotation are fixed across all folds. Regression weights,
# alpha/lambda tuning and regression fitting use training studies only.
# This is not fully nested prediction validation.
fit_subset<-function(ii,keep,a)fit_path(scores[ii,keep,drop=FALSE],meta[ii],a,cfg,'primary',pc=0L)
tuning_audit<-list();outer_audit<-list();inner_membership<-list();ti<-0L;oi<-0L;mi<-0L
oof_refit<-function(keep,label) {
 ans<-vector('list',length(studies))
 for(ff in seq_along(studies)) {
  test_study<-studies[ff];tr<-which(meta$study!=test_study);te<-which(meta$study==test_study)
  in_studies<-sort(unique(meta$study[tr]));grid<-list();gi<-0L
  for(inner in in_studies) {
   itr<-tr[meta$study[tr]!=inner];ite<-tr[meta$study[tr]==inner]
   stopifnot(!length(intersect(meta$study[itr],meta$study[ite])),!test_study%in%meta$study[c(itr,ite)])
   mi<<-mi+1L;inner_membership[[mi]]<<-data.table(variant=label,outer_test_study=test_study,inner_test_study=inner,
    train_studies=paste(sort(unique(meta$study[itr])),collapse=';'),n_train=length(itr),n_test=length(ite))
   for(a in cfg$alphas) {
    fit<-fit_subset(itr,keep,a)
    pred<-as.matrix(predict(fit,newx=model_x(scores[ite,keep,drop=FALSE])))
    gi<-gi+1L;grid[[gi]]<-data.table(alpha=a,lambda=fit$lambda,inner_study=inner,mae=colMeans(abs(pred-y[ite])))
   }
  }
  grid<-rbindlist(grid)
  summary<-grid[,.(mean_mae=mean(mae),se=sd(mae)/sqrt(.N),n_studies=.N),by=.(alpha,lambda)]
  stopifnot(all(summary$n_studies==length(in_studies)))
  choices<-summary[,{
   imin<-which.min(mean_mae);lim<-mean_mae[imin]+se[imin];j<-which.max(ifelse(mean_mae<=lim,lambda,-Inf))
   .(lambda=lambda[j],inner_MAE=mean_mae[j],inner_SE=se[j],threshold=lim)
  },by=alpha]
  setorder(choices,inner_MAE,-alpha);best<-choices[1]
  fit<-fit_subset(tr,keep,best$alpha)
  pred<-as.numeric(predict(fit,newx=model_x(scores[te,keep,drop=FALSE]),s=best$lambda))
  ans[[ff]]<-data.table(sample_index=te,sample_key=meta$sample_key[te],study=test_study,age_week=y[te],
   prediction=pred,outer_fold=ff,alpha=best$alpha,lambda=best$lambda,cv_mae=best$inner_MAE)
  ti<<-ti+1L;summary[,`:=`(variant=label,outer_test_study=test_study)];tuning_audit[[ti]]<<-summary
  oi<<-oi+1L;outer_audit[[oi]]<<-data.table(variant=label,outer_test_study=test_study,keep_pcs=paste(keep,collapse=';'),
   train_studies=paste(sort(unique(meta$study[tr])),collapse=';'),n_train=length(tr),n_test=length(te),
   alpha=best$alpha,lambda=best$lambda,inner_MAE=best$inner_MAE,
   train_weight_sum=sum(training_weights(meta[tr],cfg)),test_samples_used_in_fit=0L)
  message(label,': held out ',test_study,' (',ff,'/',length(studies),')')
 }
 rbindlist(ans)[order(sample_index)]
}
full<-oof_refit(pcs,'full_active_PC_refit');write_tsv(full,'panel_c_full_conditional_oof_predictions')
lopo<-rbindlist(lapply(pcs,function(pc){
 q<-oof_refit(setdiff(pcs,pc),paste0('omit_',pc))
 stopifnot(identical(q$sample_key,full$sample_key))
 data.table(sample_key=meta$sample_key,study=meta$study,omitted_pc=pc,age_week=y,outer_fold=q$outer_fold,
  full_prediction=full$prediction,lopo_prediction=q$prediction,full_abs_error=abs(full$prediction-y),
  lopo_abs_error=abs(q$prediction-y),paired_delta_abs_error=abs(q$prediction-y)-abs(full$prediction-y),
  alpha=q$alpha,lambda=q$lambda,cv_mae=q$cv_mae)
}))
write_tsv(lopo,'panel_c_lopo_sample_metrics')
write_tsv(rbindlist(tuning_audit),'panel_c_inner_tuning_audit')
write_tsv(rbindlist(outer_audit),'panel_c_outer_fit_audit')
write_tsv(rbindlist(inner_membership),'panel_c_inner_fold_membership')
set.seed(20260630)
boot<-function(v){b<-replicate(5000,mean(sample(v,length(v),replace=TRUE)));data.table(delta_mae=mean(v),ci_low=unname(quantile(b,.025)),ci_high=unname(quantile(b,.975)),n=length(v))}
lopo_summary<-lopo[,boot(paired_delta_abs_error),by=omitted_pc]
write_tsv(lopo_summary,'panel_c_lopo_delta_mae_bootstrap_ci')
perstudy<-lopo[,.(full_mae=mean(full_abs_error),lopo_mae=mean(lopo_abs_error),delta_mae=mean(paired_delta_abs_error)),by=.(omitted_pc,study)]
write_tsv(perstudy,'panel_c_lopo_per_study_metrics')
set.seed(20260630)
study_summary<-perstudy[,boot(delta_mae),by=omitted_pc]
write_tsv(study_summary,'panel_c_lopo_study_bootstrap_ci')
secondary<-lopo[,.(delta_rmse=sqrt(mean((lopo_prediction-age_week)^2))-sqrt(mean((full_prediction-age_week)^2)),
 delta_r_squared=cor(lopo_prediction,age_week)^2-cor(full_prediction,age_week)^2,
 delta_correlation=cor(lopo_prediction,age_week)-cor(full_prediction,age_week),full_MAE=mean(full_abs_error),lopo_MAE=mean(lopo_abs_error)),by=omitted_pc]
write_tsv(secondary,'panel_c_lopo_secondary_metrics')

# d: exact frozen deployment contributions and deterministic 20-bin summary.
contrib<-sweep(scores[,pcs,drop=FALSE],2,beta,'*');pred_parts<-m$model_intercept+rowSums(contrib)
native<-predict_model(m,inp$expression[,m$training_keys,drop=FALSE])$pred
write_tsv(data.table(sample_key=meta$sample_key,tabulated_prediction=native,reconstructed_prediction=pred_parts,
 absolute_discrepancy=abs(native-pred_parts)),'sample_contribution_sum_validation')
contrib_dt<-melt(as.data.table(contrib,keep.rownames='sample_key'),id.vars='sample_key',variable.name='pc',value.name='contribution')
contrib_dt<-merge(contrib_dt,data.table(sample_key=meta$sample_key,age_week=y,predicted_age=pred_parts),by='sample_key')
setorder(contrib_dt,age_week,sample_key)
write_tsv(contrib_dt,'panel_d_unaggregated_pc_contributions')
bins<-unique(contrib_dt[,.(sample_key,age_week,predicted_age)])[order(age_week,sample_key)]
bins[,age_bin:=sprintf('Bin %02d',pmin(20L,ceiling(seq_len(.N)/(.N/20))))]
d_plot<-merge(contrib_dt,bins[,.(sample_key,age_bin)],by='sample_key')[,.(contribution=median(contribution),
 age_week=median(age_week),predicted_age=median(predicted_age),n_samples=.N),by=.(age_bin,pc)]
setorder(d_plot,age_bin,pc);write_tsv(d_plot,'panel_d_age_binned_pc_contributions');write_tsv(bins,'panel_d_sample_bins')

# e: back-project through the saved gene scales onto actual input ln1p(CPM).
w_scaled<-sweep(rep$rotation[,pcs,drop=FALSE],2,beta,'*')
w_raw<-sweep(w_scaled,1,rep$scale,'/');w_integrated_raw<-rowSums(w_raw)
gene_long<-as.data.table(as.table(cbind(w_raw,Integrated=w_integrated_raw)));setnames(gene_long,c('gene','component','weight'))
write_tsv(gene_long,'gene_pc_and_integrated_backprojected_weights')
pick<-function(v,n)unique(c(names(sort(v,decreasing=TRUE))[seq_len(n)],names(sort(v))[seq_len(n)]))
genes<-unique(c(unlist(lapply(pcs,function(pc)pick(w_raw[,pc],2))),pick(w_integrated_raw,5)))
if(length(genes)<24)genes<-unique(c(genes,unlist(lapply(pcs,function(pc)pick(w_raw[,pc],3)))))
genes<-head(genes,40);grp<-pcs[max.col(abs(w_raw[genes,pcs,drop=FALSE]),ties.method='first')]
genes<-genes[order(grp,-abs(w_integrated_raw[genes]))]
e_dt<-gene_long[gene%in%genes];e_dt[,display_order:=match(gene,genes)];setorder(e_dt,display_order,component)
write_tsv(e_dt,'panel_e_display_gene_weights')
gene_intercept<-m$model_intercept-sum(rep$center*w_integrated_raw)
gene_pred<-gene_intercept+drop(t(x)%*%w_integrated_raw)
write_tsv(data.table(sample_key=meta$sample_key,pc_predictor=pred_parts,gene_predictor=gene_pred,
 absolute_discrepancy=abs(pred_parts-gene_pred)),'gene_backprojection_validation')
checks<-data.table(check=c('native_prediction_decomposition','gene_prediction_backprojection','gene_slopes_match_locked_export',
 'gene_intercept_matches_locked_export','selected_lambda_matches_saved_1SE_path','all_136_profiles_in_each_ablation','all_14_studies_held_out_each_variant'),
 error=c(max(abs(native-pred_parts)),max(abs(pred_parts-gene_pred)),max(abs(w_integrated_raw-m$coefficients$coefficient)),
 abs(gene_intercept-m$intercept),abs(winner$mean_mae-hp$inner_MAE),max(abs(lopo[,uniqueN(sample_key),by=omitted_pc]$V1-136)),
 max(abs(rbindlist(outer_audit)[,uniqueN(outer_test_study),by=variant]$V1-14))))
checks[,passed:=error<1e-8];stopifnot(all(checks$passed));write_tsv(checks,'numerical_validation')
saveRDS(list(pcs=pcs,candidate_pcs=candidate_pcs,lambda=deployment_lambda,alpha=hp$alpha,n_genes=hp$N,
 n_samples=nrow(meta),n_studies=length(studies),cv_dt=cv_dt,path_dt=path_dt,lopo_summary=lopo_summary,study_summary=study_summary,
 d_plot=d_plot,e_dt=e_dt,genes_show=genes,beta=beta),file.path(dd,'plot_objects.rds'))
writeLines(capture.output(sessionInfo()),file.path(dd,'sessionInfo.txt'))
message('Model selection, ablation, decomposition and gene back-projection complete.')
