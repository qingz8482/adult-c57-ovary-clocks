source(file.path(dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),fixed=TRUE))),'common.R'))
p<-rt('tables/frozen_predictions.tsv');cov<-rt('tables/model_gene_coverage.tsv');objs<-readRDS(file.path(out,'inputs/projection_objects.rds'));inp<-readRDS(file.path(out,'inputs/analysis_inputs.rds'));checks<-list()
check<-function(name,value,detail=''){checks[[name]]<<-data.table(check=name,pass=isTRUE(value),detail=as.character(detail))}
check('402_rows_201_donors_two_models',nrow(p)==402&&uniqueN(p$sample_id)==201&&uniqueN(p$model)==2)
check('GTEx_193_unique_donors',uniqueN(p[cohort=='GTEx',donor_id])==193)
check('GSE202601_8_donors_42568_nuclei',uniqueN(p[cohort=='GSE202601',donor_id])==8&&sum(unique(p[cohort=='GSE202601',.(sample_id,n_nuclei)])$n_nuclei)==42568)
check('all_predictions_finite',all(is.finite(p$raw_linear_predictor)))
check('no_mouse_calibration_applied',all(p$raw_linear_predictor==p$mouse_week_equivalent_output_audit_only)&&all(vapply(models,function(m)m$calibration=='none',logical(1))))
check('score_cohort_z_only',all(p[,.(err=max(abs(zsafe(raw_linear_predictor)-cross_species_mouse_ovarian_clock_score))),by=.(cohort,model)]$err<1e-12))
check('12000_gene_coverage_rows',nrow(cov)==12000&&all(cov[,.(n=.N),by=.(cohort,model)]$n==3000))
for(key in names(objs)){
 o<-objs[[key]];nm<-o$map$model[1];cn<-o$meta$cohort[1];m<-models[[nm]];ix<-match(o$meta$sample_id,p[cohort==cn&model==nm,sample_id]);pp<-p[cohort==cn&model==nm][ix]
 y<-drop(o$imputed%*%m$coefficients$coefficient)+m$intercept;zpred<-drop(o$z%*%m$coefficients$standardized_coefficient)+m$model_intercept
 check(paste0(key,'_independent_formula'),max(abs(y-pp$raw_linear_predictor))<1e-8&&max(abs(zpred-y))<1e-8,max(abs(y-pp$raw_linear_predictor)))
 check(paste0(key,'_missing_means'),all(o$z[,!o$map$measured,drop=FALSE]==0))
 check(paste0(key,'_frozen_predict'),max(abs(predict_model(m,o$x)$pred-y))<1e-8)
 check(paste0(key,'_no_human_train_overlap'),!length(intersect(o$meta$sample_id,m$training_keys)))
 check(paste0(key,'_one2one_keys'),!anyDuplicated(o$map$lookup_key[o$map$measured]))
 # Predictions must be invariant to input sample order and other donors' presence.
 single<-predict_model(m,o$x[,1,drop=FALSE])$pred;revp<-predict_model(m,o$x[,rev(seq_len(ncol(o$x))),drop=FALSE])$pred
 check(paste0(key,'_sample_independence'),abs(single-y[1])<1e-8&&max(abs(rev(revp)-y))<1e-8)
 co<-inp$cohorts[[cn]];xn<-log1p(sweep(co$raw,2,colSums(co$raw),'/')*1e6)
 check(paste0(key,'_normalization_replay'),max(abs(xn-co$x))<1e-12)
}
ass<-rt('tables/age_association_sensitivities.tsv')
for(nm in names(models)){
 d<-p[cohort=='GTEx'&model==nm];a<-hc3(cross_species_mouse_ovarian_clock_score~age_order+rin_z+ischemia_z+log_library_z,d,'age_order');b<-ass[model==nm&analysis=='primary_adjusted_ordered']
 check(paste0(nm,'_adjusted_statistic_replay'),max(abs(c(a$estimate-b$estimate,a$p-b$p)))<1e-10)
}
ns<-rt('tables/stageB_null_benchmarks.tsv')
for(i in seq_len(nrow(ns))){r<-ns[i];v<-rt(paste0('null_distributions/stageB_',r$model,'_',r$null_type,'.tsv'))$statistic
 check(paste0(r$model,'_',r$null_type,'_10000_and_p_replay'),length(v)==10000&&all(is.finite(v))&&abs((1+sum(v>=r$observed))/(length(v)+1)-r$empirical_one_sided_p)<1e-12)
 check(paste0(r$model,'_',r$null_type,'_same_observed_statistic'),abs(r$observed-ass[model==r$model&analysis=='primary_adjusted_ordered',estimate])<1e-12)
}
for(nm in names(models)){
 pc<-rt('tables/pca_scores_contributions.tsv');z<-pc[,.(pred=sum(contribution)+intercept[1]),by=.(cohort,sample_id)];cmp<-merge(z,p[model=='VersionB-PCA',.(cohort,sample_id,raw_linear_predictor)],by=c('cohort','sample_id'))
 check('PCA_contribution_sum_replay',max(abs(cmp$pred-cmp$raw_linear_predictor))<1e-8)
}
# Exact n=8 Spearman validation via a separate recursive donor-permutation enumerator.
cs<-rt('tables/age_correlations.tsv')
for(nm in names(models)){
 d<-p[cohort=='GSE202601'&model==nm];xr<-rank(d$age_year)-mean(rank(d$age_year));yr<-rank(d$raw_linear_predictor)-mean(rank(d$raw_linear_predictor));observed<-abs(sum(xr*yr));count<-0L;total<-0L
 enum<-function(v,remaining){if(!length(remaining)){total<<-total+1L;if(abs(sum(xr[v]*yr))>=observed-1e-10)count<<-count+1L}else for(i in seq_along(remaining))enum(c(v,remaining[i]),remaining[-i])}
 enum(integer(),seq_len(8));check(paste0(nm,'_exact_Spearman_replay'),total==40320&&abs(count/total-cs[cohort=='GSE202601'&model==nm,spearman_p])<1e-12)
}
sa<-rt('tables/stageA_direction_conservation.tsv')
for(i in seq_len(nrow(sa))){v<-rt(paste0('null_distributions/stageA_',sa$gene_set[i],'.tsv'))$sign_concordance;check(paste0('StageA_',sa$gene_set[i],'_inclusive_tail_replay'),length(v)==10000&&abs((1+sum(v>=sa$sign_concordance[i]-1e-12))/(length(v)+1)-sa$empirical_p[i])<1e-12)}
man<-rt('provenance/source_checksums.tsv');actual<-vapply(man$path,function(f)strsplit(system2('sha256sum',shQuote(f),stdout=TRUE),' +')[[1]][1],character(1));check('source_files_unchanged',identical(unname(actual),man$sha256))
files<-list.files(file.path(out,'figures'),pattern='[.]pdf$',full.names=TRUE);check('58_reference_PDF_PNG_pairs',length(files)==58&&all(file.info(files)$size>1000)&&all(file.exists(sub('[.]pdf$','.png',files))))
check('current_mouse_StageA_136_profiles',nrow(inp$mouse_metadata)==136)
check('no_missing_GTEx_adjustment_covariates',all(is.finite(as.matrix(inp$cohorts$GTEx$meta[,.(age_order,rin_z,ischemia_z,log_library_z)])) ))
rn<-rt('tables/reference_panels/reference_four_null_families.tsv')
for(i in seq_len(nrow(rn))){rr<-rn[i];vv<-rt(paste0('null_distributions/stageB_',rr$model,'_',rr$null_type,'.tsv'))$statistic
 check(paste0('reference_null_',rr$model,'_',rr$null_type),length(vv)==10000&&all(is.finite(vv))&&abs(rr$observed-ass[model==rr$model&analysis=='primary_adjusted_ordered',estimate])<1e-10&&abs(rr$empirical_one_sided_p-(1+sum(vv>=rr$observed))/10001)<1e-12)
}
gp<-rt('figures/source_data/GSE202601/standardized_score_vs_chronological_age_source_data.tsv')
check('reference_GSE_source_current_scores',all(vapply(seq_len(nrow(gp)),function(i){z<-gp[i];nm<-if(z$model_family=='elasticnet')'HVG-ElasticNet' else 'VersionB-PCA';abs(z$cross_species_mouse_ovarian_clock_score-p[cohort=='GSE202601'&model==nm&donor_id==z$sample_id,cross_species_mouse_ovarian_clock_score])<1e-10},logical(1))))
res<-rbindlist(checks);wt(res,'provenance/validation_checks.tsv');if(!all(res$pass))stop('Validation failure: ',paste(res[pass==FALSE,check],collapse=', '))
message('PASS: ',nrow(res),' validation checks.')
