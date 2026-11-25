set.seed(SEED)
palette<-c(elasticnet='#2878B5',pca_glmnet='#D65F5F',neutral='#606060',null='#B8B8B8',accent='#2A9D8F')
theme_set(theme_classic(base_size=7,base_family='sans')+theme(axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3),strip.background=element_blank(),strip.text=element_text(face='bold'),plot.title=element_text(face='bold',size=8),legend.key.height=unit(3,'mm')))
meta<-copy(analysis_current$cohorts$GTEx$meta);meta[,age_bracket:=factor(age_bracket,levels=age_levels)]
dall<-copy(cp[cohort=='GTEx']);dall[,age_bracket:=factor(age_bracket,levels=age_levels)]
recon<-data.table(count_label=c(paste0('deployment_shared_artifact_',length(shared)),'historical_all_selected_gene_sensitivity_277','historical_pca_selector_diagnostic_262'),count=c(length(shared),277L,262L),authority=c('AUTHORITATIVE_FOR_DEPLOYMENT','HISTORICAL','HISTORICAL'))
coverage<-cm[cohort=='GTEx',.(selected_gene_count=.N,mapped_one_to_one_count=sum(one2one_mapped),human_detected_count=sum(measured),raw_gene_coverage=mean(one2one_mapped),absolute_weight_coverage=sum(abs(raw_coefficient[measured]))/sum(abs(raw_coefficient))),by=model_family]
stagea_stats<-rt('tables/stageA_direction_conservation.tsv');stagea_stats[,`:=`(analysis=fcase(gene_set=='hvg_input_genes','elasticnet_hvg_3000',gene_set=='pca_input_genes','pca_versionB_3000',default=gene_set),model_family=ifelse(grepl('^hvg',gene_set),'elasticnet','pca_glmnet'),spearman_correlation=spearman_signed_t,gene_count=n_genes,matched_null_empirical_p=empirical_p)]
# Restore the reference's independent absolute-weight/sign permutation family.
objs<-readRDS(file.path(out,'inputs/projection_objects.rds'));ns<-rt('tables/stageB_null_benchmarks.tsv')
added<-list()
for(nm in names(models)){
 o<-objs[[paste('GTEx',nm,sep='__')]];d<-o$meta;X<-model.matrix(~age_order+rin_z+ischemia_z+log_library_z,d);pr<-solve(crossprod(X),t(X))[2,];obs<-drop(pr%*%o$score);mapped<-o$map$measured;Z<-o$z[,mapped,drop=FALSE];w<-models[[nm]]$coefficients$standardized_coefficient[mapped];v<-numeric(B)
 cache<-file.path(out,'null_distributions',paste0('stageB_',nm,'_independent_weight_sign_permutation.tsv'))
 if(file.exists(cache)&&all(is.finite(fread(cache)$statistic)))v<-fread(cache)$statistic else {
  for(st in seq(1L,B,by=100L)){ids<-st:min(B,st+99L);W<-replicate(length(ids),sample(abs(w))*sample(sign(w)));S<-Z%*%W;ss<-apply(S,2,sd);S<-sweep(sweep(S,2,colMeans(S),'-'),2,ifelse(ss==0,1,ss),'/');v[ids]<-drop(pr%*%S)}
  fwrite(data.table(iteration=seq_len(B),statistic=v),cache,sep='\t')
 }
 added[[nm]]<-data.table(model=nm,null_type='independent_weight_sign_permutation',observed=obs,null_q025=quantile(v,.025),null_q975=quantile(v,.975),empirical_one_sided_p=(1+sum(v>=obs))/(B+1))
}
ns<-rbindlist(list(ns,rbindlist(added)),fill=TRUE)
keepnull<-c(mapped_weight_permutation='original_matched_random_signature_primary',matched_gene_signature='matched_gene_null',independent_weight_sign_permutation='weight_sign_permutation_actual_mapped_features',age_label_permutation_sensitivity='gtex_age_bracket_label_permutation')
ns<-ns[null_type%in%names(keepnull)]
null_summary<-ns[,.(model_family=fam(model),null_type=unname(keepnull[null_type]),observed_statistic=observed,null_2.5_percentile=null_q025,null_97.5_percentile=null_q975,empirical_one_sided_p)]
null_source<-rbindlist(lapply(seq_len(nrow(ns)),function(i){d<-rt(paste0('null_distributions/stageB_',ns$model[i],'_',ns$null_type[i],'.tsv'));d[,.(model_family=fam(ns$model[i]),null_type=unname(keepnull[ns$null_type[i]]),iteration,null_statistic=statistic,observed_statistic=ns$observed[i])]}))
wt(ns,'tables/reference_panels/reference_four_null_families.tsv')
lob<-rt('tables/leave_one_bracket_batch_out.tsv')[omitted_variable=='age_bracket'];lob[,`:=`(model_family=fam(model),omitted_age_bracket=omitted_value,robust_ci_low=ci_low,robust_ci_high=ci_high)]
write_tsv(lob,prev('results/table_leave_one_age_bracket_out.tsv'))
ood_models<-rt('tables/ood_associations.tsv');ood_models[,model_family:=fam(model)]
# Matched held-out current nested-CV recipe predictions replace historical locked-model rows.
mr<-rt('tables/mouse_reference_predictions.tsv');mw<-na.omit(dcast(mr,sample_id+age_week~model,value.var='mouse_output'))
e<-fread(file.path(dev_by_model[['HVG-ElasticNet']],'HVG-ElasticNet/count_tuning/results/predictions_A.csv'));q<-fread(file.path(dev_by_model[['VersionB-PCA']],'VersionB-PCA/glmnet/count_tuning/results/predictions_A.csv'));held<-merge(e[,.(sample_key,age_week,fold,elasticnet=pred)],q[,.(sample_key,fold,pca_glmnet=pred)],by=c('sample_key','fold'))
wt(held,'tables/reference_panels/current_nested_CV_paired_mouse_predictions.tsv')
hw<-dcast(dall,donor_id+age_midpoint_sensitivity~model_family,value.var='cross_species_mouse_ovarian_clock_score')
agr<-function(label,x,y,age)data.table(dataset=label,n=length(x),pearson_r=cor(x,y),spearman_rho=cor(x,y,method='spearman'),age_residualized_correlation=cor(residuals(lm(x~age)),residuals(lm(y~age))))
agree<-rbindlist(list(agr('mouse_deployment_reference_shared_110',mw[['HVG-ElasticNet']],mw[['VersionB-PCA']],mw$age_week),agr('paired_nested_CV_mouse_predictions',held$elasticnet,held$pca_glmnet,held$age_week),agr('GTEx_193_donors',hw$elasticnet,hw$pca_glmnet,hw$age_midpoint_sensitivity)))
wt(agree,'tables/reference_panels/cross_model_reference_agreement.tsv')
# Decompose the actual adjusted clock slope in its current ln1pCPM representation.
gene_contrib<-list();pcs<-list()
for(nm in names(models)){
 o<-objs[[paste('GTEx',nm,sep='__')]];m<-models[[nm]];d<-o$meta;X<-model.matrix(~age_order+rin_z+ischemia_z+log_library_z,d);proj<-solve(crossprod(X),t(X))[2,];beta<-drop(proj%*%o$imputed);v<-beta*m$coefficients$coefficient
 z<-data.table(model_family=fam(nm),mouse_gene_symbol=m$representation$genes,gene_contribution_to_human_age_slope=v);setorder(z,-gene_contribution_to_human_age_slope)
 z<-z[order(-abs(gene_contribution_to_human_age_slope))];z[,`:=`(absolute_contribution_rank=.I,cumulative_absolute_fraction=cumsum(abs(gene_contribution_to_human_age_slope))/sum(abs(gene_contribution_to_human_age_slope)))];gene_contrib[[nm]]<-z
 stopifnot(abs(sum(v)-drop(proj%*%(o$imputed%*%m$coefficients$coefficient+m$intercept)))<1e-8)
 if(nm=='VersionB-PCA')pcs[[nm]]<-data.table(pc=colnames(o$human_low),calibrated_contribution_to_total_score=drop(proj%*%o$human_low)*m$pc_coefficients)
}
gene_contrib<-rbindlist(gene_contrib);pcs<-rbindlist(pcs);wt(gene_contrib,'tables/reference_panels/gene_contribution_concentration.tsv');wt(pcs,'tables/reference_panels/pc_age_slope_contributions.tsv')
hg<-rt('tables/HPA_gene_effects.tsv');hh<-dcast(hg,human_gene_id_stable~subset,value.var='hpa_rho');hpa_stability<-rbindlist(lapply(names(models),function(nm){z<-merge(cm[cohort=='GTEx'&model==nm&raw_coefficient!=0],hh,by='human_gene_id_stable');z<-z[is.finite(all_seven)&is.finite(adult_only)];data.table(model_family=fam(nm),all_seven_direction_concordance=mean(sign(z$raw_coefficient)==sign(z$all_seven)),adult_only_direction_concordance=mean(sign(z$raw_coefficient)==sign(z$adult_only)),direction_stable_all_vs_adult=mean(sign(z$all_seven)==sign(z$adult_only)),n_active_genes=nrow(z))}))
wt(hpa_stability,'tables/reference_panels/HPA_all_age_adult_stability.tsv')
