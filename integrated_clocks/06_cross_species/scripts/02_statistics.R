source(file.path(dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),fixed=TRUE))),'common.R'))
p<-rt('tables/frozen_predictions.tsv');objs<-readRDS(file.path(out,'inputs/projection_objects.rds'));inp<-readRDS(file.path(out,'inputs/analysis_inputs.rds'))
ass<-list();loo<-list();influence<-list();ood<-list();pcstats<-list();groups<-list();agree<-list()
form<-cross_species_mouse_ovarian_clock_score~age_order+rin_z+ischemia_z+log_library_z
for(nm in names(models)){
 d<-p[cohort=='GTEx'&model==nm]
 formulas<-list(primary_adjusted_ordered=form,age_only_ordered=cross_species_mouse_ovarian_clock_score~age_order,midpoint_decade=cross_species_mouse_ovarian_clock_score~age_decade+rin_z+ischemia_z+log_library_z,historical_RIN_only=cross_species_mouse_ovarian_clock_score~age_decade+rin_z,filtered_input_sensitivity=filtered_input_score_sensitivity~age_order+rin_z+ischemia_z+log_library_z)
 for(an in names(formulas))ass[[paste(nm,an)]]<-hc3(formulas[[an]],d,if(an%in%c('midpoint_decade','historical_RIN_only'))'age_decade' else 'age_order')[,`:=`(cohort='GTEx',model=nm,analysis=an)]
 ass[[paste(nm,'exclude_oldest')]]<-hc3(form,d[age_bracket!='70-79'],'age_order')[,`:=`(cohort='GTEx',model=nm,analysis='exclude_70_79')]
 dw<-copy(d);qq<-quantile(dw$cross_species_mouse_ovarian_clock_score,c(.025,.975));dw[,cross_species_mouse_ovarian_clock_score:=pmin(qq[2],pmax(qq[1],cross_species_mouse_ovarian_clock_score))]
 ass[[paste(nm,'winsor')]]<-hc3(form,dw,'age_order')[,`:=`(cohort='GTEx',model=nm,analysis='winsorize_2.5_97.5')]
 rr<-summary(MASS::rlm(form,data=d,maxit=100))$coefficients['age_order',];ass[[paste(nm,'robust')]]<-data.table(cohort='GTEx',model=nm,analysis='MASS_robust_sensitivity',n=nrow(d),estimate=rr[1],se=rr[2],ci_low=rr[1]-qnorm(.975)*rr[2],ci_high=rr[1]+qnorm(.975)*rr[2],p=2*pnorm(abs(rr[3]),lower.tail=FALSE),status='estimable')
 for(v in c('age_bracket','sequencing_batch'))for(omit in unique(d[[v]])){
  dd<-d[get(v)!=omit];if(nrow(dd)<20||uniqueN(dd$age_order)<3)next
  loo[[paste(nm,v,omit)]]<-hc3(form,dd,'age_order')[,`:=`(cohort='GTEx',model=nm,omitted_variable=v,omitted_value=as.character(omit),n_removed=nrow(d)-nrow(dd))]
 }
 f<-lm(form,d);influence[[nm]]<-data.table(cohort='GTEx',model=nm,sample_id=d$sample_id,cooks_distance=cooks.distance(f),leverage=hatvalues(f),studentized_residual=rstudent(f),influential=cooks.distance(f)>4/nrow(d))
 odforms<-list(ood_vs_age=projection_ood_distance~age_order+rin_z+ischemia_z+log_library_z,score_vs_ood=cross_species_mouse_ovarian_clock_score~projection_ood_distance+rin_z+ischemia_z+log_library_z,age_adjusted_for_ood=cross_species_mouse_ovarian_clock_score~age_order+projection_ood_distance+rin_z+ischemia_z+log_library_z)
 for(an in names(odforms))ood[[paste(nm,an)]]<-hc3(odforms[[an]],d,if(an=='score_vs_ood')'projection_ood_distance' else 'age_order')[,`:=`(model=nm,analysis=an)]
 for(q in c('mouse_ood_q95','mouse_ood_q99')){
  dd<-d[projection_ood_distance<=get(q)];res<-if(nrow(dd)>=20&&uniqueN(dd$age_order)>=3)hc3(form,dd,'age_order') else data.table(n=nrow(dd),estimate=NA_real_,status='not_estimable_no_common_support')
  ass[[paste(nm,q)]]<-res[,`:=`(cohort='GTEx',model=nm,analysis=paste0('common_support_',q))]
 }
}
# Exact Spearman permutation for eight donors, including ties by permuting donor labels.
perms<-function(v){if(length(v)==1)return(matrix(v,nrow=1));do.call(rbind,lapply(seq_along(v),function(i)cbind(v[i],perms(v[-i]))))}
perm8<-perms(1:8);cors<-list();gseloo<-list()
for(cn in unique(p$cohort))for(nm in names(models)){
 d<-p[cohort==cn&model==nm];x<-d$age_year;y<-d$cross_species_mouse_ovarian_clock_score
 ct<-cor.test(x,y,method='pearson');sp<-suppressWarnings(cor.test(x,y,method='spearman',exact=FALSE));ep<-NA_real_
 if(cn=='GSE202601'){
  xr<-zsafe(rank(x));yr<-zsafe(rank(y));null<-rowSums(matrix(xr[perm8],nrow=nrow(perm8))*rep(yr,each=nrow(perm8)))/(length(x)-1);ep<-mean(abs(null)>=abs(unname(sp$estimate))-1e-12)
  ass[[paste(cn,nm)]]<-hc3(cross_species_mouse_ovarian_clock_score~age_decade,d,'age_decade')[,`:=`(cohort=cn,model=nm,analysis='exact_age_per_decade_HC3')]
  for(i in seq_len(nrow(d)))gseloo[[paste(nm,i)]]<-data.table(model=nm,omitted_sample=d$sample_id[i],n=7,pearson=cor(x[-i],y[-i]),spearman=cor(x[-i],y[-i],method='spearman'))
  tt<-t.test(y[d$age_group=='reproductively_aged'],y[d$age_group=='young']);groups[[nm]]<-data.table(model=nm,n_young=sum(d$age_group=='young'),n_aged=sum(d$age_group=='reproductively_aged'),effect_aged_minus_young=diff(rev(tt$estimate)),ci_low=tt$conf.int[1],ci_high=tt$conf.int[2],p=tt$p.value,test='two-sided unpaired Welch t-test; descriptive')
 }
 cors[[paste(cn,nm)]]<-data.table(cohort=cn,model=nm,n=nrow(d),pearson=unname(ct$estimate),pearson_p=ct$p.value,pearson_ci_low=ct$conf.int[1],pearson_ci_high=ct$conf.int[2],spearman=unname(sp$estimate),spearman_p=if(cn=='GSE202601')ep else sp$p.value,spearman_method=if(cn=='GSE202601')'exact enumeration of all 40320 donor permutations; two-sided' else 'asymptotic two-sided; tied age brackets')
 w<-dcast(p[cohort==cn],sample_id~model,value.var='cross_species_mouse_ovarian_clock_score');a<-w[[2]];b<-w[[3]];r<-cor(a,b);br<-replicate(2000,{j<-sample.int(length(a),replace=TRUE);cor(a[j],b[j])})
 agree[[cn]]<-data.table(cohort=cn,n=nrow(w),pearson=r,spearman=cor(a,b,method='spearman'),bootstrap_ci_low=quantile(br,.025,na.rm=TRUE),bootstrap_ci_high=quantile(br,.975,na.rm=TRUE),bootstrap_draws=2000)
}
ass<-rbindlist(ass,fill=TRUE);ass[analysis=='primary_adjusted_ordered',q_primary_two_models:=p.adjust(p,'BH')];ass[cohort=='GSE202601',q_two_models:=p.adjust(p,'BH')]
wt(ass,'tables/age_association_sensitivities.tsv');wt(rbindlist(cors),'tables/age_correlations.tsv');wt(rbindlist(loo),'tables/leave_one_bracket_batch_out.tsv');wt(rbindlist(influence),'tables/donor_influence.tsv');wt(rbindlist(ood),'tables/ood_associations.tsv');wt(rbindlist(gseloo),'tables/GSE202601_leave_one_donor_out.tsv');wt(rbindlist(groups),'tables/GSE202601_age_group_comparison.tsv')
mp<-rt('tables/mouse_reference_predictions.tsv');mw<-na.omit(dcast(mp,sample_id~model,value.var='mouse_output'));agree[['mouse_reference']]<-data.table(cohort='mouse_reference_in_sample',n=nrow(mw),pearson=cor(mw[[2]],mw[[3]]),spearman=cor(mw[[2]],mw[[3]],method='spearman'))
wt(rbindlist(agree,fill=TRUE),'tables/cross_model_agreement.tsv')
pc<-rt('tables/pca_scores_contributions.tsv')
for(cn in unique(pc$cohort))for(pp in unique(pc$PC)){
 d<-merge(pc[cohort==cn&PC==pp],unique(p[cohort==cn,.(sample_id,age_order,rin_z,ischemia_z,log_library_z,age_decade)]),by='sample_id')
 f<-if(cn=='GTEx')contribution~age_order+rin_z+ischemia_z+log_library_z else contribution~age_decade
 pcstats[[paste(cn,pp)]]<-hc3(f,d,if(cn=='GTEx')'age_order' else 'age_decade')[,`:=`(cohort=cn,PC=pp,coefficient=d$coefficient[1])]
}
wt(rbindlist(pcstats),'tables/pca_contribution_age_trends.tsv')
# HPA is exploratory Stage A only: recompute all-seven and adult-only direction.
hpa<-readRDS(file.path(old,'derived_data/hpa/hpa_ovary_prepared.rds'));master<-inp$master;hstats<-list();hgene<-list()
for(subset in c('all_seven','adult_only')){
 ix<-if(subset=='all_seven')seq_len(nrow(hpa$metadata)) else which(hpa$metadata$age>=18)
 rho<-apply(hpa$log2_tpm_plus1[,ix,drop=FALSE],1,function(z)suppressWarnings(cor(z,hpa$metadata$age[ix],method='spearman')))
 stab<-vapply(seq_len(nrow(hpa$log2_tpm_plus1)),function(i){z<-hpa$log2_tpm_plus1[i,ix];s<-vapply(seq_along(ix),function(j)sign(suppressWarnings(cor(z[-j],hpa$metadata$age[ix][-j],method='spearman'))),numeric(1));mean(s==sign(rho[i]),na.rm=TRUE)},numeric(1))
 d<-copy(master);d[,`:=`(hpa_rho=rho[human_gene_id_stable],hpa_loo_sign_stability=stab[match(human_gene_id_stable,rownames(hpa$log2_tpm_plus1))],subset=subset)]
 hgene[[subset]]<-d[,.(mouse_gene_symbol,human_gene_id_stable,subset,mouse_age_t,hpa_rho,hpa_loo_sign_stability)]
 for(nm in names(models)){
  prefix<-if(nm=='HVG-ElasticNet')'hvg' else 'pca';dd<-d[get(paste0(prefix,'_selected'))&is.finite(hpa_rho)];w<-dd[[paste0(prefix,'_weight')]]
  hstats[[paste(subset,nm)]]<-data.table(subset=subset,model=nm,n_donors=length(ix),n_genes=nrow(dd),mouse_hpa_spearman=cor(dd$mouse_age_t,dd$hpa_rho,method='spearman'),sign_concordance=mean(sign(dd$mouse_age_t)==sign(dd$hpa_rho)),weight_direction_concordance=weighted.mean(sign(w)==sign(dd$hpa_rho),abs(w)),mean_loo_sign_stability=mean(dd$hpa_loo_sign_stability,na.rm=TRUE),exploratory=TRUE)
 }
}
wt(rbindlist(hstats),'tables/HPA_exploratory_direction.tsv');wt(rbindlist(hgene),'tables/HPA_gene_effects.tsv');wt(hpa$metadata,'inputs/HPA_metadata.tsv')
message('PASS: associations, exact n=8 tests, influence, OOD, PCA and HPA statistics.')
# Orthogonal distance complements low-dimensional Mahalanobis distance.
residuals_out<-list();coverage<-rt('tables/model_gene_coverage.tsv')
for(key in names(objs)){
 o<-objs[[key]];cn<-o$meta$cohort[1];nm<-o$map$model[1]
 hr<-sqrt(pmax(0,rowSums(o$z^2)-rowSums(o$human_low^2))/ncol(o$z));mr<-sqrt(pmax(0,rowSums(o$reference_z^2)-rowSums(o$reference_low^2))/ncol(o$reference_z));th<-quantile(mr,.99)
 residuals_out[[key]]<-data.table(cohort=cn,model=nm,sample_id=o$meta$sample_id,orthogonal_residual_rms=hr,mouse_reference_q99=th,outside_reference_q99=hr>th,all_feature_rms=sqrt(rowMeans(o$z^2)))
}
wt(rbindlist(residuals_out),'tables/orthogonal_distance_qc.tsv')
cs<-coverage[,.(n_input=.N,n_one2one=sum(one2one_mapped),n_measured=sum(measured),n_filter_pass=sum(pass_historical_expression_filter),raw_coverage=mean(measured),absolute_weight_coverage=sum(abs(raw_coefficient[measured]))/sum(abs(raw_coefficient)),positive_weight_coverage=sum(pmax(raw_coefficient[measured],0))/sum(pmax(raw_coefficient,0)),negative_weight_coverage=sum(abs(pmin(raw_coefficient[measured],0)))/sum(abs(pmin(raw_coefficient,0)))),by=.(cohort,model)]
wt(cs,'tables/coverage_summary.tsv')
wt(p[,.(n=.N,mean_score=mean(cross_species_mouse_ovarian_clock_score),sd_score=sd(cross_species_mouse_ovarian_clock_score),se_score=sd(cross_species_mouse_ovarian_clock_score)/sqrt(.N)),by=.(cohort,model,age_bracket,age_group)],'tables/age_group_score_summary.tsv')
wt(dcast(p,cohort+sample_id+donor_id+age_year~model,value.var=c('raw_linear_predictor','cross_species_mouse_ovarian_clock_score','projection_ood_distance')),'tables/frozen_predictions_wide.tsv')
