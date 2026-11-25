#!/usr/bin/env Rscript
# Reconstruct reference analyses from frozen predictions. No model fitting or calibration.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
setDTthreads(1)
a<-commandArgs(FALSE);here<-dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE)))
out<-dirname(here);dev<-file.path(dirname(out),'model_development');root<-dirname(dirname(dirname(out)))
oldroot<-file.path(root,'bulk-pseudobulk_clock/results/external_projection')
oldscript<-file.path(root,'bulk-pseudobulk_clock/scripts')
release<-file.path(dirname(out),'lockedmodel')
rd<-function(f)fread(file.path(out,f));oldrd<-function(f)fread(file.path(oldroot,f))
wr<-function(x,f){if(startsWith(f,'simplified_classes/'))return(invisible(NULL));dir.create(dirname(file.path(out,f)),recursive=TRUE,showWarnings=FALSE);fwrite(x,file.path(out,f),sep='\t',na='NA')}
models<-c('HVG-ElasticNet','VersionB-PCA/glmnet')
model_dir<-c('HVG-ElasticNet'='HVG-ElasticNet-exclude-single-age','VersionB-PCA/glmnet'='VersionB-PCA-glmnet')
p<-rd('external_projection_predictions.tsv');p[,model:=factor(model,levels=models)]
tt<-rd('tables/external_projection_contrast_statistics.tsv');mem<-rd('tables/external_projection_contrast_membership.tsv')
orig<-unique(oldrd('external_projection_predictions_all.tsv')[,.(sample_key,projection_group,projection_class,projection_role,comparison_id,natural_aging_plot_group,primary_external_projection_include)])
simple<-unique(oldrd('simplified_classes/external_projection_predictions_all_simplified.tsv')[,.(sample_key,projection_group,projection_class,projection_role)])
stopifnot(!anyDuplicated(orig$sample_key),!anyDuplicated(simple$sample_key),all(p$sample_key%in%orig$sample_key))
for(nm in setdiff(names(orig),'sample_key'))p[,(nm):=orig[[nm]][match(sample_key,orig$sample_key)]]
p[,`:=`(simple_group=simple$projection_group[match(sample_key,simple$sample_key)],simple_class=simple$projection_class[match(sample_key,simple$sample_key)],simple_role=simple$projection_role[match(sample_key,simple$sample_key)])]
newlabels<-unique(oldrd('diagnostics/newly added with pvalue/newly_added7_relative_delta_plot_source_data.tsv')[,.(sample_key,plot_group,display_projection_class)])
stopifnot(!anyDuplicated(newlabels$sample_key))
p[,`:=`(new_group=newlabels$plot_group[match(sample_key,newlabels$sample_key)],new_class=newlabels$display_projection_class[match(sample_key,newlabels$sample_key)])]
special<-rbindlist(lapply(c('PRJNA863443','PRJNA1002222'),function(ds)oldrd(paste0(ds,'/',tolower(ds),'_projection_metadata.tsv'))[,.(sample_key,projection_group,projection_class,natural_aging_plot_group,comparison_id)]))
p[,`:=`(special_group=special$projection_group[match(sample_key,special$sample_key)],special_class=special$projection_class[match(sample_key,special$sample_key)])]
wr(unique(p[,.(sample_key,dataset_cohort,projection_group,projection_class,simple_group,simple_class,new_group,new_class,special_group,special_class,natural_aging_plot_group,comparison_id)]),'provenance/reference_annotation_dictionary.tsv')
# Parse only explicit data definitions, never execute historical projection pipelines.
parse_assignment<-function(path,name){expr<-parse(path);hit<-Filter(function(e)is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name(name)),as.list(expr));stopifnot(length(hit)==1);eval(hit[[1]][[3]])}
cm<-parse_assignment(file.path(oldscript,'06e_raw_delta_age_dataset_t_tests.R'),'contrast_map')
orders<-parse_assignment(file.path(oldscript,'06e_raw_delta_age_dataset_t_tests.R'),'dataset_plot_group_levels')
newspec<-parse_assignment(file.path(oldscript,'06d_plot_newly_added_relative_delta_ttests.R'),'contrast_spec')
orders$PRJNA1002222<-c('WT_2M','WT_8M','Cd38KO_8M','WT_12M','Cd38KO_12M')
orders$CRA003645<-c('0-day developmental ovary','3-week prepubertal ovary')
orders$GSE228975<-c('Srsf1 WT','Srsf1 cKO')
# User-specified display orders: controls/references on the left.
orders$CRA017492<-c('young control','aged water control','aged NMN')
orders$GSE98497<-c('WT control','Mtor oocyte-conditional KO')
orders$GSE106339<-c('3-week WT','3-week GOAT KO')
orders$GSE128240<-c('same-study control','cyclophosphamide')
orders$GSE273724<-c('reference control','adult 27-28 week vehicle','adult 27-28 week plus HEP14','aged 73-74 week vehicle','aged 73-74 week plus HEP14')
order_for<-function(d,variant='raw'){
 ds<-unique(d$dataset_folder);present<-unique(as.character(d$plot_group))
 if(length(ds)!=1)return(sort(present))
 if(variant=='special')return(sort(present))
 if(variant=='natural')return(intersect(c('young_13w-17.4w','mid_76w','old_86w-87w','2M'),present))
 if(variant=='simplified'&&ds!='GSE235021')return(sort(present))
 # Root 04/05 figures use the same curated reference->target order as the t-test folder.
 if(variant=='new'||ds=='GSE235021'){
  z<-newspec[dataset_folder==ds];o<-unique(as.vector(rbind(z$control_label,z$experimental_label)))
 }else if(ds%in%names(orders))o<-orders[[ds]]
 else if(variant%in%c('simplified','special'))o<-sort(present)
 else {z<-cm[dataset_folder==ds];o<-if(ds=='GSE273724')as.vector(t(as.matrix(z[,.(control_group,condition_group)])))else c(z$control_group,z$condition_group)}
 unique(c(o[o%in%present],sort(setdiff(present,o))))
}
annotate_data<-function(d,variant='raw'){
 d<-copy(d);d[,`:=`(plot_group=simple_group,plot_class=simple_class,plot_role=simple_role)]
 # The historical diagnostics explicitly prefer simplified metadata when available.
 if(variant=='root')d[,`:=`(plot_group=simple_group,plot_class=simple_class)]
 if(variant=='special')d[,`:=`(plot_group=special_group,plot_class=special_class)]
 if(variant=='new'||any(d$dataset_folder=='GSE235021'))d[!is.na(new_group),`:=`(plot_group=new_group,plot_class=new_class)]
 if(variant=='natural')d[,`:=`(plot_group=natural_aging_plot_group,plot_class=natural_aging_plot_group)]
 if(variant%in%c('raw','root'))d[dataset_folder=='PRJNA1002222',`:=`(plot_group=special_group,plot_class=special_class)]
 # An old single OOD label cannot identify the two prespecified Srsf1 test arms.
 d[dataset_folder=='GSE228975',plot_group:=fifelse(grepl('cKO',genotype,ignore.case=TRUE),'Srsf1 cKO','Srsf1 WT')]
 if(!variant%in%c('simplified','special','natural'))d[plot_class=='stress',plot_class:='stress_model']
 stopifnot(!anyNA(d$plot_group),!anyNA(d$plot_class));d
}
# Restore the additional reference-only KO rescue versus WT comparison separately.
extra_keys<-p[model==models[1]&dataset_folder=='GSE279573',.(sample_key,simple_group,modality,dataset_cohort)]
extra<-rbindlist(lapply(c('WT chow fertile','KO vegetable-oil-free rescue'),function(g)extra_keys[simple_group==g,.(sample_key,modality,dataset_cohort,analysis_contrast_id=paste0('GSE279573_KO_oilfree_vs_WT_chow_reference_addition__',modality),contrast_id='GSE279573_KO_oilfree_vs_WT_chow_reference_addition',contrast_role=if(g=='WT chow fertile')'control'else'experimental')]))
memall<-rbindlist(list(mem,extra),fill=TRUE)
exstats<-rbindlist(lapply(models,function(mo){d<-p[model==mo];a<-d[sample_key%in%extra[contrast_role=='experimental',sample_key],delta_age_week];b<-d[sample_key%in%extra[contrast_role=='control',sample_key],delta_age_week];t<-t.test(a,b);data.table(dataset_folder='GSE279573',dataset_cohort=extra$dataset_cohort[1],modality=extra$modality[1],contrast_id=extra$contrast_id[1],analysis_contrast_id=extra$analysis_contrast_id[1],model=mo,n_control=length(b),n_experimental=length(a),control_mean_delta=mean(b),experimental_mean_delta=mean(a),effect_weeks=mean(a)-mean(b),ci95_low=t$conf.int[1],ci95_high=t$conf.int[2],t_statistic=unname(t$statistic),df=unname(t$parameter),welch_p=t$p.value,status='tested',placement_only=FALSE,contrast_name='GSE279573 KO oil-free rescue versus WT chow (reference-figure addition)',analysis_note='Additional contrast explicitly present in original reference figures; not in the 132-test prespecified family.')}))
ttall<-rbindlist(list(copy(tt)[,analysis_family:='prespecified'],exstats[,analysis_family:='reference_figure_addition']),fill=TRUE)
ttall[,bh_q_expanded_90:=p.adjust(welch_p,'BH')]
wr(exstats,'tables/additional_reference_contrast_statistics.tsv');wr(ttall,'tables/all_reference_contrast_statistics.tsv');wr(memall,'tables/all_reference_contrast_membership.tsv')
format_p<-function(v)ifelse(v<1e-4,'****',ifelse(v<.001,'***',ifelse(v<.01,'**',ifelse(v<.05,'*',sprintf('p = %.3f',v)))))
ttall[,`:=`(annotation=format_p(welch_p),model=factor(model,levels=models))]
test_labels<-function(d,tests=ttall){
 z<-merge(memall[,.(sample_key,analysis_contrast_id,contrast_role)],d[,.(sample_key,model,plot_group)],by='sample_key',allow.cartesian=TRUE)
 z[,role:=fifelse(contrast_role%in%c('control','placement_reference'),'Reference','Target')]
 z<-z[,.(label=paste(unique(plot_group),collapse=';'),n=uniqueN(sample_key),ng=uniqueN(plot_group)),by=.(analysis_contrast_id,model,role)]
 z<-dcast(z,analysis_contrast_id+model~role,value.var=c('label','n','ng'))
 if(!all(c('label_Reference','label_Target')%in%names(z)))return(tests[0][,`:=`(Reference=character(),Target=character())])
 z<-merge(copy(tests),z,by=c('analysis_contrast_id','model'))
 z<-z[n_Reference==n_control&n_Target==n_experimental&ng_Reference==1&ng_Target==1]
 z[,`:=`(Reference=label_Reference,Target=label_Target)]
 z
}
rawsource<-annotate_data(p);wr(rawsource[,.(sample_key,dataset_cohort,modality,model,plot_group,plot_class,age_week,predicted_age_week,delta_age_week)],'diagnostics/previous_style_plot_source.tsv')
wr(test_labels(rawsource,ttall[analysis_family=='prespecified']),'diagnostics/previous_style_annotation_statistics.tsv')

# Training projections and distance diagnostics use the same frozen representation.
source(file.path(release,'scripts/core.R'))
train<-readRDS(file.path(dev,'shared_inputs/prepared_expression.rds'));external<-readRDS(file.path(out,'inputs/model_gene_expression.rds'))$expression
refrows<-list();distrows<-list();refpc<-NULL
for(mo in models){
 model<-readRDS(file.path(release,model_dir[[mo]],'model.rds'))
 za<-project_representation(model$representation,external)$x
 zr<-project_representation(model$representation,train$expression[,model$training_keys,drop=FALSE])$x
 if(!is.null(model$pc_coefficients)){
  za<-za[,seq_along(model$pc_coefficients),drop=FALSE];zr<-zr[,seq_along(model$pc_coefficients),drop=FALSE]
  covar<-cov(zr);diag(covar)<-diag(covar)+1e-6
  da<-mahalanobis(za,colMeans(zr),covar);dr<-mahalanobis(zr,colMeans(zr),covar);method<-'squared Mahalanobis in 10 frozen PC inputs; covariance ridge 1e-6'
  refpc<-data.table(sample_key=rownames(zr),model=mo,zr)
 }else{da<-sqrt(rowSums(sweep(za,2,colMeans(zr),'-')^2));dr<-sqrt(rowSums(sweep(zr,2,colMeans(zr),'-')^2));method<-'Euclidean distance from training centroid in frozen standardized input genes'}
 percentile<-vapply(da,function(z)mean(dr<=z),numeric(1))
 distrows[[mo]]<-data.table(sample_key=rownames(za),model=mo,projection_distance=da,projection_distance_percentile=percentile,OOD_flag=ifelse(percentile>.99,'OOD_gt_99th',ifelse(percentile>.95,'OOD_gt_95th','in_reference_range')),distance_method=method)
 rr<-predict_model(model,train$expression[,model$training_keys,drop=FALSE]);rr[,`:=`(model=mo,age_week=train$metadata$age_week[match(sample_key,train$metadata$sample_key)],projection_distance=dr)]
 rr[,`:=`(predicted_age_week=pred,delta_age_week=pred-age_week)];refrows[[mo]]<-rr
}
ref<-rbindlist(refrows);ref[,model:=factor(model,levels=models)]
refstats<-ref[,.(n_reference_samples=.N,locked_mae_week=mean(abs(delta_age_week)),locked_r2=summary(lm(predicted_age_week~age_week))$r.squared,fit_intercept=coef(lm(predicted_age_week~age_week))[1],fit_slope=coef(lm(predicted_age_week~age_week))[2]),by=model]
refstats[,locked_metric_label:=sprintf('Deployment training fit (in-sample; n = %d)\nMAE = %.1f weeks; R² = %.2f',n_reference_samples,locked_mae_week,locked_r2)]
p<-merge(p,rbindlist(distrows),by=c('sample_key','model'));p[,model:=factor(model,levels=models)]
p[,`:=`(low_coverage_flag=gene_coverage_fraction<.8,OOD_flag_bool=projection_distance_percentile>.95)]
wr(p,'external_projection_predictions_all.tsv');wr(ref,'tables/frozen_deployment_reference_predictions.tsv');wr(refstats,'tables/frozen_deployment_reference_fit_stats.tsv');wr(refpc,'reference_pca_projection_scores.tsv')
wr(p,'diagnostics/projection_qc_sample_level.tsv');wr(p[,.(sample_key,model,negative_predicted_age,outside_training_age_range,low_coverage_flag,OOD_flag,locked_training_sample)],'external_projection_qc_flags.tsv')
qc<-p[,.(n=.N,min_gene_coverage=min(gene_coverage_fraction),mean_gene_coverage=mean(gene_coverage_fraction),n_OOD_95=sum(OOD_flag_bool),n_OOD_99=sum(projection_distance_percentile>.99),n_negative_age=sum(negative_predicted_age),n_age_extrapolation=sum(outside_training_age_range),n_development_overlap=sum(locked_training_sample)),by=.(dataset_folder,dataset_cohort,modality,model)]
wr(qc,'tables/projection_qc_by_cohort.tsv');wr(qc,'external_projection_gene_coverage.tsv')
wr(p[,.(n_prediction_rows=.N,n_samples=uniqueN(sample_key),mean_gene_coverage=mean(gene_coverage_fraction),min_gene_coverage=min(gene_coverage_fraction),n_low_coverage=sum(low_coverage_flag),n_OOD=sum(OOD_flag_bool),OOD_fraction=mean(OOD_flag_bool)),by=model],'diagnostics/projection_qc_coverage_ood_summary.tsv')
groups<-annotate_data(p)[,.(n=.N,mean_predicted_age=mean(predicted_age_week),mean_delta_age=mean(delta_age_week),median_delta_age=median(delta_age_week),sd_delta_age=sd(delta_age_week),mean_gene_coverage=mean(gene_coverage_fraction),OOD_fraction=mean(OOD_flag_bool),n_training_overlap=sum(locked_training_sample)),by=.(dataset_folder,dataset_cohort,modality,model,plot_group,plot_class)]
wr(groups,'external_projection_group_summary.tsv');wr(groups,'diagnostics/projection_group_summary.tsv')
wr(p[,.(n=.N,n_samples=uniqueN(sample_key),mean_delta_age=mean(delta_age_week),median_delta_age=median(delta_age_week),sd_delta_age=sd(delta_age_week),min_delta_age=min(delta_age_week),max_delta_age=max(delta_age_week)),by=model],'diagnostics/raw_delta_age_summary_by_clock.tsv')
rel<-merge(memall[,.(sample_key,analysis_contrast_id,contrast_id,contrast_role)],p,by='sample_key',allow.cartesian=TRUE)
rel[,role:=fifelse(contrast_role%in%c('control','placement_reference'),'Reference','Target')]
rel[,control_mean_delta:=mean(delta_age_week[role=='Reference']),by=.(analysis_contrast_id,model)]
rel[,relative_delta_age_week:=delta_age_week-control_mean_delta]
rel<-merge(rel,ttall[,.(analysis_contrast_id,model,placement_only,analysis_family,source_confounding_flag)],by=c('analysis_contrast_id','model'))
wr(rel,'external_projection_relative_to_controls.tsv')
no_relative<-p[!sample_key%in%rel$sample_key];no_relative[,`:=`(relative_delta_age_week=NA_real_,relative_control_status='no prespecified control contrast')]
wr(rbindlist(list(copy(rel)[,relative_control_status:=ifelse(placement_only,'placement reference','matched within contrast and modality')],no_relative),fill=TRUE),'external_projection_relative_to_controls_all.tsv')
controlaudit<-rel[,.(n_control=uniqueN(sample_key[role=='Reference']),n_target=uniqueN(sample_key[role=='Target']),control_sample_keys=paste(sort(unique(sample_key[role=='Reference'])),collapse=';'),target_sample_keys=paste(sort(unique(sample_key[role=='Target'])),collapse=';'),n_modalities=uniqueN(modality),n_datasets=uniqueN(dataset_folder),centered_control_mean=mean(relative_delta_age_week[role=='Reference']),placement_only=any(placement_only),any_training_overlap=any(locked_training_sample),age_difference=mean(age_week[role=='Target'])-mean(age_week[role=='Reference'])),by=.(dataset_folder,dataset_cohort,analysis_contrast_id,model)]
stopifnot(all(controlaudit$n_modalities==1),all(controlaudit$n_datasets==1),all(abs(controlaudit$centered_control_mean)<1e-9))
wr(controlaudit,'external_projection_control_matching_audit.tsv');wr(controlaudit,'diagnostics/matched_control_integrity_audit.tsv')
wr(controlaudit[placement_only|abs(age_difference)>1|any_training_overlap],'external_projection_unreliable_relative_controls_audit.tsv')
wr(rel[role=='Reference'],'external_projection_reference_controls.tsv');wr(unique(rel[,.(sample_key,dataset_cohort,analysis_contrast_id,role)]),'external_projection_control_pool_candidates.tsv')
status<-unique(p[,.(sample_key,dataset_folder,dataset_cohort,modality)])
status[,relative_control_status:=ifelse(sample_key%in%rel[placement_only==FALSE,sample_key],'matched control contrast',ifelse(sample_key%in%rel$sample_key,'placement reference only','no matched control'))]
availability<-status[,.(n_profiles=.N),by=.(dataset_folder,dataset_cohort,modality,relative_control_status)]
wr(status,'tables/sample_relative_control_status.tsv');wr(availability,'external_projection_dataset_control_availability.tsv');wr(availability,'diagnostics/dataset_relative_control_status.tsv')
relsummary<-annotate_data(rel)[,.(n=uniqueN(sample_key),mean_relative_delta=mean(relative_delta_age_week),sd_relative_delta=sd(relative_delta_age_week)),by=.(dataset_folder,dataset_cohort,analysis_contrast_id,model,role,plot_group,plot_class)]
wr(relsummary,'diagnostics/matched_control_relative_delta_summary.tsv')
wide<-dcast(p,sample_key+dataset_folder+dataset_cohort+modality+simple_group+simple_class~model,value.var='delta_age_week')
pairs<-combn(models,2,simplify=FALSE)
pairdata<-rbindlist(lapply(pairs,function(mm)wide[,.(sample_key,dataset_folder,dataset_cohort,modality,plot_group=simple_group,plot_class=simple_class,model_x=mm[1],model_y=mm[2],pair=paste(mm,collapse=' vs '),x=get(mm[1]),y=get(mm[2]),difference=get(mm[2])-get(mm[1]),same_direction=sign(get(mm[1]))==sign(get(mm[2])))]))
corsummary<-function(d,by)d[,.(n=.N,pearson_r=if(.N>2&&sd(x)>0&&sd(y)>0)cor(x,y)else NA_real_,spearman_r=if(.N>2&&sd(x)>0&&sd(y)>0)cor(x,y,method='spearman')else NA_real_,mean_abs_difference=mean(abs(difference)),median_abs_difference=median(abs(difference)),mean_difference=mean(difference),fraction_same_direction=mean(same_direction)),by=by]
concord<-corsummary(pairdata,c('pair','model_x','model_y'));cg<-corsummary(pairdata,c('dataset_cohort','plot_group','pair','model_x','model_y'))
wr(pairdata,'diagnostics/model_concordance_sample_level.tsv');wr(concord,'external_projection_model_concordance.tsv');wr(concord,'diagnostics/model_concordance_overall_summary.tsv');wr(cg,'diagnostics/model_concordance_group_summary.tsv')
wide[,clock_range:=do.call(pmax,.SD)-do.call(pmin,.SD),.SDcols=models]
setorder(wide,-clock_range);wr(wide,'tables/two_clock_sample_disagreement.tsv')
worst<-rbindlist(list(
 p[order(-abs(delta_age_week)),head(.SD,10),by=model][,`:=`(discrepancy_category='largest_absolute_delta_age',ranking_value=abs(delta_age_week))],
 pairdata[order(-abs(difference)),head(.SD,10),by=pair][,`:=`(model=pair,discrepancy_category='largest_pairwise_clock_disagreement',ranking_value=abs(difference))],
 p[order(gene_coverage_fraction),head(.SD,10),by=model][,`:=`(discrepancy_category='lowest_gene_coverage',ranking_value=gene_coverage_fraction)],
 p[order(-projection_distance),head(.SD,10),by=model][,`:=`(discrepancy_category='highest_projection_distance',ranking_value=projection_distance)]
),fill=TRUE)
wr(worst,'diagnostics/worst_discrepant_samples.tsv')
agreement<-ttall[,.(n_models=.N,all_clocks_agree=uniqueN(sign(effect_weeks))==1,n_nominal_significant=sum(welch_p<.05),n_min=min(n_control,n_experimental),any_placement=any(placement_only),any_training_overlap=any(n_training_controls+n_training_targets>0,na.rm=TRUE)),by=.(dataset_folder,dataset_cohort,analysis_contrast_id)]
cq<-rel[role=='Target',.(OOD_fraction=mean(OOD_flag_bool),low_coverage_fraction=mean(low_coverage_flag),any_age_extrapolation=any(outside_training_age_range)),by=analysis_contrast_id]
claims<-merge(agreement,cq,by='analysis_contrast_id')
claims[,claim_grade:=fifelse(OOD_fraction>.5|low_coverage_fraction>.5,'exclude_from_claims',fifelse(any_placement|any_training_overlap|any_age_extrapolation,'descriptive_only',fifelse(n_min>=3&all_clocks_agree&OOD_fraction<=.25,'strong_projection_signal',fifelse(n_min>=2&OOD_fraction<=.5,'moderate_projection_signal','descriptive_only'))))]
claims[,grading_note:='Heuristic projection-QC grade adapted to both locked clocks; not a test of biological rejuvenation. Development overlap, placement or age extrapolation caps interpretation at descriptive.']
wr(claims,'diagnostics/projection_claim_grading.tsv')
pc<-merge(rd('tables/pca_projection_scores.tsv'),unique(p[,.(sample_key,dataset_folder,dataset_cohort,modality,age_week,projection_class,simple_group,simple_class,treatment)]),by='sample_key')
cc<-merge(rd('tables/pca_projection_pc_contributions.tsv'),unique(p[,.(sample_key,dataset_folder,dataset_cohort,modality,age_week,simple_group,simple_class,treatment)]),by='sample_key')
active<-grep('^PC[0-9]+$',names(cc),value=TRUE);active<-active[vapply(cc[,..active],function(x)any(abs(x)>1e-12),logical(1))]
cclong<-melt(cc,id.vars=setdiff(names(cc),grep('^PC[0-9]+$',names(cc),value=TRUE)),measure.vars=active,variable.name='PC',value.name='pc_contribution_week')
heat<-cclong[,.(contribution=mean(pc_contribution_week),sd_contribution=sd(pc_contribution_week),n=.N),by=.(dataset_folder,dataset_cohort,modality,simple_group,PC)]
wr(pc,'pca_projection_scores.tsv');wr(cclong,'pca_projection_pc_contributions.tsv');wr(heat,'diagnostics/pca_pc_contribution_group_summary.tsv')

# Figure contract: frozen projections, within-cohort contrasts and deployment QC.
# R quantitative grids, original grey strips/theme_bw/palettes and exact folder labels.
# Vector PDF and 300-dpi PNG; all point and bracket sources are exported below.
# One fixed colour per projection class across every figure.
class_pal<-c(reference_control='#2A9D8F',intervention='#0072B2',genetic_model='#D55E00',diet='#E69F00',stress_model='#CC79A7',stress='#CC79A7',microbiome='#8C564B',developmental_OOD='#7F7F7F')
class_scale<-scale_colour_manual(values=class_pal)
theme_set(theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold'),strip.text=element_text(face='bold'),legend.position='bottom',axis.text.x=element_text(angle=35,hjust=1),plot.margin=margin(12,18,12,12)))
exports<-list();plot_sources<-list();ann_audit<-list();order_audit<-list();plot_objects<-list()
figure_filter<-Sys.getenv('FIGURE_FILTER',unset='')
stats_caption<-'Boxes: median, IQR, whiskers to 1.5 × IQR; points: sample profiles. Two-sided Welch t-tests, unadjusted p.\n* p < 0.05; ** < 0.01; *** < 0.001; **** < 0.0001; otherwise numerical p. All ages in weeks.'
# Delta-age figures share one physical y-axis height and one width per x group.
PANEL_H<-2.4;COL_W<-0.6;X_PAD<-1.4
# Per-class x spacing of the combined raw-delta/tAge figure, panel a (1603 px at 300 dpi over 7 classes + 2 x 0.6 padding).
COMBINED_COL_W<-1603/300/7.2
# Slightly narrower groups for the two wide multi-group plots (PRJNA1002222, GSE235021).
NARROW_COL_W<-0.55;narrow_col_ds<-c('PRJNA1002222','GSE235021')
fixed_panels<-function(g,ncat,nrow){g+theme(panel.widths=unit((ncat+X_PAD)*COL_W,'in'),panel.heights=unit(rep(PANEL_H,nrow),'in'))}
# Wrap long contrast names after underscores too, so narrow fixed-width strips do not clip.
wrap_strip<-function(width=20)function(x)vapply(x,function(v){tok<-regmatches(v,gregexpr('[^_ ]+[_ ]?',v))[[1]];out<-'';line<-'';for(t in tok){if(nchar(line)+nchar(t)>width&&nchar(line)){out<-paste0(out,sub(' $','',line),'\n');line<-''};line<-paste0(line,t)};paste0(out,line)},character(1),USE.NAMES=FALSE)
fig_size<-function(g,minw=7.5){
 pdf(NULL);on.exit(dev.off())
 gt<-if(inherits(g,'patchwork'))patchwork::patchworkGrob(g)else ggplotGrob(g)
 c(w=max(minw,grid::convertWidth(sum(gt$widths),'in',TRUE)+0.2),h=grid::convertHeight(sum(gt$heights),'in',TRUE)+0.2)
}
saveplot<-function(g,stem,w=14,h=7,source=NULL){
 if(startsWith(stem,'simplified_classes/'))return(invisible(g))
 if(is.na(w)||is.na(h)){sz<-fig_size(g);w<-sz[['w']];h<-sz[['h']]}
 dir.create(dirname(file.path(out,stem)),recursive=TRUE,showWarnings=FALSE)
 render<-!nzchar(figure_filter)||grepl(figure_filter,stem)||!file.exists(file.path(out,paste0(stem,'.pdf')))||!file.exists(file.path(out,paste0(stem,'.png')))
 if(render){ggsave(file.path(out,paste0(stem,'.pdf')),g,width=w,height=h,device=cairo_pdf,limitsize=FALSE)
 ggsave(file.path(out,paste0(stem,'.png')),g,width=w,height=h,dpi=300,bg='white',limitsize=FALSE)}
 exports[[stem]]<<-data.table(stem=stem,width_inches=w,height_inches=h)
 if(!is.null(source))plot_sources[[stem]]<<-copy(source)[,figure_stem:=stem]
 if(render)message('Rendered ',stem);invisible(g)
}
alias<-function(from,to){
 if(startsWith(to,'simplified_classes/'))return(invisible(NULL))
 dir.create(dirname(file.path(out,to)),recursive=TRUE,showWarnings=FALSE)
 if(!nzchar(figure_filter)||grepl(figure_filter,to))for(ext in c('.pdf','.png'))stopifnot(file.copy(file.path(out,paste0(from,ext)),file.path(out,paste0(to,ext)),overwrite=TRUE))
 exports[[to]]<<-copy(exports[[from]])[,stem:=to]
 if(!is.null(plot_sources[[from]]))plot_sources[[to]]<<-copy(plot_sources[[from]])[,figure_stem:=to]
 if(!is.null(order_audit[[from]]))order_audit[[to]]<<-copy(order_audit[[from]])[,stem:=to]
 if(!is.null(ann_audit[[from]]))ann_audit[[to]]<<-copy(ann_audit[[from]])[,stem:=to]
}
group_plot<-function(d,stem,title,variant='raw',tests=ttall,value='delta_age_week',relative=FALSE,forced_order=NULL,save=TRUE,min_span=10,panel_order=NULL,ylab=NULL,subtitle=NULL,class_map=NULL,panel_h=PANEL_H,col_w=NULL,box_w=.6){
 if(startsWith(stem,'simplified_classes/'))return(invisible(NULL))
 # col_w: optional exact x spacing per group (inches; the x scale pads 0.7 group per side), else the shared (n + X_PAD) * COL_W panels.
 panel_w<-function(n)if(!is.null(col_w))(n-1+2*.7)*col_w else (n+X_PAD)*COL_W
 d<-annotate_data(d,variant)
 # Optional per-figure recolouring of named groups (display only).
 if(!is.null(class_map))d[plot_group%in%names(class_map),plot_class:=class_map[plot_group]]
 if(relative){
  d[,panel:=paste(gsub('^.*?_','',contrast_id),modality,sep=' | ')]
  # Keep reference means specific to a contrast; repeated controls are explicit in separate panels.
 }else{d<-unique(d,by=c('sample_key','model'));d[,panel:=if(uniqueN(dataset_cohort)>1)modality else '']}
 ord<-if(is.null(forced_order))order_for(d,variant)else forced_order
 if(variant=='special'&&uniqueN(d$dataset_folder)==1&&d$dataset_folder[1]=='PRJNA1002222'&&grepl('/02_',stem))ord<-orders$PRJNA1002222[orders$PRJNA1002222%in%d$plot_group]
 d[,plot_group:=factor(plot_group,levels=ord)];stopifnot(!anyNA(d$plot_group))
 labs<-setNames(ord,ord)
 if(uniqueN(d$panel)==1){n<-d[,.(n=uniqueN(sample_key)),by=plot_group];labs<-setNames(paste0(n$plot_group,'\n(n = ',n$n,')'),n$plot_group)}
 ttg<-test_labels(d,tests);ttg[,`:=`(xmin=match(Reference,ord),xmax=match(Target,ord))]
 if(relative)ttg[,panel:=paste(gsub('^.*?_','',contrast_id),modality,sep=' | ')]else ttg[,panel:=if(uniqueN(d$dataset_cohort)>1)modality else '']
 if(relative&&nrow(ttg))for(pa in unique(ttg$panel)){
  local_order<-ord[ord%in%d[panel==pa,as.character(plot_group)]]
  ttg[panel==pa,`:=`(xmin=match(Reference,local_order),xmax=match(Target,local_order))]
 }
 if(nrow(ttg)){
  rr<-d[,.(lo=min(c(0,get(value))),hi=max(c(0,get(value)))),by=panel];ttg<-merge(ttg,rr,by='panel')
  pack<-function(a,b){lo<-pmin(a,b);hi<-pmax(a,b);tier<-integer(length(a));for(i in order(hi-lo)){t<-1L;while(any(tier>0&tier==t&lo<=hi[i]&hi>=lo[i]))t<-t+1L;tier[i]<-t};tier}
  ttg[,tier:=pack(xmin,xmax),by=.(model,panel)]
  ttg[,span:=pmax(hi-lo,min_span)];ttg[,`:=`(y=hi+span*(.12+.2*(tier-1)),tip=hi+span*(.095+.2*(tier-1)),ylab=hi+span*(.15+.2*(tier-1)))]
 }
 if(relative){plev<-if(is.null(panel_order))sort(unique(as.character(d$panel)))else intersect(panel_order,as.character(d$panel));d[,panel:=factor(panel,levels=plev)];if(nrow(ttg))ttg[,panel:=factor(panel,levels=plev)]}
 g<-ggplot(d,aes(plot_group,.data[[value]],colour=plot_class))+geom_hline(yintercept=0,linetype='dashed',colour='grey45')+geom_boxplot(outlier.shape=NA,alpha=.12,width=box_w)+geom_point(position=position_jitter(width=.12,height=0,seed=310),size=1.8,alpha=.85)+scale_x_discrete(labels=labs,drop=TRUE,expand=expansion(add=.7))+labs(title=title,subtitle=if(!is.null(subtitle))subtitle else if(relative)'Each contrast uses its own reference mean;\nbulk and pseudobulk are separate.'else'External projection samples | Two frozen clocks',x=if(variant=='natural')'Natural aging group'else'Projection group',y=if(!is.null(ylab))ylab else if(relative)'Relative delta age (weeks)'else'Raw delta age (weeks)',colour=if(variant=='natural')'Natural aging group'else'Projection class',caption=if(nrow(ttg))stats_caption else'Boxes: median and IQR; whiskers to 1.5 × IQR. Points: sample profiles. Descriptive placement; no prespecified test.')
 if(relative)g<-g+facet_grid(model~panel,scales='free_x',space='free_x',labeller=labeller(panel=wrap_strip(20)))else if(uniqueN(d$panel)>1)g<-g+facet_grid(panel~model,scales='free_x')else g<-g+facet_wrap(~model,ncol=if(length(ord)>8||variant=='special')1 else 3)
 if(nrow(ttg))g<-g+geom_segment(data=ttg,aes(x=xmin,xend=xmax,y=y,yend=y),inherit.aes=FALSE,linewidth=.35)+geom_segment(data=ttg,aes(x=xmin,xend=xmin,y=tip,yend=y),inherit.aes=FALSE,linewidth=.35)+geom_segment(data=ttg,aes(x=xmax,xend=xmax,y=tip,yend=y),inherit.aes=FALSE,linewidth=.35)+geom_text(data=ttg,aes(x=(xmin+xmax)/2,y=ylab,label=annotation),inherit.aes=FALSE,fontface='bold',size=3.2,vjust=0)
 g<-g+scale_y_continuous(expand=expansion(mult=c(.06,.1)))+theme(axis.text.x=element_text(angle=35,hjust=1,size=if(length(ord)>8)9 else 10),plot.caption=element_text(hjust=0,size=8),strip.text=element_text(size=10))
 if(relative){pcols<-levels(d$panel);g<-g+theme(panel.widths=unit((vapply(pcols,function(pa)uniqueN(d[panel==pa,plot_group]),numeric(1))+X_PAD)*COL_W,'in'),panel.heights=unit(rep(panel_h,uniqueN(d$model)),'in'))}
 else if(uniqueN(d$panel)>1)g<-g+theme(panel.widths=unit(rep(panel_w(length(ord)),uniqueN(d$model)),'in'),panel.heights=unit(rep(PANEL_H,uniqueN(d$panel)),'in'))
 else{nc<-if(length(ord)>8||variant=='special')1 else min(3,uniqueN(d$model));g<-g+theme(panel.widths=unit(rep(panel_w(length(ord)),nc),'in'),panel.heights=unit(rep(PANEL_H,ceiling(uniqueN(d$model)/nc)),'in'))}
 if(variant!='natural')g<-g+class_scale
 if(variant=='special')g$layers[[2]]$aes_params$colour<-'grey55'
 if(variant=='simplified'){g$layers[[3]]$mapping$shape<-aes(shape=plot_role)$shape;g<-g+labs(shape='Projection role')}
 order_audit[[stem]]<<-data.table(stem=stem,variant=variant,group_order=seq_along(ord),plot_group=ord)
 if(nrow(ttg))ann_audit[[stem]]<<-ttg[,.(stem=stem,analysis_contrast_id,model,Reference,Target,welch_p,annotation,analysis_family,panel,xmin,xmax)]
 if(save)saveplot(g,stem,w=NA,h=NA,source=d[,.(sample_key,dataset_cohort,model,plot_group=as.character(plot_group),plot_class,panel,value=get(value),value_type=value)])
 g
}
scatter<-function(d,stem,title,variant='class',training_line=FALSE,w=14){
 d<-annotate_data(d,if(variant%in%c('special','natural'))variant else 'raw')
 if(variant=='class')d[,color_group:=plot_class]else d[,color_group:=plot_group]
 if(variant=='special')d[grepl('natural_aging',comparison_id),color_group:=natural_aging_plot_group]
 if(variant=='combined')d[,color_group:=dataset_folder]
 if(variant%in%c('special','natural')){preferred<-c('young_13w-17.4w','mid_76w','old_86w-87w','2M');lev<-unique(c(preferred[preferred%in%d$color_group],sort(as.character(d$color_group))));d[,color_group:=factor(color_group,levels=lev)]}
 if(variant=='combined')d[,color_group:=factor(color_group,levels=c('PRJNA863443','GSE235021'))]
 g<-ggplot(d,aes(age_week,predicted_age_week,colour=color_group,shape=modality))+geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey45')+geom_point(size=2.4,alpha=.85)+facet_wrap(~model,nrow=1)+labs(title=title,x='Chronological age (weeks)',y='Predicted age (weeks)',colour=if(variant=='combined')'Dataset'else if(variant=='class')'Projection class'else if(variant=='natural')'Natural aging group'else'Projection group',shape='Modality')+theme(axis.text.x=element_text(angle=0))+guides(colour=guide_legend(ncol=if(uniqueN(d$color_group)>8)4 else 3))
 if(grepl('01_external_projection_predicted_vs_chronological_age|diagnostics/02_external_predicted_vs_age',stem)){
  g<-g+annotate('rect',xmin=min(ref$age_week),xmax=max(ref$age_week),ymin=-Inf,ymax=Inf,alpha=.06,fill='grey50')+labs(subtitle='Grey band marks the frozen deployment training age range.')
  lab<-d[order(-abs(delta_age_week)),head(.SD,4),by=model];lab[,outlier_label:=paste(dataset_folder,sub('^.*__','',sample_key),sep='-')]
  if(requireNamespace('ggrepel',quietly=TRUE))g<-g+ggrepel::geom_text_repel(data=lab,aes(label=outlier_label),size=2.4,max.overlaps=20,show.legend=FALSE,seed=310)else g<-g+geom_text(data=lab,aes(label=outlier_label),size=2.2,vjust=-.7,show.legend=FALSE)
 }
 if(training_line)g<-g+geom_abline(data=refstats,aes(slope=fit_slope,intercept=fit_intercept),inherit.aes=FALSE,colour='black',linewidth=.8)+geom_text(data=refstats,aes(x=-Inf,y=Inf,label=locked_metric_label),inherit.aes=FALSE,hjust=-.03,vjust=1.15,size=3)
 if(variant=='class')g<-g+class_scale
 if(variant=='combined')g<-g+scale_colour_manual(values=c(PRJNA863443='#0072B2',GSE235021='#D55E00'))+scale_shape_manual(values=c('bulk RNA-seq'=17,'scRNA-derived pseudobulk'=16),labels=c('Bulk RNA-seq','scRNA pseudobulk'))
 saveplot(g,stem,w=w,h=if(uniqueN(d$color_group)>8)9 else 7,source=d[,.(sample_key,dataset_cohort,model,plot_group=color_group,plot_class,panel=modality,value=predicted_age_week,value_type='predicted_age_week')])
}
concordance_plot<-function(d,stem,title,grouping='plot_class',w=14,h=NULL,ylim=NULL){
 s<-corsummary(d,'pair');s[,label:=sprintf('n = %d; r = %.2f\nmean |difference| = %.1f weeks',n,pearson_r,mean_abs_difference)]
 g<-ggplot(d,aes(x,y,colour=.data[[grouping]],shape=modality))+geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey45')+geom_hline(yintercept=0,linetype='dotted',colour='grey60')+geom_vline(xintercept=0,linetype='dotted',colour='grey60')+geom_point(size=2,alpha=.8)+facet_wrap(~pair,nrow=1,labeller=label_wrap_gen(25))+geom_text(data=s,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.05,vjust=1.1,size=3)+labs(title=title,x='First clock delta age (weeks)',y='Second clock delta age (weeks)',colour=if(grouping=='plot_group')'Projection group'else'Projection class',shape='Modality')+theme(axis.text.x=element_text(angle=0))+guides(colour=guide_legend(ncol=4))
 if(grouping=='plot_class')g<-g+class_scale
 if(grouping=='plot_group')g<-g+guides(shape='none')
 if(w<14)g<-g+theme(legend.box='vertical')+guides(colour=guide_legend(ncol=3))
 # Optional fixed y window (view only; summary stats use every sample).
 if(!is.null(ylim))g<-g+scale_y_continuous(breaks=seq(ylim[1],ylim[2],20))+coord_cartesian(ylim=ylim)
 saveplot(g,stem,w=w,h=if(!is.null(h))h else if(uniqueN(d[[grouping]])>8)9 else 7)
}
heat_plot<-function(hh,stem,title,groupvar='simple_group',row_order=NULL){
 # Optional top-to-bottom row order; groups not listed keep their default order below them.
 if(!is.null(row_order)){hh<-copy(hh);lv<-c(intersect(row_order,hh[[groupvar]]),setdiff(sort(unique(hh[[groupvar]])),row_order));hh[,(groupvar):=factor(get(groupvar),levels=rev(lv))]}
 g<-ggplot(hh,aes(PC,.data[[groupvar]],fill=contribution))+geom_tile(colour='white')+scale_fill_gradient2(low='#3B6EA8',mid='white',high='#C44E52',midpoint=0)+labs(title=title,x='Frozen selected PC',y='Projection group',fill='Weeks')+theme(axis.text.x=element_text(angle=0),panel.grid=element_blank())
 if(uniqueN(hh$dataset_cohort)>1&&groupvar=='simple_group')g<-g+facet_grid(dataset_cohort~.,scales='free_y',space='free_y')+theme(strip.text.y=element_text(angle=0,size=8))
 saveplot(g,stem,w=11,h=max(6,.29*uniqueN(paste(hh$dataset_cohort,hh[[groupvar]]))+2));g
}

# Generate the retained application and reference-folder plots.
for(co in unique(p$dataset_cohort)){
 d<-p[dataset_cohort==co];ds<-d$dataset_folder[1];tests<-ttall[dataset_cohort==co]
 stem<-paste0('diagnostics/',co,'_all_groups_delta_age')
 group_plot(d,stem,paste('Raw delta age by projection group:',ds),tests=tests)
 ids<-unique(tt[dataset_cohort==co,analysis_contrast_id])
 if(length(ids))for(page in seq_len(ceiling(length(ids)/3))){wanted<-ids[seq.int((page-1)*3+1,min(page*3,length(ids)))];keys<-memall[analysis_contrast_id%in%wanted,sample_key];group_plot(d[sample_key%in%keys],paste0('diagnostics/',co,'_delta_age_ttests_page',page),paste('Delta-age comparisons:',ds,'| Set',page),tests=tests[analysis_contrast_id%in%wanted])}
}
# 05 relative plots use the 04 layout: one value per sample, groups in the 04 order.
# A sample in several contrasts normally shares one reference mean; where a group is the
# target in one contrast and the reference in another, the whole dataset is centred on one
# common reference group instead. A constant shift leaves every Welch test unchanged.
common_ref<-c(CRA017492='young control',GSE273724='reference control',GSE279573='WT chow fertile')
relative_04_style<-function(ds){
 r<-rel[dataset_folder==ds];v<-r[,.(relative_delta_age_week=relative_delta_age_week[1],nv=uniqueN(round(relative_delta_age_week,8))),by=.(sample_key,model)]
 note<-'Delta age minus the mean of each contrast\'s matched reference group | Two frozen clocks'
 if(any(v$nv>1)){
  stopifnot(ds%in%names(common_ref));dd<-annotate_data(p[dataset_folder==ds],'root');m<-dd[plot_group==common_ref[[ds]],.(m=mean(delta_age_week)),by=model]
  stopifnot(nrow(m)==length(models));v<-dd[sample_key%in%r$sample_key,.(sample_key,model,relative_delta_age_week=delta_age_week-m$m[match(model,m$model)])]
  note<-paste0('Delta age minus the mean of the ',common_ref[[ds]],' group (common reference) | Two frozen clocks')
 }
 d<-merge(p[dataset_folder==ds&sample_key%in%v$sample_key],v[,.(sample_key,model=as.character(model),relative_delta_age_week)],by.x=c('sample_key',as.character('model')),by.y=c('sample_key','model'))
 d[,model:=factor(model,levels=models)];list(d=d,note=note)
}
for(ds in unique(p$dataset_folder)){
 d<-p[dataset_folder==ds]
 group_plot(if(ds=='GSE235021')d[modality=='bulk RNA-seq']else d,paste0('diagnostics/04_raw_delta_age_by_dataset_treatment_',ds),paste0('Raw delta age: ',ds,if(ds=='GSE235021')' (bulk RNA-seq)'else''),variant='root')
 # Group (column) width and box width match panel a of diagnostics/03_raw_delta_age_and_tage_shift_combined
 # (5.343 in facets / 7.2 x units = 0.742 in per class; boxes width .7).
 group_plot(d,paste0('diagnostics/raw delta_age by dataset_with t test/04_raw_delta_age_by_dataset_treatment_',ds),paste('Raw delta age by projection group:',ds),col_w=if(ds%in%narrow_col_ds)NARROW_COL_W else COMBINED_COL_W,box_w=.7)
 if(ds%in%rel$dataset_folder){d5<-relative_04_style(ds)
  group_plot(if(ds=='GSE235021')d5$d[modality=='bulk RNA-seq']else d5$d,paste0('diagnostics/05_relative_delta_age_matched_controls_',ds),paste0('Matched-control relative delta age: ',ds,if(ds=='GSE235021')' (bulk RNA-seq)'else''),variant='root',value='relative_delta_age_week',ylab='Relative delta age (weeks)',subtitle=d5$note,class_map=if(ds=='GSE228975')c('Srsf1 WT'='reference_control','Srsf1 cKO'='intervention'))}
}
for(stem in c('plots/01_external_projection_predicted_vs_chronological_age','plots/external_projection_predicted_vs_age','diagnostics/02_external_predicted_vs_age'))scatter(p[primary_external_projection_include==TRUE],stem,'External projection samples: primary display cohort')
allscatter<-copy(p);allscatter[,interpretation_frame:=ifelse(primary_external_projection_include,'primary','OOD/supplementary')]
g<-ggplot(allscatter,aes(age_week,predicted_age_week,colour=dataset_folder,shape=interpretation_frame))+geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey45')+geom_point(size=2.4,alpha=.85)+facet_wrap(~model,nrow=1)+labs(title='All external projections including OOD/supplementary samples',x='Chronological age (weeks)',y='Projected transcriptomic age (weeks)',colour='Dataset',shape='Interpretation frame')+theme(axis.text.x=element_text(angle=0))+guides(colour=guide_legend(ncol=5))
saveplot(g,'plots/10_all_projection_predicted_vs_age_including_ood',w=16,h=8)
alias('plots/10_all_projection_predicted_vs_age_including_ood','diagnostics/13_all_projection_predicted_vs_age_including_ood')
g<-ggplot(qc,aes(model,dataset_cohort,fill=min_gene_coverage))+geom_tile(colour='white')+geom_text(aes(label=sprintf('%.1f%%',100*min_gene_coverage)),size=3.2)+scale_fill_gradient(low='white',high='#00BFC4',limits=c(0,1))+labs(title='External projection gene coverage',x=NULL,y=NULL,fill='Measured fraction')+theme(axis.text.x=element_text(angle=0),axis.text.y=element_text(size=9))
saveplot(g,'plots/03_input_gene_coverage',w=11,h=8)

# Layout revision: fixed physical widths; all PC scores and full ranges retained.
pca_layout<-function(g,panel_width,panel_height,title_width,subtitle_width){
 g$labels$title<-paste(strwrap(g$labels$title,width=title_width),collapse='\n')
 g$labels$subtitle<-paste(strwrap(g$labels$subtitle,width=subtitle_width),collapse='\n')
 g+theme(panel.widths=unit(panel_width,'in'),panel.heights=unit(panel_height,'in'),
         plot.title=element_text(size=12),legend.box='vertical',legend.text=element_text(size=9),
         legend.title=element_text(size=10),legend.spacing.y=unit(3,'pt'))+
   guides(colour=guide_legend(ncol=2,title.position='top',order=2),
          shape=guide_legend(ncol=1,title.position='top',order=1))
}
g<-ggplot(pc,aes(PC1,PC2,colour=ifelse(simple_class=='stress','stress_model',simple_class),shape=modality))+geom_point(data=refpc,aes(PC1,PC2),inherit.aes=FALSE,colour='grey75',size=1.6)+geom_point(size=2,alpha=.8)+labs(title='Frozen PCA: external projection and training reference',subtitle='Grey points: 136 development profiles; frozen centering, scaling and loadings',colour='Projection class',shape='Modality')+theme(axis.text.x=element_text(angle=0));g<-g+class_scale
g<-pca_layout(g,4.05,4.85,45,52)
sz<-fig_size(g,minw=0);saveplot(g,'plots/04_frozen_pca_projection',w=sz[['w']],h=sz[['h']])
g<-ggplot(pc,aes(PC1,PC2,colour=ifelse(simple_class=='stress','stress_model',simple_class),shape=modality))+geom_hline(yintercept=0,colour='grey85')+geom_vline(xintercept=0,colour='grey85')+geom_point(size=2.5,alpha=.88)+labs(title='PCA PC1/PC2 projection map',subtitle='External samples only',x='Frozen deployment PC1 score',y='Frozen deployment PC2 score',colour='Projection class',shape='Modality')+theme(axis.text.x=element_text(angle=0));g<-g+class_scale
g<-pca_layout(g,8.27,5.85,80,85)
sz<-fig_size(g,minw=0);saveplot(g,'plots/04_pca_pc1_pc2_projection_map',w=sz[['w']],h=sz[['h']])
alias('plots/04_pca_pc1_pc2_projection_map','diagnostics/09_pca_projection_map')
heat_plot(heat,'plots/05_pca_pc_contribution_heatmap','Frozen PCA: mean PC contributions')
# Supplementary Fig. 4b format: one row per cohort x projection class (mean contribution),
# single panel, rows labelled 'DATASET | class', four active PCs.
hc<-cclong[,.(contribution=mean(pc_contribution_week),n_samples=uniqueN(sample_key)),by=.(dataset_cohort,simple_class,PC)]
hc[,plot_class:=ifelse(simple_class=='stress','stress_model',simple_class)]
hc[,PC:=factor(as.character(PC),levels=active)]
hc[,cohort_label:=sub('__.*$','',dataset_cohort)]
hc[dataset_cohort=='GSE235021__bulk_rnaseq',cohort_label:='GSE235021 (bulk)']
hc[dataset_cohort=='GSE235021__scrna_pseudobulk',cohort_label:='GSE235021 (pseudobulk)']
hc[,row_label:=paste(cohort_label,gsub('_',' ',plot_class),sep=' | ')]
# Datasets ascend from the bottom; within each dataset the reference control is the top row,
# followed by the other classes (alphabetical from the bottom).
rows<-unique(hc[,.(cohort_label,plot_class,row_label)])[order(cohort_label,plot_class=='reference_control',plot_class)]
hc[,row_label:=factor(row_label,levels=rows$row_label)]
g<-ggplot(hc,aes(PC,row_label,fill=contribution))+geom_tile()+
 scale_fill_gradient2(low='#2b6cb0',mid='white',high='#b8322a',midpoint=0)+
 labs(x='Selected PC',y=NULL,fill='Weeks')+
 theme_bw(base_size=11)+theme(panel.grid=element_blank(),legend.position='right',axis.text.x=element_text(angle=0,hjust=.5),
       axis.text.y=element_text(size=8),axis.title.x=element_text(face='bold'),legend.title=element_text(face='bold'),
       panel.widths=unit(length(active)*0.9,'in'),panel.heights=unit(0.17*nrow(rows),'in'))
sz<-fig_size(g,minw=0)
saveplot(g,'plots/05_pca_selected_pc_contribution_heatmap',w=sz[['w']],h=sz[['h']])
wr(hc[,.(dataset_cohort,projection_class=plot_class,row_label,PC,mean_contribution_week=contribution,n_samples)],'tables/selected_pc_contribution_heatmap_source.tsv')
stopifnot(length(active)==4,nrow(hc)==nrow(rows)*length(active),!anyDuplicated(hc[,.(row_label,PC)]),!anyNA(hc$contribution))
for(s in c('plots/pca_projection_pc_contribution_heatmap','diagnostics/10_pca_pc_contribution_heatmap'))alias('plots/05_pca_selected_pc_contribution_heatmap',s)
pairdata[plot_class=='stress',plot_class:='stress_model']
concordance_plot(pairdata,'plots/external_projection_model_concordance','Delta-age concordance: the two locked clocks',w=9,ylim=c(-40,40))
alias('plots/external_projection_model_concordance','diagnostics/06_clock_concordance_delta_age')
allraw<-annotate_data(p)
# Class order matches compare_tAge figure 03: alphabetical, developmental_OOD last.
class_order<-c('diet','genetic_model','intervention','microbiome','reference_control','stress_model','developmental_OOD')
stopifnot(all(allraw$plot_class%in%class_order))
g<-ggplot(copy(allraw)[,plot_class:=factor(plot_class,levels=class_order)],aes(plot_class,delta_age_week,colour=plot_class))+geom_hline(yintercept=0,linetype='dashed',colour='grey45')+geom_boxplot(outlier.shape=NA,alpha=.12,width=.7)+geom_point(position=position_jitter(width=.13,height=0,seed=310),size=2.2,alpha=.65)+facet_wrap(~model,nrow=1)+labs(x='Projection class',y='Delta age (weeks)',colour='Projection class',title='Raw delta age by projection class');g<-g+class_scale
# Size and type scale match compare_tAge figure 03 (11.5 x 7 in, base size 10).
g<-g+theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold'),strip.text=element_text(face='bold'),legend.position='bottom',axis.text.x=element_text(angle=35,hjust=1),legend.text=element_text(size=8))
saveplot(g,'plots/02_delta_age_by_projection_class',w=11.5,h=7)
# Plot object kept for the side-by-side figure with compare_tAge figure 03.
saveRDS(g,file.path(out,'plots/02_delta_age_by_projection_class.rds'))
alias('plots/02_delta_age_by_projection_class','diagnostics/03_raw_delta_age_by_projection_class')
# Retain the reference's sample-level, dataset-faceted overview rather than replacing it with means.
g<-ggplot(allraw,aes(plot_group,delta_age_week,colour=plot_class,shape=plot_role))+geom_hline(yintercept=0,linetype='dashed',colour='grey45')+geom_point(position=position_jitter(width=.16,height=0,seed=310),size=1.6,alpha=.82)+facet_grid(model~dataset_cohort,scales='free_x',space='free_x')+labs(title='Raw delta age by dataset and treatment',x='Projection group',y='Raw delta age (weeks)',colour='Projection class',shape='Projection role')+theme(axis.text.x=element_text(angle=90,hjust=1,size=8),strip.text.x=element_text(size=7));g<-g+class_scale
saveplot(g,'plots/03_delta_age_by_treatment_within_dataset',w=30,h=11)
alias('plots/03_delta_age_by_treatment_within_dataset','plots/external_projection_raw_delta_by_dataset_treatment')
allraw[,interpretation_frame:=ifelse(primary_external_projection_include,'primary_projection','OOD_supplement')]
g<-ggplot(allraw,aes(dataset_cohort,delta_age_week,colour=interpretation_frame))+geom_hline(yintercept=0,linetype='dashed',colour='grey45')+geom_boxplot(outlier.shape=NA,alpha=.12)+geom_point(position=position_jitter(width=.16,seed=310),size=1.4,alpha=.8)+facet_wrap(~model,ncol=1)+labs(title='All external raw delta age by dataset',x='Dataset cohort',y='Delta age (projected - chronological weeks)',colour='Interpretation frame')+theme(axis.text.x=element_text(angle=35,hjust=1,size=9))
saveplot(g,'diagnostics/14_all_projection_raw_delta_by_dataset_including_unmatched',w=16,h=12)
alias('diagnostics/14_all_projection_raw_delta_by_dataset_including_unmatched','plots/11_all_projection_raw_delta_by_dataset_including_unmatched')
for(kind in c('dataset','class')){
 counts<-unique(allraw[,.(sample_key,dataset_cohort,plot_class,modality,primary_external_projection_include)])[,.(n=.N),by=c(if(kind=='dataset')'dataset_cohort'else'plot_class','modality','primary_external_projection_include')]
 counts[,cohort_status:=ifelse(primary_external_projection_include,'primary display','developmental/OOD supplement')]
 g<-ggplot(counts,aes(.data[[if(kind=='dataset')'dataset_cohort'else'plot_class']],n,fill=modality))+geom_col()+facet_wrap(~cohort_status,ncol=1)+labs(x=if(kind=='dataset')'Dataset cohort'else'Projection class',y='Sample profiles',fill='Modality',title='External projection cohort accounting')+theme(axis.text.x=element_text(angle=35,hjust=1,size=9))
 # Keep the overview copies; the duplicate diagnostics copies were removed by the user.
 stem<-if(kind=='dataset')'plots/08_all_projection_cohort_flow_dataset'else'plots/09_all_projection_class_modality_primary_vs_ood'
 saveplot(g,stem,w=15,h=9)
}
# Reinstated on user request (2026-10-13): reference diagnostics/01_projection_cohort_flow_class_modality,
# ported from bulk-pseudobulk_clock/scripts/06_external_projection_diagnostics.R (same geometry, theme and size).
cc<-unique(allraw[,.(sample_key,projection_class=plot_class,modality,primary_external_projection_include)])
stopifnot(!anyDuplicated(cc$sample_key))
cc<-cc[,.(n=.N,n_non_primary=sum(primary_external_projection_include==FALSE,na.rm=TRUE)),by=.(projection_class,modality)]
cc<-merge(CJ(projection_class=sort(unique(cc$projection_class)),modality=sort(unique(cc$modality))),cc,by=c('projection_class','modality'),all.x=TRUE)
cc[is.na(n),n:=0L];cc[is.na(n_non_primary),n_non_primary:=0L]
cl<-cc[,.(n=max(n,na.rm=TRUE),n_non_primary=sum(n_non_primary,na.rm=TRUE)),by=projection_class][n_non_primary>0]
cl[,label:=paste0('OOD/not primary: ',n_non_primary)]
g<-ggplot(cc,aes(projection_class,n,fill=modality))+geom_col(position=position_dodge(width=.78),width=.68,color='grey25')+
 geom_text(data=cl,aes(projection_class,n,label=label),inherit.aes=FALSE,hjust=-.05,size=3.2)+
 expand_limits(y=max(cc$n,na.rm=TRUE)*1.15)+coord_flip()+
 labs(title='Projected samples by projection class and modality',subtitle=paste0(uniqueN(allraw$sample_key),' all projected samples; non-primary/OOD classes are labeled.'),
  x='Projection class',y='Projected sample count',fill='Modality')+
 theme_bw(base_size=14)+theme(plot.title=element_text(face='bold'),plot.subtitle=element_text(size=11),strip.text=element_text(face='bold'),legend.position='bottom',
  axis.text.x=element_text(angle=0,hjust=.5),axis.text.y=element_text(angle=0,hjust=1))
saveplot(g,'diagnostics/01_projection_cohort_flow_class_modality',w=10,h=7,source=cc)
g<-ggplot(availability,aes(dataset_cohort,n_profiles,fill=relative_control_status))+geom_col()+labs(title='Relative-control availability',x='Dataset cohort',y='Sample profiles',fill='Control status')+theme(axis.text.x=element_text(angle=35,hjust=1,size=9))
saveplot(g,'diagnostics/15_relative_control_status_by_dataset',w=15,h=8);alias('diagnostics/15_relative_control_status_by_dataset','plots/12_relative_control_status_by_dataset')
# Sample-level relative-delta overviews preserve the original point/role encodings.
rr<-annotate_data(rel)
g<-ggplot(rr,aes(dataset_cohort,relative_delta_age_week,colour=plot_class,shape=plot_role))+geom_hline(yintercept=0,linetype='dashed',colour='grey45')+geom_point(position=position_jitter(width=.16,seed=310),size=1.6,alpha=.82)+facet_wrap(~model,ncol=1)+labs(x='Dataset cohort',y='Matched-control relative delta age (weeks)',colour='Projection class',shape='Projection role',subtitle='Reference means remain specific to each contrast and modality')+theme(axis.text.x=element_text(angle=35,hjust=1,size=9));g<-g+class_scale
saveplot(g,'plots/external_projection_relative_delta_by_dataset',w=16,h=12)
g<-ggplot(rr,aes(plot_class,relative_delta_age_week,colour=plot_class,shape=plot_role))+geom_hline(yintercept=0,linetype='dashed',colour='grey45')+geom_point(position=position_jitter(width=.16,seed=310),size=1.7,alpha=.82)+facet_wrap(~model,nrow=1)+labs(x='Projection class',y='Matched-control relative delta age (weeks)',colour='Projection class',shape='Projection role',subtitle='Reference means remain specific to each contrast and modality');g<-g+class_scale
saveplot(g,'plots/external_projection_relative_delta_by_projection_class',w=15,h=8)
alias('plots/external_projection_relative_delta_by_projection_class','plots/07_matched_control_relative_delta')
for(ds in c('GSE240547','GSE279446','GSE279573','PRJNA1076509')){
 nm<-switch(ds,GSE240547='GSE240547_HFHS_relative_delta',GSE279446='GSE279446_cohort_relative_delta',GSE279573='GSE279573_relative_delta',PRJNA1076509='PRJNA1076509_FMT_relative_delta')
 group_plot(rel[dataset_folder==ds],paste0('plots/',nm),paste('Relative delta age:',ds),relative=TRUE,value='relative_delta_age_week')
}
gg<-lapply(c('GSE194196','GSE252791'),function(ds)group_plot(rel[dataset_folder==ds],paste0('internal/',ds),paste('Relative delta age:',ds),relative=TRUE,value='relative_delta_age_week',save=FALSE))
gg<-lapply(gg,function(g)g+labs(subtitle=NULL,caption=NULL));sz<-sapply(gg,fig_size,minw=0)
saveplot(wrap_plots(gg,ncol=2)+plot_layout(guides='collect')+plot_annotation(subtitle='Each contrast uses its own reference mean; bulk and pseudobulk are separate.',caption=stats_caption,theme=theme(plot.caption=element_text(hjust=0,size=8)))&theme(legend.position='bottom'),'plots/GSE194196_GSE252791_Bnc1_relative_delta',w=max(9,sum(sz['w',])+0.4),h=max(sz['h',])+1)
pairdata[,absolute_disagreement:=abs(difference)]
g<-ggplot(pairdata,aes(plot_class,absolute_disagreement,colour=plot_class))+geom_boxplot(outlier.shape=NA,alpha=.15)+geom_point(position=position_jitter(width=.15,seed=310),size=1.2,alpha=.8)+facet_grid(dataset_cohort~pair,scales='free_x',labeller=label_wrap_gen(26))+labs(title='Model disagreement by group',x='Projection class',y='Absolute pairwise delta-age difference (weeks)',colour='Projection class')+theme(axis.text.x=element_text(angle=35,hjust=1,size=8),strip.text.y=element_text(size=8));g<-g+class_scale
saveplot(g,'diagnostics/07_model_disagreement_by_group',w=16,h=27)
g<-ggplot(allraw,aes(plot_class,gene_coverage_fraction,colour=model))+geom_hline(yintercept=.8,linetype='dashed',colour='grey45')+geom_boxplot(outlier.shape=NA,alpha=.12)+geom_point(position=position_jitter(width=.16,seed=310),size=1.4,alpha=.8)+facet_wrap(~dataset_cohort,scales='free_x',ncol=4)+labs(title='Projection gene coverage',x='Projection class',y='Gene coverage fraction',colour='Model family')+theme(axis.text.x=element_text(angle=35,hjust=1,size=8),strip.text=element_text(size=8))
saveplot(g,'diagnostics/08_projection_qc_coverage_ood_gene_coverage',w=16,h=15)
alias('diagnostics/08_projection_qc_coverage_ood_gene_coverage','plots/external_projection_gene_coverage')
g<-ggplot(allraw,aes(plot_class,projection_distance,colour=OOD_flag))+geom_boxplot(outlier.shape=NA,alpha=.12)+geom_point(position=position_jitter(width=.16,seed=310),size=1.8,alpha=.8)+facet_wrap(~model,scales='free_y',nrow=1)+labs(title='Projection OOD distance',x='Projection class',y='Projection distance',colour='OOD flag')
saveplot(g,'diagnostics/08_projection_qc_coverage_ood_distance',w=15,h=8)
g<-ggplot(allraw,aes(gene_coverage_fraction,projection_distance,colour=OOD_flag,shape=low_coverage_flag))+geom_vline(xintercept=.8,linetype='dashed',colour='grey45')+geom_point(size=2.2,alpha=.85)+facet_wrap(~model,scales='free_y',nrow=1)+labs(title='Projection coverage and OOD diagnostic scatter',x='Gene coverage fraction',y='Projection distance',colour='OOD flag',shape='Low coverage')+theme(axis.text.x=element_text(angle=0))
saveplot(g,'diagnostics/08_projection_qc_coverage_ood_scatter',w=15,h=8)
alias('diagnostics/08_projection_qc_coverage_ood_scatter','plots/06_ood_gene_coverage')
grade_counts<-claims[,.(n_contrasts=.N),by=claim_grade]
g<-ggplot(grade_counts,aes(claim_grade,n_contrasts,fill=claim_grade))+geom_col(colour='grey25')+geom_text(aes(label=n_contrasts),vjust=-.25,size=4)+labs(title='Projection claim grading',subtitle='Heuristic QC assessment of modality-specific contrasts across two clocks',x='Claim grade',y='Number of contrasts',fill='Claim grade')+theme(legend.position='none')
saveplot(g,'diagnostics/11_projection_claim_grading',w=12,h=7)
w<-worst[order(discrepancy_category,-ranking_value)][seq_len(min(.N,30))];w[,row_label:=paste(discrepancy_category,model,sample_key,sep=' | ')];w[,row_label:=factor(make.unique(row_label),levels=rev(make.unique(row_label))) ]
g<-ggplot(w,aes(ranking_value,row_label,colour=discrepancy_category))+geom_point(size=2.5)+labs(title='Worst/discrepant projected samples',subtitle='Full table includes extreme delta ages, all clock pairs, coverage and distance; ranking units are category-specific.',x='Ranking value',y=NULL,colour='Category')+theme(axis.text.x=element_text(angle=0),axis.text.y=element_text(size=7))
saveplot(g,'diagnostics/12_worst_discrepant_samples',w=18,h=10)

for(ds in c('PRJNA863443','PRJNA1002222')){
 prefix<-paste0(ds,'/plots/');low<-tolower(ds);d<-p[dataset_folder==ds];r<-rel[dataset_folder==ds]
 scatter(d,paste0(prefix,'01_',low,'_predicted_vs_age'),paste(ds,': predicted versus chronological age'),'special',TRUE)
 alias(paste0(prefix,'01_',low,'_predicted_vs_age'),paste0(prefix,'predicted_vs_chronological_age'))
 group_plot(d,paste0(prefix,'02_',low,'_raw_delta_by_group'),paste('Raw delta age:',ds),'special')
 n<-d[grepl('natural_aging',comparison_id)]
 group_plot(n,paste0(prefix,'08_',low,'_natural_aging_raw_delta_by_group'),paste(ds,': natural aging'),'natural',tests=ttall[0])
 scatter(n,paste0(prefix,'09_',low,'_predicted_vs_age_natural_aging_bins'),paste(ds,': natural aging'),'natural')
 contrastsets<-if(ds=='PRJNA863443')list('03_prjna863443_vcd_relative_delta'='VCD','04_prjna863443_foxl2_relative_delta'='Foxl2')else list('03_prjna1002222_cd38ko_8m_relative_delta'='8M','04_prjna1002222_cd38ko_12m_relative_delta'='12M')
 # Foxl2: wild type left of heterozygote; young contrast before mid-age.
 foxl2_order<-c('Foxl2+/+_young','Foxl2+/-_young','Foxl2+/+_mid_age','Foxl2+/-_mid_age')
 for(nm in names(contrastsets)){fx<-grepl('foxl2',nm)
  group_plot(r[grepl(contrastsets[[nm]],contrast_id)],paste0(prefix,nm),paste(ds,': matched-control relative delta age'),'special',relative=TRUE,value='relative_delta_age_week',forced_order=if(fx)foxl2_order,panel_order=if(fx)c('Foxl2_het_vs_WT_young | scRNA-derived pseudobulk','Foxl2_het_vs_WT_mid_age | scRNA-derived pseudobulk'),panel_h=if(ds=='PRJNA863443')1.7 else PANEL_H)}
 pd<-copy(pairdata[dataset_folder==ds]);pd[,`:=`(plot_class=special$projection_class[match(sample_key,special$sample_key)],plot_group=special$projection_group[match(sample_key,special$sample_key)])]
 concordance_plot(pd,paste0(prefix,'05_',low,'_clock_concordance'),paste(ds,': clock concordance'),'plot_group')
 dd<-annotate_data(d,'special');dd[,plot_group:=factor(plot_group,levels=order_for(dd,'special'))]
 g<-ggplot(dd,aes(plot_group,gene_coverage_fraction,colour=model))+geom_hline(yintercept=.8,linetype='dashed',colour='grey45')+geom_boxplot(outlier.shape=NA,colour='grey55')+geom_point(position=position_jitter(width=.12,seed=310),alpha=.85,size=1.8)+labs(x='Projection group',y='Model gene coverage fraction',colour='Model family',title=paste(ds,': gene coverage'))
 saveplot(g,paste0(prefix,'06_',low,'_gene_coverage'),w=16,h=8)
 # Dedicated heatmaps use the dedicated dataset metadata labels.
 hm<-unique(p[dataset_folder==ds,.(simple_group,special_group)])
 special_heat<-merge(heat[dataset_folder==ds],hm,by='simple_group');special_heat[,simple_group:=special_group]
 # PRJNA863443: each VCD group directly below its matched control; Foxl2 WT above het per age.
 heat_order<-if(ds=='PRJNA863443')c(as.vector(outer(c('CTL_','VCD_'),c('3m_30d','3m_90d','10m_30d','10m_90d'),paste0)),foxl2_order,'natural ageing young adult reference','natural ageing old reference')
 heat_plot(special_heat,paste0(prefix,'07_',low,'_pca_pc_contribution_heatmap'),paste(ds,': frozen PC contributions'),row_order=heat_order)
 alias(paste0(prefix,'07_',low,'_pca_pc_contribution_heatmap'),paste0(prefix,'07_pca_pc_contribution_heatmap'))
 for(co in unique(d$dataset_cohort))for(s in names(exports)[startsWith(names(exports),paste0('diagnostics/',co,'_'))])alias(s,paste0(prefix,basename(s)))
}
ds<-'GSE235021';prefix<-'GSE235021/plots/'
scatter(p[dataset_folder==ds],paste0(prefix,'predicted_vs_chronological_age'),'GSE235021: predicted versus chronological age','class',TRUE)
scatter(p[dataset_folder%in%c('GSE235021','PRJNA863443')],paste0(prefix,'01_prjna863443_gse235021_combined_predicted_vs_age'),'PRJNA863443 and GSE235021: external projection','combined',TRUE,w=11)
z<-copy(cclong[dataset_folder==ds]);z[,`:=`(age_month=ifelse(age_week<32,6,9),diet_label=fcase(treatment=='polished rice','Polished rice',treatment=='black rice','Black rice',default='White/brown rice'),modality_label=ifelse(modality=='bulk RNA-seq','Bulk RNA-seq','scRNA pseudobulk'))]
z[,display_group:=paste0(age_month,' months · ',diet_label,' · ',modality_label)]
lev<-CJ(modality_label=c('Bulk RNA-seq','scRNA pseudobulk'),age_month=c(6,9),diet_label=c('Polished rice','Black rice','White/brown rice'),sorted=FALSE)[,paste0(age_month,' months · ',diet_label,' · ',modality_label)]
zh<-z[,.(contribution=mean(pc_contribution_week),sd_contribution=sd(pc_contribution_week),n=.N),by=.(dataset_cohort,modality,display_group,PC)]
interleaved<-CJ(age_month=c(6,9),diet_label=c('Polished rice','Black rice','White/brown rice'),modality_label=c('Bulk RNA-seq','scRNA pseudobulk'),sorted=FALSE)[,paste0(age_month,' months · ',diet_label,' · ',modality_label)]
zh[,display_group:=factor(display_group,levels=rev(interleaved))]
heat_plot(zh,paste0(prefix,'07_gse235021_pca_pc_contribution_heatmap'),'GSE235021: mean frozen PC contributions','display_group')
alias(paste0(prefix,'07_gse235021_pca_pc_contribution_heatmap'),paste0(prefix,'07_pca_pc_contribution_heatmap'))
zh[,display_group:=factor(display_group,levels=rev(lev))]
heat_plot(zh,paste0(prefix,'08_gse235021_pca_pc_contribution_heatmap_modality_blocked'),'GSE235021: mean frozen PC contributions','display_group')
heat_plot(zh[modality=='bulk RNA-seq'],paste0(prefix,'09_gse235021_pca_pc_contribution_heatmap_bulk_only'),'GSE235021 bulk RNA-seq: mean frozen PC contributions','display_group')
wr(zh,'GSE235021/tables/gse235021_pca_pc_contribution_group_summary.tsv')
for(co in unique(p[dataset_folder==ds,dataset_cohort]))for(s in names(exports)[startsWith(names(exports),paste0('diagnostics/',co,'_'))])alias(s,paste0(prefix,basename(s)))

# Corresponding diagnostics and dataset source tables.
wr(test_labels(rawsource),'diagnostics/raw delta_age by dataset_with t test/raw_delta_age_condition_vs_control_tests_by_dataset_model.tsv');wr(groups,'diagnostics/raw delta_age by dataset_with t test/raw_delta_age_group_summary_by_dataset_model.tsv');wr(cm[dataset_folder!='GSE154890'],'diagnostics/raw delta_age by dataset_with t test/raw_delta_age_condition_control_contrast_map.tsv')
for(ds in c('GSE235021','PRJNA863443','PRJNA1002222')){
 low<-tolower(ds);d<-p[dataset_folder==ds]
 wr(d,paste0(ds,'/',low,'_projection_predictions.tsv'));wr(unique(d[,setdiff(names(d),c('model','predicted_age_week','delta_age_week','projection_distance','projection_distance_percentile','OOD_flag','distance_method','gene_coverage_fraction','n_model_genes','missing_gene_fraction')),with=FALSE]),paste0(ds,'/',low,'_projection_metadata.tsv'))
 wr(groups[dataset_folder==ds],paste0(ds,'/',low,'_group_summary.tsv'));wr(qc[dataset_folder==ds],paste0(ds,'/',low,'_gene_coverage.tsv'));wr(rel[dataset_folder==ds],paste0(ds,'/',low,'_relative_to_controls.tsv'));wr(corsummary(pairdata[dataset_folder==ds],c('pair','model_x','model_y')),paste0(ds,'/',low,'_model_concordance.tsv'))
 wr(pc[dataset_folder==ds],paste0(ds,'/',low,'_pca_scores.tsv'));wr(cclong[dataset_folder==ds],paste0(ds,'/',low,'_pca_pc_contributions.tsv'))
 wr(test_labels(annotate_data(d,if(ds=='GSE235021')'new'else'special')),paste0(ds,'/02_',low,'_raw_delta_by_group_ttest.tsv'))
 wr(refstats,paste0(ds,'/tables/frozen_deployment_reference_fit_stats.tsv'))
}
wr(p[dataset_folder%in%c('PRJNA863443','GSE235021')],'GSE235021/tables/PRJNA863443_GSE235021_combined_predicted_vs_age_source.tsv');wr(ref,'GSE235021/tables/PRJNA863443_GSE235021_frozen_deployment_reference_fit.tsv');wr(refstats,'GSE235021/tables/PRJNA863443_GSE235021_frozen_deployment_reference_fit_stats.tsv')
wr(rbindlist(exports),'provenance/complete_figure_exports.tsv');wr(rbindlist(order_audit),'provenance/complete_group_order_audit.tsv');wr(rbindlist(ann_audit),'provenance/complete_plot_annotation_audit.tsv');wr(rbindlist(plot_sources,fill=TRUE),'provenance/complete_plot_source_data.tsv')
# Preserve earlier validation interface while limiting it to the 88 prespecified tests.
wr(rbindlist(exports),'provenance/plot_revision_exports.tsv');wr(rbindlist(ann_audit)[analysis_family=='prespecified'],'provenance/plot_annotation_audit.tsv')
wr(data.table(check=c('281 samples / 20 cohorts / 19 accessions','562 predictions','90 tests including 2 reference additions','control matching separate by modality','relative controls centered','246 in-sample reference predictions (110 HVG + 136 PCA)'),pass=c(uniqueN(p$sample_key)==281&&uniqueN(p$dataset_cohort)==20&&uniqueN(p$dataset_folder)==19,nrow(p)==562,nrow(ttall)==90,all(controlaudit$n_modalities==1),all(abs(controlaudit$centered_control_mean)<1e-9),nrow(ref)==246)),'provenance/comprehensive_analysis_validation.tsv')
message('PASS: comprehensive reference analyses and ',length(exports),' PDF/PNG pairs regenerated.')
