#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
setDTthreads(1)
a<-commandArgs(FALSE);here<-dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE)));out<-dirname(here)
rd<-function(p)fread(file.path(out,p));wr<-function(x,p)fwrite(x,file.path(out,p),sep='\t',na='NA')
p<-rd('external_projection_predictions.tsv');mem<-rd('tables/external_projection_contrast_membership.tsv')
p[,`:=`(negative_predicted_age=predicted_age_week<0,high_missing_gene_fraction=missing_gene_fraction>.1)]
wr(p,'external_projection_predictions.tsv');wr(p,'tables/external_projection_predictions.tsv')
wr(p[,.(sample_key,dataset_folder,dataset_cohort,modality,model,predicted_age_week,delta_age_week,negative_predicted_age,outside_training_age_range,missing_gene_fraction,high_missing_gene_fraction,locked_training_sample,locked_training_study)],'tables/external_projection_qc_flags.tsv')
model_levels<-c('HVG-ElasticNet','VersionB-PCA/glmnet');p[,model:=factor(model,levels=model_levels)]
# Figure contract: quantify observed intervention-associated delta-age differences,
# without interpreting projection shifts as independently validated rejuvenation.
# Quantitative grids; R-only PDF and 300-dpi PNG; individual biological profiles,
# mean +/- sample SD for distributions; Welch effect CIs for contrast summaries.
# All plots use weeks and separate modalities. Exact p, t, df, CI and BH q are saved.
writeLines(c('Figure contract: intervention-associated differences in frozen-clock delta age.',
 'Archetype: quantitative grid, two model columns, one contrast per row.',
 'R/ggplot2 exclusively; PDF with editable text and 300-dpi PNG previews.',
 'Distributions: each dot is one input sample-level profile; bar is mean; whiskers are sample SD.',
 'Contrasts: unpaired two-sided Welch t-test on delta age; target minus reference; 95% Welch confidence intervals.',
 'Asterisks use unadjusted P: * <0.05, ** <0.01, *** <0.001, **** <0.0001. Otherwise show numeric P.',
 'BH q-values across all prespecified modality-specific model tests and separately within each model are supplied.',
 'Natural-age placement and developmental comparisons remain descriptive, explicitly labelled.',
 'No fitting, new calibration, batch correction, or filtering on projection outcomes.'),file.path(out,'provenance/figure_contract.txt'))
theme_set(theme_classic(base_size=9,base_family='sans')+theme(strip.background=element_blank(),strip.text=element_text(face='bold',size=8),plot.title=element_text(face='bold',size=11),plot.subtitle=element_text(size=8),axis.text=element_text(colour='black'),legend.position='bottom'))
saveplot<-function(g,stem,w=9,h=4.5){dir.create(dirname(file.path(out,stem)),recursive=TRUE,showWarnings=FALSE);ggsave(file.path(out,paste0(stem,'.pdf')),g,width=w,height=h,device=cairo_pdf,limitsize=FALSE);ggsave(file.path(out,paste0(stem,'.png')),g,width=w,height=h,dpi=300,limitsize=FALSE,bg='white')}
labelp<-function(v)ifelse(!is.finite(v),'not estimable',ifelse(v<1e-4,'****',ifelse(v<.001,'***',ifelse(v<.01,'**',ifelse(v<.05,'*',paste0('P = ',formatC(v,digits=3,format='g')))))))
dat<-merge(mem,p,by=c('sample_key','modality','dataset_cohort'),allow.cartesian=TRUE,sort=FALSE)
dat[,role:=fifelse(contrast_role%in%c('control','placement_reference'),'Reference','Target')]
dat[,role:=factor(role,levels=c('Reference','Target'))]
tests<-dat[,{
 ctl<-delta_age_week[role=='Reference'];trt<-delta_age_week[role=='Target'];tt<-tryCatch(t.test(trt,ctl,var.equal=FALSE,alternative='two.sided'),error=function(e)e)
 ok<-inherits(tt,'htest')
 list(n_control=length(ctl),n_experimental=length(trt),control_mean_delta=mean(ctl),experimental_mean_delta=mean(trt),effect_weeks=mean(trt)-mean(ctl),
  ci95_low=if(ok)tt$conf.int[1] else NA_real_,ci95_high=if(ok)tt$conf.int[2] else NA_real_,t_statistic=if(ok)unname(tt$statistic) else NA_real_,df=if(ok)unname(tt$parameter) else NA_real_,welch_p=if(ok)tt$p.value else NA_real_,
  status=if(ok)'estimated' else conditionMessage(tt),n_training_controls=sum(locked_training_sample&role=='Reference'),n_training_targets=sum(locked_training_sample&role=='Target'),
  outside_training_age_range=any(outside_training_age_range),contrast_type=contrast_type[1],contrast_name=contrast_name[1],placement_only=tage_placement_only[1],
  source_contrast_warning=contrast_warning[1],source_confounding_flag=tage_confounding_flag[1],
  reference_sample_keys=paste(sample_key[role=='Reference'],collapse=';'),target_sample_keys=paste(sample_key[role=='Target'],collapse=';'))
},by=.(dataset_folder,dataset_cohort,modality,contrast_id,analysis_contrast_id,model)]
tests[,`:=`(bh_q_all_tests=p.adjust(welch_p,'BH'),plot_label=labelp(welch_p),test='two-sided unpaired Welch t-test on delta_age_week')]
tests[,bh_q_within_model:=p.adjust(welch_p,'BH'),by=model]
tests[,analysis_note:=fifelse(dataset_folder=='GSE235021','Historical pooled contrast split by modality for this application','Prespecified contrast membership retained')]
wr(tests,'tables/external_projection_contrast_statistics.tsv');wr(tests,'external_projection_contrast_statistics.tsv')
wr(tests,'diagnostics/delta_age_ttest_statistics.tsv')
wide_tests<-dcast(tests,dataset_folder+dataset_cohort+modality+contrast_id+contrast_type+contrast_name+n_control+n_experimental~model,value.var=c('effect_weeks','ci95_low','ci95_high','welch_p','bh_q_all_tests'))
wr(wide_tests,'tables/two_locked_models_contrast_summary.tsv')
wide_predictions<-dcast(p,sample_key+dataset_folder+dataset_cohort+modality+age_week+locked_training_sample~model,value.var=c('predicted_age_week','delta_age_week'))
wr(wide_predictions,'tables/two_locked_models_sample_predictions.tsv')
dat[,control_mean_delta:=mean(delta_age_week[role=='Reference']),by=.(analysis_contrast_id,model)]
dat[,relative_delta_age_week:=delta_age_week-control_mean_delta]
wr(dat,'tables/external_projection_relative_to_controls.tsv')
wr(dat[,.(sample_key,dataset_cohort,modality,model,analysis_contrast_id,contrast_name,contrast_type,role,age_week,predicted_age_week,delta_age_week,relative_delta_age_week,locked_training_sample)],'diagnostics/delta_age_contrast_plot_source.tsv')
gs<-dat[,.(n=.N,mean_delta_age_week=mean(delta_age_week),sd_delta_age_week=sd(delta_age_week),mean_predicted_age_week=mean(predicted_age_week)),by=.(dataset_cohort,analysis_contrast_id,model,role)]
wr(gs,'tables/external_projection_group_summary.tsv')
# All samples, including profiles without a prespecified contrast, remain visible.
p[,group_label:=paste0(treatment,' | ',genotype,' | ',format(round(age_week,2),trim=TRUE),' wk')]
allgs<-p[,.(n=.N,mean_delta=mean(delta_age_week),sd_delta=sd(delta_age_week)),by=.(dataset_cohort,model,group_label)]
wr(p,'diagnostics/delta_age_all_groups_plot_source.tsv');wr(allgs,'tables/all_sample_group_summary.tsv')
palette<-c(Reference='#8A8F98',Target='#6B94B8')
for(co in unique(p$dataset_cohort)){
 pp<-p[dataset_cohort==co];ss<-allgs[dataset_cohort==co];ds<-pp$dataset_folder[1];mo<-pp$modality[1]
 # Horizontal layout accommodates descriptive treatment/genotype names.
 labs<-unique(pp[order(age_week,treatment,genotype),group_label]);pp[,group_label:=factor(group_label,levels=rev(labs))];ss[,group_label:=factor(group_label,levels=rev(labs))]
 g<-ggplot(pp,aes(delta_age_week,group_label))+geom_vline(xintercept=0,linetype=2,colour='grey60',linewidth=.3)+geom_errorbar(data=ss,aes(x=mean_delta,xmin=mean_delta-sd_delta,xmax=mean_delta+sd_delta),orientation='y',width=.15,colour='#527596',na.rm=TRUE)+geom_point(position=position_jitter(width=0,height=.12,seed=310),size=1.6,alpha=.8)+geom_point(data=ss,aes(x=mean_delta),shape=18,size=2.5,colour='#527596')+facet_wrap(~model,nrow=1)+scale_y_discrete(labels=function(x)vapply(x,function(s)paste(strwrap(s,width=36),collapse='\n'),character(1)))+labs(title=paste(ds,'— delta age by group'),subtitle=paste(mo,'| All profiles; points = samples, diamond = mean, whiskers = SD'),x='Delta age (predicted − chronological, weeks)',y=NULL)+theme(axis.text.y=element_text(size=7))
 saveplot(g,paste0('diagnostics/',co,'_all_groups_delta_age'),w=12,h=max(4,1.6+length(labs)*.58))
 dd<-dat[dataset_cohort==co];tt<-tests[dataset_cohort==co]
 if(nrow(dd)){
  ids<-unique(dd$analysis_contrast_id)
  for(page in seq_len(ceiling(length(ids)/3))){
   wanted<-ids[seq.int((page-1)*3+1,min(page*3,length(ids)))];d<-dd[analysis_contrast_id%in%wanted];t<-tt[analysis_contrast_id%in%wanted]
   d[,panel:=factor(paste0(sub(paste0('^',ds,'_'),'',contrast_id),'\n[',contrast_type,']'),levels=unique(paste0(sub(paste0('^',ds,'_'),'',contrast_id),'\n[',contrast_type,']')))]
   t[,panel:=d$panel[match(analysis_contrast_id,d$analysis_contrast_id)]]
   sy<-d[,.(mean=mean(delta_age_week),sd=sd(delta_age_week),n=.N),by=.(panel,model,role)]
   bounds<-d[,.(lo=min(delta_age_week),hi=max(delta_age_week)),by=panel]
   bounds<-merge(unique(d[,.(panel,model)]),bounds,by='panel');bounds[,span:=pmax(hi-lo,3)]
   t<-merge(t,bounds,by=c('panel','model'));t[,y:=hi+.23*span]
   nn<-merge(sy,bounds,by=c('panel','model'));nn[,y:=lo-.23*span]
   g<-ggplot(d,aes(role,delta_age_week,colour=role))+geom_hline(yintercept=0,linetype=2,colour='grey65',linewidth=.3)+geom_errorbar(data=sy,aes(y=mean,ymin=mean-sd,ymax=mean+sd),width=.14,linewidth=.5,na.rm=TRUE)+geom_point(position=position_jitter(width=.09,height=0,seed=311),size=1.7)+geom_point(data=sy,aes(y=mean),shape=18,size=2.8)+geom_segment(data=t,aes(x=1,xend=2,y=y,yend=y),inherit.aes=FALSE,colour='black',linewidth=.3)+geom_text(data=t,aes(x=1.5,y=y+.07*span,label=plot_label),inherit.aes=FALSE,size=3)+geom_text(data=nn,aes(x=role,y=y,label=paste0('n = ',n)),inherit.aes=FALSE,size=2.5)+facet_grid(panel~model,scales='free_y')+scale_colour_manual(values=palette)+scale_y_continuous(expand=expansion(mult=c(.12,.18)))+labs(title=paste(ds,'— delta-age contrasts'),subtitle=paste(mo,'| Two-sided Welch t-test; asterisks use unadjusted P'),x=NULL,y='Delta age (weeks)',caption='* P < 0.05; ** < 0.01; *** < 0.001; **** < 0.0001. Mean ± SD; each point is one sample.\nNatural-age placement and developmental/OOD comparisons are descriptive.')+theme(legend.position='none',strip.text.y=element_text(angle=0,size=7),plot.caption=element_text(size=7))
   saveplot(g,paste0('diagnostics/',co,'_delta_age_ttests_page',page),w=11,h=1.8+2.35*length(wanted))
  }
 }
}
# Retained overview plots live outside diagnostics.
coverage<-p[,.(n_samples=.N,min_gene_coverage=min(gene_coverage_fraction),mean_gene_coverage=mean(gene_coverage_fraction),n_above_10pct_missing=sum(missing_gene_fraction>.1)),by=.(dataset_cohort,model)]
wr(coverage,'tables/external_projection_gene_coverage.tsv')
g<-ggplot(coverage,aes(model,dataset_cohort,fill=min_gene_coverage))+geom_tile()+geom_text(aes(label=sprintf('%.1f%%',100*min_gene_coverage)),size=2.7)+scale_fill_gradient(low='#D9C6A8',high='#7292B3',limits=c(0,1))+labs(title='Measured input-gene coverage',x=NULL,y=NULL,fill='Fraction')+theme(axis.text.x=element_text(angle=20,hjust=1),axis.text.y=element_text(size=7))
saveplot(g,'plots/03_input_gene_coverage',w=8,h=7)
pc<-merge(rd('tables/pca_projection_scores.tsv'),p[model=='VersionB-PCA/glmnet',.(sample_key,dataset_folder,dataset_cohort,modality,age_week)],by='sample_key')
g<-ggplot(pc,aes(PC1,PC2,colour=age_week,shape=modality))+geom_point(size=1.8)+scale_colour_viridis_c()+labs(title='Frozen PCA projection',subtitle='Training centering, scaling and loadings; no external PCA fitting',colour='Age (weeks)',shape=NULL)
saveplot(g,'plots/04_frozen_pca_projection',w=5,h=4.5)
for(ds in c('GSE235021','PRJNA863443','PRJNA1002222')){
 dir.create(file.path(out,ds,'plots'),recursive=TRUE,showWarnings=FALSE)
 wr(p[dataset_folder==ds],paste0(ds,'/projection_predictions.tsv'));wr(tests[dataset_folder==ds],paste0(ds,'/contrast_statistics.tsv'));wr(dat[dataset_folder==ds],paste0(ds,'/relative_to_controls.tsv'));wr(coverage[grepl(ds,dataset_cohort)],paste0(ds,'/gene_coverage.tsv'));wr(pc[dataset_folder==ds],paste0(ds,'/pca_scores.tsv'))
 cc<-rd('tables/pca_projection_pc_contributions.tsv');wr(cc[sample_key%in%p[dataset_folder==ds,sample_key]],paste0(ds,'/pca_pc_contributions.tsv'))
 pp<-p[dataset_folder==ds]
 g<-ggplot(pp,aes(age_week,predicted_age_week,colour=treatment))+geom_abline(slope=1,intercept=0,linetype=2,colour='grey60')+geom_point(size=2)+facet_grid(modality~model)+labs(title=paste(ds,'— locked model projections'),x='Chronological age (weeks)',y='Predicted age (weeks)',colour='Treatment')+theme(legend.text=element_text(size=7))
 saveplot(g,paste0(ds,'/plots/predicted_vs_chronological_age'),w=10,h=3+2*uniqueN(pp$modality))
 ff<-list.files(file.path(out,'diagnostics'),pattern=paste0('^',ds,'.*\\.(pdf|png)$'),full.names=TRUE);stopifnot(all(file.copy(ff,file.path(out,ds,'plots'),overwrite=TRUE)))
}
# Export pairwise model agreement on exactly matched profiles.
wide<-dcast(p,sample_key+dataset_cohort~model,value.var='delta_age_week')
agreement<-rbindlist(lapply(combn(model_levels,2,simplify=FALSE),function(ab)data.table(model_a=ab[1],model_b=ab[2],n_samples=nrow(wide),pearson_cor=cor(wide[[ab[1]]],wide[[ab[2]]]),mean_absolute_difference=mean(abs(wide[[ab[1]]]-wide[[ab[2]]])))))
wr(agreement,'tables/external_projection_model_concordance.tsv')
checks<-data.table(check=c('281 profiles','20 modality-specific datasets','19 study accessions','562 predictions','no GSE154890 projections','44 modality-specific contrasts','88 model-by-contrast tests','each contrast has both groups','all tests finite','no cross-modality tests'),pass=c(uniqueN(p$sample_key)==281,uniqueN(p$dataset_cohort)==20,uniqueN(p$dataset_folder)==19,nrow(p)==562,!any(p$dataset_folder=='GSE154890'),uniqueN(tests$analysis_contrast_id)==44,nrow(tests)==88,all(tests$n_control>=2&tests$n_experimental>=2),all(is.finite(tests$welch_p)),all(dat[,uniqueN(modality),by=analysis_contrast_id]$V1==1)))
wr(checks,'provenance/validation_checks.tsv');stopifnot(all(checks$pass))
message('PASS: ',nrow(tests),' Welch tests; plots and dataset outputs exported.')
