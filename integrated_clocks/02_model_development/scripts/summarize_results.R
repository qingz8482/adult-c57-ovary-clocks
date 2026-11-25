#!/usr/bin/env Rscript
a<-commandArgs(FALSE);here<-dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE)));dev<-dirname(here)
source(file.path(here,'core.R'));suppressPackageStartupMessages(library(ggplot2))
args<-commandArgs(TRUE);cp<-args[startsWith(args,'--config=')];cfg<-fromJSON(if(length(cp))sub('^--config=','',cp[1])else file.path(dev,'config/production.json'))
base<-if(grepl('^/',cfg$output_dir))cfg$output_dir else file.path(dev,cfg$output_dir)
meta<-fread(file.path(dev,'shared_inputs/sample_metadata.csv'))
metric<-function(d,study_scope=FALSE) {
 err<-d$pred-d$age_week;den<-sum((d$age_week-mean(d$age_week))^2); required_ages<-if(study_scope)3L else 2L
 data.table(n=nrow(d),n_studies=if(study_scope)1L else uniqueN(d$study),n_distinct_ages=uniqueN(d$age_week),
 Pearson_r=if(uniqueN(d$age_week)>=required_ages && is.finite(sd(d$pred)) && sd(d$pred)>0)cor(d$age_week,d$pred)else NA_real_,
 R2=if(den>0)1-sum(err^2)/den else NA_real_,MAE=mean(abs(err)),median_AE=median(abs(err)),bias=mean(err),RMSE=sqrt(mean(err^2)))
}
all_pred<-list();status<-list();all_metrics<-list();study_metrics<-list()
for(recipe in c('HVG-ElasticNet','VersionB-ElasticNet','VersionB-PCA'))for(mode in c('fixed_2000','count_tuning')) {
 path<-file.path(base,if(recipe=='VersionB-PCA')'VersionB-PCA/glmnet'else recipe,mode,'results')
 checkpoints<-list.files(file.path(path,'checkpoints'),pattern='\\.rds$',full.names=TRUE)
 if(!length(checkpoints))next
 obj<-lapply(checkpoints,readRDS);aud<-rbindlist(lapply(obj,'[[','audit'));pred<-rbindlist(lapply(obj,'[[','predictions'),fill=TRUE)
 fwrite(aud,file.path(path,'fit_audit.csv'));if(!nrow(pred))next
 group<-pred[,.(n=.N,n_keys=uniqueN(sample_key),n_folds=uniqueN(fold)),by=.(analysis,repeat_id)]
 group[,expected_n:=fifelse(analysis=='B',18L,fifelse(analysis=='C',29L,fifelse(analysis=='A_verified_inputs',78L,136L)))]
 group[,complete:=n==expected_n&n_keys==expected_n]
 assert(all(group$n==group$n_keys),'Duplicate outer predictions within analysis/repeat')
 fwrite(group,file.path(path,'analysis_completion.csv'));status[[paste(recipe,mode)]]<-copy(group)[,`:=`(recipe=recipe,mode=mode)]
 # Descriptive LOWESS residualization only on complete main/repeated A predictions.
 pred[,descriptive_LOWESS_age_gap:=NA_real_]
 for(i in which(group$complete & group$analysis%in%c('A','A_repeated'))) {
  an<-group$analysis[i];ri<-group$repeat_id[i];ix<-which(pred$analysis==an&pred$repeat_id==ri)
  lo<-lowess(pred$age_week[ix],pred$pred[ix],f=2/3)
  pred$descriptive_LOWESS_age_gap[ix]<-pred$pred[ix]-approx(lo$x,lo$y,xout=pred$age_week[ix],ties=mean,rule=2)$y
 }
 for(an in unique(pred$analysis))fwrite(pred[analysis==an],file.path(path,paste0('predictions_',an,'.csv')))
 # Partial folds remain available as predictions; cohort metrics are emitted only for complete groups.
 keep<-merge(pred,group[complete==TRUE,.(analysis,repeat_id)],by=c('analysis','repeat_id'),sort=FALSE)
 metrics<-list();ms<-list()
 for(an in unique(keep$analysis))for(ri in unique(keep[analysis==an,repeat_id])) {
  d<-keep[analysis==an&repeat_id==ri];st<-d[,metric(.SD,TRUE),by=study]
  ss<-copy(st)[,`:=`(scope='study',label=study)];ss[,study:=NULL]
  pooled<-metric(d)[,`:=`(scope='pooled',label='all')]
  bins<-d[,metric(.SD),by=age_bin][,`:=`(scope='age_bin',label=age_bin)];bins[,age_bin:=NULL]
  macro<-data.table(n=nrow(d),n_studies=nrow(st),scope='macro_study',label='all',MAE=mean(st$MAE),bias=mean(st$bias))
  v<-rbindlist(list(pooled,bins,ss,macro),fill=TRUE)[,`:=`(analysis=an,repeat_id=ri,recipe=recipe,mode=mode)]
  metrics[[paste(an,ri)]]<-v;ms[[paste(an,ri)]]<-copy(st)[,`:=`(analysis=an,repeat_id=ri,recipe=recipe,mode=mode)]
 }
 if(length(metrics)) {
  met<-rbindlist(metrics,fill=TRUE);fwrite(met,file.path(path,'metrics_summary.csv'));all_metrics[[paste(recipe,mode)]]<-met
  study_metrics[[paste(recipe,mode)]]<-rbindlist(ms)
  rr<-met[analysis=='A_repeated'];if(nrow(rr)) {
   long<-melt(rr,id.vars=c('scope','label','repeat_id'),measure.vars=c('Pearson_r','R2','MAE','median_AE','bias','RMSE'),variable.name='metric')
   summary<-long[is.finite(value),.(n_completed_repeats=uniqueN(rr$repeat_id),n_defined_repeats=.N,median=median(value),q25=quantile(value,.25),q75=quantile(value,.75)),by=.(scope,label,metric)]
   fwrite(summary,file.path(path,'repeat_sensitivity_summary.csv'))
  }
 }
 # Stability denominators explicitly use completed evaluation folds, not final model fits.
 for(block in c('A','A_repeated')) {
  use<-obj[vapply(obj,function(z)z$audit$analysis==block,logical(1))]
  if(!length(use))next
  cf<-rbindlist(lapply(use,function(z)copy(z$model$coefficients)[,.(gene,coefficient)]))
  st<-cf[,.(n_selected=.N,n_nonzero=sum(coefficient!=0),mean_coef=sum(coefficient)/length(use)),by=gene]
  st[,`:=`(n_completed_models=length(use),selection_frequency=n_selected/length(use),nonzero_frequency=n_nonzero/length(use))]
  setorder(st,-n_nonzero,-n_selected,gene);fwrite(st,file.path(path,paste0('feature_stability_',block,'.csv')))
 }
 fwrite(aud[,.(n_models=.N),by=.(analysis,requested_N,PCs)],file.path(path,'selected_count_frequencies.csv'))
 dir.create(file.path(path,'plots'),showWarnings=FALSE)
 if(mode=='count_tuning') {
  curves<-list()
  for(o in obj) {
   au<-o$audit;id<-sprintf('%s_r%02d_f%d',au$analysis,au$repeat_id,au$fold)
   cp<-file.path(dirname(dirname(path)),'tuning_audit',id,'count_tuning_candidate_summary.csv')
   if(!file.exists(cp))next
   ca<-fread(cp);ca<-ca[, .SD[which.min(inner_MAE)],by=N]
   ca[,`:=`(analysis=au$analysis,repeat_id=au$repeat_id,fold=au$fold,selected_N=au$requested_N)]
   curves[[id]]<-ca
  }
  if(length(curves)) {
   cc<-rbindlist(curves);fwrite(cc,file.path(path,'inner_feature_count_sensitivity.csv'))
   cc<-cc[analysis=='A'];if(nrow(cc)) {
    p<-ggplot(cc,aes(N,inner_MAE,color=factor(fold),group=fold))+geom_line()+geom_point()+
      scale_x_continuous(breaks=sort(unique(cc$N)))+theme_bw(base_size=10)+
      labs(x='Candidate input gene count',y='Inner study-balanced MAE (weeks)',color='Outer fold',
       title=paste(recipe,'inner count sensitivity'),subtitle='Training-study validation only; lambda selected by the one-SE rule')
    ggsave(file.path(path,'plots/inner_count_sensitivity_A.pdf'),p,width=7,height=4.5)
   }
  }
 }
 for(an in setdiff(unique(keep$analysis),'A_repeated')) {
  d<-keep[analysis==an]
  p<-ggplot(d,aes(age_week,pred,color=study))+geom_abline(slope=1,intercept=0,linetype=2,color='grey50')+
    geom_point(size=1.8,alpha=.8)+coord_equal()+theme_bw(base_size=10)+theme(legend.position='bottom',legend.text=element_text(size=6))+
    guides(color=guide_legend(ncol=2))+labs(x='Chronological age (weeks)',y='Predicted age (weeks)',title=paste(recipe,mode,an),color='Held-out study')
  ggsave(file.path(path,'plots',paste0('predicted_vs_age_',an,'.pdf')),p,width=9,height=7)
  b<-d[,.(bias=mean(pred-age_week)),by=study]
  p<-ggplot(b,aes(reorder(study,bias),bias))+geom_col(fill='#35618F')+coord_flip()+theme_bw(base_size=9)+labs(x=NULL,y='Mean prediction bias (weeks)',title=paste(recipe,mode,an))
  ggsave(file.path(path,'plots',paste0('study_bias_',an,'.pdf')),p,width=10,height=6)
 }
 if(length(metrics)) {
  rrplot<-rbindlist(metrics,fill=TRUE)[analysis=='A_repeated'&scope=='pooled']
  if(nrow(rrplot)) {
   zz<-melt(rrplot,id.vars='repeat_id',measure.vars=c('Pearson_r','R2','MAE','median_AE','bias','RMSE'),variable.name='metric')
   zz<-zz[is.finite(value)]
   p<-ggplot(zz,aes(x='',y=value))+geom_boxplot(width=.45,outlier.shape=NA,fill='#B7CEE4')+
     geom_point(position=position_jitter(width=.08,height=0,seed=cfg$seed),size=1.2,alpha=.7)+
     facet_wrap(~metric,scales='free_y',ncol=3)+theme_bw(base_size=10)+
     labs(x=NULL,y=NULL,title=paste(recipe,mode,'repeated study partitions'),
       subtitle='One point per complete partition; box spans IQR, not an independent-sample confidence interval')
   ggsave(file.path(path,'plots/repeated_partition_metrics.pdf'),p,width=9,height=5.5)
  }
 }
 all_pred[[paste(recipe,mode)]]<-keep
}
sumdir<-file.path(base,'comparison');dir.create(sumdir,recursive=TRUE,showWarnings=FALSE)
if(length(status))fwrite(rbindlist(status),file.path(sumdir,'run_status.csv'))
if(length(all_metrics))fwrite(rbindlist(all_metrics,fill=TRUE),file.path(sumdir,'all_metrics.csv'))
if(length(study_metrics)) {
 st<-rbindlist(study_metrics);ref<-st[recipe=='HVG-ElasticNet'&mode=='fixed_2000',.(analysis,repeat_id,study,baseline_MAE=MAE,baseline_bias=bias)]
 pair<-merge(st,ref,by=c('analysis','repeat_id','study'));pair[,`:=`(delta_MAE_vs_HVG_fixed=MAE-baseline_MAE,delta_bias_vs_HVG_fixed=bias-baseline_bias)]
 fwrite(pair,file.path(sumdir,'paired_per_study_comparisons.csv'))
 # Count tuning contrasts within the same family and outer studies.
 fixed<-st[mode=='fixed_2000',.(analysis,repeat_id,study,recipe,fixed_MAE=MAE,fixed_bias=bias)]
 tune<-merge(st[mode=='count_tuning'],fixed,by=c('analysis','repeat_id','study','recipe'))
 tune[,`:=`(delta_MAE_tuned_minus_fixed=MAE-fixed_MAE,delta_bias_tuned_minus_fixed=bias-fixed_bias)]
 fwrite(tune,file.path(sumdir,'paired_count_tuning_comparisons.csv'))
}
writeLines(c('Repeated-fold IQR describes dependent split sensitivity, not an independent-sample confidence interval.',
'Metrics are reported only for complete analysis/repeat groups. Check run_status.csv and fit_audit.csv for progress.',
'LOWESS age gaps use observed held-out ages for descriptive residualization; they are not deployment calibration.',
'No final-model in-sample performance is reported. Choosing among recipes from these results remains model development.'),file.path(sumdir,'INTERPRETATION.txt'))
message('Summaries updated in ',sumdir)
