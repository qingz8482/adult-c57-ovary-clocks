#!/usr/bin/env Rscript
a<-commandArgs(FALSE);script<-gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE)
here<-dirname(normalizePath(script));dev<-dirname(here);source(file.path(here,'core.R'))
args<-commandArgs(TRUE); opt<-function(key,default=NULL){v<-args[startsWith(args,paste0('--',key,'='))];if(length(v))sub(paste0('^--',key,'='),'',v[1])else default}
config_path<-normalizePath(opt('config',file.path(dev,'config/production.json')));cfg<-fromJSON(config_path)
assert(cfg$training_gene_measurement_fraction==1,'Only complete training measurement policy is implemented')
assert(cfg$missing_gene_policy=='training_mean'&&cfg$calibration=='none'&&cfg$batch_correction=='none','Unsupported preprocessing configuration')
assert(cfg$fixed_count%in%cfg$counts,'Fixed count must belong to common count grid')
recipe<-opt('recipe');assert(recipe%in%c('HVG-ElasticNet','VersionB-ElasticNet','VersionB-PCA'),'Specify valid --recipe=')
base<-if(grepl('^/',cfg$output_dir))cfg$output_dir else file.path(dev,cfg$output_dir)
folder<-file.path(base,if(recipe=='VersionB-PCA')'VersionB-PCA/glmnet'else recipe)
input<-file.path(dev,'shared_inputs/prepared_expression.rds');dat<-readRDS(input);x<-dat$expression;m<-dat$metadata;folds<-dat$folds
assert(identical(colnames(x),m$sample_key),'Prepared input alignment failed')
if(any(!m$primary_reconstructed_training_ready)) warning('Mixed-provenance expression benchmark: readiness FALSE; retaining intended cohort explicitly.')
source_files<-c(file.path(here,c('core.R','run_models.R')),file.path(dev,'shared_inputs/source_checksums.tsv'),input)
signature<-object_hash(list(config=cfg,source_hashes=tools::md5sum(source_files),R=R.version.string,glmnet=as.character(packageVersion('glmnet'))))
cache_dir<-file.path(base,'cache',signature);dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
for(mode in c('fixed_2000','count_tuning')) {
 d<-file.path(folder,mode,'results');dir.create(d,recursive=TRUE,showWarnings=FALSE)
 stamp<-file.path(d,'run_signature.txt')
 if(file.exists(stamp))assert(identical(readLines(stamp),signature),'Config, inputs or code changed: use a new output_dir to avoid mixing runs') else writeLines(signature,stamp)
 write_json(cfg,file.path(d,'config_used.json'),auto_unbox=TRUE,pretty=TRUE)
 fwrite(folds,file.path(d,'fold_assignments.csv'));capture.output(sessionInfo(),file=file.path(d,'sessionInfo.txt'))
}
tasks<-list();add<-function(analysis,repeat_id,fold,train,test,variant='primary') {
 tasks[[length(tasks)+1L]]<<-list(analysis=analysis,repeat_id=repeat_id,fold=fold,train=train,test=test,variant=variant)
}
allids<-m$sample_key
if('A'%in%cfg$analyses) for(k in 1:4)add('A',0,k,folds[cv_fold!=k,sample_key],folds[cv_fold==k,sample_key])
for(an in intersect(c('B','C'),cfg$analyses)){
 s<-if(an=='B')'NIHMS1554107'else'2024_Cell_Winkler_FRT_scRNA_spatial';add(an,0,1,m[study!=s,sample_key],m[study==s,sample_key])
}
if('final'%in%cfg$analyses)add('final',0,0,allids,character())
if('A_repeated'%in%cfg$analyses)for(r in seq_len(cfg$repeats))for(k in 1:4){c<-paste0('repeat_',sprintf('%02d',r));add('A_repeated',r,k,folds[get(c)!=k,sample_key],folds[get(c)==k,sample_key])}
if('sensitivities'%in%cfg$analyses) for(v in cfg$sensitivities) for(k in 1:4) {
 tr<-folds[cv_fold!=k,sample_key];te<-folds[cv_fold==k,sample_key]
 if(v=='exclude_single_age'){valid<-m[sample_key%in%tr,.(n=uniqueN(age_week)),by=study][n>1,study];tr<-m[sample_key%in%tr&study%in%valid,sample_key]}
 if(v=='verified_inputs'){valid<-m[count_provenance_gate_pass==TRUE,sample_key];tr<-intersect(tr,valid);te<-intersect(te,valid)}
 if(length(te))add(paste0('A_',v),0,k,tr,te,v)
}
# Rebutal4: optional training-only exclusion, applied identically to every outer task; held-out sets unchanged.
if(identical(cfg$training_exclusion,'single_age_studies')) for(i in seq_along(tasks)) {
 tr<-tasks[[i]]$train;valid<-m[sample_key%in%tr,.(n=uniqueN(age_week)),by=study][n>1,study]
 tasks[[i]]$train<-m[sample_key%in%tr&study%in%valid,sample_key]
} else assert(is.null(cfg$training_exclusion),'Unknown training_exclusion')
only<-opt('only');if(!is.null(only))tasks<-Filter(function(t)t$analysis%in%strsplit(only,',',fixed=TRUE)[[1]],tasks)
selected_fold<-opt('fold');if(!is.null(selected_fold))tasks<-Filter(function(t)t$fold==as.integer(selected_fold),tasks)
selected_repeat<-opt('repeat');if(!is.null(selected_repeat))tasks<-Filter(function(t)t$repeat_id==as.integer(selected_repeat),tasks)
limit<-as.integer(opt('max-tasks',length(tasks)));tasks<-head(tasks,limit)
for(task in tasks) {
 id<-sprintf('%s_r%02d_f%d',task$analysis,task$repeat_id,task$fold)
 paths<-setNames(file.path(folder,c('fixed_2000','count_tuning'),'results','checkpoints',paste0(id,'.rds')),c('fixed_2000','count_tuning'))
 if(all(file.exists(paths))){message('RESUME skip ',recipe,' ',id);next}
 message(format(Sys.time()),' START ',recipe,' ',id)
 tr<-m[match(task$train,sample_key)];te<-m[match(task$test,sample_key)];assert_split(tr,te)
 assert(nrow(tr)>0,'Empty training partition')
 selector<-if(recipe=='HVG-ElasticNet')'HVG'else'VersionB';rank<-fit_ranking(x[,tr$sample_key,drop=FALSE],tr,selector,cfg)
 inner<-inner_grid(x,tr,recipe,cfg,task$variant,cache_dir,signature)
 auditdir<-file.path(folder,'tuning_audit',id);dir.create(auditdir,recursive=TRUE,showWarnings=FALSE)
 fwrite(inner$audit,file.path(auditdir,'inner_split_audit.csv'));fwrite(inner$feasibility,file.path(auditdir,'candidate_feasibility.csv'))
 for(mode in names(paths)) {
  if(file.exists(paths[[mode]]))next
  counts<-if(mode=='fixed_2000')cfg$fixed_count else cfg$counts
  sel<-select_parameters(inner$grid,uniqueN(tr$study),counts,nrow(rank$table),if(recipe=='VersionB-PCA')min(nrow(tr)-1,max(cfg$pcs))else Inf)
  fwrite(sel$candidate_summary,file.path(auditdir,paste0(mode,'_candidate_summary.csv')))
  fwrite(sel$path_summary,file.path(auditdir,paste0(mode,'_lambda_summary.csv')))
  model<-fit_selected(x,tr,rank,recipe,sel$best,cfg,task$variant)
  assert(identical(model$representation$training_keys,tr$sample_key),'Scaler fit IDs invalid')
  assert(!length(intersect(model$representation$training_keys,te$sample_key)),'Outer preprocessing leakage')
  model$run_signature<-signature
  pred<-NULL
  if(nrow(te)) {
   pred<-merge(te[,.(sample_key,study,modality,age_week,age_bin,expression_provenance_status)],predict_model(model,x[,te$sample_key,drop=FALSE]),by='sample_key',sort=FALSE)
   pred<-merge(pred,age_support(tr,te),by='sample_key',sort=FALSE)
   pred[,`:=`(analysis=task$analysis,repeat_id=task$repeat_id,fold=task$fold,recipe=recipe,mode=mode,selected_N=sel$best$N,selected_pc=sel$best$pc)]
   assert(nrow(pred)==nrow(te)&&all(is.finite(pred$pred)),'Invalid outer predictions')
  }
  audit<-data.table(analysis=task$analysis,repeat_id=task$repeat_id,fold=task$fold,variant=task$variant,
   n_train=nrow(tr),n_test=nrow(te),train_studies=paste(unique(tr$study),collapse=';'),test_studies=paste(unique(te$study),collapse=';'),
   train_age_min=min(tr$age_week),train_age_max=max(tr$age_week),test_age_min=if(nrow(te))min(te$age_week)else NA_real_,
   test_age_max=if(nrow(te))max(te$age_week)else NA_real_,requested_N=sel$best$N,actual_N=length(model$representation$genes),
   alpha=sel$best$alpha,lambda=sel$best$lambda,PCs=sel$best$pc,inner_macro_MAE=sel$best$inner_MAE,
   n_nonzero_effective_genes=sum(model$coefficients$coefficient!=0),
   lambda_at_grid_boundary=sel$best$lambda%in%range(inner$grid$lambda),
   actual_count_at_grid_boundary=sel$best$N%in%range(cfg$counts),
   test_max_missing_gene_fraction=if(nrow(te))max(pred$missing_gene_fraction)else NA_real_)
  # Checkpoint contains fitted artifact, predictions, chosen settings and exact training membership.
  atomic_rds(list(model=model,predictions=pred,audit=audit,signature=signature),paths[[mode]])
  if(task$analysis=='final')export_model(model,file.path(folder,mode,'results/final_model'))
  message('  SAVED ',mode,': N=',sel$best$N,' PC=',sel$best$pc,' alpha=',sel$best$alpha,' lambda=',signif(sel$best$lambda,4))
 }
 message(format(Sys.time()),' DONE ',recipe,' ',id)
}
message('REQUESTED TASKS COMPLETE: ',recipe)
