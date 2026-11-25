source(file.path(dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),fixed=TRUE))),'common.R'))
.libPaths(c(file.path(tage,'environment/R_library'),.libPaths()))
Sys.setenv(RETICULATE_PYTHON=file.path(tage,'environment/python_env_locked/bin/python'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
suppressPackageStartupMessages({library(Biobase);library(tAge)})
pp <- file.path(app,'external_projection_predictions.tsv')
sp <- file.path(app,'tables/all_reference_contrast_statistics.tsv')
ap <- file.path(app,'provenance/sample_expression_alignment.tsv')
ip <- read_tsv(pp)[model %chin% models[1:2]]
samples <- ip[,.(dataset_folder=dataset_folder[1],dataset_cohort=dataset_cohort[1],modality=modality[1],age_week=age_week[1],sample_projection_class=sample_projection_class[1],locked_training_sample=any(locked_training_sample),locked_training_study=any(locked_training_study)),by=sample_key]
stopifnot(nrow(samples)==281,!anyDuplicated(samples$sample_key),!anyDuplicated(ip[,.(sample_key,model)]),all(is.finite(ip$predicted_age_week)),!any(samples$dataset_folder=='GSE154890'))
ct <- read_tsv(sp)[model==models[1]]
ct[,model:=NULL]
stopifnot(!anyDuplicated(ct$analysis_contrast_id))
extra_id<-'GSE279573_KO_oilfree_vs_WT_chow_reference_addition'
ct[,membership_source:='explicit application reference/target sample keys']
ct[contrast_id==extra_id,`:=`(
  reference_sample_keys=ct[contrast_id=='GSE279573_KO_chow_vs_WT_chow',reference_sample_keys],
  target_sample_keys=ct[contrast_id=='GSE279573_KO_oilfree_vs_KO_chow',target_sample_keys],
  contrast_type='genetic_model',
  membership_source='WT controls and oil-free targets recovered from the two existing GSE279573 prespecified contrasts; verified against application effect')]
membership <- rbindlist(lapply(seq_len(nrow(ct)),function(i){
  r<-ct[i]; keys<-c(strsplit(r$reference_sample_keys,';',fixed=TRUE)[[1]],strsplit(r$target_sample_keys,';',fixed=TRUE)[[1]])
  roles<-c(rep('control',r$n_control),rep('experimental',r$n_experimental))
  stopifnot(length(keys)==length(roles),!anyDuplicated(keys),all(keys%chin%samples$sample_key))
  data.table(analysis_contrast_id=r$analysis_contrast_id,sample_key=keys,role=roles)
}))
membership<-merge(membership,samples,by='sample_key',sort=FALSE)
ct[,developmental:=vapply(analysis_contrast_id,function(cid) any(membership[analysis_contrast_id==cid,age_week]<7)||any(membership[analysis_contrast_id==cid,sample_projection_class]=='developmental OOD'),logical(1))|contrast_type%chin%'developmental_OOD']
ct[,primary:=analysis_family=='prespecified' & !placement_only & !developmental]
ct[,exclusion_reason:=fcase(analysis_family!='prespecified','additional reference-figure contrast',placement_only,'natural-age placement',developmental,'developmental/prepubertal',default='included primary adult contrast')]
write_tsv(ct,'tables/contrast_definitions.tsv');write_tsv(membership,'tables/contrast_membership.tsv')
samples[,n_contrasts:=vapply(sample_key,function(k)sum(membership$sample_key==k),integer(1))]
samples[,comparison_status:=fifelse(n_contrasts>0,'included in declared contrasts','no prespecified matched-control contrast; no reference invented')]
write_tsv(samples,'tables/sample_inclusion_audit.tsv')
align<-read_tsv(ap)
registry<-read_tsv(file.path(tage,'audit/tage_model_registry.tsv'))
modelpaths<-file.path(tage,'environment/models',c('EN_Chronoage_Mouse_Multitissue_scaleddiff.pkl','EN_Chronoage_Mouse_Multitissue_yugenediff.pkl'))
stopifnot(identical(unname(tools::md5sum(modelpaths)),c('c31d3b5ea8af7aa0891aee12d95109b5','7dff7ed53329448476836fa052ed6097')))
ovpaths<-file.path(release,c('HVG-ElasticNet-exclude-single-age/model.rds','VersionB-PCA-glmnet/model.rds'))
cachepaths<-file.path(tage,'cache',paste0('tage_preprocessed_',sort(unique(membership$dataset_folder)),'.rds'))
sourcepaths<-c(pp,sp,ap,ovpaths,modelpaths,cachepaths,file.path(app,'inputs/model_gene_expression.rds'),file.path(tage,'audit/tage_official_mouse_gene_mapping.tsv'),file.path(tage,'tables/tage_sample_contrast_predictions.tsv'))
before<-data.table(path=normalizePath(sourcepaths),md5_before=unname(tools::md5sum(sourcepaths)))
write_tsv(before,'provenance/input_checksums_before.tsv')
allrows<-list();alignment_rows<-list();coverage_rows<-list();integrity<-list()
for(ds in sort(unique(membership$dataset_folder))){
  message('tAge projection: ',ds)
  cache<-readRDS(file.path(tage,'cache',paste0('tage_preprocessed_',ds,'.rds')))
  aa<-as.data.table(cache$alignment_audit)
  cur<-align[sample_key%chin%samples[dataset_folder==ds,sample_key]]
  # Match by physical source column within modality; correct historical replicate labels without altering values.
  aa[,modality_id:=sub('__.*','',sample_key)]
  cur[,modality_id:=sub('__.*','',sample_key)]
  idx<-match(paste(cur$modality_id,cur$expression_column),paste(aa$modality_id,aa$raw_count_matrix_column))
  stopifnot(!anyNA(idx),!anyDuplicated(idx))
  cur[,cache_sample_key:=aa$sample_key[idx]]
  cur[,sample_label_corrected:=sample_key!=cache_sample_key]
  alignment_rows[[ds]]<-cur
  for(cid in ct[dataset_folder==ds,analysis_contrast_id]){
    cm<-membership[analysis_contrast_id==cid]; cr<-ct[analysis_contrast_id==cid]
    ctl<-cm$role=='control';target<-!ctl
    stopifnot(sum(ctl)==cr$n_control,sum(target)==cr$n_experimental,uniqueN(cm$dataset_cohort)==1)
    integrity[[cid]]<-data.table(analysis_contrast_id=cid,n_control=sum(ctl),n_experimental=sum(target),same_modality=uniqueN(cm$modality)==1,no_control_target_overlap=!anyDuplicated(cm$sample_key),fallback_used=FALSE)
    ck<-cur$cache_sample_key[match(cm$sample_key,cur$sample_key)]
    for(j in 1:2){
      variant<-c('scaleddiff','yugenediff')[j]
      source_eset<-if(j==1)cache$scaled else cache$yugene
      vals<-exprs(source_eset)[,ck,drop=FALSE];colnames(vals)<-cm$sample_key
      pd<-data.frame(sample_key=cm$sample_key,control_flag=ifelse(ctl,'declared_control','target'),row.names=cm$sample_key)
      eset<-ExpressionSet(assayData=vals,phenoData=AnnotatedDataFrame(pd))
      diff<-tAge::control_subtraction(eset,'control_flag','declared_control',verbose=FALSE)
      features<-strsplit(registry[model_filename==basename(modelpaths[j]),required_features],';',fixed=TRUE)[[1]]
      x<-matrix(NA_real_,length(features),nrow(cm),dimnames=list(features,cm$sample_key))
      hit<-intersect(features,rownames(diff));x[hit,]<-exprs(diff)[hit,,drop=FALSE]
      stopifnot(all(is.finite(x[!is.na(x)])))
      aligned<-ExpressionSet(assayData=x,phenoData=AnnotatedDataFrame(pd))
      res<-tAge::predict_tAge_one(aligned,modelpaths[j],species='mouse',mode='EN')
      # Official wrapper preserves the input sample row order.
      stopifnot(nrow(res)==nrow(cm))
      if('sample_id'%in%names(res)) stopifnot(identical(as.character(res$sample_id),cm$sample_key))
      months<-as.numeric(res$EN_tAge);weeks<-months*52.142857/12
      stopifnot(length(weeks)==nrow(cm),all(is.finite(weeks)))
      z<-copy(cm);z[,`:=`(model=models[j+2],model_variant=variant,native_tage_output_months=months,value_weeks=weeks,relative_shift_weeks=weeks-mean(weeks[ctl]),input_gene_coverage=colMeans(is.finite(x)))]
      allrows[[length(allrows)+1L]]<-z
      coverage_rows[[length(coverage_rows)+1L]]<-data.table(analysis_contrast_id=cid,model_variant=variant,n_input_genes=length(features),min_sample_coverage=min(colMeans(is.finite(x))),max_sample_coverage=max(colMeans(is.finite(x))))
    }
    for(nm in models[1:2]){
      p<-ip[model==nm];ii<-match(cm$sample_key,p$sample_key);stopifnot(!anyNA(ii))
      vals<-p$predicted_age_week[ii]
      z<-copy(cm);z[,`:=`(model=nm,value_weeks=vals,relative_shift_weeks=vals-mean(vals[ctl]),delta_age_weeks=p$delta_age_week[ii])]
      allrows[[length(allrows)+1L]]<-z
    }
  }
}
pred<-rbindlist(allrows,fill=TRUE)
pred<-merge(pred,ct[,.(analysis_contrast_id,contrast_id,contrast_type,analysis_family,placement_only,developmental,primary,exclusion_reason)],by='analysis_contrast_id')
pred[,projection_class:=gsub(' ','_',sample_projection_class,fixed=TRUE)]
pred[projection_class=='stress',projection_class:='stress_model']
pred[tolower(projection_class)=='developmental_ood',projection_class:='developmental_OOD']
stopifnot(!anyDuplicated(pred[,.(analysis_contrast_id,sample_key,model)]))
center_error<-pred[role=='control',.(mean_shift=mean(relative_shift_weeks)),by=.(analysis_contrast_id,model)]
stopifnot(max(abs(center_error$mean_shift))<1e-9)
write_tsv(center_error,'provenance/control_centering_checks.tsv')
write_tsv(pred,'tables/sample_contrast_predictions.tsv')
write_tsv(rbindlist(alignment_rows),'provenance/tage_cache_sample_alignment.tsv')
write_tsv(rbindlist(coverage_rows),'tables/tage_gene_coverage.tsv')
write_tsv(rbindlist(integrity),'provenance/control_integrity.tsv')
effects<-pred[,.(n_control=sum(role=='control'),n_experimental=sum(role=='experimental'),effect_weeks=mean(value_weeks[role=='experimental'])-mean(value_weeks[role=='control']),delta_age_effect_weeks=if(all(is.na(delta_age_weeks)))NA_real_ else mean(delta_age_weeks[role=='experimental'])-mean(delta_age_weeks[role=='control']),chronological_age_difference=mean(age_week[role=='experimental'])-mean(age_week[role=='control']),n_training_profiles=sum(locked_training_sample),development_study=any(locked_training_study)),by=.(analysis_contrast_id,dataset_folder,dataset_cohort,model,primary,analysis_family,placement_only,developmental)]
write_tsv(effects,'tables/contrast_effects.tsv')
wide<-dcast(effects,analysis_contrast_id+dataset_folder+primary+analysis_family+placement_only+developmental~model,value.var='effect_weeks')
write_tsv(wide,'tables/four_clock_contrast_effects.tsv')
pair_stats<-function(d,label)rbindlist(lapply(combn(models,2,simplify=FALSE),function(mm)cbind(data.table(analysis_set=label,model_1=mm[1],model_2=mm[2]),safe_stats(d[[mm[1]]],d[[mm[2]]]))))
pairs<-rbindlist(list(pair_stats(wide[primary==TRUE],'primary_adult_prespecified'),pair_stats(wide[analysis_family=='prespecified' & !placement_only],'all_prespecified_nonplacement'),pair_stats(wide[analysis_family=='prespecified'],'all_prespecified'),pair_stats(wide,'all_including_reference_addition')))
pairs[,`:=`(pearson_bh_q=p.adjust(pearson_p,'BH'),spearman_bh_q=p.adjust(spearman_p,'BH')),by=analysis_set]
write_tsv(pairs,'tables/clock_pairwise_agreement.tsv')
tp<-pred[model%chin%models[3:4] & analysis_family=='prespecified']
tw<-dcast(tp,analysis_contrast_id+sample_key~model,value.var='relative_shift_weeks')
te<-wide[analysis_family=='prespecified']
ag<-rbind(cbind(data.table(level='Sample within contrast'),safe_stats(tw[[models[3]]],tw[[models[4]]])),cbind(data.table(level='Contrast effect'),safe_stats(te[[models[3]]],te[[models[4]]])))
write_tsv(ag,'tables/tage_variant_agreement.tsv')
# Check the ovary delta-age estimand against the supplied application, without conflating it with predicted-age differences.
current<-read_tsv(sp)[model%chin%models[1:2],.(analysis_contrast_id,model,application_effect=effect_weeks)]
replay<-merge(effects,current,by=c('analysis_contrast_id','model'))
replay[,absolute_error:=abs(delta_age_effect_weeks-application_effect)]
stopifnot(nrow(replay)==2*nrow(ct),max(replay$absolute_error)<1e-8)
write_tsv(replay,'provenance/application_effect_replay.tsv')
# Recompute gene weights from current artifacts, using standardized feature space as in the historical figure.
om<-lapply(ovpaths,readRDS)
stopifnot(om[[1]]$recipe=='HVG-ElasticNet',om[[2]]$recipe=='VersionB-PCA',length(om[[1]]$representation$genes)==3000,length(om[[2]]$representation$genes)==3000)
en<-copy(om[[1]]$coefficients);en[,model:=models[1]]
pc<-copy(om[[2]]$coefficients);pc[,model:=models[2]]
back<-drop(om[[2]]$representation$rotation[,seq_along(om[[2]]$pc_coefficients),drop=FALSE]%*%om[[2]]$pc_coefficients)
stopifnot(max(abs(pc$standardized_coefficient-back))<1e-10)
write_tsv(en,'tables/ovary_elasticnet_features.tsv');write_tsv(pc,'tables/ovary_pca_backprojected_features.tsv')
# The two locked clocks were filtered on different training sets (110 vs 136 profiles); the shared background is
# the union of their eligible genes so that every gene either clock can use is represented (user decision 2026-10-12).
universe<-sort(union(om[[1]]$representation$filtered_genes,om[[2]]$representation$filtered_genes))
write_tsv(data.table(gene=universe),'tables/current_ovary_gene_universe.tsv')
stopifnot(all(c(en$gene,pc$gene)%chin%universe))
ex<-readRDS(file.path(app,'inputs/model_gene_expression.rds'))$expression
modelchecks<-rbindlist(lapply(1:2,function(i){
  m<-om[[i]];x<-ex[match(m$coefficients$gene,rownames(ex)),,drop=FALSE];missing<-which(!is.finite(x),arr.ind=TRUE);x[missing]<-m$coefficients$center[missing[,1]]
  val<-drop(crossprod(m$coefficients$coefficient,x))+m$intercept
  err<-max(abs(val-ip[model==models[i]]$predicted_age_week[match(colnames(ex),ip[model==models[i]]$sample_key)]))
  stopifnot(err<1e-8);data.table(model=models[i],input_genes=nrow(m$coefficients),active_genes=sum(m$coefficients$coefficient!=0),prediction_replay_max_error=err)
}))
write_tsv(modelchecks,'provenance/locked_model_validation.tsv')
before[,md5_after:=unname(tools::md5sum(path))];before[,unchanged:=md5_before==md5_after];stopifnot(all(before$unchanged))
write_tsv(before,'provenance/input_integrity.tsv')
capture.output(sessionInfo(),file=file.path(out,'provenance/R_session_info.txt'))
message('PASS: ',nrow(ct),' contrasts; ',sum(ct$primary),' primary adult contrasts; ',uniqueN(pred$sample_key),' distinct contrasted profiles.')
