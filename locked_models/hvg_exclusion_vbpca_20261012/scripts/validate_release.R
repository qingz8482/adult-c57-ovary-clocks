a<-commandArgs(FALSE);script<-gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE)
release<-dirname(dirname(normalizePath(script)))
source(file.path(release,'scripts/core.R'))
manifest<-fromJSON(file.path(release,'manifest.json'))
rows<-list();preds<-list()
for(i in seq_len(nrow(manifest$models))) {
 entry<-manifest$models[i,];dev<-entry$source_dev
 dat<-readRDS(file.path(dev,'shared_inputs/prepared_expression.rds'))
 m<-readRDS(file.path(release,entry$artifact));src<-readRDS(file.path(dev,entry$source_model))
 excl<-unlist(entry$excluded_studies);expected<-dat$metadata[!study%in%excl,sample_key]
 stopifnot(identical(m,src),length(m$training_keys)==entry$training_profiles,setequal(m$training_keys,expected),
           uniqueN(dat$metadata[sample_key%in%m$training_keys,study])==entry$training_studies,m$target=='age_week',m$calibration=='none')
 stopifnot(length(m$representation$genes)==entry$input_genes,!anyDuplicated(m$representation$genes))
 x<-dat$expression[,m$training_keys,drop=FALSE];p<-predict_model(m,x);ref<-predict_model(src,x)
 co<-m$coefficients;v<-t(x[co$gene,,drop=FALSE]);for(j in seq_len(ncol(v)))v[is.na(v[,j]),j]<-co$center[j]
 expanded<-drop(v%*%co$coefficient+m$intercept)
 err<-max(abs(p$pred-ref$pred));coeferr<-max(abs(p$pred-expanded))
 stopifnot(err<1e-10,coeferr<1e-8,all(is.finite(p$pred)))
 rows[[i]]<-data.table(recipe=entry$recipe,training_profiles=length(m$training_keys),training_studies=entry$training_studies,selected_genes=length(m$representation$genes),
  nonzero_gene_coefficients=sum(co$coefficient!=0),PC_inputs=m$hyperparameters$pc,nonzero_PC_count=sum(m$pc_coefficients!=0),
  source_prediction_max_error=err,gene_coefficient_prediction_max_error=coeferr,status='PASS')
 p[,`:=`(recipe=entry$recipe,scope='all_training_profiles_replay_not_validation')];preds[[i]]<-p
}
dir.create(file.path(release,'validation'),showWarnings=FALSE)
fwrite(rbindlist(rows),file.path(release,'validation/replay_checks.csv'))
fwrite(rbindlist(preds),file.path(release,'validation/training_prediction_replay.csv'))
print(rbindlist(rows))
