suppressPackageStartupMessages({library(data.table); library(glmnet); library(jsonlite)})
setDTthreads(1)
glmnet.control(fdev=0,devmax=1)
assert <- function(ok,msg) if(!isTRUE(ok)) stop(msg,call.=FALSE)
atomic_rds <- function(x,p) {dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);tmp<-paste0(p,'.',Sys.getpid(),'.tmp');saveRDS(x,tmp);assert(file.rename(tmp,p),'Atomic output rename failed')}
object_hash <- function(x) {p<-tempfile();on.exit(unlink(p));saveRDS(x,p,version=3);unname(tools::md5sum(p))}
row_variance <- function(x) {a<-rowMeans(x);pmax(0,rowSums((x-a)^2)/(ncol(x)-1))}
assert_split <- function(tr,te) {
 assert(!anyDuplicated(tr$sample_key)&&!anyDuplicated(te$sample_key),'Duplicate sample key')
 assert(!length(intersect(tr$sample_key,te$sample_key)),'Sample leakage')
 assert(!length(intersect(tr$study,te$study)),'Study leakage')
}
training_weights <- function(m,cfg,variant='primary') {
 mode<-if(variant=='inverse_stratum') 'inverse_stratum' else cfg$weighting
 w<-switch(mode,study_balanced=1/as.numeric(table(m$study)[m$study]),
  inverse_stratum=1/as.numeric(table(m$split_age_stratum)[m$split_age_stratum]),stop('Unknown weighting'))
 assert(all(is.finite(w)&w>0),'Invalid regression weights');w/mean(w)
}
fit_ranking <- function(x,m,selector,cfg) {
 assert(identical(colnames(x),m$sample_key),'Training expression and metadata misaligned')
 measured<-rowMeans(is.finite(x))>=cfg$training_gene_measurement_fraction
 detected<-rowSums(x>=log1p(cfg$detection_cpm),na.rm=TRUE)>=ceiling(cfg$detection_fraction*ncol(x))
 x<-x[measured&detected,,drop=FALSE]
 # Default requires complete training measurements; no outcome-dependent imputation.
 assert(!anyNA(x),'Training matrix contains missing values after measurement filter')
 vv<-row_variance(x);x<-x[is.finite(vv)&vv>0,,drop=FALSE];vv<-vv[is.finite(vv)&vv>0]
 assert(nrow(x)>1,'Too few eligible measured genes')
 if(selector=='HVG') tab<-data.table(gene=rownames(x),score=vv,beta=NA_real_) else {
  d<-data.frame(age=m$age_week,batch=factor(m$source_batch_for_feature_selection))
  assert(nlevels(d$batch)>=2,'Version B requires multiple source batches')
  full<-model.matrix(~age+batch,d);q<-qr(full);rdf<-nrow(d)-q$rank
  assert(q$rank==ncol(full)&&rdf>=3,'Version B age-plus-batch design not identifiable')
  y<-t(x);rss<-colSums(qr.resid(q,y)^2)
  ra<-colSums(qr.resid(qr(model.matrix(~batch,d)),y)^2)
  rd<-colSums(qr.resid(qr(model.matrix(~age,d)),y)^2)
  beta<-qr.coef(q,y)[2,];se<-sqrt(rss/rdf*chol2inv(qr.R(q))[2,2])
  nv<-rowSums(vapply(split(seq_len(nrow(m)),d$batch),function(ii) {
   if(length(ii)<2) rep(FALSE,nrow(x)) else row_variance(x[,ii,drop=FALSE])>0
  },logical(nrow(x))))
  pa<-pmax(0,pmin(1,(ra-rss)/ra));pd<-pmax(0,pmin(1,(rd-rss)/rd))
  score<-abs(beta/se)*sqrt(pa)*(nv/nlevels(d$batch))/(1+pd/(pa+1e-8))
  score[nv<2|pd>0.9|!is.finite(score)]<-NA_real_
  tab<-data.table(gene=rownames(x),score=score,beta=beta,age_t=beta/se,
    partial_r2_age=pa,partial_r2_dataset=pd,n_variable_batches=nv)
 }
 tab<-tab[is.finite(score)];setorder(tab,-score,gene)
 list(table=tab,training_keys=m$sample_key,filtered_genes=rownames(x),selector=selector)
}
choose_genes <- function(rank,n) {
 tab<-rank$table;assert(nrow(tab)>=n,paste('Infeasible feature count',n,'eligible',nrow(tab)))
 if(rank$selector=='HVG') g<-head(tab$gene,n) else {
  g<-c(head(tab[beta>=0,gene],ceiling(n/2)),head(tab[beta<0,gene],floor(n/2)))
  g<-head(unique(c(g,tab$gene)),n)
 }
 assert(length(g)==n&&!anyDuplicated(g),'Feature selection count mismatch');g
}
fit_representation <- function(x,m,rank,n,family,cfg) {
 assert(identical(rank$training_keys,m$sample_key),'Ranking reused across training partitions')
 g<-choose_genes(rank,n);z<-t(x[g,m$sample_key,drop=FALSE])
 mu<-colMeans(z);sdv<-apply(z,2,sd);assert(all(is.finite(sdv)&sdv>0),'Zero variance selected feature')
 z<-sweep(sweep(z,2,mu,'-'),2,sdv,'/')
 rotation<-NULL;maxpc<-0L
 if(family=='PCA') {
  p<-prcomp(z,center=FALSE,scale.=FALSE,rank.=min(max(cfg$pcs),nrow(z)-1,ncol(z)),tol=1e-8)
  rotation<-p$rotation;z<-z%*%rotation;maxpc<-ncol(z)
 }
 list(genes=g,center=mu,scale=sdv,rotation=rotation,maxpc=maxpc,train=z,
  training_keys=m$sample_key,filtered_genes=rank$filtered_genes)
}
project_representation <- function(rep,x) {
 j<-match(rep$genes,rownames(x));z<-t(x[j,,drop=FALSE]);colnames(z)<-rep$genes
 missing<-rowMeans(!is.finite(z));z<-sweep(sweep(z,2,rep$center,'-'),2,rep$scale,'/')
 z[!is.finite(z)]<-0
 if(!is.null(rep$rotation)) z<-z%*%rep$rotation
 list(x=z,missing_fraction=missing)
}
model_x <- function(z,pc=0L) {
 if(pc>0) z<-z[,seq_len(pc),drop=FALSE]
 if(ncol(z)==1L) z<-cbind(z,zero_padding=0)
 z
}
fit_path <- function(z,m,alpha,cfg,variant,pc=0L) {
 target<-if(variant=='log_age'||cfg$target=='log_age') log(m$age_week) else m$age_week
 lambdas<-10^seq(cfg$lambda_log10_max,cfg$lambda_log10_min,length.out=cfg$lambda_values)
 fit<-glmnet(model_x(z,pc),target,alpha=alpha,lambda=lambdas,weights=training_weights(m,cfg,variant),
  standardize=FALSE,intercept=TRUE,control=list(thresh=cfg$glmnet_thresh,maxit=cfg$glmnet_maxit))
 assert(fit$jerr==0L&&length(fit$lambda)==length(lambdas),'glmnet failed or truncated requested path')
 assert(max(abs(fit$lambda-lambdas)/lambdas)<1e-10,'glmnet path differs from requested lambda grid')
 # glmnet rescales lambda through the response SD, causing fold-dependent last-bit drift.
 # Use exact prespecified grid keys so all 100 settings pool across every inner study.
 fit$lambda<-lambdas
 fit
}
backtransform <- function(y,cfg,variant) if(variant=='log_age'||cfg$target=='log_age') exp(y) else y
inner_splits <- function(m,cfg) {
 studies<-sort(unique(m$study));assert(length(studies)>=3,'Too few training studies')
 if(cfg$inner_scheme=='LOSO') folds<-seq_along(studies) else {
  assert(cfg$inner_scheme=='GroupKFold','Unknown inner CV scheme')
  # Deterministic study assignment balanced by profile count, never age/expression.
  k<-min(cfg$inner_group_k,length(studies));folds<-integer(length(studies));sizes<-numeric(k)
  counts<-as.numeric(table(m$study)[studies]);for(j in order(-counts,studies)){v<-which.min(sizes);folds[j]<-v;sizes[v]<-sizes[v]+counts[j]}
 }
 lapply(sort(unique(folds)),function(f) studies[folds==f])
}
inner_grid <- function(x,m,recipe,cfg,variant,cache_dir,signature) {
 selector<-if(recipe=='HVG-ElasticNet') 'HVG' else 'VersionB';family<-if(recipe=='VersionB-PCA') 'PCA' else 'EN'
 grids<-list();audits<-list();feasibility<-list()
 splits<-inner_splits(m,cfg)
 for(h in seq_along(splits)) {
  te<-m[study%in%splits[[h]]];tr<-m[!study%in%splits[[h]]];assert_split(tr,te)
  key<-object_hash(list(signature,recipe,variant,tr$sample_key,te$sample_key))
  path<-file.path(cache_dir,paste0(key,'.rds'))
  if(file.exists(path)) val<-readRDS(path) else {
   rank<-fit_ranking(x[,tr$sample_key,drop=FALSE],tr,selector,cfg)
   res<-list();fr<-list();i<-0L
   for(n in cfg$counts) {
    if(nrow(rank$table)<n) {fr[[as.character(n)]]<-data.table(N=n,pc=NA_integer_,feasible=FALSE,reason='insufficient_eligible_genes',eligible=nrow(rank$table));next}
    rep<-fit_representation(x,tr,rank,n,family,cfg);test<-project_representation(rep,x[,te$sample_key,drop=FALSE])
    pcs<-if(family=='PCA') cfg$pcs[cfg$pcs<=rep$maxpc] else 0L
    fr[[as.character(n)]]<-data.table(N=n,pc=if(family=='PCA')cfg$pcs else 0L,
      feasible=if(family=='PCA')cfg$pcs<=rep$maxpc else TRUE,reason='rank_check',eligible=nrow(rank$table))
    for(pc in pcs) for(alpha in cfg$alphas) {
     fit<-fit_path(rep$train,tr,alpha,cfg,variant,pc)
     pred<-backtransform(as.matrix(predict(fit,model_x(test$x,pc))),cfg,variant)
     assert(all(is.finite(pred)),'Nonfinite inner prediction')
     for(study in unique(te$study)) {jj<-which(te$study==study);i<-i+1L
      res[[i]]<-data.table(study=study,N=n,pc=pc,alpha=alpha,lambda=fit$lambda,
        MAE=colMeans(abs(pred[jj,,drop=FALSE]-te$age_week[jj])))
     }
    }
   }
   val<-list(grid=rbindlist(res),feasibility=rbindlist(fr,fill=TRUE),audit=data.table(inner_fold=h,
    train_keys=paste(tr$sample_key,collapse=';'),test_keys=paste(te$sample_key,collapse=';'),
    train_studies=paste(unique(tr$study),collapse=';'),test_studies=paste(unique(te$study),collapse=';'),
    n_train=nrow(tr),n_test=nrow(te),eligible_genes=nrow(rank$table)))
   atomic_rds(val,path)
  }
  grids[[h]]<-val$grid;audits[[h]]<-copy(val$audit)[,inner_fold:=h]
  feasibility[[h]]<-copy(val$feasibility)[,inner_fold:=h]
  message('  ',recipe,' ',variant,' inner ',h,'/',length(splits),' complete')
 }
 list(grid=rbindlist(grids),audit=rbindlist(audits),feasibility=rbindlist(feasibility))
}
select_parameters <- function(grid,nstudies,counts,outer_eligible,outer_maxpc=Inf) {
 z<-grid[N%in%counts&N<=outer_eligible&pc<=outer_maxpc,
   .(mean_mae=mean(MAE),se=sd(MAE)/sqrt(.N),n_studies=.N),by=.(N,pc,alpha,lambda)]
 z<-z[n_studies==nstudies & is.finite(mean_mae)&is.finite(se)]
 assert(nrow(z)>0,'No consistently feasible candidate across inner studies and outer training')
 assert(all(z[, .(path_size=.N), by=.(N,pc,alpha)]$path_size == uniqueN(grid$lambda)),
   'Incomplete shared lambda path after cross-study aggregation')
 # Each path must be complete over every validation study.
 winners<-z[,{
   imin<-which.min(mean_mae);allow<-which(mean_mae<=mean_mae[imin]+se[imin]);j<-allow[which.max(lambda[allow])]
   .(lambda=lambda[j],inner_MAE=mean_mae[j],inner_SE=se[j],minimum_path_MAE=mean_mae[imin],threshold=mean_mae[imin]+se[imin])
 },by=.(N,pc,alpha)]
 setorder(winners,inner_MAE,N,pc,-alpha)
 list(best=winners[1],path_summary=z,candidate_summary=winners)
}
fit_selected <- function(x,m,rank,recipe,hp,cfg,variant) {
 family<-if(recipe=='VersionB-PCA')'PCA' else 'EN'
 rep<-fit_representation(x,m,rank,hp$N,family,cfg);fit<-fit_path(rep$train,m,hp$alpha,cfg,variant,hp$pc)
 cf<-as.numeric(coef(fit,s=hp$lambda));b<-cf[-1]
 if(hp$pc>0) {b<-b[seq_len(hp$pc)];eff<-drop(rep$rotation[,seq_len(hp$pc),drop=FALSE]%*%b)} else eff<-b[seq_along(rep$genes)]
 rawcoef<-eff/rep$scale;intercept<-cf[1]-sum(rawcoef*rep$center)
 rep$train<-NULL
 list(representation=rep,glmnet=fit,hyperparameters=hp,recipe=recipe,variant=variant,
  target=if(variant=='log_age'||cfg$target=='log_age')'log_age' else 'age_week',
  coefficients=data.table(gene=rep$genes,coefficient=rawcoef,standardized_coefficient=eff,
   center=rep$center,scale=rep$scale),intercept=intercept,model_intercept=cf[1],pc_coefficients=if(hp$pc>0)b else NULL,
  config=cfg,training_keys=m$sample_key,training_studies=unique(m$study),training_age_range=range(m$age_week),
  normalization=cfg$normalization,calibration='none')
}
predict_model <- function(model,x) {
 z<-project_representation(model$representation,x);hp<-model$hyperparameters
 y<-as.numeric(predict(model$glmnet,model_x(z$x,hp$pc),s=hp$lambda))
 if(model$target=='log_age') y<-exp(y)
 data.table(sample_key=colnames(x),pred=y,missing_gene_fraction=z$missing_fraction)
}
export_model <- function(model,dir) {
 dir.create(dir,recursive=TRUE,showWarnings=FALSE);atomic_rds(model,file.path(dir,'model.rds'))
 fwrite(model$coefficients,file.path(dir,'coefficients.csv'))
 fwrite(model$coefficients[,.(gene,center,scale)],file.path(dir,'scaler.csv'))
 writeLines(model$representation$genes,file.path(dir,'genes.txt'))
 writeLines(format(model$intercept,digits=17),file.path(dir,'intercept.txt'))
 write_json(as.list(model$hyperparameters),file.path(dir,'hyperparameters.json'),auto_unbox=TRUE,pretty=TRUE)
 if(!is.null(model$representation$rotation)) {
  fwrite(data.table(gene=model$representation$genes,model$representation$rotation),file.path(dir,'pca_loadings.csv'))
  fwrite(data.table(PC=seq_along(model$pc_coefficients),coefficient=model$pc_coefficients),file.path(dir,'pc_coefficients.csv'))
 }
}
age_support <- function(train,test) {
 data.table(sample_key=test$sample_key,train_age_min=min(train$age_week),train_age_max=max(train$age_week),
  outside_training_range=test$age_week<min(train$age_week)|test$age_week>max(train$age_week),
  nearest_training_age_gap=vapply(test$age_week,function(a)min(abs(train$age_week-a)),numeric(1)),
  nearest_same_modality_age_gap=vapply(seq_len(nrow(test)),function(i){a<-train[modality==test$modality[i],age_week];if(length(a))min(abs(a-test$age_week[i])) else NA_real_},numeric(1)))
}
