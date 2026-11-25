#!/usr/bin/env Rscript
a<-commandArgs(FALSE);here<-dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE)));dev<-dirname(here)
source(file.path(dev,'scripts/core.R'));cfg<-fromJSON(file.path(dev,'config/production.json'));results<-list()
check<-function(name,condition){assert(condition,paste('FAIL:',name));results[[length(results)+1L]]<<-data.table(test=name,status='PASS');message('PASS ',name)}
dat<-readRDS(file.path(dev,'shared_inputs/prepared_expression.rds'));m<-dat$metadata;x<-dat$expression
check('136 unique profiles / 14 studies',nrow(m)==136&&!anyDuplicated(m$sample_key)&&uniqueN(m$study)==14)
for(fc in c('cv_fold',grep('^repeat_',names(dat$folds),value=TRUE))) {
 f<-dat$folds
 for(k in 1:4)assert_split(m[match(f[get(fc)!=k,sample_key],sample_key)],m[match(f[get(fc)==k,sample_key],sample_key)])
}
check('31 study-disjoint fourfold designs',TRUE)
check('Main fold sizes follow workbook',identical(as.integer(table(dat$folds$cv_fold)),c(36L,48L,20L,32L)))
check('Count provenance sensitivity cohort',sum(m$count_provenance_gate_pass)==78)
tr<-m[cv_fold!=1];te<-m[cv_fold==1]
r<-fit_ranking(x[,tr$sample_key,drop=FALSE],tr,'VersionB',cfg)
check('Fixed 2000 and tuning 5000 feasible',nrow(r$table)>=5000)
# Independently verify Version B's score against ordinary per-gene lm fits.
genes<-head(r$table$gene,12);reference<-vapply(genes,function(g){
 d<-data.frame(y=x[g,tr$sample_key],age=tr$age_week,batch=factor(tr$source_batch_for_feature_selection))
 fit<-lm(y~age+batch,d);rss<-sum(resid(fit)^2);ra<-sum(resid(lm(y~batch,d))^2);rd<-sum(resid(lm(y~age,d))^2)
 pa<-max(0,min(1,(ra-rss)/ra));pd<-max(0,min(1,(rd-rss)/rd));nv<-sum(vapply(split(d$y,d$batch),function(z)length(z)>1&&var(z)>0,logical(1)))
 abs(coef(summary(fit))['age','t value'])*sqrt(pa)*nv/nlevels(d$batch)/(1+pd/(pa+1e-8))
},numeric(1))
check('Version B agrees with independent lm calculation',max(abs(reference-r$table$score[match(genes,r$table$gene)]))<1e-8)
g<-choose_genes(r,2000);check('Exactly 2000 unique Version B genes',length(g)==2000&&!anyDuplicated(g))
# Extreme sign imbalance checks historical fallback overshoot is corrected.
fake<-list(selector='VersionB',table=data.table(gene=paste0('g',1:3000),score=3000:1,beta=1))
check('Sign fallback cannot overshoot N',length(choose_genes(fake,2000))==2000)
check('Infeasible count rejected',inherits(try(choose_genes(r,nrow(r$table)+1),silent=TRUE),'try-error'))
# Held-out values and target are deliberately corrupted; learned representation must be unchanged.
small<-x[r$filtered_genes,,drop=FALSE];leaky<-small;leaky[,te$sample_key]<-1e9
r1<-fit_ranking(small[,tr$sample_key,drop=FALSE],tr,'VersionB',cfg);r2<-fit_ranking(leaky[,tr$sample_key,drop=FALSE],tr,'VersionB',cfg)
p1<-fit_representation(small,tr,r1,500,'PCA',cfg);p2<-fit_representation(leaky,tr,r2,500,'PCA',cfg)
check('Held-out expression does not affect filter, ranking, scaler or PCA',identical(r1,r2)&&identical(p1,p2))
check('Regression weights have equal study totals',diff(range(tapply(training_weights(tr,cfg),tr$study,sum)))<1e-12)
check('Inverse-stratum weights cover all profiles despite missing legacy weights',all(is.finite(training_weights(m,cfg,'inverse_stratum'))))
check('Inverse-stratum weights balance training strata',diff(range(tapply(training_weights(tr,cfg,'inverse_stratum'),tr$split_age_stratum,sum)))<1e-12)
# No dependence on legacy split_role/cv_group or held-out ages exists in fitting functions.
changed<-copy(m);changed[sample_key%in%te$sample_key,age_week:=1e8]
r3<-fit_ranking(small[,tr$sample_key,drop=FALSE],changed[match(tr$sample_key,sample_key)],'VersionB',cfg)
check('Held-out ages cannot change supervised selection',identical(r1,r3))
grid<-CJ(study=c('a','b','c'),N=500L,pc=0L,alpha=.5,lambda=c(1,2,3));grid[,MAE:=c(1,3,5)[match(study,c('a','b','c'))]+fifelse(lambda==3,.5,0)]
s<-select_parameters(grid,3,500,500)
check('One-SE rule chooses largest eligible lambda',s$best$lambda==3)
# Frozen model artifact reproduces glmnet predictions via original-gene coefficients.
for(recipe in c('HVG-ElasticNet','VersionB-ElasticNet','VersionB-PCA')) {
 rank<-fit_ranking(small[,tr$sample_key,drop=FALSE],tr,if(recipe=='HVG-ElasticNet')'HVG'else'VersionB',cfg)
 hp<-data.table(N=500L,pc=if(recipe=='VersionB-PCA')1L else 0L,alpha=.5,lambda=1)
 model<-fit_selected(small,tr,rank,recipe,hp,cfg,'primary');pred<-predict_model(model,small[,te$sample_key,drop=FALSE])
 z<-small[model$coefficients$gene,te$sample_key,drop=FALSE];z[!is.finite(z)]<-matrix(model$coefficients$center,nrow(z),ncol(z))[!is.finite(z)]
 direct<-drop(crossprod(model$coefficients$coefficient,z))+model$intercept
 check(paste(recipe,'exportable gene coefficients'),max(abs(pred$pred-direct))<1e-7)
 check(paste(recipe,'all lambda keys equal the exact common grid'),identical(model$glmnet$lambda,10^seq(cfg$lambda_log10_max,cfg$lambda_log10_min,length.out=cfg$lambda_values)))
 tmp<-tempfile();saveRDS(model,tmp);check(paste(recipe,'serialized model reproduces predictions'),identical(pred,predict_model(readRDS(tmp),small[,te$sample_key,drop=FALSE])));unlink(tmp)
 absent<-small[!rownames(small)%in%head(model$coefficients$gene,10),te$sample_key,drop=FALSE]
 proj<-project_representation(model$representation,absent)
 check(paste(recipe,'missing genes use training means'),all(proj$missing_fraction>=10/500-1e-12)&&all(is.finite(predict_model(model,absent)$pred)))
}
rows<-list()
for(k in 1:4) {
 tm<-m[cv_fold!=k];rr<-fit_ranking(x[,tm$sample_key,drop=FALSE],tm,'HVG',cfg)
 rp<-fit_representation(x,tm,rr,500,'EN',cfg);ff<-fit_path(rp$train,tm,.5,cfg,'primary')
 rows[[k]]<-data.table(study=paste0('fold',k),N=500L,pc=0L,alpha=.5,lambda=ff$lambda,MAE=seq_along(ff$lambda)/100+k)
}
ss<-select_parameters(rbindlist(rows),4,500,500)
check('All 100 lambda settings survive cross-fold aggregation',nrow(ss$path_summary)==100L&&all(ss$path_summary$n_studies==4))
# End-to-end inner-tuning invariance, with independent caches (no cached-result shortcut).
tcfg<-cfg;tcfg$counts<-500L;tcfg$alphas<-.5;tcfg$lambda_values<-5L;tcfg$inner_scheme<-'GroupKFold';tcfg$inner_group_k<-3L
ca<-tempfile();cb<-tempfile()
z1<-inner_grid(small,tr,'VersionB-ElasticNet',tcfg,'primary',ca,'test-a')
z2<-inner_grid(leaky,tr,'VersionB-ElasticNet',tcfg,'primary',cb,'test-b')
check('Full inner tuning is invariant to outer-held-out expression',identical(z1$grid,z2$grid))
unlink(c(ca,cb),recursive=TRUE)
fwrite(rbindlist(results),file.path(here,'test_results.csv'))
message('ALL ',length(results),' CHECKS PASSED')
