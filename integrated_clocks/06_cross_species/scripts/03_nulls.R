source(file.path(dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),fixed=TRUE))),'common.R'))
inp<-readRDS(file.path(out,'inputs/analysis_inputs.rds'));objs<-readRDS(file.path(out,'inputs/projection_objects.rds'));master<-inp$master
bg<-copy(master[is.finite(human_age_t)&is.finite(mouse_age_t)])
stratify<-function(d){d[,`:=`(expression_bin=pmin(10L,pmax(1L,ceiling(frank(human_mean_expression)/.N*10))),detection_bin=pmin(5L,pmax(1L,ceiling(frank(human_detection_fraction)/.N*5))),variance_bin=pmin(10L,pmax(1L,ceiling(frank(human_expression_variance)/.N*10))))];d[,stratum:=paste(expression_bin,detection_bin,variance_bin,sep='_')];d}
bg<-stratify(bg);pools<-split(seq_len(nrow(bg)),bg$stratum);sets<-list(testable_mouse_universe=bg,shared_model_genes=bg[hvg_selected&pca_selected])
for(prefix in c('hvg','pca')){
 d<-bg[get(paste0(prefix,'_selected'))];sets[[paste0(prefix,'_input_genes')]]<-d
 d<-d[abs(get(paste0(prefix,'_weight')))>0];sets[[paste0(prefix,'_active_genes')]]<-d
 for(n in c(30,100,250,500))sets[[paste0(prefix,'_top',n)]]<-d[order(-abs(get(paste0(prefix,'_weight'))))][seq_len(min(n,nrow(d)))]
}
rows<-list();enrich<-list();adeq<-list()
for(sn in names(sets)){
 message('Stage A matched null: ',sn)
 d<-sets[[sn]];wi<-rep(1,nrow(d));wn<-wi/sum(wi);obs<-mean(sign(d$mouse_age_t)==sign(d$human_age_t));null<-numeric(B);nabs<-numeric(B)
 for(j in seq_len(nrow(d))){pool<-pools[[d$stratum[j]]];draw<-bg$human_age_t[pool[sample.int(length(pool),B,replace=TRUE)]];null<-null+wn[j]*(sign(d$mouse_age_t[j])==sign(draw));nabs<-nabs+abs(draw)/nrow(d)}
 rows[[sn]]<-data.table(gene_set=sn,n_genes=nrow(d),spearman_signed_t=cor(d$mouse_age_t,d$human_age_t,method='spearman'),spearman_standardized_beta=cor(d$mouse_standardized_age_beta,d$human_standardized_age_beta,method='spearman'),sign_concordance=obs,matched_null_mean=mean(null),matched_null_q025=quantile(null,.025),matched_null_q975=quantile(null,.975),empirical_p=(1+sum(null>=obs-1e-12))/(B+1),n_draws=B)
 adeq[[sn]]<-data.table(gene_set=sn,n_genes=nrow(d),mean_abs_human_age_t=mean(abs(d$human_age_t)),matched_null_mean=mean(nabs),empirical_p=(1+sum(nabs>=mean(abs(d$human_age_t))))/(B+1),n_draws=B)
 wt(data.table(iteration=seq_len(B),sign_concordance=null,mean_abs_human_age_t=nabs),paste0('null_distributions/stageA_',sn,'.tsv'))
}
for(prefix in c('hvg','pca'))for(direction in c(-1,1)){
 d<-bg[get(paste0(prefix,'_selected'))&sign(get(paste0(prefix,'_weight')))==direction];w<-abs(d[[paste0(prefix,'_weight')]]);w<-w/sum(w);obs<-sum(w*direction*d$human_age_t);null<-numeric(B)
 for(j in seq_len(nrow(d)))null<-null+w[j]*direction*bg$human_age_t[pools[[d$stratum[j]]][sample.int(length(pools[[d$stratum[j]]]),B,replace=TRUE)]]
 sn<-paste0(prefix,if(direction==1)'_older' else '_younger');enrich[[sn]]<-data.table(signature=sn,n_genes=nrow(d),observed_weighted_signed_t=obs,weighted_sign_compatibility=sum(w*(sign(d$human_age_t)==direction)),null_mean=mean(null),null_q025=quantile(null,.025),null_q975=quantile(null,.975),empirical_p=(1+sum(null>=obs-1e-12))/(B+1),predicted_direction_pass=obs>0,n_draws=B)
 wt(data.table(iteration=seq_len(B),statistic=null),paste0('null_distributions/stageA_weight_',sn,'.tsv'))
}
rows<-rbindlist(rows);rows[,q:=p.adjust(empirical_p,'BH')];enrich<-rbindlist(enrich);enrich[,q:=p.adjust(empirical_p,'BH')];adeq<-rbindlist(adeq);adeq[,q:=p.adjust(empirical_p,'BH')]
wt(rows,'tables/stageA_direction_conservation.tsv');wt(enrich,'tables/stageA_weight_direction_enrichment.tsv');wt(adeq,'tables/current_features_human_age_information.tsv');wt(bg[,.(mouse_gene_symbol,human_gene_id_stable,stratum,expression_bin,detection_bin,variance_bin)],'provenance/stageA_null_matching.tsv')
# Stage B: use the SAME normalized score and technical-adjusted ordered-age slope
# for every observed and randomized signature; preserve inter-gene covariance.
gt<-inp$cohorts$GTEx;orth<-rt('inputs/orthology_one2one.tsv');cand<-orth[human_gene_id_stable%in%rownames(gt$x)]
cx<-gt$x[cand$human_gene_id_stable,,drop=FALSE]
cand[,`:=`(human_mean_expression=rowMeans(cx),human_expression_variance=apply(cx,1,var),human_detection_fraction=rowMeans(gt$raw[human_gene_id_stable,,drop=FALSE]>=10))]
cand<-stratify(cand);cpools<-split(seq_len(nrow(cand)),cand$stratum);Xc<-t(cx);summaryrows<-list();matchaudit<-list()
nullsummary<-function(nm,type,obs,v){p1<-(1+sum(v>=obs))/(B+1);data.table(model=nm,null_type=type,observed=obs,expected_direction_pass=obs>0,null_mean=mean(v),null_q025=quantile(v,.025),null_q975=quantile(v,.975),empirical_one_sided_p=p1,empirical_two_sided_p=(1+sum(abs(v-median(v))>=abs(obs-median(v))))/(B+1),monte_carlo_se=sqrt(p1*(1-p1)/(B+1)),n_draws=B,seed=SEED)}
for(nm in names(models)){
 message('Stage B nulls: ',nm);o<-objs[[paste('GTEx',nm,sep='__')]];m<-models[[nm]];d<-o$meta;score<-o$score
 design<-model.matrix(~age_order+rin_z+ischemia_z+log_library_z,d);take<-which(colnames(design)=='age_order');projector<-solve(crossprod(design),t(design))[take,];obs<-drop(projector%*%score)
 mapped<-o$map$measured;Xs<-o$z[,mapped,drop=FALSE];w<-m$coefficients$standardized_coefficient[mapped]
 stopifnot(max(abs(zsafe(drop(Xs%*%w))-score))<1e-10)
 nullw<-numeric(B);nullg<-numeric(B);active<-which(mapped&m$coefficients$coefficient!=0);wg<-m$coefficients$coefficient[active]
 ci<-match(o$map$human_gene_id_stable[active],cand$human_gene_id_stable);stopifnot(!anyNA(ci));strata<-cand$stratum[ci]
 matchaudit[[nm]]<-data.table(model=nm,mouse_gene_symbol=o$map$mouse_gene_symbol[active],human_gene_id_stable=o$map$human_gene_id_stable[active],stratum=strata,pool_size=lengths(cpools[strata]),raw_weight=wg)
 for(st in seq(1L,B,by=100L)){
  ids<-st:min(st+99L,B);k<-length(ids);W<-replicate(k,sample(w,length(w),replace=FALSE));S<-Xs%*%W;S<-sweep(sweep(S,2,colMeans(S),'-'),2,apply(S,2,sd),'/');nullw[ids]<-drop(projector%*%S)
  # Slot-wise within-stratum sampling with replacement, fixed raw weights.
  # Duplicate candidate draws add; original donor covariance is retained.
  idx<-vapply(strata,function(s)cpools[[s]][sample.int(length(cpools[[s]]),k,replace=TRUE)],integer(k));idx<-t(idx)
  Wg<-sparseMatrix(i=as.vector(idx),j=rep(seq_len(k),each=length(wg)),x=rep(wg,k),dims=c(nrow(cand),k))
  S<-as.matrix(Xc%*%Wg);S<-sweep(sweep(S,2,colMeans(S),'-'),2,apply(S,2,sd),'/');nullg[ids]<-drop(projector%*%S)
 }
 # Freedman-Lane residual permutation preserves the technical-covariate fit.
 reduced<-lm(score~rin_z+ischemia_z+log_library_z,d);res<-residuals(reduced);fitted0<-fitted(reduced)
 nullfl<-replicate(B,drop(projector%*%zsafe(fitted0+sample(res))))
 # Historical age-label permutation is a sensitivity, because age and quality may covary.
 nulllab<-replicate(B,{dp<-copy(d);dp[,age_order:=sample(age_order)];coef(lm.fit(model.matrix(~age_order+rin_z+ischemia_z+log_library_z,dp),score))[2]})
 for(type in c('mapped_weight_permutation','matched_gene_signature','Freedman_Lane_residual_permutation','age_label_permutation_sensitivity')){
  v<-switch(type,mapped_weight_permutation=nullw,matched_gene_signature=nullg,Freedman_Lane_residual_permutation=nullfl,age_label_permutation_sensitivity=nulllab)
  stopifnot(all(is.finite(v)));summaryrows[[paste(nm,type)]]<-nullsummary(nm,type,obs,v);wt(data.table(iteration=seq_len(B),statistic=v),paste0('null_distributions/stageB_',nm,'_',type,'.tsv'))
 }
}
z<-rbindlist(summaryrows);z[,q_within_null_across_two_models:=p.adjust(empirical_one_sided_p,'BH'),by=null_type];wt(z,'tables/stageB_null_benchmarks.tsv');wt(rbindlist(matchaudit),'provenance/stageB_null_matching.tsv')
message('PASS: 10000-draw Stage A and Stage B nulls completed.')
