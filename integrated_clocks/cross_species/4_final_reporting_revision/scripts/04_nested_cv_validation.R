#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]))), "common_final.R"))

message("[4/5] Leakage-controlled nested-CV age-information audit")
x<-load_inputs();gtex<-x$gtex;audit<-x$audit
counts<-gtex$counts_filtered
lib_total<-audit$raw_library_size
X<-log2(t(t(counts)/lib_total*1e6)+2) # sample-local transform; no cohort scaling
y<-audit$age_bracket_midpoint_sensitivity_only
mids<-c(24.5,34.5,44.5,54.5,64.5,74.5);names(mids)<-age_levels
lambda_grid<-10^seq(2,-3,length.out=30)
MAX_FEATURES<-2000L

make_folds<-function(strata,k=5L,seed=SEED){set.seed(seed);out<-integer(length(strata));for(s in unique(strata)){ii<-which(strata==s);out[ii]<-sample(rep(seq_len(k),length.out=length(ii)))};out}
nearest_bin<-function(z)age_levels[max.col(-abs(outer(z,mids,"-")),ties.method="first")]

select_features<-function(train_idx){
  min_n<-ceiling(.20*length(train_idx));keep<-rowSums(counts[,train_idx,drop=FALSE]>=10)>=min_n
  vv<-rowVars(X[keep,train_idx,drop=FALSE]);idx<-which(keep)[order(vv,decreasing=TRUE)[seq_len(min(MAX_FEATURES,length(vv)))]]
  idx
}

inner_lambda<-function(train_idx,ytrain,seed){
  inner<-make_folds(nearest_bin(ytrain),5L,seed);err<-matrix(NA_real_,5,length(lambda_grid))
  for(fi in 1:5){itr<-train_idx[inner!=fi];iva<-train_idx[inner==fi];sel<-select_features(itr)
    fit<-glmnet(t(X[sel,itr,drop=FALSE]),y[itr],alpha=0,family="gaussian",lambda=lambda_grid,standardize=TRUE)
    pr<-predict(fit,newx=t(X[sel,iva,drop=FALSE]),s=lambda_grid)
    err[fi,]<-colMeans((pr-y[iva])^2)
  }
  lambda_grid[which.min(colMeans(err,na.rm=TRUE))]
}

nested_predict_folds<-function(y_target,outer,seed_offset=0L,eligible=rep(TRUE,length(y_target))){
  old_y<-y;assign("y",y_target,envir=.GlobalEnv);on.exit(assign("y",old_y,envir=.GlobalEnv),add=TRUE)
  pred<-rep(NA_real_,length(y_target));chosen<-numeric()
  for(fo in sort(unique(outer[eligible]))){tr<-eligible&outer!=fo;te<-eligible&outer==fo
    lam<-inner_lambda(which(tr),y_target[tr],SEED+seed_offset+fo);sel<-select_features(which(tr))
    fit<-glmnet(t(X[sel,tr,drop=FALSE]),y_target[tr],alpha=0,family="gaussian",lambda=lam,standardize=TRUE)
    pred[te]<-as.numeric(predict(fit,newx=t(X[sel,te,drop=FALSE]),s=lam));chosen<-c(chosen,lam)
  }
  list(pred=pred,lambda=chosen)
}

metric_fun<-function(pred,truth){
  ok<-is.finite(pred)&is.finite(truth);pred<-pred[ok];truth<-truth[ok];tb<-nearest_bin(truth);pb<-nearest_bin(pred)
  recalls<-vapply(age_levels,function(a){ii<-tb==a;if(any(ii))mean(pb[ii]==a)else NA_real_},numeric(1))
  cal<-lm(truth~pred)
  c(pearson=cor(pred,truth),spearman=cor(pred,truth,method="spearman"),mae_midpoint_years=mean(abs(pred-truth)),
    balanced_age_bracket_accuracy=mean(recalls,na.rm=TRUE),ordinal_concordance_kendall=cor(pred,truth,method="kendall"),
    calibration_intercept=unname(coef(cal)[1]),calibration_slope=unname(coef(cal)[2]),exact_bracket_accuracy=mean(pb==tb),n=length(truth))
}

outer_baseline<-make_folds(audit$age_bracket,5L,SEED)
outer_batch_strat<-make_folds(interaction(audit$age_bracket,audit$sequencing_center,drop=TRUE),5L,SEED+10L)
hardy_strata<-interaction(audit$age_bracket,ifelse(is.finite(audit$donor_death_hardy_scale),audit$donor_death_hardy_scale,"missing"),drop=TRUE)
outer_hardy<-make_folds(hardy_strata,5L,SEED+20L)

base<-nested_predict_folds(y,outer_baseline,0L)
batch_strat<-nested_predict_folds(y,outer_batch_strat,100L)
hardy_strat<-nested_predict_folds(y,outer_hardy,200L)

## Two-direction sequencing-center blocked sensitivity; singleton D1 excluded.
eligible_block<-audit$sequencing_center%in%c("B1","C1")
outer_block<-ifelse(audit$sequencing_center=="B1",1L,ifelse(audit$sequencing_center=="C1",2L,NA_integer_))
batch_block<-nested_predict_folds(y,outer_block,300L,eligible_block)

predictions<-rbindlist(list(
  data.table(sample_id=audit$sample_id,specification="baseline_age_stratified_5fold",observed=y,predicted=base$pred,outer_fold=outer_baseline),
  data.table(sample_id=audit$sample_id,specification="batch_stratified_5fold",observed=y,predicted=batch_strat$pred,outer_fold=outer_batch_strat),
  data.table(sample_id=audit$sample_id,specification="hardy_age_stratified_5fold",observed=y,predicted=hardy_strat$pred,outer_fold=outer_hardy),
  data.table(sample_id=audit$sample_id,specification="sequencing_center_blocked",observed=y,predicted=batch_block$pred,outer_fold=outer_block)))
predictions<-merge(predictions,audit[,.(sample_id,donor_id,age_bracket,sequencing_center,donor_death_hardy_scale)],by="sample_id")
write_tsv(predictions,pfinal("results/table_nested_cv_predictions_validated.tsv"))

perf<-predictions[,as.list(metric_fun(predicted,observed)),by=specification]
write_tsv(perf,pfinal("results/table_nested_cv_performance_validated.tsv"))
bybin<-predictions[specification=="baseline_age_stratified_5fold",.(n=.N,
  mae_midpoint_years=mean(abs(predicted-observed)),mean_prediction=mean(predicted),prediction_sd=sd(predicted),bias=mean(predicted-observed),
  exact_bracket_accuracy=mean(nearest_bin(predicted)==as.character(age_bracket))),by=age_bracket]
write_tsv(bybin,pfinal("results/table_nested_cv_performance_by_age_bracket.tsv"))

## Stratified donor-label permutation null; the entire nested workflow is rerun.
N_PERM<-100L;set.seed(SEED+30L);perm_metrics<-list()
for(b in seq_len(N_PERM)){
  yp<-y
  for(s in unique(audit$sequencing_center)){ii<-which(audit$sequencing_center==s);yp[ii]<-sample(yp[ii])}
  outerp<-make_folds(nearest_bin(yp),5L,SEED+1000L+b)
  z<-nested_predict_folds(yp,outerp,2000L+b*10L)
  perm_metrics[[b]]<-as.data.table(c(list(iteration=b),as.list(metric_fun(z$pred,yp))))
  if(b%%10L==0L)message("  nested permutation ",b,"/",N_PERM)
}
perm<-rbindlist(perm_metrics);write_tsv(perm,pfinal("results/table_nested_cv_stratified_permutation_null.tsv"))
obs<-perf[specification=="baseline_age_stratified_5fold"]
metric_names<-c("pearson","spearman","mae_midpoint_years","balanced_age_bracket_accuracy","ordinal_concordance_kendall","exact_bracket_accuracy")
ps<-rbindlist(lapply(metric_names,function(m){higher<-!grepl("mae",m);o<-obs[[m]];nv<-perm[[m]];p<-(1+if(higher)sum(nv>=o)else sum(nv<=o))/(N_PERM+1);data.table(metric=m,observed=o,null_median=median(nv),null_q025=quantile(nv,.025),null_q975=quantile(nv,.975),empirical_p=p,n_permutations=N_PERM,seed=SEED)}))
write_tsv(ps,pfinal("results/table_nested_cv_permutation_summary.tsv"))

audit_steps<-data.table(step=c("prespecified gene universe","fold-specific count filtering","sample normalization","variance feature selection","scaling","dimensional reduction","ridge tuning","outer evaluation"),
  location=c("fixed before adequacy audit","inside every inner and outer training split","sample-local library-size CPM; no cross-sample fit","inside every inner and outer training split","glmnet estimates center/scale from each training split","not used","manual inner five-fold CV with fold-specific filtering/selection/scaling","held-out outer fold only"),
  leakage_status=c("fixed upstream scope","PASS","PASS","PASS","PASS","NOT_APPLICABLE","PASS","PASS"),
  detail=c("25,927 prespecified testable genes; no outcome-guided filtering","counts>=10 in >=20% of the current training donors","log2(CPM+2) uses only each sample's library size","top 2,000 training-fold variances","test donors transformed by training glmnet fit","ridge acts directly on selected genes","fixed lambda grid; validation labels used only for error calculation","no donor contributes to its own preprocessing or tuning"))
write_tsv(audit_steps,pfinal("results/table_nested_cv_leakage_audit.tsv"))

## Figures
pd<-predictions[specification%in%c("baseline_age_stratified_5fold","batch_stratified_5fold","hardy_age_stratified_5fold")]
p_pred<-ggplot(pd,aes(observed,predicted,colour=age_bracket))+geom_abline(slope=1,intercept=0,colour="#BDBDBD")+geom_point(size=.8,alpha=.65)+
  facet_wrap(~specification)+scale_colour_manual(values=palette_age)+labs(x="Bracket midpoint (descriptive only)",y="Nested-CV prediction",colour="Age",title="Leakage-controlled age-information prediction",
    subtitle="Dataset-adequacy analysis; not a validated human ovarian clock")+theme(legend.position="top")
save_pub_r(p_pred,pfinal("figures/14_nested_cv_validation_predictions"),183,82)

pm<-melt(perf,id.vars="specification",measure.vars=c("pearson","spearman","balanced_age_bracket_accuracy","ordinal_concordance_kendall","calibration_slope"),variable.name="metric",value.name="value")
p_perf<-ggplot(pm,aes(specification,value,fill=specification))+geom_col()+facet_wrap(~metric,scales="free_y")+coord_flip()+
  scale_fill_brewer(palette="Set2",guide="none")+labs(x=NULL,y=NULL,title="Nested-CV robustness metrics")+theme(axis.text.y=element_text(size=5.5))
save_pub_r(p_perf,pfinal("figures/15_nested_cv_validation_metrics"),155,105)

p_null<-ggplot(perm,aes(pearson))+geom_histogram(bins=30,fill="#D8D8D8",colour="white")+geom_vline(xintercept=obs$pearson,colour=palette_contract["accent_red"],linewidth=.7)+
  labs(x="Pearson correlation",y="Stratified permutations",title=sprintf("Nested-CV age-label null (B=%d)",N_PERM))
save_pub_r(p_null,pfinal("figures/16_nested_cv_stratified_permutation"),100,72)
write_tsv(pd,pfinal("figures/source_data/nested_cv_validation_predictions.tsv"));write_tsv(pm,pfinal("figures/source_data/nested_cv_validation_metrics.tsv"))

saveRDS(list(predictions=predictions,performance=perf,by_age_bin=bybin,permutation=perm,permutation_summary=ps,leakage_audit=audit_steps),pfinal("derived_data/nested_cv_validation_objects.rds"))
message("Nested-CV validation complete; baseline Pearson=",fmt(obs$pearson),", Spearman=",fmt(obs$spearman))
