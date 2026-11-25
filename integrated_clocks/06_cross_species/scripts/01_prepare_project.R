source(file.path(dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),fixed=TRUE))),'common.R'))
message('Loading current mouse reference and verified human inputs')
refpath<-file.path(dev,'shared_inputs/prepared_expression.rds');ref<-readRDS(refpath)
stopifnot(length(models[['HVG-ElasticNet']]$training_keys)==110,identical(models[['VersionB-PCA']]$training_keys,ref$metadata$sample_key),all(vapply(models,function(m)all(m$training_keys%in%ref$metadata$sample_key)&&m$target=='age_week'&&m$calibration=='none'&&length(m$representation$genes)==3000,logical(1))))
orthpath<-file.path(old,'derived_data/orthology/mouse_human_orthology_primary_one2one.tsv');orth<-fread(orthpath)
# Exclude ambiguous SYMBOL keys instead of silently choosing the first row.
orth<-unique(orth);orth[,mouse_n:=.N,by=mouse_gene_symbol];orth[,human_n:=.N,by=human_gene_id_stable]
wt(orth[mouse_n>1|human_n>1],'provenance/ambiguous_orthology_excluded.tsv');orth<-orth[mouse_n==1&human_n==1]
stopifnot(all(orth$orthology_confidence==1),all(orth$orthology_type=='ortholog_one2one'))
wt(orth,'inputs/orthology_one2one.tsv')
gtexpath<-file.path(old,'derived_data/gtex/gtex_ovary_prepared.rds');gt<-readRDS(gtexpath)
rawpath<-file.path(old,'derived_data/gtex/raw/gene_reads_adult_gtex_v11_ovary.gct.gz')
gct<-fread(cmd=paste('gzip -cd',shQuote(rawpath)),skip=2);ids<-gt$metadata$sample_id
raw<-as.matrix(gct[,..ids]);rownames(raw)<-sub('\\.[0-9]+$','',gct$Name);raw<-rowsum(raw,rownames(raw),reorder=FALSE)
stopifnot(max(abs(raw[rownames(gt$counts_filtered),ids]-gt$counts_filtered))==0)
den<-colSums(raw);stopifnot(max(abs(den-gt$metadata$library_size))<1e-6)
gm<-copy(gt$metadata);gm[,`:=`(cohort='GTEx',age_year=age_midpoint_sensitivity,age_order=as.integer(age_bracket),rin_z=zsafe(rin),ischemia_z=zsafe(ischemic_time),log_library_z=zsafe(log10(library_size)),age_value_type='public_age_bracket_midpoint_for_sensitivity')]
gm[,age_decade:=age_year/10]
gsexpath<-file.path(gseold,'input/GSE202601_sample_pseudobulk_counts_before_filter.tsv');gr<-readmatrix(gsexpath)
gsempath<-file.path(gseold,'input/GSE202601_projection_sample_metadata.tsv');sm<-fread(gsempath)
stopifnot(nrow(sm)==8,sum(sm$n_nuclei)==42568,setequal(colnames(gr),sm$gsm));gr<-gr[,sm$gsm,drop=FALSE]
stopifnot(all(colSums(gr)==sm$pseudobulk_library_size_raw),all(gr>=0),all(gr==round(gr)))
sm[,`:=`(cohort='GSE202601',donor_id=sample_id,sample_id=projection_sample_id,age_decade=age_year/10,age_value_type='exact_years',library_size=as.numeric(colSums(gr)),normalization='ln1p_CPM_full_available_prefilter_column_sum')]
colnames(gr)<-sm$sample_id
cohorts<-list(GTEx=list(raw=raw,x=log1p(sweep(raw,2,den,'/')*1e6),meta=gm,key='human_gene_id_stable',filter=rownames(gt$counts_filtered)),GSE202601=list(raw=gr,x=log1p(sweep(gr,2,colSums(gr),'/')*1e6),meta=sm,key='human_gene_symbol',filter=rownames(gr)[rowSums(gr>=10)>=2]))
wt(rbindlist(lapply(cohorts,function(o)o$meta),fill=TRUE),'inputs/human_sample_metadata.tsv')
wt(rbindlist(lapply(cohorts,function(o)data.table(cohort=o$meta$cohort,sample_id=o$meta$sample_id,prefilter_total=colSums(o$raw),n_source_genes=nrow(o$raw),n_filter_pass=length(o$filter)))),'provenance/normalization_audit.tsv')
manifest<-list();preds<-list();maps<-list();objs<-list();checks<-list();pcs<-list();supports<-list();mousepred<-list()
for(nm in names(models)){
 m<-models[[nm]];r<-m$representation;g<-r$genes;w<-m$coefficients$coefficient;ws<-m$coefficients$standardized_coefficient
 stopifnot(identical(g,m$coefficients$gene),all(is.finite(ref$expression[g,m$training_keys])))
 tr<-t(ref$expression[g,m$training_keys]);trz<-sweep(sweep(tr,2,r$center,'-'),2,r$scale,'/')
 mousepred[[nm]]<-data.table(sample_id=m$training_keys,model=nm,age_week=ref$metadata$age_week[match(m$training_keys,ref$metadata$sample_key)],predict_model(m,ref$expression[,m$training_keys,drop=FALSE])$pred)
 setnames(mousepred[[nm]],'V4','mouse_output')
 if(nm=='VersionB-PCA'){trlow<-trz%*%r$rotation;diagfit<-NULL}else{diagfit<-prcomp(trz,center=FALSE,scale.=FALSE,rank.=10);trlow<-trz%*%diagfit$rotation}
 covr<-cov(trlow);diag(covr)<-diag(covr)+1e-6;cen<-colMeans(trlow);rd<-mahalanobis(trlow,cen,covr);qs<-quantile(rd,c(.95,.99))
 manifest[[nm]]<-data.table(model=nm,n_input_genes=length(g),n_nonzero_gene_weights=sum(w!=0),active_PCs=if(is.null(m$pc_coefficients))NA_character_ else paste(paste0('PC',which(m$pc_coefficients!=0)),collapse=';'),n_mouse_profiles=nrow(tr),normalization=m$normalization,calibration=m$calibration)
 for(cn in names(cohorts)){
  o<-cohorts[[cn]];mp<-orth[match(g,mouse_gene_symbol)];key<-mp[[o$key]]
  if(o$key=='human_gene_symbol'){dups<-orth[, .N,by=human_gene_symbol][N>1|is.na(human_gene_symbol)|human_gene_symbol=='',human_gene_symbol];key[key%in%dups]<-NA_character_}
  measured<-!is.na(key)&key%in%rownames(o$x);passed<-measured&key%in%o$filter
  xx<-matrix(NA_real_,length(g),ncol(o$x),dimnames=list(g,colnames(o$x)));xx[measured,]<-o$x[key[measured],,drop=FALSE]
  pr<-predict_model(m,xx);stopifnot(identical(pr$sample_key,o$meta$sample_id))
  xi<-t(xx);ii<-which(!is.finite(xi),arr.ind=TRUE);xi[ii]<-r$center[ii[,2]]
  zz<-sweep(sweep(xi,2,r$center,'-'),2,r$scale,'/');yy<-drop(xi%*%w)+m$intercept
  err<-max(abs(yy-pr$pred));stopifnot(err<1e-8,all(zz[,!measured,drop=FALSE]==0))
  low<-zz%*%if(nm=='VersionB-PCA')r$rotation else diagfit$rotation
  distance<-mahalanobis(low,cen,covr)
  fxx<-xx;fxx[!passed,]<-NA_real_;filtered<-predict_model(m,fxx)$pred
  cf<-sum(abs(w[measured]))/sum(abs(w))
  p<-copy(o$meta);p[,`:=`(model=nm,raw_linear_predictor=pr$pred,mouse_week_equivalent_output_audit_only=pr$pred,cross_species_mouse_ovarian_clock_score=zsafe(pr$pred),filtered_input_score_sensitivity=zsafe(filtered),n_model_genes=length(g),n_measured_genes=sum(measured),raw_gene_coverage=mean(measured),absolute_weight_coverage=cf,projection_ood_distance=distance,mouse_ood_q95=qs[1],mouse_ood_q99=qs[2],fraction_features_abs_z_gt3=rowMeans(abs(zz[,measured,drop=FALSE])>3),fraction_features_abs_z_gt5=rowMeans(abs(zz[,measured,drop=FALSE])>5),fraction_features_outside_mouse_range=rowMeans(sweep(zz[,measured,drop=FALSE],2,apply(trz[,measured,drop=FALSE],2,min),'<')|sweep(zz[,measured,drop=FALSE],2,apply(trz[,measured,drop=FALSE],2,max),'>')),nearest_mouse_reference_distance=vapply(seq_len(nrow(low)),function(i)min(sqrt(rowSums(sweep(trlow,2,low[i,],'-')^2))),numeric(1)))]
  preds[[paste(cn,nm)]]<-p
  maps[[paste(cn,nm)]]<-data.table(cohort=cn,model=nm,model_gene_order=seq_along(g),mouse_gene_symbol=g,human_gene_id_stable=mp$human_gene_id_stable,human_gene_symbol=mp$human_gene_symbol,lookup_key=key,one2one_mapped=!is.na(key),measured=measured,pass_historical_expression_filter=passed,raw_coefficient=w,standardized_coefficient=ws,training_center=r$center,training_scale=r$scale,missing_policy=ifelse(measured,'observed_ln1p_CPM_including_measured_zero','frozen_training_mean'))
  checks[[paste(cn,nm)]]<-data.table(cohort=cn,model=nm,max_formula_error=err,missing_to_zero=all(zz[,!measured,drop=FALSE]==0),n=length(yy))
  objs[[paste(cn,nm,sep='__')]]<-list(x=xx,imputed=xi,z=zz,reference_z=trz,human_low=low,reference_low=trlow,reference_ood=rd,map=maps[[paste(cn,nm)]],meta=o$meta,score=p$cross_species_mouse_ovarian_clock_score)
  for(q in seq_along(qs))supports[[paste(cn,nm,q)]]<-data.table(cohort=cn,model=nm,mouse_quantile=names(qs)[q],threshold=qs[q],n_within=sum(distance<=qs[q]),n_total=length(distance))
  if(nm=='VersionB-PCA'){
   cc<-sweep(low,2,m$pc_coefficients,'*');stopifnot(max(abs(rowSums(cc)+m$model_intercept-pr$pred))<1e-8)
   pcs[[cn]]<-rbindlist(lapply(seq_len(ncol(low)),function(j)data.table(cohort=cn,sample_id=o$meta$sample_id,PC=colnames(low)[j],pc_score=low[,j],coefficient=m$pc_coefficients[j],contribution=cc[,j],intercept=m$model_intercept,age_year=o$meta$age_year)))
  }
 }
}
wt(rbindlist(manifest),'provenance/model_manifest.tsv');wt(rbindlist(preds,fill=TRUE),'tables/frozen_predictions.tsv');wt(rbindlist(maps),'tables/model_gene_coverage.tsv');wt(rbindlist(checks),'provenance/prediction_formula_checks.tsv');wt(rbindlist(pcs),'tables/pca_scores_contributions.tsv');wt(rbindlist(supports),'tables/ood_common_support.tsv');wt(rbindlist(mousepred),'tables/mouse_reference_predictions.tsv')
saveRDS(objs,file.path(out,'inputs/projection_objects.rds'))
# Stage A: recompute mouse effects using the CURRENT development samples and batch design.
eligible<-Reduce(union,lapply(models,function(m)m$representation$filtered_genes));eligible<-eligible[rowSums(is.na(ref$expression[eligible,,drop=FALSE]))==0];mx<-ref$expression[eligible,,drop=FALSE]
stopifnot(!anyNA(mx));md<-ref$metadata;design<-model.matrix(~age_week+factor(source_batch_for_feature_selection),data=md);stopifnot(qr(design)$rank==ncol(design))
fit<-eBayes(lmFit(mx,design));k<-match('age_week',colnames(design))
mt<-data.table(mouse_gene_symbol=rownames(mx),mouse_age_beta=fit$coefficients[,k],mouse_age_t=fit$t[,k],mouse_age_p=fit$p.value[,k],mouse_standardized_age_beta=fit$coefficients[,k]*sd(md$age_week)/apply(mx,1,sd))
mt[,mouse_age_fdr:=p.adjust(mouse_age_p,'BH')]
master<-merge(mt,orth,by='mouse_gene_symbol');ae<-copy(gt$age_effects);ae[,human_gene_symbol:=NULL];master<-merge(master,ae,by='human_gene_id_stable',all.x=TRUE)
master[,human_standardized_age_beta:=human_age_beta_per_decade*sd(gm$age_decade)/sqrt(human_expression_variance)]
for(nm in names(models)){m<-models[[nm]];prefix<-if(nm=='HVG-ElasticNet')'hvg' else 'pca';master[,(paste0(prefix,'_selected')):=mouse_gene_symbol%in%m$representation$genes];master[,(paste0(prefix,'_weight')):=m$coefficients$coefficient[match(mouse_gene_symbol,m$coefficients$gene)]]}
wt(master,'tables/stageA_gene_effect_master.tsv');wt(data.table(mouse_gene_symbol=intersect(models[[1]]$representation$genes,models[[2]]$representation$genes)),'tables/shared_model_genes.tsv')
saveRDS(list(cohorts=cohorts,master=master,mouse_metadata=md),file.path(out,'inputs/analysis_inputs.rds'))
paths<-c(refpath,orthpath,gtexpath,rawpath,gsexpath,gsempath,file.path(old,'derived_data/hpa/hpa_ovary_prepared.rds'),file.path(release,c('HVG-ElasticNet-exclude-single-age/model.rds','VersionB-PCA-glmnet/model.rds','scripts/core.R','checksums.tsv')))
sha<-vapply(paths,function(p)strsplit(system2('sha256sum',shQuote(p),stdout=TRUE),' +')[[1]][1],character(1));wt(data.table(path=normalizePath(paths),sha256=sha),'provenance/source_checksums.tsv')
capture.output(sessionInfo(),file=file.path(out,'provenance/R_sessionInfo.txt'))
message('PASS: 402 human projections, frozen formula replay, current Stage A effects.')
