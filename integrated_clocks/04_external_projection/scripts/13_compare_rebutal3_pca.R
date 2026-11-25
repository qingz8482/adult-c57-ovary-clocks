#!/usr/bin/env Rscript
# Version B-PCA is unchanged between releases: its projections and tests must reproduce rebutal3 exactly.
suppressPackageStartupMessages(library(data.table))
a<-commandArgs(FALSE);here<-dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE)));out<-dirname(here)
r3<-file.path(dirname(dirname(out)),'rebutal3_batcheffects removal_CV_LOO/lockedmodel_application')
mo<-'VersionB-PCA/glmnet'
p4<-fread(file.path(out,'external_projection_predictions.tsv'))[model==mo];p3<-fread(file.path(r3,'external_projection_predictions.tsv'))[model==mo]
pp<-merge(p3[,.(sample_key,old=predicted_age_week,old_cov=missing_gene_fraction,old_train=locked_training_sample)],p4[,.(sample_key,new=predicted_age_week,new_cov=missing_gene_fraction,new_train=locked_training_sample)],by='sample_key')
t4<-fread(file.path(out,'external_projection_contrast_statistics.tsv'))[model==mo];t3<-fread(file.path(r3,'external_projection_contrast_statistics.tsv'))[model==mo]
tt<-merge(t3[,.(analysis_contrast_id,e3=effect_weeks,p3=welch_p,l3=ci95_low,h3=ci95_high)],t4[,.(analysis_contrast_id,e4=effect_weeks,p4=welch_p,l4=ci95_low,h4=ci95_high)],by='analysis_contrast_id')
c4<-fread(file.path(out,'tables/pca_projection_pc_contributions.tsv'));c3<-fread(file.path(r3,'tables/pca_projection_pc_contributions.tsv'))
pcs<-grep('^PC[0-9]+$',names(c3),value=TRUE);cc<-merge(c3,c4,by=c('sample_key','model'),suffixes=c('.3','.4'))
pcerr<-max(sapply(pcs,function(k)max(abs(cc[[paste0(k,'.3')]]-cc[[paste0(k,'.4')]]))))
res<-data.table(check=c('281 PCA predictions matched','max |prediction difference| (weeks)','missing-gene fractions identical','training-overlap flags identical',
 '44 PCA contrasts matched','max |effect difference| (weeks)','max |Welch P difference|','max |CI bound difference|','max |PC contribution difference| (weeks)'),
 value=c(nrow(pp),max(abs(pp$new-pp$old)),all(pp$old_cov==pp$new_cov),all(pp$old_train==pp$new_train),nrow(tt),max(abs(tt$e4-tt$e3)),max(abs(tt$p4-tt$p3)),max(abs(c(tt$l4-tt$l3,tt$h4-tt$h3))),pcerr))
fwrite(res,file.path(out,'provenance/rebutal3_pca_reuse_check.tsv'),sep='\t')
print(res)
stopifnot(nrow(pp)==281,max(abs(pp$new-pp$old))<1e-10,all(pp$old_cov==pp$new_cov),all(pp$old_train==pp$new_train),nrow(tt)==44,max(abs(tt$e4-tt$e3))<1e-10,max(abs(tt$p4-tt$p3))<1e-12,pcerr<1e-10)
message('PASS: Version B-PCA application results reproduce rebutal3 exactly.')
