#!/usr/bin/env Rscript
# Version B-PCA is unchanged between releases. Deterministic human projections and age associations must
# reproduce rebutal3 exactly; permutation p-values and Stage A (union gene background) are reported, not required equal.
suppressPackageStartupMessages(library(data.table))
a<-commandArgs(FALSE);here<-dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE)));out<-dirname(here)
r3<-file.path(dirname(dirname(out)),'rebutal3_batcheffects removal_CV_LOO/generalization_tohuman')
rd<-function(base,f)fread(file.path(base,f));mo<-'VersionB-PCA';rows<-list()
add<-function(check,value,exact)rows[[length(rows)+1L]]<<-data.table(check=check,value=value,must_be_exact=exact,passed=!exact||value<1e-9)
p3<-rd(r3,'tables/frozen_predictions.tsv')[model==mo];p4<-rd(out,'tables/frozen_predictions.tsv')[model==mo]
pp<-merge(p3,p4,by=c('cohort','sample_id'),suffixes=c('.3','.4'));stopifnot(nrow(pp)==nrow(p3),nrow(pp)==201)
add('PCA human raw outputs max |difference|',max(abs(pp$raw_linear_predictor.3-pp$raw_linear_predictor.4)),TRUE)
add('PCA standardized scores max |difference|',max(abs(pp$cross_species_mouse_ovarian_clock_score.3-pp$cross_species_mouse_ovarian_clock_score.4)),TRUE)
add('PCA gene coverage max |difference|',max(abs(pp$absolute_weight_coverage.3-pp$absolute_weight_coverage.4)),TRUE)
c3<-rd(r3,'tables/age_correlations.tsv')[model==mo];c4<-rd(out,'tables/age_correlations.tsv')[model==mo]
num<-intersect(names(c3)[vapply(c3,is.numeric,logical(1))],names(c4));keys<-setdiff(intersect(names(c3),names(c4)),num)
cc<-merge(c3,c4,by=keys,suffixes=c('.3','.4'));stopifnot(nrow(cc)==nrow(c3))
add('PCA age correlations max |difference|',max(abs(as.matrix(cc[,paste0(num,'.3'),with=FALSE])-as.matrix(cc[,paste0(num,'.4'),with=FALSE])),na.rm=TRUE),TRUE)
s3<-rd(r3,'tables/age_association_sensitivities.tsv')[model==mo];s4<-rd(out,'tables/age_association_sensitivities.tsv')[model==mo]
ss<-merge(s3[,.(cohort,analysis,e3=estimate,p3=p)],s4[,.(cohort,analysis,e4=estimate,p4=p)],by=c('cohort','analysis'));stopifnot(nrow(ss)==nrow(s3))
add('PCA adjusted age-association estimates max |difference|',max(abs(ss$e3-ss$e4),na.rm=TRUE),TRUE)
add('PCA adjusted age-association p max |difference|',max(abs(ss$p3-ss$p4),na.rm=TRUE),TRUE)
n3<-rd(r3,'tables/stageB_null_benchmarks.tsv')[model==mo];n4<-rd(out,'tables/stageB_null_benchmarks.tsv')[model==mo]
nn<-merge(n3[,.(null_type,o3=observed,p3=empirical_one_sided_p)],n4[,.(null_type,o4=observed,p4=empirical_one_sided_p)],by='null_type')
add('PCA Stage B null observed statistics max |difference|',max(abs(nn$o3-nn$o4)),TRUE)
add('PCA Stage B permutation p max |difference| (Monte Carlo; shared RNG stream)',max(abs(nn$p3-nn$p4)),FALSE)
a3<-rd(r3,'tables/stageA_direction_conservation.tsv');a4<-rd(out,'tables/stageA_direction_conservation.tsv')
aa<-merge(a3[grepl('pca',gene_set),.(gene_set,n3=n_genes,r3=spearman_signed_t,p3=empirical_p)],a4[,.(gene_set,n4=n_genes,r4=spearman_signed_t,p4=empirical_p)],by='gene_set')
add('PCA Stage A sets: max |gene-count difference| (union background)',max(abs(aa$n3-aa$n4)),FALSE)
add('PCA Stage A sets: max |Spearman difference| (union background)',max(abs(aa$r3-aa$r4)),FALSE)
add('PCA Stage A sets: max |empirical p difference|',max(abs(aa$p3-aa$p4)),FALSE)
res<-rbindlist(rows);fwrite(res,file.path(out,'provenance/rebutal3_pca_reuse_check.tsv'),sep='\t');print(res)
stopifnot(all(res$passed));message('PASS: Version B-PCA human projections and age associations reproduce rebutal3; Monte Carlo/Stage A shifts reported.')
