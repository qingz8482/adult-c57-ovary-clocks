#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(jsonlite); library(readxl)})
setDTthreads(1)
a <- commandArgs(FALSE); here <- dirname(normalizePath(gsub('~+~',' ',sub('^--file=', '', a[grepl('^--file=', a)][1]),fixed=TRUE)))
dev <- dirname(here); project <- dirname(dev); root <- dirname(dirname(project))
out <- file.path(dev, 'shared_inputs'); dir.create(out, recursive=TRUE, showWarnings=FALSE)
cfg <- fromJSON(file.path(dev,'config/production.json'))
meta_path <- file.path(project,'batch-effect assessment/metadata/sample_expression_metadata.csv')
m <- fread(meta_path); stopifnot(nrow(m)==136L, !anyDuplicated(m$sample_key), uniqueN(m$study)==14L)
wb <- file.path(project,'data_plan/ovary_elasticnetclock_cv_fold_plan.xlsx')
w <- as.data.table(read_excel(wb, sheet='Samples'))
stopifnot(setequal(m$sample_key,w$sample_key)); w <- w[match(m$sample_key,sample_key)]
stopifnot(max(abs(m$age_week-w$age_week))<0.001)
m[,cv_fold:=as.integer(w$cv_fold)]
m[,age_bin:=as.character(cut(age_week,c(-Inf,12,26,52,Inf),labels=c('0-12','12-26','26-52','52+'),right=TRUE))]
dict_path <- file.path(project,'batch-effect assessment/current_input_assessment_2026-10-10/provenance/reconstructed_inputs/identifier_dictionary.tsv')
dict <- fread(dict_path); stopifnot(!anyDuplicated(dict$gene_id))
fwrite(dict,file.path(out,'identifier_dictionary.tsv'),sep='\t')
map_ids <- function(g) { b <- sub('\\.[0-9]+$','',g); hit <- match(b,dict$gene_id); change <- grepl('^ENSMUSG',b)&!is.na(hit); g[change] <- dict$gene_symbol[hit[change]]; g }
read_gene <- function(path,cols) { z<-fread(path); stopifnot(all(cols %in% names(z))); v<-as.matrix(z[,..cols]); storage.mode(v)<-'double'; rownames(v)<-z[[1]]; v }
matrices <- list(); provenance <- list(); per_sample <- list(); mapping <- list()
legacy_root <- file.path(Sys.getenv("OVARY_CLOCK_UPSTREAM_ROOT"), 'outputs/model_ready_counts')
for(ds in unique(m$dataset_folder)) {
 s<-m[dataset_folder==ds]; cols<-s$expression_column; pre<-s$expression_counts_before_filter_path[1]
 if(!is.na(pre)&&pre!='NA'&&file.exists(pre)) {path<-pre; raw<-read_gene(path,cols); kind<-'organized_prefilter_gene_matrix'}
 else if(ds=='GSE154890') {
  path<-file.path(root,'data_input/bulkseq_dataset/GSE154890/raw/GSE154890_rawcount.txt.gz')
  z<-fread(cmd=paste('gzip -cd',shQuote(path)),skip=4)
  rc<-paste0(sub('_.*','',cols),'(raw)'); stopifnot(all(rc%in%names(z)))
  raw<-as.matrix(z[,..rc]); colnames(raw)<-cols
  g<-as.character(z[['Gene Symbol']]); bad<-is.na(g)|g==''|g=='---';g[bad]<-as.character(z[['Gene ID']])[bad]
  stopifnot(!anyNA(g),all(nzchar(g))); rownames(raw)<-g;kind<-'deposited_prefilter_abundance_not_verified_raw'
 } else {
  path<-file.path(legacy_root,paste0(ds,'.model_ready_counts.csv'));z<-fread(path)
  stopifnot(all(cols%in%names(z)));raw<-as.matrix(z[,..cols]);g<-as.character(z$gene_symbol)
  bad<-is.na(g)|tolower(trimws(g))%in%c('','na','nan','none','null');g[bad]<-as.character(z$gene_id)[bad]
  stopifnot(!anyNA(g),all(nzchar(g)));rownames(raw)<-g;kind<-'upstream_prefilter_model_ready_quantities'
 }
 stopifnot(is.numeric(raw),all(is.finite(raw)),all(raw>=0),all(colSums(raw)>0))
 # Verify recovered source agrees with the current filtered export BEFORE identifier remapping.
 agg<-rowsum(raw,rownames(raw),reorder=FALSE)
 current<-read_gene(s$expression_counts_path[1],cols)
 stopifnot(all(rownames(current)%in%rownames(agg)))
 err<-max(abs(agg[rownames(current),,drop=FALSE]-current)/pmax(1,current))
 stopifnot(err<1e-6)
 keys<-map_ids(rownames(agg)); mapping[[ds]]<-data.table(study=ds,source_gene=rownames(agg),model_gene=keys)
 den<-colSums(agg);counts<-rowsum(agg,keys,reorder=FALSE)
 lx<-log1p(sweep(counts,2,den,'/')*1e6);colnames(lx)<-s$sample_key;matrices[[ds]]<-lx
 oldden<-s$normalization_denominator_value
 per_sample[[ds]]<-data.table(sample_key=s$sample_key,study=ds,prefilter_denominator=den,legacy_denominator=oldden,
   denominator_ratio=den/oldden,source_column=cols,n_measured_genes=nrow(counts))
 provenance[[ds]]<-data.table(study=ds,n_profiles=nrow(s),source_path=normalizePath(path),source_kind=kind,
  n_source_rows=nrow(raw),n_gene_rows=nrow(counts),filtered_source_max_relative_error=err,
  provenance_status=s$expression_provenance_status[1],normalization=cfg$normalization)
 message(ds,': ',nrow(counts),' genes; ',nrow(s),' profiles; source agreement ',signif(err,3))
}
genes<-sort(unique(unlist(lapply(matrices,rownames))));x<-matrix(NA_real_,length(genes),nrow(m),dimnames=list(genes,m$sample_key))
for(z in matrices) x[rownames(z),colnames(z)]<-z
stopifnot(all(colSums(is.finite(x))>10000),!anyDuplicated(rownames(x)),identical(colnames(x),m$sample_key))
fwrite(rbindlist(provenance),file.path(out,'input_provenance.tsv'),sep='\t')
fwrite(rbindlist(per_sample),file.path(out,'normalization_comparison.tsv'),sep='\t')
fwrite(rbindlist(mapping),file.path(out,'gene_identifier_mapping.tsv'),sep='\t')
# Generate unique metadata-only constrained study assignments; fold numbers are canonicalized.
set.seed(cfg$seed); studies<-sort(unique(m$study));anchor<-c('2024_Cell_Winkler_FRT_scRNA_spatial','NIHMS1554107','2022_eLife_Morris_cycling_murine_ovary_atlas')
sm<-sapply(studies,function(s) as.numeric(table(factor(m[study==s,age_bin],levels=c('0-12','12-26','26-52','52+')))))
mm<-sapply(studies,function(s) as.numeric(table(factor(m[study==s,modality],levels=sort(unique(m$modality))))))
features<-rbind(sm,mm);target<-rowSums(features)/4
pool<-list();scores<-numeric(); seen<-new.env(hash=TRUE)
main<-m$cv_fold[match(studies,m$study)];canon<-function(v) match(v,unique(v))
assign(paste(canon(main),collapse=''),TRUE,envir=seen)
for(i in seq_len(cfg$repeat_random_candidates)) {
 f<-sample(rep(1:4,length.out=length(studies))); if(uniqueN(f[match(anchor,studies)])!=3) next
 key<-paste(canon(f),collapse='');if(exists(key,envir=seen,inherits=FALSE)) next
 assign(key,TRUE,envir=seen);bal<-sapply(1:4,function(k)rowSums(features[,f==k,drop=FALSE]))
 score<-sum((bal-target)^2/pmax(target,1));pool[[length(pool)+1L]]<-canon(f);scores<-c(scores,score)
}
stopifnot(length(pool)>=cfg$repeats)
folds<-m[,.(sample_key,study,modality,age_week,age_bin,cv_fold)]
for(j in seq_len(cfg$repeats)) folds[,(paste0('repeat_',sprintf('%02d',j))):=pool[[order(scores)[j]]][match(study,studies)]]
# Audit all folds and preserve independent read-only metadata snapshot.
for(col in c('cv_fold',grep('^repeat_',names(folds),value=TRUE))) {
 stopifnot(all(folds[,uniqueN(get(col)),by=study]$V1==1),setequal(folds[[col]],1:4))
 if(col!='cv_fold') stopifnot(uniqueN(folds[match(anchor,study),get(col)])==3)
}
fwrite(folds,file.path(out,'fold_assignments.csv'));fwrite(m,file.path(out,'sample_metadata.csv'))
modelmeta<-copy(m);ns<-rbindlist(per_sample);pr<-rbindlist(provenance)
modelmeta[,`:=`(model_expression_path=file.path(out,'prepared_expression.rds'),model_expression_object='expression',
 model_expression_column=sample_key,model_expression_normalization=cfg$normalization,
 model_CPM_denominator=ns$prefilter_denominator[match(sample_key,ns$sample_key)],
 model_count_source_path=pr$source_path[match(study,pr$study)],model_count_source_kind=pr$source_kind[match(study,pr$study)],
 model_missing_gene_policy='training mean for missing held-out genes; training requires measurement in every profile',
 model_inverse_stratum_weight_definition=cfg$inverse_stratum_definition)]
fwrite(modelmeta,file.path(dev,'metadata_allsamples/sample_metadata.csv'))
saveRDS(list(expression=x,metadata=m,folds=folds,normalization=cfg$normalization),file.path(out,'prepared_expression.rds'),compress=FALSE)
paths<-c(meta_path,wb,dict_path,unique(rbindlist(provenance)$source_path))
sha<-vapply(paths,function(p)strsplit(system2('sha256sum',shQuote(p),stdout=TRUE),' +')[[1]][1],character(1))
fwrite(data.table(path=paths,sha256=sha),file.path(out,'source_checksums.tsv'),sep='\t')
writeLines(c('Prepared expression: ln(1 + 1e6 * gene quantity / full available prefilter column sum).',
'No TMM, DESeq2 or batch correction was fitted here. Input provenance remains mixed.',
'Full available upstream gene universe is not a guarantee that authors supplied every genomic feature.',
'Missing gene rows are NA (unmeasured/absent upstream), not biological zero.',
'Identifier union stores inputs only. Each training split requires genes measured in all its profiles,',
'CPM >=1 in >=50%, and positive variance. Held-out missing genes use frozen training means.',
'Legacy denominator differences are recorded in normalization_comparison.tsv; source inputs are unchanged.',
'All expression-dependent filtering, gene selection, scaling and PCA occur inside nested CV.',
'All 136 profiles retained despite primary_reconstructed_training_ready=FALSE.'),file.path(out,'README.txt'))
message('PREPARED: ',nrow(x),' union genes x ',ncol(x),' profiles; ',cfg$repeats,' repeated assignments.')
