#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
setDTthreads(1)
a<-commandArgs(FALSE);here<-dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE)))
out<-dirname(here); project<-dirname(out);root<-dirname(dirname(project))
dev<-file.path(project,'model_development');release<-file.path(project,'lockedmodel')
old<-file.path(root,'bulk-pseudobulk_clock');bulk<-file.path(root,'data_input/bulkseq_dataset')
legacy<-file.path(Sys.getenv("OVARY_CLOCK_UPSTREAM_ROOT"), 'outputs/model_ready_counts')
source(file.path(release,'scripts/core.R'))
source(file.path(bulk,'prepare_external_bulk_counts_common.R'))
for(d in c('inputs','inputs/cohort_cache','provenance','tables','diagnostics','plots','logs'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
write_tsv<-function(x,p)fwrite(x,file.path(out,p),sep='\t',na='NA')
sup<-file.path(old,'methods_tables/publication_supplementary_tables')
sp<-file.path(sup,'Supplementary_Table_3b_external_projection_samples.tsv')
cp<-file.path(sup,'Supplementary_Table_4_external_projection_contrasts.tsv')
mp<-file.path(old,'methods_tables/19_tage_sample_contrast_membership.tsv')
m0<-fread(sp);m<-m0[dataset_folder!='GSE154890']
stopifnot(nrow(m)==281,!anyDuplicated(m$sample_key),uniqueN(m$dataset_name)==20)
m[,modality_id:=fifelse(modality=='bulk RNA-seq','bulk_rnaseq','scrna_pseudobulk')]
m[,dataset_cohort:=paste(dataset_folder,modality_id,sep='__')]
# Rebutal4 release: HVG-ElasticNet trained without single-age studies (110 profiles); Version B-PCA unchanged (136).
models<-list('HVG-ElasticNet'=readRDS(file.path(release,'HVG-ElasticNet-exclude-single-age/model.rds')),
 'VersionB-PCA/glmnet'=readRDS(file.path(release,'VersionB-PCA-glmnet/model.rds')))
training<-fread(file.path(dev,'shared_inputs/sample_metadata.csv'))
stopifnot(all(vapply(models,function(z)all(z$training_keys%in%training$sample_key),logical(1))),length(models[['HVG-ElasticNet']]$training_keys)==110,setequal(models[['VersionB-PCA/glmnet']]$training_keys,training$sample_key))
m[,`:=`(locked_training_sample=sample_key%in%training$sample_key,locked_training_study=dataset_folder%in%training$dataset_folder)]
m[,current_development_relationship:=fifelse(locked_training_sample,'sample used in locked model development',fifelse(locked_training_study,'new sample from a development study','external study'))]
dictpath<-file.path(dev,'shared_inputs/identifier_dictionary.tsv');dict<-fread(dictpath)
map_ids<-function(g){b<-sub('\\.[0-9]+$','',g);i<-match(b,dict$gene_id);use<-grepl('^ENSMUSG',b)&!is.na(i);g[use]<-dict$gene_symbol[i[use]];g}
base_meta<-fread(file.path(old,'metadata/bulk_pseudobulk_sample_metadata.csv'))
oldalign<-fread(file.path(old,'results/integrated_bulk_pseudobulk_sample_expression_alignment.tsv'))
runpath<-file.path(Sys.getenv("OVARY_CLOCK_UPSTREAM_ROOT"), 'outputs/fastq_raw_counts/fastq_to_counts_metadata.csv')
runmeta<-fread(runpath)[dataset_id=='GSE228975']
runmeta[,sample_id:=paste0(fifelse(grepl('cKO',sample_title),'cKO','ctl'),sub('.*rep','',sample_title))]
paths<-c(sp,cp,mp,dictpath,runpath,file.path(project,'MODEL_SELECTION.json'),file.path(release,'checksums.tsv'),file.path(old,'metadata/bulk_pseudobulk_sample_metadata.csv'),file.path(old,'results/integrated_bulk_pseudobulk_sample_expression_alignment.tsv'))
annotation_map<-function(ds){p<-file.path(bulk,ds,'input',paste0(ds,'_gene_annotation.csv'));paths<<-c(paths,p);z<-fread(p);z[,.(gene_id=unlist(strsplit(gene_id,';',fixed=TRUE))),by=gene_id_or_symbol]}
read_matrix<-function(p){z<-fread(p);x<-as.matrix(z[,-1,with=FALSE]);storage.mode(x)<-'double';rownames(x)<-z[[1]];x}
all_genes<-sort(unique(unlist(lapply(models,function(z)z$representation$genes))))
x<-matrix(NA_real_,length(all_genes),nrow(m),dimnames=list(all_genes,m$sample_key))
alignments<-list();prov<-list();norms<-list();mapping<-list()
for(cohort in unique(m$dataset_cohort)){
 s<-m[dataset_cohort==cohort];ds<-s$dataset_folder[1];mo<-s$modality_id[1];kind<-''
 message('Preparing ',cohort,' (',nrow(s),')')
 cachepath<-file.path(out,'inputs/cohort_cache',paste0(cohort,'.rds'))
 if(file.exists(cachepath)){
  cache<-readRDS(cachepath)
  valid<-identical(cache$sample_keys,s$sample_key)&&identical(cache$dictionary_md5,unname(tools::md5sum(dictpath)))&&identical(cache$source_info,file.info(cache$paths)[,c('size','mtime')])
  if(valid){x[,s$sample_key]<-cache$x;alignments[[cohort]]<-cache$alignment;norms[[cohort]]<-cache$norm;prov[[cohort]]<-cache$provenance;mapping[[cohort]]<-cache$mapping;paths<-c(paths,cache$paths);message('  verified preparation cache');next}
 }
 path_start<-length(paths)
 if(mo=='scrna_pseudobulk'){
  p<-file.path(root,'data_input/aggregated_pseudobulk',ds,'input_genelevel/sample_gene_counts_before_expression_filter.tsv');raw<-read_matrix(p);kind<-'prefilter sample-level pseudobulk counts'
 }else if(ds%in%c('GSE106339','GSE249151')){
  cfg<-dataset_config(ds);p<-file.path(bulk,ds,'raw',cfg$raw_file);z<-read_source_counts(p,cfg$parser,cfg$sample)
  g<-z$gene_symbol;bad<-is.na(g)|g==''|g=='---';g[bad]<-z$gene_id[bad];raw<-z$counts;rownames(raw)<-g;kind<-'deposited prefilter gene counts'
 }else if(ds%in%c('GSE128240','GSE98497')){
  q<-fread(file.path(bulk,ds,'input',paste0(ds,'_salmon_mapping_qc.tsv')));am<-annotation_map(ds)
  p<-file.path(bulk,ds,'processed/salmon',q$run,'quant.sf')
  ll<-lapply(seq_len(nrow(q)),function(i){v<-fread(p[i],select=c('Name','NumReads'));v[,gene_id:=sub('\\.[0-9]+$','',tstrsplit(Name,'|',fixed=TRUE)[[2]])];v<-v[,.(count=sum(NumReads)),by=gene_id];v[,sample:=q$sample_name[i]];v})
  z<-dcast(rbindlist(ll),gene_id~sample,value.var='count');raw<-as.matrix(z[,-1,with=FALSE]);g<-am$gene_id_or_symbol[match(z$gene_id,am$gene_id)];g[is.na(g)]<-z$gene_id[is.na(g)];rownames(raw)<-g;kind<-'prefilter Salmon counts aggregated from all transcripts'
 }else if(ds=='GSE235021'){
  p<-list.files(file.path(bulk,ds,'raw'),pattern='txt.gz$',full.names=TRUE);am<-annotation_map(ds)
  ll<-lapply(p,function(f){z<-fread(cmd=paste('gzip -cd',shQuote(f)));sn<-toupper(gsub('_','-',names(z)[2]));list(g=z[[1]],v=z[[2]],sn=sn)})
  stopifnot(all(vapply(ll,function(z)identical(z$g,ll[[1]]$g),logical(1))))
  raw<-do.call(cbind,lapply(ll,`[[`,'v'));colnames(raw)<-vapply(ll,`[[`,character(1),'sn');g<-am$gene_id_or_symbol[match(ll[[1]]$g,am$gene_id)];g[is.na(g)]<-ll[[1]]$g[is.na(g)];rownames(raw)<-g;kind<-'deposited prefilter featureCounts'
 }else if(ds=='PRJNA1002222'){
  p<-file.path(bulk,ds,'input/PRJNA1002222_GENCODE_M39_gene_id_raw_estimated_counts.csv');raw<-read_matrix(p);am<-annotation_map(ds);g<-am$gene_id_or_symbol[match(rownames(raw),am$gene_id)];g[is.na(g)]<-rownames(raw)[is.na(g)];rownames(raw)<-g;kind<-'prefilter Salmon gene counts'
 }else if(ds=='PRJNA787028'){
  p<-file.path(bulk,ds,'processed/PRJNA787028_STAR_column3_stranded_raw_counts.csv');raw<-read_matrix(p);kind<-'prefilter STAR column-3 strand-specific gene counts'
 }else if(ds=='GSE329024'){
  p<-file.path(bulk,ds,'input/GSE329024_gene_counts_before_expression_filter.csv');raw<-read_matrix(p);kind<-'organized prefilter gene counts'
 }else{
  p<-file.path(legacy,paste0(ds,'.model_ready_counts.csv'));z<-fread(p);g<-as.character(z$gene_symbol);bad<-is.na(g)|tolower(trimws(g))%in%c('','na','nan','none','null');g[bad]<-z$gene_id[bad];cols<-setdiff(names(z),c('gene_id','gene_symbol'));raw<-as.matrix(z[,..cols]);rownames(raw)<-g;kind<-'upstream prefilter model-ready quantities'
 }
 paths<-c(paths,p);stopifnot(all(is.finite(raw)),all(raw>=0),!anyNA(rownames(raw)),all(nzchar(rownames(raw))))
 raw<-rowsum(raw,rownames(raw),reorder=FALSE)
 ec<-character(nrow(s));method<-character(nrow(s))
 for(i in seq_len(nrow(s))){
  hits<-intersect(c(s$sample_id[i],s$accession[i],s$sample_key[i]),colnames(raw))
  if(length(hits)==1){ec[i]<-hits;method[i]<-'exact sample identifier'}
  else if(ds=='GSE228975'){ec[i]<-runmeta$run_accession[match(s$sample_id[i],runmeta$sample_id)];method[i]<-'SRA title and replicate number verified'}
  else if(ds=='CRA003645'){ec[i]<-oldalign$expression_column[match(s$sample_key[i],oldalign$sample_key)];method[i]<-'preserved curated alignment (historical dataset_position)'}
  else stop('Ambiguous or missing alignment: ',s$sample_key[i])
 }
 stopifnot(!anyNA(ec),all(ec%in%colnames(raw)),!anyDuplicated(ec));raw<-raw[,ec,drop=FALSE]
 # Check the recovered quantity scale against the existing filtered matrix, before frozen dictionary mapping.
 filtpath<-if(ds=='PRJNA863443')file.path(root,'data_input/aggregated_pseudobulk',ds,'input_genelevel/sample_gene_counts_after_expression_filter.tsv')else unique(base_meta[sample_key%in%s$sample_key,expression_counts_path])
 stopifnot(length(filtpath)==1);paths<-c(paths,filtpath);cur<-read_matrix(filtpath)[,ec,drop=FALSE]
 stopifnot(all(rownames(cur)%in%rownames(raw)));err<-max(abs(raw[rownames(cur),,drop=FALSE]-cur)/pmax(1,cur));stopifnot(err<1e-6)
 den<-colSums(raw);stopifnot(all(den>0));g<-map_ids(rownames(raw));mapping[[cohort]]<-data.table(dataset_cohort=cohort,source_gene=rownames(raw),model_gene=g)
 counts<-rowsum(raw,g,reorder=FALSE);lx<-log1p(sweep(counts,2,den,'/')*1e6);colnames(lx)<-s$sample_key
 common<-intersect(all_genes,rownames(lx));x[common,s$sample_key]<-lx[common,,drop=FALSE]
 alignments[[cohort]]<-data.table(sample_key=s$sample_key,dataset_cohort=cohort,expression_column=ec,alignment_method=method)
 norms[[cohort]]<-data.table(sample_key=s$sample_key,dataset_cohort=cohort,prefilter_total=den,filtered_total=colSums(cur),n_measured_genes=nrow(counts))
 prov[[cohort]]<-data.table(dataset_cohort=cohort,n_samples=nrow(s),source_kind=kind,source_paths=paste(p,collapse=';'),n_source_genes=nrow(raw),n_mapped_genes=nrow(counts),filtered_agreement_max_relative_error=err)
 usedpaths<-unique(paths[seq.int(path_start+1,length(paths))]);saveRDS(list(sample_keys=s$sample_key,dictionary_md5=unname(tools::md5sum(dictpath)),paths=usedpaths,source_info=file.info(usedpaths)[,c('size','mtime')],x=x[,s$sample_key,drop=FALSE],alignment=alignments[[cohort]],norm=norms[[cohort]],provenance=prov[[cohort]],mapping=mapping[[cohort]]),cachepath)
}
write_tsv(m,'tables/external_projection_samples.tsv');write_tsv(m0[dataset_folder=='GSE154890'],'provenance/excluded_GSE154890_samples.tsv')
write_tsv(rbindlist(alignments),'provenance/sample_expression_alignment.tsv');write_tsv(rbindlist(norms),'provenance/normalization_audit.tsv');write_tsv(rbindlist(prov),'provenance/input_provenance.tsv');write_tsv(rbindlist(mapping),'provenance/gene_identifier_mapping.tsv')
saveRDS(list(expression=x,metadata=m,normalization='ln(1 + 1e6 * gene quantity / full prefilter sample total)'),file.path(out,'inputs/model_gene_expression.rds'),compress=FALSE)
preds<-list();pcs<-list();contrib<-list();audits<-list()
for(nm in names(models)){
 model<-models[[nm]];p<-predict_model(model,x);p<-merge(m,p,by='sample_key',sort=FALSE)
 # Development overlap is model-specific: the HVG model did not train on single-age studies.
 tk<-model$training_keys;ts<-unique(training[sample_key%in%tk,dataset_folder])
 p[,`:=`(locked_training_sample=sample_key%in%tk,locked_training_study=dataset_folder%in%ts)]
 p[,current_development_relationship:=fifelse(locked_training_sample,'sample used in this locked model',fifelse(locked_training_study,'new sample from a study used by this locked model','external to this locked model'))]
 p[,`:=`(model=nm,predicted_age_week=pred,delta_age_week=pred-age_week,gene_coverage_fraction=1-missing_gene_fraction,n_model_genes=length(model$representation$genes),outside_training_age_range=age_week<model$training_age_range[1]|age_week>model$training_age_range[2])];p[,pred:=NULL]
 # Independent exported-coefficient check, including training-mean imputation.
 xx<-x[match(model$coefficients$gene,rownames(x)),,drop=FALSE];ii<-which(!is.finite(xx),arr.ind=TRUE);xx[ii]<-model$coefficients$center[ii[,1]]
 yp<-drop(crossprod(model$coefficients$coefficient,xx))+model$intercept
 err<-max(abs(yp-p$predicted_age_week[match(colnames(x),p$sample_key)]));stopifnot(err<1e-8)
 preds[[nm]]<-p;audits[[nm]]<-data.table(model=nm,n_samples=nrow(p),coefficient_replay_max_error=err)
 if(!is.null(model$pc_coefficients)){
  z<-project_representation(model$representation,x)$x[,seq_along(model$pc_coefficients),drop=FALSE]
  pcs[[nm]]<-data.table(sample_key=rownames(z),model=nm,z)
  cc<-sweep(z,2,model$pc_coefficients,'*');contrib[[nm]]<-data.table(sample_key=rownames(cc),model=nm,model_intercept=model$model_intercept,cc)
 }
}
p<-rbindlist(preds);stopifnot(nrow(p)==281*length(models),all(is.finite(p$predicted_age_week)))
write_tsv(p,'external_projection_predictions.tsv');write_tsv(p,'tables/external_projection_predictions.tsv');write_tsv(rbindlist(audits),'provenance/prediction_replay_checks.tsv')
write_tsv(rbindlist(pcs),'tables/pca_projection_scores.tsv');write_tsv(rbindlist(contrib),'tables/pca_projection_pc_contributions.tsv')
membership<-fread(mp);membership<-membership[sample_key%in%m$sample_key]
setnames(membership,c('dataset','tage_contrast_id','tage_contrast_name','tage_contrast_type','tage_contrast_role'),c('dataset_id','contrast_id','contrast_name','contrast_type','contrast_role'))
membership<-merge(membership,m[,.(sample_key,modality,dataset_cohort)],by='sample_key',sort=FALSE)
membership[,analysis_contrast_id:=paste(contrast_id,modality,sep='__')]
write_tsv(membership,'tables/external_projection_contrast_membership.tsv');write_tsv(fread(cp)[dataset_id!='GSE154890'],'inputs/prespecified_contrasts.tsv')
write_tsv(m[,.(n_samples=.N,n_training_samples=sum(locked_training_sample)),by=.(dataset_folder,dataset_cohort,modality)],'tables/cohort_counts.tsv')
paths<-unique(c(paths,unlist(lapply(c('HVG-ElasticNet-exclude-single-age','VersionB-PCA-glmnet'),function(d)file.path(release,d,'model.rds')))))
sha<-vapply(paths,function(p)strsplit(system2('sha256sum',shQuote(p),stdout=TRUE),' +')[[1]][1],character(1))
write_tsv(data.table(path=paths,sha256=sha),'provenance/source_checksums.tsv')
capture.output(sessionInfo(),file=file.path(out,'provenance/R_sessionInfo.txt'))
message('PASS: 281 samples, 20 modality-specific datasets, 562 frozen predictions (two locked models).')
