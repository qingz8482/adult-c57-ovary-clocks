#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
arg<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
arg<-gsub('~+~',' ',arg,fixed=TRUE)
root<-normalizePath(file.path(dirname(arg),'..'))
gdir<-file.path(root,'12_PCA deployment model GSEA')
ddir<-file.path(root,'10_pca_diagnostics')
getg<-function(n)fread(file.path(gdir,paste0(n,'.tsv')))
m<-readRDS(file.path(root,'../../locked_models/count_tuning_20261011/VersionB-PCA-glmnet/model.rds'))
r<-getg('final_pca_pc_gsea_rank_vectors')
ir<-getg('final_pca_integrated_gsea_rank_vector')
a<-getg('final_pca_pc_gsea_all_results')
ia<-getg('final_pca_integrated_gsea_all_results')
sets<-getg('gene_set_universe_audit')[primary_eligible==TRUE]
active<-colnames(m$representation$rotation)[m$pc_coefficients!=0]
checks<-list()
check<-function(name,ok,error=NA_real_) {
 checks[[length(checks)+1L]]<<-data.table(check=name,passed=isTRUE(ok),error=error)
 if(!isTRUE(ok))stop('Validation failed: ',name)
}
check('active_PC_identity',setequal(unique(r$pc),active))
check('complete_rank_vectors',all(r[,uniqueN(gene_id),by=pc]$V1==3000)&&nrow(r)==12000&&nrow(ir)==3000)
e<-max(abs(ir$rank_score-m$coefficients$standardized_coefficient[match(ir$gene_id,m$coefficients$gene)]))
check('integrated_rank_matches_frozen_gene_coefficients',e<1e-12,e)
sumr<-r[,.(expected=sum(rank_score)),by=gene_id]
e<-max(abs(ir$rank_score-sumr$expected[match(ir$gene_id,sumr$gene_id)]))
check('integrated_rank_is_sum_of_active_PC_contributions',e<1e-12,e)
check('every_eligible_set_tested_once_per_rank',nrow(a)==length(active)*nrow(sets)&&nrow(ia)==nrow(sets)&&
 !anyDuplicated(a[,.(pc,gene_set_id)])&&!anyDuplicated(ia$gene_set_id))
both<-rbindlist(list(a,ia),fill=TRUE)
check('finite_enrichment_statistics',all(is.finite(both$ES))&&all(is.finite(both$NES))&&all(is.finite(both$empirical_p)))
check('valid_empirical_P_and_BH_ranges',all(both$empirical_p>0&both$empirical_p<=1&both$fdr_within_pc>=0&both$fdr_within_pc<=1))
qcheck<-both[,.(error=max(abs(p.adjust(empirical_p,'BH')-fdr_within_pc))),by=pc]
check('BH_within_each_rank_recomputed',max(qcheck$error)<1e-12,max(qcheck$error))
qcheck<-both[,.(error=max(abs(p.adjust(empirical_p,'BH')-fdr_within_pc_collection))),by=.(pc,collection_label)]
check('BH_within_each_rank_collection_recomputed',max(qcheck$error)<1e-12,max(qcheck$error))
check('size_filters',all(both$gene_set_size_in_universe>=10&both$gene_set_size_in_universe<=300))
check('direction_uses_model_coefficient_orientation',all((both$NES>0)==(both$direction=='older_aligned')))
check('permutation_count',all(both$n_permutations==500))
check('nominal_P_resolution',max(abs(both$empirical_p*501-round(both$empirical_p*501)))<1e-9)
for(prefix in c('final_pca_pc','final_pca_integrated')) {
 results<-getg(paste0(prefix,'_gsea_all_results'))
 sig<-getg(paste0(prefix,'_gsea_significant_terms'))
 check(paste0(prefix,'_significance_filter'),setequal(paste(sig$pc,sig$gene_set_id),paste(results[empirical_p<.05&fdr_within_pc<.1]$pc,results[empirical_p<.05&fdr_within_pc<.1]$gene_set_id)))
}
# Independent dense running-sum replay (rather than the sparse hit-position
# implementation): 100 systematically distributed pathways per rank.
replay<-rbindlist(lapply(unique(both$pc),function(pcname) {
 rr<-if(pcname=='integrated_PCA_clock')copy(ir) else copy(r[pc==pcname])
 setorder(rr,-rank_score)
 rows<-both[pc==pcname][unique(round(seq(1,.N,length.out=min(.N,100))))]
 rbindlist(lapply(seq_len(nrow(rows)),function(i) {
  row<-rows[i];genes<-strsplit(sets$genes[match(row$gene_set_id,sets$gene_set_id)],';',fixed=TRUE)[[1]]
  hit<-rr$gene_id %in% genes
  steps<-rep(-1/sum(!hit),nrow(rr));steps[hit]<-abs(rr$rank_score[hit])/sum(abs(rr$rank_score[hit]))
  curve<-c(0,cumsum(steps));j<-which.max(abs(curve));es<-curve[j];peak<-j-1L
  le<-if(es>=0)rr$gene_id[hit&seq_len(nrow(rr))<=peak] else rr$gene_id[hit&seq_len(nrow(rr))>peak]
  data.table(pc=pcname,gene_set_name=row$gene_set_name,ES_error=abs(es-row$ES),
   leading_edge_matches=setequal(le,strsplit(row$leading_edge_genes,';',fixed=TRUE)[[1]]))
 }))
}))
fwrite(replay,file.path(root,'provenance/gsea_running_sum_replay.tsv'),sep='\t')
check('500_independent_running_sum_ES_replays',nrow(replay)==500&&max(replay$ES_error)<1e-10,max(replay$ES_error))
check('500_leading_edge_replays',all(replay$leading_edge_matches))
pcsum<-fread(file.path(ddir,'final_pca_clock_all_input_pc_summary.tsv'))
check('all_10_PC_inputs_documented',nrow(pcsum)==10&&sum(pcsum$active_in_clock)==4)
check('cumulative_variance_includes_intervening_zero_coefficient_PCs',max(abs(cumsum(pcsum$variance_explained)-pcsum$cumulative_variance_through_pc))<1e-12)
qc<-rbindlist(checks)
fwrite(qc,file.path(root,'provenance/gsea_validation.tsv'),sep='\t')
cat(nrow(qc),'GSEA and scope checks passed.\n')
