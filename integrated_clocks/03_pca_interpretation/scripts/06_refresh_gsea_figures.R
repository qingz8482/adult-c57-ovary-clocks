#!/usr/bin/env Rscript
# Reuse the canonical plotting expressions with saved results; no permutations.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(grid)})
setDTthreads(1)
arg<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
arg<-gsub('~+~',' ',arg,fixed=TRUE)
root<-normalizePath(file.path(dirname(arg),'..'))
out_dir<-file.path(root,'12_PCA deployment model GSEA')
getg<-function(n)fread(file.path(out_dir,paste0(n,'.tsv')))
pc_summary<-fread(file.path(root,'10_pca_diagnostics/final_pca_clock_selected_pc_summary.tsv'))
expected_pcs<-pc_summary$pc
all_results<-getg('final_pca_pc_gsea_all_results')
integrated_results<-getg('final_pca_integrated_gsea_all_results')
significant_terms<-getg('final_pca_pc_gsea_significant_terms')
rank_vectors<-getg('final_pca_pc_gsea_rank_vectors')
collection_summary<-getg('final_pca_pc_gsea_collection_summary')
gene_sets_primary<-getg('gene_set_universe_audit')[primary_eligible==TRUE]
gene_sets_primary[,genes:=strsplit(genes,';',fixed=TRUE)]
score_type<-1;top_n_plot<-8
exprs<-parse(file.path(root,'scripts/02_deployment_gsea.R'))
heatmaps_only <- '--heatmaps-only' %in% commandArgs(trailingOnly = TRUE)
started<-FALSE
for(e in exprs) {
 lhs<-if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]]))as.character(e[[2]]) else ''
 if (heatmaps_only) {
   if (lhs %in% c('theme_gsea','pretty_term','save_plot','empty_plot','biology_rules',
                  'term_biology_group','biology_group_levels')) eval(e)
   if (lhs == 'heat_top_n') {started <- TRUE; plot_files <- character()}
   if (lhs == 'running_curve_dt') break
   if (started) eval(e)
   next
 }
 if(lhs=='theme_gsea')started<-TRUE
 if(lhs=='format_top_terms_md')break
 if(started)eval(e)
}
cat(if (heatmaps_only) 'Refreshed heatmaps 05 and 13 from saved numerical results.\n' else
  'Refreshed all 13 GSEA figures from saved numerical results.\n')
