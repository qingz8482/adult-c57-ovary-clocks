#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table));setDTthreads(1)
arg<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),fixed=TRUE)
out<-normalizePath(file.path(dirname(arg),'..'));diagnostic<-dirname(out)
dd<-file.path(out,'derived_data');gd<-file.path(diagnostic,'12_PCA deployment model GSEA')
write_tsv<-function(x,n)fwrite(x,file.path(dd,paste0(n,'.tsv')),sep='\t')
gsea<-fread(file.path(gd,'final_pca_integrated_gsea_all_results.tsv'))
rank_dt<-fread(file.path(gd,'final_pca_integrated_gsea_rank_vector.tsv'))
rank_vec<-setNames(rank_dt$rank_score,rank_dt$gene_id);rank_vec<-sort(rank_vec,decreasing=TRUE)
stopifnot(length(rank_vec)==3000,!anyDuplicated(names(rank_vec)),all(is.finite(rank_vec)))
sets<-fread(file.path(gd,'gene_set_universe_audit.tsv'))
split_genes<-function(x)if(is.na(x)||!nzchar(x))character() else unique(strsplit(x,';',fixed=TRUE)[[1]])
leading<-setNames(lapply(gsea$leading_edge_genes,split_genes),gsea$gene_set_name)
jaccard<-function(a,b){u<-union(a,b);if(length(u))length(intersect(a,b))/length(u) else 0}
audit<-list();slots<-list();ai<-0L;si<-0L
select_terms<-function(panel,category,pattern='.',direction=NULL,n=2,exclude=character(),compare_prior=FALSE) {
 candidates<-gsea[grepl(pattern,gene_set_name)]
 direction_filter<-direction
 if(!is.null(direction_filter))candidates<-candidates[direction==direction_filter]
 candidates<-candidates[order(fdr_within_pc,-abs(NES),gene_set_name)]
 selected<-character()
 for(i in seq_len(nrow(candidates))) {
  r<-candidates[i];term<-r$gene_set_name
  comp<-unique(c(if(compare_prior)exclude else character(),selected))
  overlap<-if(length(comp))max(vapply(comp,function(t)jaccard(leading[[term]],leading[[t]]),numeric(1))) else 0
  decision<-if(r$fdr_within_pc>=.05)'not_significant_at_FDR_0.05' else if(term%in%exclude)'already_displayed' else if(length(selected)>=n)'quota_filled' else if(overlap>=.5)'leading_edge_Jaccard_at_least_0.5' else 'selected'
  if(decision=='selected')selected<-c(selected,term)
  ai<<-ai+1L;audit[[ai]]<<-data.table(panel=panel,category=category,gene_set_name=term,NES=r$NES,fdr_within_pc=r$fdr_within_pc,
   max_jaccard_at_evaluation=overlap,selected=decision=='selected',decision=decision,category_rank=i)
 }
 for(k in seq_len(n)) {
  si<<-si+1L;slots[[si]]<<-data.table(panel=panel,category=category,category_slot=k,
   gene_set_name=if(k<=length(selected))selected[k] else NA_character_,
   status=if(k<=length(selected))'selected' else if(!nrow(candidates[fdr_within_pc<.05]))'no_terms_pass_FDR_0.05' else 'no_further_nonredundant_terms')
 }
 selected
}
immune_pattern<-'IMMUNE|INFLAM|INTERFERON|CYTOKINE|T_CELL|B_CELL|NK_CELL|LEUKOCYTE|MACROPHAGE|COMPLEMENT'
f_immune<-select_terms('f','older immune',immune_pattern,'older_aligned',2)
f_mito<-select_terms('f','younger energy','MITOCHON|OXIDATIVE_PHOSPHORYLATION|RESPIRATORY|ENERGY|TCA|ATP','younger_aligned',1,exclude=f_immune)
f_lipid<-select_terms('f','younger lipid','STEROID|LIPID|FATTY_ACID|CHOLESTEROL','younger_aligned',1,exclude=c(f_immune,f_mito))
f<-c(f_immune,f_mito,f_lipid)
f2_older<-select_terms('f2','additional older','.', 'older_aligned',2,exclude=f,compare_prior=TRUE)
f2_younger<-select_terms('f2','additional younger','.', 'younger_aligned',2,exclude=c(f,f2_older),compare_prior=TRUE)
f2<-c(f2_older,f2_younger)
f3_steroid<-select_terms('f3','steroid / cholesterol','STEROID|CHOLESTEROL',n=2,exclude=c(f,f2))
f3_cycle<-select_terms('f3','cell cycle','CELL_CYCLE|MITOTIC|DNA_REPLICATION|G2M|E2F',n=2,exclude=c(f,f2,f3_steroid))
# f4: literal top two integrated results among the rows displayed in heatmap 13.
# No exclusion of earlier panels or redundancy filter: those would change "top two".
heatmap_file<-file.path(gd,'heatmap_term_selection_audit.tsv')
heatmap_terms<-fread(heatmap_file)
f4_candidates<-merge(heatmap_terms[,.(gene_set_id,biology_group)],gsea,by='gene_set_id',sort=FALSE)
f4_candidates[,sort_key:=fdr_within_pc+empirical_p/1000]
f4_candidates<-f4_candidates[order(sort_key,-abs(NES),gene_set_name)]
f4_candidates[,category_rank:=seq_len(.N),by=biology_group]
f4_candidates[,selected:=category_rank<=2]
stopifnot(uniqueN(f4_candidates$biology_group)==9L,
 all(f4_candidates[selected==TRUE,fdr_within_pc]<.05),
 all(f4_candidates[selected==TRUE,empirical_p]<.05))
for(group in unique(heatmap_terms$biology_group)) {
 chosen<-f4_candidates[biology_group==group&selected==TRUE][order(category_rank)]
 stopifnot(nrow(chosen)==2L)
 for(k in 1:2){si<-si+1L;slots[[si]]<-data.table(panel='f4',category=group,category_slot=k,
   gene_set_name=chosen$gene_set_name[k],status='selected')}
}
ai<-ai+1L;audit[[ai]]<-f4_candidates[,.(panel='f4',category=biology_group,gene_set_name,NES,fdr_within_pc,
 empirical_p,selected,decision=ifelse(selected,'selected','below_top_two'),category_rank)]
# f5: two user-specified terms, plotted as requested (no automated selection or redundancy filter).
f5<-c('HALLMARK_TNFA_SIGNALING_VIA_NFKB','GOBP_CYTOKINE_PRODUCTION')
stopifnot(all(f5%in%gsea$gene_set_name),all(gsea[match(f5,gene_set_name),fdr_within_pc]<.05))
for(k in seq_along(f5)){si<-si+1L;slots[[si]]<-data.table(panel='f5',category='requested older immune',category_slot=k,
  gene_set_name=f5[k],status='selected')}
ai<-ai+1L;audit[[ai]]<-gsea[match(f5,gene_set_name),.(panel='f5',category='requested older immune',gene_set_name,NES,fdr_within_pc,
 empirical_p,selected=TRUE,decision='user_specified',category_rank=seq_along(f5))]
slot_dt<-rbindlist(slots);slot_dt[,display_order:=seq_len(.N),by=panel]
write_tsv(slot_dt,'gsea_panel_slots');audit_dt<-rbindlist(audit,fill=TRUE)
selected_all<-unique(na.omit(slot_dt$gene_set_name))
stopifnot(all(gsea[gene_set_name%in%selected_all,fdr_within_pc]<.05))
write_tsv(data.table(input_file=normalizePath(heatmap_file),
 sha256=strsplit(system2('sha256sum',shQuote(heatmap_file),stdout=TRUE),' ')[[1]][1]),'panel_f4_input_files_used')
cache_files<-list.files(file.path(diagnostic,'provenance/msigdb'),pattern='Mm.*rds$',full.names=TRUE)
membership<-rbindlist(lapply(cache_files,function(p){z<-as.data.table(readRDS(p));unique(z[gs_name%in%selected_all,.(gene_set_name=gs_name,gene_id=db_gene_symbol)])}))
membership<-unique(membership);membership[,in_universe:=gene_id%in%names(rank_vec)]
curve<-function(term) {
 genes<-membership[gene_set_name==term&in_universe==TRUE,gene_id]
 hit<-names(rank_vec)%in%genes;pos<-which(hit);nh<-sum(hit);nn<-length(rank_vec)
 weights<-abs(rank_vec[pos]);weights<-weights/sum(weights)
 increments<-rep(-1/(nn-nh),nn);increments[pos]<-weights
 running<-cumsum(increments);cum<-cumsum(weights);idx<-seq_along(pos)
 after<-cum-(pos-idx)/(nn-nh);before<-c(0,head(cum,-1))-(pos-idx)/(nn-nh)
 positive<-abs(max(after))>=abs(min(before))
 peak<-if(positive)pos[which.max(after)] else pos[which.min(before)]
 es<-if(positive)max(after) else min(before)
 le<-names(rank_vec)[if(positive)pos[pos<=peak] else pos[pos>=peak]]
 saved<-gsea[gene_set_name==term]
 stopifnot(nrow(saved)==1,peak==saved$peak_rank,abs(es-saved$ES)<1e-12,setequal(le,leading[[term]]),nh==saved$gene_set_size_in_universe)
 data.table(gene_set_name=term,rank=seq_len(nn),gene_id=names(rank_vec),rank_score=as.numeric(rank_vec),hit=hit,
  running_es=running,peak_rank_curve_extremum=which.max(abs(running)),peak_rank_pipeline=peak,es_pipeline=es)
}
curves<-rbindlist(lapply(selected_all,curve))
for(pn in c('f','f2','f3','f4','f5')) {
 ss<-slot_dt[panel==pn&status=='selected']
 stats<-merge(ss,gsea,by='gene_set_name',sort=FALSE);setorder(stats,display_order)
 write_tsv(stats,paste0('panel_',pn,'_selected_terms'))
 write_tsv(audit_dt[panel==pn],paste0('panel_',pn,'_selection_audit'))
 write_tsv(membership[gene_set_name%in%ss$gene_set_name],paste0('panel_',pn,'_selected_gene_set_membership'))
 cc<-curves[gene_set_name%in%ss$gene_set_name]
 write_tsv(cc,paste0('panel_',pn,'_enrichment_curves'))
 qc<-merge(unique(cc[,.(gene_set_name,peak_rank_curve_extremum,peak_rank_pipeline,es_pipeline)]),
  stats[,.(gene_set_name,peak_rank_saved=peak_rank,es_saved=ES)],by='gene_set_name')
 qc[,`:=`(peak_exact_match=peak_rank_pipeline==peak_rank_saved,es_absolute_discrepancy=abs(es_pipeline-es_saved))]
 write_tsv(qc,paste0('panel_',pn,'_peak_rank_validation'))
}
write_tsv(gsea[grepl('STEROID|CHOLESTEROL',gene_set_name)],'panel_f3_all_tested_steroid_terms')
message('GSEA curves complete; selected slots: ',sum(slot_dt$status=='selected'),'/',nrow(slot_dt),'.')
print(slot_dt[,.(panel,category,gene_set_name,status)])
