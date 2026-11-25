source(file.path(dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),fixed=TRUE))),'common.R'))
suppressPackageStartupMessages(library(ggrepel))
inputfiles<-file.path(out,'tables',c('feature_overlap_gene_universe.tsv','pathway_library_membership.tsv','ovary_elasticnet_features.tsv','ovary_pca_backprojected_features.tsv','tage_scaleddiff_features.tsv','tage_yugenediff_features.tsv'))
protected<-c(inputfiles,list.files(file.path(out,'figures'),pattern='^(shared_ageing_pathways|03_|10_|feature_coefficient).*',full.names=TRUE))
integrity<-data.table(path=normalizePath(protected),md5_before=unname(tools::md5sum(protected)))
universe<-sort(read_tsv(inputfiles[1])$gene);N<-length(universe)
member<-read_tsv(inputfiles[2]);member[,pathway_id:=paste(database,pathway,sep=' | ')]
sets<-split(member$gene,member$pathway_id)
sets<-lapply(sets,unique);indices<-lapply(sets,match,table=universe)
meta<-unique(member[,.(pathway_id,database,pathway)]);meta<-meta[match(names(sets),pathway_id)]
meta[,pathway_size:=lengths(sets)]
weights<-matrix(0,N,4,dimnames=list(universe,models))
for(i in 1:4){
  z<-read_tsv(inputfiles[i+2]);col<-if(i<=2)'standardized_coefficient' else 'coefficient'
  ii<-match(universe,z$gene);weights[!is.na(ii),i]<-z[[col]][ii[!is.na(ii)]]
}
stopifnot(all(is.finite(weights)),all(meta$pathway_size>=5 & meta$pathway_size<=500))
write_tsv(data.table(gene=universe,weights,check.names=FALSE),'tables/model_resolved_gene_weights.tsv')
stats_rows<-list()
for(i in 1:4)for(analysis in c('absolute','signed')){
  w<-weights[,i];v<-if(analysis=='absolute')abs(w) else w
  ranks<-rank(v,ties.method='average');ties<-as.numeric(table(v));tie_term<-sum(ties^3-ties)/(N*(N-1))
  m<-meta$pathway_size;n<-N-m
  u<-vapply(indices,function(ii)sum(ranks[ii]),numeric(1))-m*(m+1)/2
  effect<-2*u/(m*n)-1
  sdU<-sqrt(m*n/12*(N+1-tie_term))
  z<-if(analysis=='absolute')(u-m*n/2-.5)/sdU else (abs(u-m*n/2)-.5)/sdU
  p<-pnorm(z,lower.tail=FALSE)*if(analysis=='absolute')1 else 2
  p<-pmin(1,p);q<-p.adjust(p,'BH')
  d<-copy(meta);d[,`:=`(model=models[i],analysis=analysis,rank_biserial=effect,U=u,p_value=p,FDR=q,n_background=N,active_genes_in_model=sum(w!=0),active_genes_in_pathway=vapply(indices,function(ii)sum(w[ii]!=0),integer(1)),positive_genes_in_pathway=vapply(indices,function(ii)sum(w[ii]>0),integer(1)),negative_genes_in_pathway=vapply(indices,function(ii)sum(w[ii]<0),integer(1)))]
  stats_rows[[length(stats_rows)+1L]]<-d
}
st<-rbindlist(stats_rows)
# Exchange gene labels jointly across models. This handles the sparse ElasticNet ties
# without relying on a normal approximation for small pathway/active-gene overlaps.
st[,`:=`(asymptotic_p_value=p_value,asymptotic_FDR=FDR)]
rankmatrix<-do.call(cbind,lapply(1:4,function(i)cbind(rank(abs(weights[,i]))-(N+1)/2,rank(weights[,i])-(N+1)/2)))
cache_dir<-file.path(out,'cache');dir.create(cache_dir,showWarnings=FALSE)
cache_path<-file.path(cache_dir,'pathway_rank_permutations.rds')
signature<-unname(tools::md5sum(c(inputfiles,file.path(dirname(script),'pathway_permutation.cpp'))))
B<-100000L
cached<-if(file.exists(cache_path))readRDS(cache_path) else NULL
if(is.null(cached)||!identical(cached$signature,signature)||cached$B!=B){
  source(file.path(dirname(script),'permutation_build.R'))
  Rcpp::sourceCpp(file.path(dirname(script),'pathway_permutation.cpp'),cacheDir=file.path(cache_dir,'Rcpp'))
  counts<-pathway_permutation(rankmatrix,indices,B,20261012L)
  cached<-list(signature=signature,B=B,seed=20261012L,counts=counts,pathway_ids=meta$pathway_id)
  saveRDS(cached,cache_path)
}
stopifnot(identical(cached$pathway_ids,meta$pathway_id))
for(i in 1:4)for(a in 1:2){
  analysis_name<-c('absolute','signed')[a];k<-cached$counts[,2*(i-1)+a];pp<-(k+1)/(B+1)
  st[model==models[i] & analysis==analysis_name,`:=`(p_value=pp,FDR=p.adjust(pp,'BH'),permutation_exceedances=k,permutations=B)]
}
write_tsv(st,'tables/model_resolved_pathway_statistics.tsv')
absolute<-st[analysis=='absolute'];signed<-st[analysis=='signed']
effects<-dcast(absolute,pathway_id~model,value.var='rank_biserial')[match(meta$pathway_id,pathway_id)]
fdr<-dcast(absolute,pathway_id~model,value.var='FDR')[match(meta$pathway_id,pathway_id)]
em<-as.matrix(effects[,..models]);qm<-as.matrix(fdr[,..models])
sig<-qm<.05 & em>0
programs<-copy(meta)
programs[,`:=`(ovary_mean_effect=rowMeans(em[,1:2]),tage_mean_effect=rowMeans(em[,3:4]),ovary_model_support=rowSums(sig[,1:2]),tage_model_support=rowSums(sig[,3:4]),shared_strength=pmin(apply(em[,1:2],1,max),apply(em[,3:4],1,max)),all_four_supported=rowSums(sig)==4)]
programs[,difference:=ovary_mean_effect-tage_mean_effect]
programs[,support_class:=fcase(ovary_model_support>0 & tage_model_support>0,'Shared support',ovary_model_support>0,'Ovary support only',tage_model_support>0,'tAge support only',default='No FDR support')]
sm<-as.matrix(dcast(signed,pathway_id~model,value.var='rank_biserial')[match(meta$pathway_id,pathway_id),..models])
sq<-as.matrix(dcast(signed,pathway_id~model,value.var='FDR')[match(meta$pathway_id,pathway_id),..models])
programs[,`:=`(all_four_signed_same_direction=apply(sign(sm),1,function(z)all(z>0)||all(z<0)),all_four_signed_supported=rowSums(sq<.05)==4)]
write_tsv(programs,'tables/model_resolved_pathway_comparison.tsv')
groups<-c('Shared pathway support','Ovary-preferential candidates','tAge-preferential candidates')
chosen<-list();audit<-list();selected_ids<-character()
for(gr in groups){
  cand<-if(gr==groups[1])copy(programs[support_class=='Shared support']) else if(gr==groups[2])copy(programs[support_class=='Ovary support only' & difference>0]) else copy(programs[support_class=='tAge support only' & difference<0])
  cand[,priority:=if(gr==groups[1])shared_strength else if(gr==groups[2])difference else -difference]
  setorder(cand,-priority,pathway,database)
  nkeep<-0L
  for(i in seq_len(nrow(cand))){
    id<-cand$pathway_id[i];sim<-if(length(selected_ids))vapply(selected_ids,function(k)length(intersect(sets[[id]],sets[[k]]))/length(union(sets[[id]],sets[[k]])),numeric(1)) else numeric()
    maxj<-if(length(sim))max(sim) else 0;near<-if(length(sim))selected_ids[which.max(sim)] else NA_character_
    keep<-cand$pathway_size[i]>=15 && nkeep<5 && maxj<.35
    audit[[length(audit)+1L]]<-data.table(display_group=gr,pathway_id=id,candidate_rank=i,priority=cand$priority[i],selected=keep,max_Jaccard_to_already_selected=maxj,nearest_already_selected=near,reason=if(keep)'selected' else if(cand$pathway_size[i]<15)'fewer than 15 mapped genes; tested but excluded from representative display' else if(nkeep>=5)'five representatives already selected' else 'overlapping representative already selected')
    if(keep){row<-copy(cand[i]);row[,display_group:=gr];chosen[[length(chosen)+1L]]<-row;selected_ids<-c(selected_ids,id);nkeep<-nkeep+1L}
  }
}
sel<-rbindlist(chosen);sel[,row_index:=seq_len(.N)]
stopifnot(nrow(sel)>0,!anyDuplicated(sel$pathway_id))
write_tsv(sel,'tables/model_resolved_selected_pathways.tsv');write_tsv(rbindlist(audit),'tables/model_resolved_representative_selection.tsv')
source(file.path(dirname(script),'model_resolved_figures.R'))
integrity[,md5_after:=unname(tools::md5sum(path))];integrity[,unchanged:=md5_before==md5_after];stopifnot(all(integrity$unchanged));write_tsv(integrity,'provenance/model_resolved_input_integrity.tsv')
counts<-programs[,.N,by=support_class];write_tsv(counts,'tables/model_resolved_support_counts.tsv')
write_tsv(data.table(model=models,active_genes_in_background=colSums(weights!=0),background_genes=N),'tables/model_resolved_gene_coverage.tsv')
message('PASS: revised pathway comparisons generated.');print(counts);cat('All four absolute supported:',sum(programs$all_four_supported),'\nAll four signed supported and same direction:',sum(programs$all_four_signed_supported & programs$all_four_signed_same_direction),'\n')
