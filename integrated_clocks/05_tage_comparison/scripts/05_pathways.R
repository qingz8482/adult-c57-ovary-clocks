source(file.path(dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),fixed=TRUE))),'common.R'))
ref<-file.path(archive,'external_projection_before_reference_control_patch_20260727_164947/tAge/feature_overlap_analysis')
inputdir<-file.path(out,'inputs/gene_sets');dir.create(inputdir,recursive=TRUE,showWarnings=FALSE)
dbnames<-c('GO Biological Process 2025','Reactome Pathways 2024','KEGG Mouse 2019')
gmtfiles<-c('GO_Biological_Process_2025.gmt','Reactome_Pathways_2024.gmt','KEGG_2019_Mouse.gmt')
gmtpaths<-file.path(inputdir,gmtfiles)
for(i in seq_along(gmtfiles)){
  original<-file.path(ref,'environment/gene_sets',gmtfiles[i])
  if(!file.exists(gmtpaths[i]))stopifnot(file.copy(original,gmtpaths[i]))
  stopifnot(unname(tools::md5sum(original))==unname(tools::md5sum(gmtpaths[i])))
}
protected<-c(file.path(out,'tables',c('ovary_elasticnet_features.tsv','ovary_pca_backprojected_features.tsv','tage_scaleddiff_features.tsv','tage_yugenediff_features.tsv','feature_overlap_gene_universe.tsv')),file.path(ref,'figures/shared_ageing_pathways.pdf'),file.path(ref,'scripts/run_feature_overlap_analysis.py'),file.path(ref,'tables',c('pathway_overlap_enrichment.tsv','pathway_program_classification.tsv')),list.files(file.path(out,'figures'),pattern='^(03_|10_|feature_coefficient).*',full.names=TRUE),gmtpaths)
integrity<-data.table(path=normalizePath(protected),md5_before=unname(tools::md5sum(protected)))
read_sets<-function(universe){
  canonical<-setNames(universe,toupper(universe));stopifnot(!anyDuplicated(names(canonical)))
  setNames(lapply(gmtpaths,function(path){
    sets<-list()
    for(line in readLines(path,warn=FALSE)){
      z<-strsplit(line,'\t',fixed=TRUE)[[1]];if(length(z)<3)next
      genes<-sort(unique(unname(canonical[toupper(z[-c(1,2)])])));genes<-genes[!is.na(genes)]
      if(length(genes)>=5&&length(genes)<=500)sets[[z[1]]]<-genes
    };sets
  }),dbnames)
}
make_weights<-function(en,pc,scaled,yugene){
  setorder(pc,-abs_weight,gene)
  list(queries=list(ovary_elasticnet_active=en[weight!=0,gene],ovary_pca_top100=head(pc$gene,100),ovary_pca_top500=head(pc$gene,500),tage_scaleddiff_active=scaled[weight!=0,gene],tage_yugenediff_active=yugene[weight!=0,gene]),ranks=list(ovary_elasticnet_ranked=en,ovary_pca_backprojected_ranked=pc,tage_scaleddiff_ranked=scaled,tage_yugenediff_ranked=yugene))
}
en<-read_tsv(file.path(out,'tables/ovary_elasticnet_features.tsv'))[,.(gene,weight=standardized_coefficient)]
pc<-read_tsv(file.path(out,'tables/ovary_pca_backprojected_features.tsv'))[,.(gene,weight=standardized_coefficient,abs_weight=abs(standardized_coefficient))]
scaled<-read_tsv(file.path(out,'tables/tage_scaleddiff_features.tsv'))[,.(gene,weight=coefficient)]
yugene<-read_tsv(file.path(out,'tables/tage_yugenediff_features.tsv'))[,.(gene,weight=coefficient)]
universe<-sort(read_tsv(file.path(out,'tables/feature_overlap_gene_universe.tsv'))$gene)
weights<-make_weights(en,pc,scaled,yugene);libraries<-read_sets(universe)
analyse<-function(weights,libraries,universe){
  N<-length(universe);rows<-list()
  for(db in names(libraries)){
    sets<-libraries[[db]];M<-length(sets);K<-lengths(sets);indices<-lapply(sets,match,table=universe)
    for(nm in names(weights$queries)){
      query<-intersect(weights$queries[[nm]],universe);n<-length(query)
      shared<-lapply(sets,intersect,y=query);k<-lengths(shared)
      den<-(n-k)*(K-k);odds<-k*(N-n-K+k)/den;odds[den==0 & k*(N-n-K+k)>0]<-Inf
      p<-phyper(k-1,K,N-K,n,lower.tail=FALSE)
      z<-data.table(gene_set=nm,pathway=names(sets),database=db,method='ORA_Fisher',overlap_gene_count=k,query_size=n,pathway_size=K,background_size=N,NES_or_OR=odds,p_value=p,direction='not_applicable',shared_genes=vapply(shared,function(g)paste(sort(g),collapse=';'),character(1)))
      z<-z[overlap_gene_count>0]
      z[,`:=`(FDR=p.adjust(p_value,'BH'),FDR_full_pathway_family=p.adjust(p_value,'BH',n=M),n_eligible_pathways=M,n_tested_legacy=.N)]
      rows[[length(rows)+1L]]<-z
    }
    for(nm in names(weights$ranks)){
      tab<-weights$ranks[[nm]];v<-setNames(rep(0,N),universe);hit<-match(universe,tab$gene);v[!is.na(hit)]<-tab$weight[hit[!is.na(hit)]]
      ranks<-rank(abs(v),ties.method='average')
      rank_sum<-vapply(indices,function(ii)sum(ranks[ii]),numeric(1))
      zscore<-(rank_sum-K*(K+1)/2-K*(N-K)/2)/sqrt(K*(N-K)*(N+1)/12)
      med<-vapply(indices,function(ii)median(v[ii]),numeric(1))
      z<-data.table(gene_set=nm,pathway=names(sets),database=db,method='rank_based_absolute_score',overlap_gene_count=K,query_size=sum(v!=0),pathway_size=K,background_size=N,NES_or_OR=zscore,p_value=pnorm(zscore,lower.tail=FALSE),direction=ifelse(med>0,'positive_older_aligned',ifelse(med<0,'negative_younger_aligned','zero')),shared_genes=vapply(sets,paste,collapse=';',FUN.VALUE=character(1)))
      z<-z[pathway_size>=5 & background_size-pathway_size>=5]
      z[,`:=`(FDR=p.adjust(p_value,'BH'),FDR_full_pathway_family=p.adjust(p_value,'BH'),n_eligible_pathways=M,n_tested_legacy=.N)]
      rows[[length(rows)+1L]]<-z
    }
  };rbindlist(rows)
}
classify<-function(e,fdr='FDR'){
  p<-e[,{
    ovarian<-grepl('^ovary_',gene_set);q<-get(fdr);ov<-min(c(q[ovarian],1));ta<-min(c(q[!ovarian],1))
    list(ovary_min_FDR=ov,tage_min_FDR=ta,ovary_support=paste(sort(unique(gene_set[ovarian&q<.05])),collapse=';'),tage_support=paste(sort(unique(gene_set[!ovarian&q<.05])),collapse=';'))
  },by=.(database,pathway)]
  p[,program_class:=fcase(ovary_min_FDR<.05&tage_min_FDR<.05,'conserved',ovary_min_FDR<.05,'ovary_specific',tage_min_FDR<.05,'tage_specific',default='not_FDR_supported')]
  p[,support_score:=fcase(program_class=='conserved',-log10(pmax(ovary_min_FDR,tage_min_FDR,1e-300)),program_class=='ovary_specific',-log10(pmax(ovary_min_FDR,1e-300)),program_class=='tage_specific',-log10(pmax(tage_min_FDR,1e-300)),default=0)]
  setorder(p,-support_score,pathway,database);p
}
message('Computing current pathway analysis')
enrich<-analyse(weights,libraries,universe);programs<-classify(enrich);sensitivity<-classify(enrich,'FDR_full_pathway_family')
enrich<-merge(enrich,programs[,.(database,pathway,program_class)],by=c('database','pathway'))
write_tsv(enrich,'tables/pathway_overlap_enrichment.tsv');write_tsv(programs,'tables/pathway_program_classification.tsv');write_tsv(sensitivity,'tables/pathway_program_classification_full_family_sensitivity.tsv')
query_audit<-rbindlist(lapply(names(weights$queries),function(nm){g<-weights$queries[[nm]];data.table(gene_set=nm,gene=g,in_common_universe=g%chin%universe)}))
write_tsv(query_audit,'tables/pathway_query_gene_membership.tsv')
library_members<-rbindlist(lapply(names(libraries),function(db)rbindlist(lapply(names(libraries[[db]]),function(pathway)data.table(database=db,pathway=pathway,gene=libraries[[db]][[pathway]])))))
write_tsv(library_members,'tables/pathway_library_membership.tsv')
write_tsv(data.table(database=dbnames,gmt_file=gmtfiles,n_eligible_pathways=lengths(libraries),background_size=length(universe)),'tables/pathway_library_summary.tsv')

# Replay every historical statistic before accepting the port of the original method.
message('Replaying historical pathway analysis')
olduniv<-sort(read_tsv(file.path(ref,'audit/feature_overlap_gene_universe.tsv'))$gene)
oe<-read_tsv(file.path(ref,'tables/ovary_elasticnet_features.tsv'))[,.(gene,weight=coefficient)]
op<-read_tsv(file.path(ref,'tables/ovary_pca_backprojected_features.tsv'))[,.(gene,weight=gene_score,abs_weight=abs(gene_score))]
os<-read_tsv(file.path(ref,'tables/tage_scaleddiff_features.tsv'))[,.(gene,weight=coefficient)]
oy<-read_tsv(file.path(ref,'tables/tage_yugenediff_features.tsv'))[,.(gene,weight=coefficient)]
replay<-analyse(make_weights(oe,op,os,oy),read_sets(olduniv),olduniv)
original<-read_tsv(file.path(ref,'tables/pathway_overlap_enrichment.tsv'))
joined<-merge(replay,original,by=c('gene_set','database','pathway','method'),suffixes=c('_replay','_original'))
stopifnot(nrow(joined)==nrow(original),nrow(joined)==nrow(replay),max(abs(joined$p_value_replay-joined$p_value_original))<1e-9,max(abs(joined$FDR_replay-joined$FDR_original))<1e-9)
oldprograms<-read_tsv(file.path(ref,'tables/pathway_program_classification.tsv'))
cp<-merge(classify(replay),oldprograms,by=c('database','pathway'),suffixes=c('_replay','_original'))
stopifnot(nrow(cp)==nrow(oldprograms),all(cp$program_class_replay==cp$program_class_original),max(abs(cp$support_score_replay-cp$support_score_original))<1e-8)
write_tsv(data.table(check=c('historical_enrichment_row_count','historical_max_p_error','historical_max_FDR_error','historical_classification_identical','historical_max_score_error'),value=c(nrow(joined),max(abs(joined$p_value_replay-joined$p_value_original)),max(abs(joined$FDR_replay-joined$FDR_original)),as.numeric(all(cp$program_class_replay==cp$program_class_original)),max(abs(cp$support_score_replay-cp$support_score_original))),passed=TRUE),'provenance/pathway_historical_replay.tsv')

chosen<-rbindlist(lapply(c('conserved','ovary_specific'),function(cls)head(programs[program_class==cls],10)))
write_tsv(chosen,'tables/shared_ageing_pathways_source.tsv')
plot_one<-function(cls,title,col,xlab){
  d<-copy(chosen[program_class==cls]);nall<-programs[program_class==cls,.N]
  if(!nrow(d))return(ggplot()+annotate('text',x=.5,y=.5,label='No pathways met FDR < 0.05')+xlim(0,1)+ylim(0,1)+theme_void()+labs(title=title))
  d[,label:=vapply(sub(' \\(GO:[0-9]+\\)$','',pathway),function(z)paste(strwrap(z,width=43),collapse='\n'),character(1))]
  d[,key:=paste(database,pathway,sep=' | ')]
  d[,key:=factor(key,levels=rev(key))]
  ggplot(d,aes(support_score,key))+geom_point(color=col,size=3)+scale_y_discrete(labels=setNames(d$label,d$key))+scale_x_continuous(expand=expansion(mult=c(.08,.12)))+labs(title=title,subtitle=sprintf('Showing %d of %d supported pathways',nrow(d),nall),x=xlab,y=NULL)+theme_classic(base_size=10)+theme(plot.title=element_text(face='bold',size=12),plot.subtitle=element_text(size=9),axis.text.y=element_text(color='black',size=9),axis.text.x=element_text(color='black'),plot.margin=margin(8,15,8,8))
}
p1<-plot_one('conserved','Conserved programs','#5E8C61','Evidence score: -log10(max group minimum FDR)')
p2<-plot_one('ovary_specific','Ovary-supported only','#B45A64','Evidence score: -log10(ovary minimum FDR)')
plot<-(p1+p2)+plot_annotation(tag_levels='a',caption='Conserved = enrichment support in both model groups; ovary-only support does not establish tissue specificity.\nHistorical exploratory rule: minimum within-analysis BH FDR across feature-set/rank tests; group threshold <0.05.',theme=theme(plot.tag=element_text(face='bold',size=13),plot.caption=element_text(size=9,hjust=0)))
ggsave(file.path(out,'figures/shared_ageing_pathways.pdf'),plot,device=cairo_pdf,width=15,height=8)
ggsave(file.path(out,'figures/shared_ageing_pathways.png'),plot,device=ragg::agg_png,width=15,height=8,dpi=300,bg='white')
integrity[,md5_after:=unname(tools::md5sum(path))];integrity[,unchanged:=md5_before==md5_after];stopifnot(all(integrity$unchanged))
write_tsv(integrity,'provenance/pathway_input_integrity.tsv')
counts<-programs[,.N,by=program_class];write_tsv(counts,'tables/pathway_program_counts.tsv')
section<-c('## Shared ageing pathways','','[Updated two-panel plot](figures/shared_ageing_pathways.pdf) · [Detailed methods and results](PATHWAY_ANALYSIS.md)','','The historical pathway analysis has been rerun with the current HVG-ElasticNet and VersionB-PCA/glmnet features. See the detailed report for supported pathway counts, the complete statistics, and the full-pathway multiplicity sensitivity. The original “ovary-specific” class is labelled “Ovary-supported only” because absence of tAge significance does not establish tissue specificity.','')
rp<-file.path(out,'README.md');txt<-readLines(rp,warn=FALSE);marker<-which(txt=='## Shared ageing pathways');if(length(marker))txt<-head(txt,marker[1]-1)
writeLines(c(txt,'',section),rp)
report<-c('# Shared ageing pathway analysis','',
 sprintf('Common gene universe: **%d genes**. Tested **%d enrichment rows** across the archived GO Biological Process 2025, Reactome Pathways 2024 and KEGG Mouse 2019 libraries.',length(universe),nrow(enrich)),'',
 paste0('- ',counts$program_class,': ',counts$N),'',
 '## Method','',
 'The query sets are active HVG-ElasticNet genes, the top 100 and top 500 VersionB-PCA genes ranked by absolute standardized backprojected weight, and active genes from each tAge variant. PCA top-N selection precedes intersection with the common gene universe, as in the original. The background is the current ovary eligible gene universe intersected with both mapped tAge feature tables. GMT symbols are matched case-insensitively; eligible pathways contain 5–500 genes after background intersection.','',
 'Overrepresentation uses a one-sided Fisher exact test (hypergeometric upper tail). The historical procedure tests pathways having at least one query overlap and applies BH within each query and library. `FDR_full_pathway_family` additionally includes all eligible pathways in the multiplicity family, with zero-overlap pathways assigned P=1. The corresponding classification sensitivity is saved separately.','',
 'Four rank analyses use absolute weights of all four clocks, with genes outside each model assigned weight zero in the common background. The historical rank-sum normal approximation is reproduced exactly: average ranks for ties, variance n1*n0*(N+1)/12 without a tie correction or continuity correction, one-sided upper tail, and BH within each analysis/library. This is an exploratory absolute-weight rank test, not GSEA/NES or a directional pathway activity test. `NES_or_OR` retains the historical column name but contains a sample odds ratio for ORA and a z score for rank tests.','',
 'Within each pathway/library, the minimum FDR across ovary analyses and the minimum across tAge analyses define the two model groups. Both below 0.05 gives conserved/shared support; only the ovary group below 0.05 gives the historical ovary_specific class. These minima are descriptive and are not adjusted again across the multiple feature-set/rank analyses. Conserved scores are -log10(max of the two minima); ovary-only scores are -log10(ovary minimum FDR). Up to ten pathways per panel are selected by decreasing score, then pathway name.','',
 'The panel-b x-axis label is corrected to reflect its actual score. The original category name is retained in source tables for compatibility, but the figure says “Ovary-supported only”. Neither significance in one group but not another nor this classification is a formal between-clock difference test. Shared support does not establish agreement in gene direction or biologically conserved activity. This analysis is based on frozen model features, not differential expression from external projection samples.','',
 '## Results shown','',
 '| Class | Database | Pathway | Ovary minimum FDR | tAge minimum FDR | Evidence score |','|---|---|---|---:|---:|---:|',
 vapply(seq_len(nrow(chosen)),function(i){z<-chosen[i];sprintf('| %s | %s | %s | %.4g | %.4g | %.3f |',z$program_class,z$database,gsub('|','/',z$pathway,fixed=TRUE),z$ovary_min_FDR,z$tage_min_FDR,z$support_score)},character(1)),'',
 '## Files and reproduction','',
 '- `figures/shared_ageing_pathways.pdf` and `.png`: updated plot.',
 '- `tables/pathway_overlap_enrichment.tsv`: all tested rows, effect statistics, nominal P, BH FDR, sensitivity FDR and gene overlaps.',
 '- `tables/pathway_program_classification.tsv`: all pathway classifications and the individual supporting analyses.',
 '- `tables/pathway_program_classification_full_family_sensitivity.tsv`: ORA multiplicity sensitivity.',
 '- `tables/shared_ageing_pathways_source.tsv`: exact plotted rows.',
 '- `tables/pathway_query_gene_membership.tsv` and `pathway_library_membership.tsv`: auditable gene sets.',
 '- `inputs/gene_sets/`: unchanged archived library snapshots.',
 '- `provenance/pathway_historical_replay.tsv`: every historical P/FDR and class reproduced with the ported implementation.',
 '- `provenance/pathway_input_integrity.tsv`: input and pre-existing figure integrity.',
 '', 'Run the recorded R runtime on `scripts/05_pathways.R`, followed by the pinned Python runtime on `scripts/06_validate_pathways.py`. Both are included in `scripts/run_all.sh`. No locked model, external projection or historical result is modified.')
writeLines(report,file.path(out,'PATHWAY_ANALYSIS.md'))
message('PASS: pathway analysis and figure generated.');print(counts)
