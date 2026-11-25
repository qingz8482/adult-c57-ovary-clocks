suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg);library(ggrepel);library(matrixStats)})
figdir<-file.path(out,'figures');compat<-file.path(out,'inputs/reference_adapter/cross_species');gcompat<-file.path(out,'inputs/reference_adapter/gse');ref_tables<-file.path(out,'tables/reference_panels')
fam<-function(x)ifelse(x=='HVG-ElasticNet','elasticnet','pca_glmnet')
modelname<-c(elasticnet='HVG-ElasticNet',pca_glmnet='VersionB-PCA')
write_tsv<-function(x,path){stopifnot(startsWith(normalizePath(dirname(path),mustWork=FALSE),normalizePath(out)));dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);fwrite(x,path,sep='\t',na='NA')}
cwrite<-function(x,rel){p<-file.path(compat,rel);dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);fwrite(x,p,sep='\t',na='NA')}
path_module<-function(...){r<-file.path(...);if(r=='1_publication/panels')return(figdir);if(grepl('^1_publication/(panels|assembled)/',r))return(file.path(figdir,basename(r)));if(grepl('^1_publication/source_data/',r))return(file.path(figdir,'source_data',basename(r)));if(grepl('^1_publication/captions/',r))return(file.path(out,'reports',basename(r)));if(grepl('^1_publication/',r))return(file.path(ref_tables,basename(r)));file.path(compat,r)}
map_ref<-function(r,part){if(grepl('^figures/source_data/',r))return(file.path(figdir,'source_data',part,basename(r)));if(grepl('^figures/',r))return(file.path(figdir,basename(r)));if(part=='2_revision_defensibility')return(file.path(ref_tables,basename(r)));p<-file.path(out,'inputs/reference_adapter',part,r);if(file.exists(p))p else file.path(old,part,r)}
prev<-function(...)map_ref(file.path(...),'2_revision_defensibility')
paudit<-function(...)map_ref(file.path(...),'3_gtex_age_signal_adequacy')
pfinal<-function(...)map_ref(file.path(...),'4_final_reporting_revision')
pbase<-paudit
for(d in c('2_revision_defensibility','3_gtex_age_signal_adequacy','4_final_reporting_revision'))dir.create(file.path(figdir,'source_data',d),recursive=TRUE,showWarnings=FALSE)
save_bundle<-function(plot,stem,width_mm=183,height_mm=125){w<-width_mm/25.4;h<-height_mm/25.4;grDevices::svg(paste0(stem,'.svg'),width=w,height=h,pointsize=7,family='sans');print(plot);dev.off();grDevices::cairo_pdf(paste0(stem,'.pdf'),width=w,height=h,pointsize=7,family='sans');print(plot);dev.off();ragg::agg_png(paste0(stem,'.png'),width=w,height=h,units='in',res=600);print(plot);dev.off();ragg::agg_tiff(paste0(stem,'.tiff'),width=w,height=h,units='in',res=600,compression='lzw');print(plot);dev.off()}
save_pub_r<-save_bundle
palette_age<-c('20-29'='#4C78A8','30-39'='#72A0C1','40-49'='#9ECAC3','50-59'='#E3B65F','60-69'='#E07A5F','70-79'='#A8434C')
palette_contract<-c(neutral_dark='#272727',neutral_mid='#767676',neutral_light='#D8D8D8',signal_blue='#3182BD',signal_teal='#33B5A5',accent_red='#D24B40',accent_orange='#E28E2C')
age_levels<-names(palette_age)
theme_context<-function()theme_classic(base_size=7,base_family='sans')+theme(axis.line=element_line(linewidth=.35),axis.ticks=element_line(linewidth=.35),strip.background=element_blank(),strip.text=element_text(face='bold'),plot.title=element_text(face='bold',size=7.7),panel.grid=element_blank())
# Adapter columns preserve the reference plotting interfaces; values come only from current results.
cp<-rt('tables/frozen_predictions.tsv');cp[,model_family:=fam(model)]
cm<-rt('tables/model_gene_coverage.tsv');cm[,model_family:=fam(model)]
master0<-rt('tables/stageA_gene_effect_master.tsv');master0[,`:=`(elasticnet_hvg_3000=hvg_selected,pca_versionB_3000=pca_selected,elasticnet_nonzero=is.finite(hvg_weight)&hvg_weight!=0,elasticnet_expression_space_weight=hvg_weight,pca_expression_space_weight=pca_weight)]
cwrite(master0,'derived_data/stageA/mouse_human_gene_effect_master.tsv')
genes_ref<-Reduce(union,lapply(models,function(m)m$representation$filtered_genes));shared<-intersect(models[[1]]$representation$genes,models[[2]]$representation$genes)
orth<-rt('inputs/orthology_one2one.tsv')
analysis_current<-readRDS(file.path(out,'inputs/analysis_inputs.rds'))
coverage_fun<-function(genes,label,weights=rep(1,length(genes))){d<-orth[match(genes,mouse_gene_symbol)];meas<-d$human_gene_id_stable%in%rownames(analysis_current$cohorts$GTEx$x);mapped<-!is.na(d$human_gene_id_stable);data.table(gene_set=label,n_total=length(genes),n_one2one=sum(mapped),n_human_detected=sum(meas),raw_gene_coverage=mean(mapped),human_detected_coverage=mean(meas),absolute_model_weight_coverage=sum(abs(weights[meas]))/sum(abs(weights)))}
ccov<-rbindlist(list(coverage_fun(genes_ref,'integrated_universe'),coverage_fun(models[[1]]$representation$genes,'elasticnet_hvg_3000',models[[1]]$coefficients$coefficient),coverage_fun(models[[2]]$representation$genes,'pca_versionB_3000',models[[2]]$coefficients$coefficient),coverage_fun(shared,paste0('deployment_shared_artifact_',length(shared))),coverage_fun(models[[1]]$coefficients[coefficient!=0,gene],'elasticnet_nonzero',models[[1]]$coefficients[coefficient!=0,coefficient])))
cwrite(ccov,'results/tables/table_orthology_coverage.tsv')
enr0<-rt('tables/stageA_weight_direction_enrichment.tsv');enr0[,signature:=sub('^hvg','elasticnet',signature)]
cwrite(enr0[,.(signature,observed_weighted_signed_t)],'results/tables/table_stageA_enrichment.tsv')
enr0[,null_q95:=vapply(signature,function(sn)quantile(rt(paste0('null_distributions/stageA_weight_',sub('^elasticnet','hvg',sn),'.tsv'))$statistic,.95),numeric(1))]
cwrite(enr0[,.(signature,null_mean,null_q95,empirical_p)],'results/tables/table_random_signature_benchmarks.tsv')
cwrite(cp[cohort=='GTEx'],'derived_data/stageB/gtex_frozen_zero_shot_scores.tsv')
gp<-copy(cp[cohort=='GSE202601']);gp[,`:=`(projection_sample_id=sample_id,sample_id=donor_id,mapped_gene_count=n_measured_genes,absolute_model_weight_coverage=absolute_weight_coverage,projection_qc_class=ifelse(absolute_weight_coverage>=.8,'HIGH_WEIGHT_COVERAGE','MODERATE_WEIGHT_COVERAGE'))]
fwrite(gp,file.path(gcompat,'projection/GSE202601_two_integrated_clock_projections.tsv'),sep='\t')
gcm<-cm[cohort=='GSE202601'];gcm[,`:=`(high_confidence_one_to_one_mapping=one2one_mapped,human_gene_detected_after_filter=measured,absolute_expression_space_weight=abs(raw_coefficient))]
fwrite(gcm,file.path(gcompat,'projection/GSE202601_model_gene_coverage.tsv'),sep='\t')
# Reference label 'detected' now means measured in the current frozen input, including zero.
write_tsv(data.table(reference_field='human_gene_detected_after_filter',current_definition='Measured in the full prefilter human input; historical filtering is a separate sensitivity'),file.path(ref_tables,'reference_field_semantics.tsv'))
