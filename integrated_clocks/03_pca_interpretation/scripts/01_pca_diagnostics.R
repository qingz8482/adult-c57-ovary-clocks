#!/usr/bin/env Rscript
# Figure contract: quantify the frozen model's variance, age associations and
# signed gene contributions. Quantitative panels; PDF and 300-dpi PNG with TSV
# source data. Training associations are descriptive, not validation accuracy.
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
arg <- sub('^--file=', '', grep('^--file=',commandArgs(),value=TRUE)[1])
arg <- gsub("~+~", " ", arg, fixed=TRUE)
root <- normalizePath(file.path(dirname(arg),'..'))
development <- normalizePath(file.path(root,'../..'))
source(file.path(development,'scripts/core.R'))
out <- file.path(root,'10_pca_diagnostics')
model_path <- file.path(development,'locked_models/count_tuning_20261011/VersionB-PCA-glmnet/model.rds')
m <- readRDS(model_path)
inp <- readRDS(file.path(development,'shared_inputs/prepared_expression.rds'))
stopifnot(m$recipe=='VersionB-PCA', m$calibration=='none', m$target=='age_week')
meta <- as.data.table(inp$metadata)[match(m$training_keys,sample_key)]
x <- inp$expression[,m$training_keys,drop=FALSE]
rep <- m$representation
stopifnot(identical(meta$sample_key,m$training_keys),identical(rownames(rep$rotation),rep$genes))
z <- sweep(sweep(t(x[rep$genes,,drop=FALSE]),2,rep$center,'-'),2,rep$scale,'/')
stopifnot(all(is.finite(z)))
scores <- z %*% rep$rotation
pcs <- colnames(scores)
coef <- m$pc_coefficients
active <- pcs[coef!=0]
scope <- 'locked_count_tuning_all_136'
score_dt <- cbind(meta[,.(sample_key,age_week,study,modality,dataset_folder)],as.data.table(scores))
total_variance <- sum(apply(z,2,var))
variance <- apply(scores,2,var)/total_variance
all_summary <- rbindlist(lapply(seq_along(pcs),function(i) data.table(
 artifact_scope=scope,pc=pcs[i],pc_number=i,n_training_samples_for_pca=nrow(scores),
 variance_explained=variance[i],cumulative_variance_through_pc=sum(variance[seq_len(i)]),
 pc_age_pearson=cor(scores[,i],meta$age_week),pc_age_spearman=cor(scores[,i],meta$age_week,method='spearman'),
 raw_pc_coefficient=coef[i],pc_coefficient=coef[i],active_in_clock=coef[i]!=0,
 clock_coefficient_direction=if(coef[i]>0)'higher_PC_score_increases_predicted_age' else if(coef[i]<0)'higher_PC_score_decreases_predicted_age' else 'zero_coefficient',
 pc_score_age_axis=if(cor(scores[,i],meta$age_week)>=0)'higher_PC_score_older' else 'higher_PC_score_younger'
)))
summary <- all_summary[active_in_clock==TRUE]
# Reconstruct only the existing selector statistics; assert frozen gene membership.
ranking <- fit_ranking(x,meta,'VersionB',m$config)
stopifnot(identical(choose_genes(ranking,length(rep$genes)),rep$genes))
features <- copy(ranking$table)[gene %in% rep$genes]
setnames(features,'gene','gene_id')
features[,version_b_age_direction:=ifelse(beta>=0,'age-increasing','age-decreasing')]
loadings <- melt(as.data.table(rep$rotation[,active,drop=FALSE],keep.rownames='gene_id'),
 id.vars='gene_id',variable.name='pc',value.name='pc_loading',variable.factor=FALSE)
loadings <- merge(loadings,summary[,.(artifact_scope,pc,pc_number,pc_age_pearson,pc_coefficient)],by='pc')
loadings <- merge(loadings,features,by='gene_id')
loadings[,`:=`(abs_pc_loading=abs(pc_loading),derived_gene_pc_coefficient=pc_loading*pc_coefficient)]
loadings[,`:=`(abs_derived_gene_pc_coefficient=abs(derived_gene_pc_coefficient),
 pc_loading_age_direction=ifelse(pc_loading*pc_age_pearson>=0,'aged-aligned-loading','young-aligned-loading'),
 clock_contribution_direction=ifelse(derived_gene_pc_coefficient>0,'positive-derived-coefficient','negative-derived-coefficient'))]
setorder(loadings,pc_number,-abs_pc_loading,gene_id)
loadings[,pc_loading_rank:=seq_len(.N),by=pc]
counts <- loadings[,.(n_genes_in_pc=.N,n_genes_nonzero_loading=sum(pc_loading!=0),
 n_aged_predictor_by_clock_direction=sum(derived_gene_pc_coefficient>0),
 n_young_predictor_by_clock_direction=sum(derived_gene_pc_coefficient<0),
 n_age_increasing_by_version_b_beta=sum(beta>=0),n_age_decreasing_by_version_b_beta=sum(beta<0),
 total_unique_genes_across_selected_pcs=length(rep$genes)),by=.(artifact_scope,pc,pc_number)]
write_tsv <- function(d,name) fwrite(d,file.path(out,paste0(name,'.tsv')),sep='\t')
write_tsv(summary,'final_pca_clock_selected_pc_summary')
write_tsv(all_summary,'final_pca_clock_all_input_pc_summary')
write_tsv(counts,'final_pca_clock_gene_counts_by_pc')
write_tsv(loadings,'final_pca_clock_gene_pc_loadings_with_feature_scores')
write_tsv(loadings[pc_loading_rank<=50],'final_pca_clock_top50_gene_pc_loadings')
write_tsv(all_summary[,.(artifact_scope,pc,raw_pc_coefficient,pc_coefficient,active_in_clock)],'pca_glmnet_pc_coefficients')
write_tsv(score_dt,'pca_score_summary_from_rds')
write_tsv(all_summary[,.(artifact_scope,pc,pc_number,variance_explained,cumulative_variance_through_pc,active_in_clock)],'pca_variance_summary')
write_tsv(as.data.table(rep$rotation,keep.rownames='gene_id'),'all_input_pc_loadings')

pred <- drop(scores %*% coef)+m$model_intercept
gene_pred <- drop(t(x[rep$genes,,drop=FALSE]) %*% m$coefficients$coefficient)+m$intercept
native <- predict_model(m,x)$pred
eff <- drop(rep$rotation %*% coef)
checks <- data.table(check=c('PC_vs_native_prediction','gene_vs_native_prediction','integrated_rank_vs_standardized_coefficient',
 'rotation_orthonormal','training_center','training_scale','variance_sum_at_most_one'),
 error=c(max(abs(pred-native)),max(abs(gene_pred-native)),max(abs(eff-m$coefficients$standardized_coefficient)),
 max(abs(crossprod(rep$rotation)-diag(ncol(rep$rotation)))),max(abs(colMeans(z))),max(abs(apply(z,2,sd)-1)),max(0,sum(variance)-1)))
checks[,passed:=error<1e-8]
fwrite(checks,file.path(root,'provenance/pca_validation.tsv'),sep='\t')
stopifnot(all(checks$passed),length(active)==4,nrow(meta)==136,length(rep$genes)==3000)
write_tsv(data.table(sample_key=meta$sample_key,native_prediction=native,pc_prediction=pred,gene_prediction=gene_pred),'prediction_replay')
contrib <- sweep(scores,2,coef,'*')
write_tsv(cbind(meta[,.(sample_key,age_week,study,modality)],as.data.table(contrib)),'training_pc_contributions')
theme_set(theme_bw(base_size=12)+theme(plot.title=element_text(face='bold'),legend.position='bottom',strip.text=element_text(face='bold')))
save_plot <- function(p,stem,w=9,h=5.8) {
 ggsave(file.path(out,paste0(stem,'.pdf')),p,width=w,height=h,device=cairo_pdf)
 ggsave(file.path(out,paste0(stem,'.png')),p,width=w,height=h,dpi=300,bg='white')
}
all_summary[,pc:=factor(pc,levels=pcs)];summary[,pc:=factor(pc,levels=pcs)]
p<-ggplot(all_summary,aes(pc,100*variance_explained,fill=active_in_clock))+geom_col(width=.7)+
 geom_line(aes(y=100*cumulative_variance_through_pc,group=1),colour='grey30')+geom_point(aes(y=100*cumulative_variance_through_pc),size=1.5)+
 scale_fill_manual(values=c('FALSE'='grey80','TRUE'='#4C78A8'),labels=c('Zero coefficient','Active coefficient'))+
 labs(title='Locked VersionB-PCA: variance explained',subtitle='Bars: individual variance; line: cumulative variance across the 10 PC inputs',x='PC input',y='Variance explained (%)',fill=NULL)
save_plot(p,'13_pca_diagnostics_variance_explained')
p<-ggplot(all_summary,aes(pc,pc_coefficient,fill=active_in_clock))+geom_hline(yintercept=0,colour='grey50')+geom_col(width=.7)+
 scale_fill_manual(values=c('FALSE'='grey80','TRUE'='#4C78A8'))+labs(title='Locked VersionB-PCA: PC coefficients',subtitle='Frozen coefficients in weeks per score unit; no calibration',x='PC input',y='PC coefficient (weeks / score unit)')+theme(legend.position='none')
save_plot(p,'13_pca_diagnostics_pc_coefficients')
p<-ggplot(score_dt,aes(PC1,PC2,colour=age_week,shape=modality))+geom_point(size=2.8,alpha=.85)+scale_colour_viridis_c()+
 labs(title='Locked VersionB-PCA: training PC scores',subtitle='136 training profiles; descriptive projection in the frozen basis',
 x=sprintf('PC1 (%.1f%%)',100*variance[1]),y=sprintf('PC2 (%.1f%%)',100*variance[2]),colour='Age (weeks)',shape='Modality')+theme(aspect.ratio=1)
save_plot(p,'13_pca_diagnostics_pc1_pc2_scores',6.8,7.4)
p<-ggplot(summary,aes(pc,100*variance_explained))+geom_col(fill='#4C78A8',width=.68)+
 geom_text(aes(label=sprintf('%.1f%%',100*variance_explained)),vjust=-.4)+scale_y_continuous(expand=expansion(mult=c(0,.12)))+
 labs(title='Locked VersionB-PCA: active PC variance',subtitle='PCs with nonzero frozen regression coefficients',x='Active PC',y='Variance explained (%)')
save_plot(p,'final_pca_clock_selected_pc_variance')
p<-ggplot(all_summary,aes(pc_age_pearson,pc_coefficient,label=pc,colour=active_in_clock))+geom_hline(yintercept=0,linetype=2,colour='grey60')+
 geom_vline(xintercept=0,linetype=2,colour='grey60')+geom_point(size=3)+geom_text(data=all_summary[active_in_clock==TRUE],nudge_y=.015)+
 scale_colour_manual(values=c('FALSE'='grey55','TRUE'='#4C78A8'))+labs(title='PC age association and frozen clock coefficient',
 subtitle='Training associations are descriptive; coefficient direction defines GSEA orientation',x='Pearson correlation with chronological age',y='PC coefficient (weeks / score unit)')+theme(legend.position='none')
save_plot(p,'final_pca_clock_pc_age_correlation_vs_coefficient')
top<-copy(loadings[pc_loading_rank<=20]);top[,label:=paste(pc,gene_id,sep=' | ')];top[,label:=factor(label,levels=rev(label))]
p<-ggplot(top,aes(abs_derived_gene_pc_coefficient,label,fill=clock_contribution_direction))+geom_col(width=.72)+
 facet_wrap(~pc,scales='free_y',ncol=2)+scale_y_discrete(labels=function(x)sub('^.* \\| ','',x))+
 scale_fill_manual(values=c('positive-derived-coefficient'='#4C78A8','negative-derived-coefficient'='#B5524A'))+
 labs(title='Locked VersionB-PCA: top 20 genes per active PC',subtitle='Absolute PC loading × frozen PC coefficient; gene inputs are standardized',x='Absolute coefficient contribution (weeks / gene SD)',y=NULL,fill=NULL)+theme(axis.text.y=element_text(size=9))
save_plot(p,'final_pca_clock_deployment_top20_genes_per_pc',12,11)
capture.output(sessionInfo(),file=file.path(root,'provenance/sessionInfo_pca.txt'))
cat('PCA diagnostics complete:',length(rep$genes),'genes;',nrow(meta),'profiles; active PCs:',paste(active,collapse=', '),'\n')
