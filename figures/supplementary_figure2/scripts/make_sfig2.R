#!/usr/bin/env Rscript
# Supplementary Figure 2 (rebutal4 release hvg_exclusion_vbpca_20261012): interpretation of the
# two locked clocks (a-d Version B-PCA, e HVG-EN). Re-plots the saved rebutal3 supple2_PCA analysis outputs (read only);
# nothing is recomputed or refitted. Changes from the rebutal3 figure:
#  - the grouped-CV lambda curve (old a) is dropped; it duplicates Supplementary Fig. 1;
#  - coefficient paths (new a) use log10 lambda, as in Supplementary Fig. 1, and carry the
#    nonzero-PC axis; the deployed lambda is the red line, the minimum-MAE lambda dashed;
#  - the contribution heatmap (new c) labels each age bin and shows chronological and predicted
#    bin-median ages on one shared absolute scale (weeks) instead of within-track z-scores;
#  - panel e is the HVG-EN nonzero-gene ORA (younger above older) in place of the PCA GSEA curves,
#    so the figure now covers gene-level interpretation of both locked clocks.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(scales)})
arg<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),fixed=TRUE)
out<-normalizePath(file.path(dirname(arg),'..'));rb<-dirname(dirname(out))
src<-file.path(rb,'rebutal3_batcheffects removal_CV_LOO/model_development/VersionB-PCA/diagnostic/supple2_PCA/derived_data')
lock<-file.path(rb,'rebutal4/lockedmodel/VersionB-PCA-glmnet')
r3lock<-file.path(rb,'rebutal3_batcheffects removal_CV_LOO/model_development/locked_models/count_tuning_20261011/VersionB-PCA-glmnet')
for(d in c('panels','final','source_data'))dir.create(file.path(out,d),showWarnings=FALSE)
sd_out<-function(x,n)fwrite(x,file.path(out,'source_data',paste0(n,'.tsv')),sep='\t')

# The saved analysis was run on the rebutal3 lock; it must be the rebutal4 locked model.
for(f in c('model.rds','coefficients.csv','pc_coefficients.csv','pca_loadings.csv','scaler.csv','intercept.txt'))
  stopifnot(unname(tools::md5sum(file.path(lock,f)))==unname(tools::md5sum(file.path(r3lock,f))))
obj<-readRDS(file.path(src,'plot_objects.rds'));invisible(list2env(obj,envir=environment()))
stopifnot(abs(lambda-jsonlite::fromJSON(file.path(lock,'hyperparameters.json'))$lambda)<1e-3)

theme_pub<-theme_classic(base_size=6.5,base_family='Arial')+theme(
 axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3),
 plot.margin=margin(3,4,3,4),legend.title=element_text(size=5.5),legend.text=element_text(size=5),
 plot.title=element_text(size=6.5,face='bold'),plot.subtitle=element_text(size=5.5),
 plot.tag=element_text(size=8,face='bold'),panel.grid=element_blank())
theme_set(theme_pub)

## a | coefficient paths, log10 lambda, nonzero-PC top axis
path_dt[,l10:=log10(lambda)];cv_dt[,l10:=log10(lambda)]
dep<-lambda;min_lambda<-cv_dt$lambda[which.min(cv_dt$mean_mae)]
top<-cv_dt[sapply(c(-4,-2,0,1,log10(dep),2.5),function(v)which.min(abs(l10-v)))]
mark<-path_dt[selected==TRUE,.SD[which.min(abs(lambda-dep))],by=pc]
sd_out(path_dt[,.(pc,lambda,log10_lambda=l10,coefficient,active_in_deployed_fit=selected)],'SFig2a_coefficient_paths')
p_a<-ggplot(path_dt,aes(l10,coefficient,group=pc))+
 geom_hline(yintercept=0,colour='grey70',linewidth=.25)+
 geom_line(data=path_dt[selected==FALSE],colour='grey78',linewidth=.3)+
 geom_line(data=path_dt[selected==TRUE],aes(colour=coefficient),linewidth=.5)+
 geom_vline(xintercept=log10(min_lambda),linetype=2,colour='grey45',linewidth=.3)+
 geom_vline(xintercept=log10(lambda),colour='#B2182B',linewidth=.4)+
 geom_text(data=mark,aes(label=pc),nudge_x=.25,hjust=0,size=1.8,show.legend=FALSE)+
 scale_colour_gradient2(low='#2166AC',mid='white',high='#B2182B',guide='none')+
 scale_x_continuous(expression(log[10]~lambda),breaks=-4:3,
   sec.axis=dup_axis(name='Nonzero PCs in deployment fit',breaks=top$l10,labels=top$n_nonzero))+
 labs(title='Version B–PCA: coefficient paths of the 10 PC inputs',
   subtitle=sprintf('3,000 genes; alpha = %g; deployed lambda = %.1f (red); no calibration',alpha,lambda),
   y='PC coefficient (weeks / score)',tag='a')

## b | conditional PC ablation (unchanged)
make_b<-function(d,subtitle,tag='b')ggplot(d,aes(omitted_pc,delta_mae))+geom_hline(yintercept=0,colour='grey60',linewidth=.3)+
 geom_errorbar(aes(ymin=ci_low,ymax=ci_high),width=.15,linewidth=.4,colour='grey30')+
 geom_point(aes(colour=delta_mae),size=1.8)+scale_colour_gradient2(low='#2166AC',mid='grey50',high='#B2182B',guide='none')+
 labs(title='Version B–PCA: conditional PC ablation',subtitle=subtitle,x='Omitted active PC',y='Change in MAE (weeks)',tag=tag)
p_b<-make_b(lopo_summary,'Paired sample bootstrap; fixed PCA space')
p_b_study<-make_b(study_summary,'Equal study weighting; study bootstrap','')
sd_out(lopo_summary,'SFig2b_pc_ablation_sample_bootstrap');sd_out(study_summary,'SFig2b_pc_ablation_study_bootstrap')

## c | age-binned PC contributions, with bin ages on the axis and absolute age tracks
bins<-fread(file.path(src,'panel_d_sample_bins.tsv'))
br<-bins[,.(n=.N,age_min=min(age_week),age_max=max(age_week),age_median=median(age_week),pred_median=median(predicted_age)),by=age_bin][order(age_bin)]
stopifnot(nrow(br)==20,sum(br$n)==136,all(abs(br$age_median-unique(d_plot[,.(age_bin,age_week)])[order(age_bin)]$age_week)<1e-9))
br[,lab:=ifelse(round(age_min)==round(age_max),sprintf('%d',round(age_median)),sprintf('%d–%d',round(age_min),round(age_max)))]
d_plot[,age_bin:=factor(age_bin,levels=br$age_bin)]
lim_d<-max(abs(d_plot$contribution))
p_c_main<-ggplot(d_plot,aes(age_bin,factor(pc,levels=rev(pcs)),fill=contribution))+geom_tile()+
 scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-lim_d,lim_d),name='Median PC\ncontribution\n(weeks)')+
 guides(fill=guide_colourbar(barheight=unit(13,'mm'),barwidth=unit(2.5,'mm')))+
 scale_x_discrete(labels=setNames(br$lab,br$age_bin))+
 labs(x='Age bin (weeks; 6–7 profiles per bin, 136 training profiles ordered by chronological age)',y=NULL)+
 theme(axis.ticks=element_blank(),axis.line=element_blank(),axis.text.x=element_text(angle=45,hjust=1,size=4.8))
ann<-melt(br[,.(age_bin,`Chronological age`=age_median,`Predicted age`=pred_median)],id.vars='age_bin',variable.name='track',value.name='weeks')
ann[,age_bin:=factor(age_bin,levels=br$age_bin)]
sd_out(merge(br[,.(age_bin,n,age_min,age_max,age_median,predicted_median=pred_median)],
  dcast(d_plot,age_bin~pc,value.var='contribution'),by='age_bin'),'SFig2c_age_bin_contributions')
p_c_ann<-ggplot(ann,aes(age_bin,factor(track,levels=c('Predicted age','Chronological age')),fill=weeks))+geom_tile()+
 geom_text(aes(label=round(weeks),colour=weeks>45),size=1.5)+
 scale_colour_manual(values=c(`TRUE`='white',`FALSE`='grey15'),guide='none')+
 scale_fill_gradient(low='grey96',high='grey20',limits=c(0,80),breaks=c(0,40,80),name='Bin median\nage (weeks)')+
 guides(fill=guide_colourbar(barheight=unit(8,'mm'),barwidth=unit(2.5,'mm')))+
 scale_x_discrete(labels=NULL,breaks=NULL)+labs(x=NULL,y=NULL)+theme(axis.ticks=element_blank(),axis.line=element_blank())
p_c<-wrap_elements(full=p_c_ann/p_c_main+plot_layout(heights=c(.55,1)))+labs(tag='c')

## d | back-projected gene weights (unchanged)
e_dt[,gene:=factor(gene,levels=rev(genes_show))];e_dt[,component:=factor(component,levels=c(pcs,'Integrated'))]
lim_e<-max(abs(e_dt$weight))
sd_out(e_dt[order(display_order,component),.(gene,component,weight_weeks_per_ln1pCPM=weight)],'SFig2d_gene_weights')
p_d<-ggplot(e_dt,aes(component,gene,fill=weight))+geom_tile(colour='white',linewidth=.12)+
 scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-lim_e,lim_e),name='Gene weight\n(weeks / ln1p CPM)')+
 guides(fill=guide_colourbar(barheight=unit(18,'mm'),barwidth=unit(3,'mm')))+
 labs(title='Version B–PCA: back-projected gene weights',x=NULL,y=NULL,tag='d')+
 theme(axis.text.y=element_text(size=5.6),axis.text.x=element_text(angle=45,hjust=1),axis.line=element_blank(),axis.ticks=element_blank())

## e | HVG-EN nonzero-gene pathway over-representation (younger above older, shared x scale)
# Reads the saved ORA results of figure1_model_development/scripts/enrich_elasticnet_nonzero_genes.R; nothing is re-tested.
ora_dir<-file.path(rb,'rebutal4/figure1_model_development/diagnostics/09_selected_gene_diagnostics/elasticnet_hvg_coefficients/pathway_enrichment')
q_cut<-0.05
sig<-fread(file.path(ora_dir,'ora_significant_BH0.05.tsv'));stopifnot(nrow(sig)>0,all(sig$p.adjust<q_cut))
mapped<-fread(file.path(ora_dir,'gene_id_mapping_summary.tsv'))[,setNames(n_entrez,set)]
db_lab<-c(GO_BP='BP',GO_CC='CC',GO_MF='MF',KEGG='KEGG',Reactome='Reactome',Hallmark='Hallmark')
sig[,`:=`(neg_log_q=-log10(p.adjust),db=factor(db_lab[database],levels=db_lab))]
setorder(sig,direction,db,p.adjust,Description);sig[,term:=sprintf('%s (%s)',Description,db)]
n_dir<-sig[,.N,by=direction][,setNames(N,direction)];stopifnot(all(c('older','younger')%in%names(n_dir)))
x_max<-max(sig$neg_log_q)*1.18
dir_col<-c(older='#F04E37',younger='#4DB6CC')
sd_out(sig[,!c('db','term')],'SFig2e_hvg_en_ora_significant')
ora_panel<-function(d,title,show_x){
 dt<-sig[direction==d];dt[,term:=factor(term,levels=rev(term))]
 ggplot(dt,aes(neg_log_q,term))+
  geom_col(width=.72,fill=dir_col[[d]],colour='black',linewidth=.2)+
  geom_text(aes(label=GeneRatio),hjust=0,nudge_x=.08,size=1.6)+
  geom_vline(xintercept=-log10(q_cut),linetype='dashed',linewidth=.3,colour='grey40')+
  scale_x_continuous(limits=c(0,x_max),expand=expansion(mult=c(0,.02)))+
  labs(title=title,x=if(show_x)expression(-log[10]~'BH-adjusted p') else NULL,y=NULL)+
  theme(axis.text.y=element_text(size=5.2,colour='black'),
   plot.title=element_text(size=6,face='bold',colour=dir_col[[d]]),plot.title.position='plot',
   axis.text.x=if(show_x)element_text() else element_blank(),axis.ticks.length.x=unit(if(show_x)1.5 else 0,'pt'))
}
p_e<-wrap_elements(full=(ora_panel('younger',sprintf('Toward younger predicted age (%d terms)',n_dir[['younger']]),FALSE)/
   ora_panel('older',sprintf('Toward older predicted age (%d terms)',n_dir[['older']]),TRUE)+
   plot_layout(heights=c(n_dir[['younger']],n_dir[['older']])))+
  plot_annotation(title='HVG–EN: pathway over-representation of nonzero genes',
   subtitle=sprintf('Hypergeometric ORA, BH q < %.2f; %d younger / %d older genes;\nuniverse = %d HVG candidates; labels = genes in term / genes in list',
     q_cut,mapped[['younger']],mapped[['older']],mapped[['universe']]),
   theme=theme(plot.title=element_text(size=6.5,face='bold'),plot.subtitle=element_text(size=5.2),plot.margin=margin(4,4,3,12))))+labs(tag='e')

## export
save_all<-function(p,stem,w,h){
 ggsave(paste0(stem,'.pdf'),p,width=w,height=h,units='mm',device=cairo_pdf,bg='white')
 ragg::agg_png(paste0(stem,'.png'),width=w,height=h,units='mm',res=600,background='white');print(p);invisible(dev.off())
}
panels<-list(a=list(p_a,90,58),b=list(p_b,90,58),b_study_weighted=list(p_b_study,75,58),c=list(p_c,183,52),d=list(p_d,105,108),e_ora_vertical=list(p_e,90,118))
for(n in names(panels))save_all(panels[[n]][[1]],file.path(out,'panels',paste0('SFig2',n)),panels[[n]][[2]],panels[[n]][[3]])
final_plot<-(p_a+p_b+plot_layout(widths=c(1.25,1)))/p_c/(p_d+p_e+plot_layout(widths=c(1,1.05)))+plot_layout(heights=c(1,.95,2.1))
save_all(final_plot,file.path(out,'final','Supplementary_Figure2_clock_gene_interpretation'),183,232)
cat('done\n')
