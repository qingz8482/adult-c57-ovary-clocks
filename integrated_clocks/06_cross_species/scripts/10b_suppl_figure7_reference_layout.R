# Supplementary Fig. 7, reference layout: same panel structure as the 117-sample-clock version
# (design, coverage, slope scatter, GTEx brackets, GSE202601, cross-model agreement), with current-clock values
# (GTEx V11 n=193, GSE202601 n=8). The message-led redesign is 10_suppl_figure7.R. No human refitting or calibration.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(ggrepel);library(ragg)})
a<-commandArgs(FALSE);here<-dirname(normalizePath(sub('^--file=','',a[grepl('^--file=',a)][1])));out<-dirname(here)
src<-file.path(out,'figures/source_data');gsrc<-file.path(src,'GSE202601')
rd<-function(p)fread(file.path(out,p))

pal<-c(elasticnet='#2878B5',pca_glmnet='#D17A22')
pal_grp<-c('Young'='#3E78B2','Reproductively aged'='#D98C3F')
mlab<-c(elasticnet='ElasticNet (HVG)',pca_glmnet='PCA/glmnet')
elab<-setNames(paste('GSE202601 ·',mlab),mlab)
fmt_p<-function(p)ifelse(p<1e-3,sprintf('%.1e',p),sprintf('%.3f',p))
theme_set(theme_classic(base_size=6.5,base_family='sans')+theme(
  axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3),axis.text=element_text(colour='black'),
  strip.background=element_rect(fill='#EFEFEF',colour='black',linewidth=.3),strip.text=element_text(face='bold',size=6.5),
  plot.tag=element_text(face='bold',size=9),legend.key.size=unit(3,'mm'),plot.margin=margin(3,4,3,4)))

# a: design
box<-function(x0,x1,y0,y1,lab,fill,col,size=1.6)list(annotate('rect',xmin=x0,xmax=x1,ymin=y0,ymax=y1,fill=fill,colour=col,linewidth=.4),annotate('text',x=(x0+x1)/2,y=(y0+y1)/2,label=lab,size=size,lineheight=.9))
arr<-function(x,xe,y,ye)annotate('segment',x=x,xend=xe,y=y,yend=ye,arrow=arrow(length=unit(1.2,'mm')),linewidth=.3)
pA<-ggplot()+coord_cartesian(xlim=c(0,10),ylim=c(0,10),expand=FALSE)+
  box(.1,4.5,8.5,9.9,'2 frozen mouse ovarian clocks\nElasticNet (HVG) + PCA/glmnet','#E7F1F8','#2878B5')+
  box(5.5,9.9,8.5,9.9,'Ensembl 116 one-to-one\northologs (high confidence)','#F1F2F3','#6F7882')+
  arr(4.6,5.4,9.2,9.2)+arr(6.2,2.6,8.4,7.65)+arr(7.7,7.7,8.4,7.65)+
  box(.1,5,6.3,7.55,'GTEx V11 ovary\n193 donors (age brackets)','#FFF0DF','#D17A22')+
  box(5.5,9.9,6.0,7.55,'GSE202601 ovary (n = 8)\n4 young (23–29 years)\n4 reproductively aged (49–54 years)','#FFF0DF','#D17A22')+
  arr(2.55,2.55,6.2,5.55)+
  box(.1,5,4.2,5.45,'Cross-gene Spearman of\nmouse vs human age slopes','white','black')+
  annotate('path',x=c(2.55,2.55,7.7,7.7),y=c(4.1,3.5,3.5,5.9),linewidth=.3)+arr(5,5,3.5,3.05)+
  box(.6,9.4,2.0,2.95,'Frozen projection of both clocks (no human refitting)','white','black')+arr(5,5,1.9,1.35)+
  box(.6,9.4,.3,1.25,'Agreement between ElasticNet and PCA/glmnet','white','black')+
  theme_void()+theme(plot.margin=margin(3,4,3,4),plot.tag=element_text(face='bold',size=9))

# b: orthology and detectability per cohort x model
cov<-fread(file.path(gsrc,'GTEx_GSE202601_gene_coverage_source_data.tsv'))
cov[,row_label:=factor(row_label,levels=rev(row_label))]
covl<-melt(cov,id.vars='row_label',measure.vars=c('raw_gene_coverage','human_detected_coverage','absolute_model_weight_coverage'),variable.name='metric',value.name='coverage')
covl[,metric:=factor(metric,levels=c('raw_gene_coverage','human_detected_coverage','absolute_model_weight_coverage'),labels=c('raw gene coverage','human detected coverage','absolute model weight coverage'))]
pB<-ggplot(covl,aes(coverage,row_label,shape=metric,colour=metric))+
  geom_vline(xintercept=.6,linetype=3,colour='#AAAAAA',linewidth=.3)+geom_vline(xintercept=.8,linetype=2,colour='#AAAAAA',linewidth=.3)+
  geom_hline(yintercept=2.5,colour='#CCCCCC',linewidth=.3)+
  geom_point(size=1.6,position=position_dodge(width=.35))+
  scale_x_continuous(labels=scales::percent_format(),limits=c(0,1),breaks=seq(0,1,.2))+
  scale_colour_manual(values=c('#777777','#4E9B88','#222222'))+scale_shape_manual(values=c(16,17,15))+
  labs(x='Coverage',y=NULL,colour=NULL,shape=NULL)+
  guides(colour=guide_legend(ncol=1),shape=guide_legend(ncol=1))+theme(legend.position='bottom',axis.text.y=element_text(size=6),legend.margin=margin(-4,0,0,0),legend.text=element_text(size=5.8))

# c: gene-wise mouse vs GTEx expression-age slopes
cs<-fread(file.path(src,'panel_C_age_direction.tsv'))
cs[,model_family:=factor(model_family,levels=c('ElasticNet','PCA/glmnet'))]
cst<-cs[,{t<-suppressWarnings(cor.test(mouse_age_beta,human_age_beta_per_decade,method='spearman',exact=FALSE));.(rho=unname(t$estimate),p=t$p.value,n=.N)},by=model_family]
cst[,label:=sprintf('ρ = %.3f, P = %s\nn = %s genes',rho,fmt_p(p),format(n,big.mark=','))]
pC<-ggplot(cs,aes(mouse_age_beta,human_age_beta_per_decade))+
  geom_hline(yintercept=0,colour='#BBBBBB',linewidth=.25)+geom_vline(xintercept=0,colour='#BBBBBB',linewidth=.25)+
  geom_point(colour='#6B6B6B',alpha=.35,size=.35,stroke=0)+
  geom_text(data=cst,aes(x=-Inf,y=Inf,label=label),hjust=-.04,vjust=1.15,size=1.8,lineheight=.9)+
  facet_wrap(~model_family,ncol=1)+
  labs(x='Mouse expression–age slope\n(ln1pCPM per week)',y='GTEx expression–age slope (logCPM per decade)')

# d: GTEx standardized predicted age by public bracket
gt<-fread(file.path(src,'panel_E_frozen_scores.tsv'))
gt[,fam:=ifelse(model=='HVG-ElasticNet','elasticnet','pca_glmnet')][,model_label:=factor(mlab[fam],levels=mlab)]
br<-gt[order(age_order),unique(age_bracket)]
gs<-gt[,.(mean=mean(cross_species_mouse_ovarian_clock_score),se=sd(cross_species_mouse_ovarian_clock_score)/sqrt(.N)),by=.(model_label,age_order)]
ac<-rd('tables/age_correlations.tsv')[cohort=='GTEx'];ac[,model_label:=factor(mlab[ifelse(model=='HVG-ElasticNet','elasticnet','pca_glmnet')],levels=mlab)]
ac[,label:=sprintf('Spearman ρ = %.2f\nP = %s (n = %d)',spearman,fmt_p(spearman_p),n)]
set.seed(20261013)
pD<-ggplot(gt,aes(age_order,cross_species_mouse_ovarian_clock_score,colour=model_label,fill=model_label))+
  geom_smooth(method='lm',formula=y~x,se=TRUE,linewidth=.55,alpha=.12)+
  geom_jitter(width=.13,height=0,alpha=.35,size=.5,stroke=0)+
  geom_errorbar(data=gs,aes(y=mean,ymin=mean-1.96*se,ymax=mean+1.96*se),width=.15,linewidth=.4)+
  geom_point(data=gs,aes(y=mean),size=1.4)+
  geom_text(data=ac,aes(x=.6,y=Inf,label=label),inherit.aes=FALSE,hjust=0,vjust=1.1,size=1.75,lineheight=.9)+
  facet_wrap(~model_label,ncol=1)+scale_x_continuous(breaks=seq_along(br),labels=br)+
  scale_y_continuous(expand=expansion(mult=c(.05,.3)))+
  scale_colour_manual(values=setNames(unname(pal),mlab))+scale_fill_manual(values=setNames(unname(pal),mlab))+
  labs(x='GTEx public age bracket',y=expression('Standardized predicted age ('*italic(z)[i]*')'))+
  theme(legend.position='none',axis.text.x=element_text(angle=40,hjust=1))

# e: GSE202601 standardized predicted age vs exact chronological age
ge<-fread(file.path(gsrc,'standardized_score_vs_chronological_age_source_data.tsv'))
ge[,model_label:=factor(mlab[model_family],levels=mlab)][,group:=factor(ifelse(age_group=='young','Young','Reproductively aged'),levels=names(pal_grp))]
gst<-fread(file.path(gsrc,'standardized_score_vs_chronological_age_statistics.tsv'))
gst[,model_label:=factor(mlab[ifelse(grepl('^ElasticNet',model_label),'elasticnet','pca_glmnet')],levels=mlab)]
gst[,label:=sprintf('Spearman ρ = %.2f\nExact P = %.3f\nn = %d',spearman_rho,spearman_exact_p,n_donors)]
pE<-ggplot(ge,aes(age_year,cross_species_mouse_ovarian_clock_score))+
  geom_smooth(method='lm',formula=y~x,se=TRUE,colour='#4A4A4A',fill='#BFC3C7',alpha=.35,linewidth=.45)+
  geom_point(aes(colour=group,shape=group),size=1.6)+
  geom_text_repel(aes(label=sample_id,colour=group),size=1.9,seed=1,box.padding=.25,min.segment.length=.2,segment.colour='#999999',segment.size=.2,show.legend=FALSE)+
  geom_text(data=gst,aes(x=-Inf,y=Inf,label=label),hjust=-.08,vjust=1.1,size=1.9,lineheight=.9,fontface='italic')+
  facet_wrap(~model_label,nrow=1,labeller=as_labeller(elab))+scale_x_continuous(expand=expansion(mult=c(.08,.12)))+
  scale_colour_manual(values=pal_grp)+scale_shape_manual(values=c(16,17))+
  scale_y_continuous(expand=expansion(mult=c(.05,.25)))+
  labs(x='Human chronological age (years)',y=expression('Standardized predicted age ('*italic(z)[i]*')'),colour=NULL,shape=NULL)+
  theme(legend.position='bottom',legend.margin=margin(-4,0,0,0))

# f: cross-model agreement, GTEx + GSE202601, ElasticNet OOD distance
ag<-fread(file.path(gsrc,'GTEx_GSE202601_cross_model_agreement_source_data.tsv'))
ag[,cohort:=factor(cohort,levels=c('GTEx','GSE202601'))]
pt<-cor.test(ag$elasticnet_score,ag$pca_glmnet_score)
ags<-fread(file.path(gsrc,'GTEx_GSE202601_cross_model_agreement_statistics.tsv'))[cohort%in%c('GTEx','GSE202601')]
flab<-paste0(sprintf('Pooled r = %.2f, P = %s\n',pt$estimate,fmt_p(pt$p.value)),
  paste(sprintf('%s r = %.2f, P = %s',ags$cohort,ags$pearson_r,fmt_p(ags$pearson_p)),collapse='\n'))
pF<-ggplot(ag,aes(elasticnet_score,pca_glmnet_score,colour=elasticnet_ood_distance,shape=cohort,size=cohort))+
  geom_hline(yintercept=0,colour='#CCCCCC',linewidth=.25)+geom_vline(xintercept=0,colour='#CCCCCC',linewidth=.25)+
  geom_point(alpha=.85)+
  annotate('text',x=-Inf,y=Inf,label=flab,hjust=-.03,vjust=1.1,size=1.7,lineheight=.9,fontface='italic')+
  scale_y_continuous(expand=expansion(mult=c(.05,.32)))+
  scale_colour_viridis_c(option='C',name='ElasticNet\nOOD distance')+scale_shape_manual(values=c(GTEx=16,GSE202601=17),name='Cohort')+
  scale_size_manual(values=c(GTEx=.8,GSE202601=1.8),guide='none')+
  guides(shape=guide_legend(override.aes=list(size=1.6,colour='#555555'),order=1,title.position='top',ncol=1),
         colour=guide_colourbar(barwidth=unit(14,'mm'),barheight=unit(2,'mm'),title.position='top'))+
  labs(x='ElasticNet standardized\npredicted age (zᵢ)',y='PCA/glmnet standardized\npredicted age (zᵢ)')+
  theme(legend.position='bottom',legend.box='horizontal',legend.margin=margin(-4,0,0,0),legend.title=element_text(size=5.8),legend.text=element_text(size=5.6))

fig<-wrap_plots(A=pA+labs(tag='a'),B=pB+labs(tag='b'),C=pC+labs(tag='c'),D=pD+labs(tag='d'),E=pE+labs(tag='e'),F=pF+labs(tag='f'),
  design='AAAAAABBBBCCC\nDDDEEEEEEFFFF',heights=c(1.05,1.2))+
  plot_annotation(title='Supplementary Fig. 7',theme=theme(plot.title=element_text(face='bold',size=9)))

stem<-file.path(out,'figures/Suppl_Figure_7_human_generalization_reference_layout');w<-200/25.4;h<-145/25.4
grDevices::cairo_pdf(paste0(stem,'.pdf'),width=w,height=h,family='sans');print(fig);invisible(dev.off())
grDevices::svg(paste0(stem,'.svg'),width=w,height=h,family='sans');print(fig);invisible(dev.off())
ragg::agg_png(paste0(stem,'.png'),width=w,height=h,units='in',res=600);print(fig);invisible(dev.off())
ragg::agg_tiff(paste0(stem,'.tiff'),width=w,height=h,units='in',res=600,compression='lzw');print(fig);invisible(dev.off())

# panel d alone, wider x axis
dstem<-file.path(out,'figures/Suppl_Figure_7d_GTEx_age_bracket');dw<-90/25.4;dh<-95/25.4
pDs<-pD;pDs$layers[[2]]$aes_params$size<-1.6;pDs$layers[[2]]$aes_params$alpha<-.45  # donor dots
pDs$layers[[3]]$aes_params$linewidth<-.6;pDs$layers[[4]]$aes_params$size<-3            # bracket means
grDevices::cairo_pdf(paste0(dstem,'.pdf'),width=dw,height=dh,family='sans');print(pDs);invisible(dev.off())
grDevices::svg(paste0(dstem,'.svg'),width=dw,height=dh,family='sans');print(pDs);invisible(dev.off())
ragg::agg_png(paste0(dstem,'.png'),width=dw,height=dh,units='in',res=600);print(pDs);invisible(dev.off())
ragg::agg_tiff(paste0(dstem,'.tiff'),width=dw,height=dh,units='in',res=600,compression='lzw');print(pDs);invisible(dev.off())

sd7<-file.path(src,'Suppl_Figure_7_reference_layout');dir.create(sd7,showWarnings=FALSE)
fwrite(cov,file.path(sd7,'panel_b_coverage.tsv'),sep='\t');fwrite(cst[,!'label'],file.path(sd7,'panel_c_slope_spearman.tsv'),sep='\t')
fwrite(gs,file.path(sd7,'panel_d_bracket_means.tsv'),sep='\t');fwrite(ag,file.path(sd7,'panel_f_agreement.tsv'),sep='\t')
fwrite(data.table(scope=c('pooled',ags$cohort),n=c(nrow(ag),ags$n_samples),pearson_r=c(unname(pt$estimate),ags$pearson_r),pearson_p=c(pt$p.value,ags$pearson_p)),file.path(sd7,'panel_f_agreement_statistics.tsv'),sep='\t')
message('Supplementary Fig. 7 (reference layout) written: ',stem)
