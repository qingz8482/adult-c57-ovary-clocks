#!/usr/bin/env Rscript
# Variant of diagnostics/06_clock_concordance_delta_age: developmental_OOD samples removed,
# both axes fixed to -40..40 weeks, narrower x panel, axes named by clock.
# Reads the sample-level table written by 07_complete_reference_outputs.R; no refitting.
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
a<-commandArgs(FALSE);here<-dirname(normalizePath(sub('^--file=','',a[grepl('^--file=',a)][1])))
out<-dirname(here)
d<-fread(file.path(out,'diagnostics/model_concordance_sample_level.tsv'))
d[plot_class=='stress',plot_class:='stress_model']
d<-d[plot_class!='developmental_OOD']
stopifnot(uniqueN(d$model_x)==1,uniqueN(d$model_y)==1)
class_pal<-c(reference_control='#2A9D8F',intervention='#0072B2',genetic_model='#D55E00',diet='#E69F00',stress_model='#CC79A7',microbiome='#8C564B')
# Summary stats use every retained sample, including any outside the -40..40 view.
s<-d[,.(n=.N,pearson_r=cor(x,y),mean_abs_difference=mean(abs(difference))),by=pair]
s[,label:=sprintf('n = %d; r = %.2f\nmean |difference| = %.1f weeks',n,pearson_r,mean_abs_difference)]
lim<-c(-40,40)
g<-ggplot(d,aes(x,y,colour=plot_class,shape=modality))+
 geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey45')+
 geom_hline(yintercept=0,linetype='dotted',colour='grey60')+geom_vline(xintercept=0,linetype='dotted',colour='grey60')+
 geom_point(size=2,alpha=.8)+
 facet_wrap(~pair,nrow=1,labeller=label_wrap_gen(25))+
 geom_text(data=s,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.05,vjust=1.1,size=3)+
 scale_colour_manual(values=class_pal)+
 scale_x_continuous(breaks=seq(lim[1],lim[2],20))+scale_y_continuous(breaks=seq(lim[1],lim[2],20))+
 coord_cartesian(xlim=lim,ylim=lim)+
 labs(title='Delta-age concordance:\nthe two locked clocks',
      x=paste0(d$model_x[1],' delta age (weeks)'),y=paste0(d$model_y[1],' delta age (weeks)'),
      colour='Projection class',shape='Modality')+
 theme_bw(base_size=12)+
 theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold'),strip.text=element_text(face='bold'),
       legend.position='bottom',legend.box='vertical',plot.margin=margin(12,18,12,12),
       panel.widths=unit(6,'in'),panel.heights=unit(5,'in'))+
 guides(colour=guide_legend(ncol=2,order=1),shape=guide_legend(ncol=1,order=2))
pdf(NULL);gt<-ggplotGrob(g);w<-grid::convertWidth(sum(gt$widths),'in',TRUE)+.2;h<-grid::convertHeight(sum(gt$heights),'in',TRUE)+.2;invisible(dev.off())
stem<-file.path(out,'diagnostics/06_clock_concordance_delta_age_no_developmental_OOD')
ggsave(paste0(stem,'.pdf'),g,width=w,height=h,device=cairo_pdf)
ggsave(paste0(stem,'.png'),g,width=w,height=h,dpi=300,bg='white')
fwrite(s,paste0(stem,'_summary.tsv'),sep='\t')
print(s)
