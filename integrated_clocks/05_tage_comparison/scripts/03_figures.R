source(file.path(dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),fixed=TRUE))),'common.R'))
pred<-read_tsv(file.path(out,'tables/sample_contrast_predictions.tsv'))
pairs<-read_tsv(file.path(out,'tables/clock_pairwise_agreement.tsv'))[analysis_set=='primary_adult_prespecified']
agreement<-read_tsv(file.path(out,'tables/tage_variant_agreement.tsv'))
theme_set(theme_bw(base_size=10,base_family='sans')+theme(panel.grid.minor=element_blank(),strip.background=element_rect(fill='grey95'),legend.position='bottom'))
export<-function(p,stem,w,h){
  ggsave(file.path(out,'figures',paste0(stem,'.pdf')),p,device=cairo_pdf,width=w,height=h,limitsize=FALSE)
  ggsave(file.path(out,'figures',paste0(stem,'.png')),p,device=ragg::agg_png,width=w,height=h,dpi=300,limitsize=FALSE,bg='white')
}
stars<-function(p)fcase(is.na(p),'',p<1e-4,'****',p<1e-3,'***',p<.01,'**',p<.05,'*',default='ns')
pairs[,comparison:=paste(model_1,'vs',model_2)]
long<-melt(pairs,id.vars=c('comparison','n','pearson_p','spearman_p'),measure.vars=c('pearson','spearman','direction_agreement'),variable.name='metric',value.name='agreement')
long[,p_value:=fcase(metric=='pearson',pearson_p,metric=='spearman',spearman_p,default=NA_real_)]
long[,label:=paste0(sprintf('%.2f',agreement),ifelse(metric=='direction_agreement','',paste0('\n',stars(p_value))))]
long[,comparison:=factor(comparison,levels=rev(pairs$comparison))]
metric_levels<-c('pearson','spearman','direction_agreement');metric_labels<-c('Pearson r','Spearman rho','Direction agreement')
long[,metric:=factor(metric,levels=metric_levels,labels=metric_labels)]
ag<-copy(agreement);ag[,comparison:=sprintf('%s (n=%d)',level,n)]
al<-melt(ag,id.vars='comparison',measure.vars=metric_levels,variable.name='metric',value.name='agreement')
al[,`:=`(metric=factor(metric,levels=metric_levels,labels=metric_labels),comparison=factor(comparison,levels=rev(ag$comparison)),label=sprintf('%.2f',agreement))]
# Use the full observed range; never clip lower correlations to the historical 0.5 floor.
lo<-min(c(long$agreement,al$agreement),na.rm=TRUE)
lim<-c(if(lo<0)-1 else 0,1)
heat<-function(d,title,subtitle){
  p<-ggplot(d,aes(metric,comparison,fill=agreement))+geom_tile(color='white',linewidth=.8)+geom_text(aes(label=label),size=3.2,fontface='bold',lineheight=.85)+labs(title=title,subtitle=subtitle,x=NULL,y=NULL,fill='Agreement')+theme(panel.grid=element_blank(),legend.position='right',axis.text.x=element_text(face='bold'),axis.text.y=element_text(size=9),plot.title=element_text(size=12),plot.subtitle=element_text(size=9))
  if(lo<0)p+scale_fill_gradient2(low='#BB6A57',mid='white',high='#087E8B',midpoint=0,limits=lim) else p+scale_fill_gradient(low='#F1F4F5',high='#087E8B',limits=lim)
}
p1<-heat(long,'a  Four-clock agreement on matched-control effects',sprintf('%d adult prespecified non-placement contrasts; nominal two-sided correlation tests',unique(pairs$n)))
p2<-heat(al,'b  Agreement between the two tAge variants','All prespecified contrasts; samples can occur in multiple contrasts')
p<-(p1/p2)+plot_layout(heights=c(2.2,1))+plot_annotation(title='Clock agreement for matched-control ovarian stress-test effects',subtitle='HVG-ElasticNet and VersionB-PCA/glmnet versus tAge scaled-diff and YuGene-diff',caption='*P<0.05, **P<0.01, ***P<0.001, ****P<0.0001; ns, not significant. Contrasts may share controls/studies.\ntAge shifts are relative to declared controls, not conventional biological age acceleration.',theme=theme(plot.title=element_text(size=14),plot.subtitle=element_text(size=10),plot.caption=element_text(size=8,hjust=0)))
export(p,'10_tage_clock_agreement_readable',12,8.5)
write_tsv(long,'tables/figure10_panel_a_source.tsv');write_tsv(al,'tables/figure10_panel_b_source.tsv')

pal<-c(diet='#F28E2B',genetic_model='#E15759',intervention='#4E79A7',microbiome='#59A14F',reference_control='#7F7F7F',stress_model='#9C755F',developmental_OOD='#B07AA1')
d<-pred[model%chin%models[3:4] & analysis_family=='prespecified']
stopifnot(all(d$projection_class%chin%names(pal)))
d[,projection_class:=factor(projection_class,levels=names(pal))]
d[,model_variant:=factor(model_variant,levels=c('scaleddiff','yugenediff'))]
set.seed(20261012)
p<-ggplot(d,aes(projection_class,relative_shift_weeks,color=projection_class))+geom_hline(yintercept=0,linetype='dashed',color='grey55')+geom_boxplot(outlier.shape=NA,linewidth=.5,width=.7)+geom_point(position=position_jitter(width=.13,height=0,seed=20261012),size=2.2,alpha=.65)+facet_wrap(~model_variant,nrow=1)+scale_color_manual(values=pal,drop=FALSE)+labs(title='Reference-centred tAge shift distributions',subtitle='Positive values indicate an older transcriptomic state relative to the matched reference; negative values indicate a younger state',x='Projection class',y='tAge shift relative to matched reference (weeks)',color='Projection class',caption=sprintf('%d sample-within-contrast observations per variant (%d unique profiles); reused controls/targets are shown in each declared contrast.\nBoxes: median and IQR; whiskers: 1.5 IQR. tAge shifts are not conventional biological age acceleration.',nrow(d)/2,uniqueN(d$sample_key)))+theme(axis.text.x=element_text(angle=35,hjust=1),plot.caption=element_text(size=8,hjust=0),legend.text=element_text(size=8),plot.subtitle=element_text(size=9))
export(p,'03_tage_shift_relative_to_reference_distributions',11.5,7)
# Plot object kept for the side-by-side figure with the locked-clock raw delta-age panel.
# The session theme (theme_set) is not stored in the plot object, so keep it alongside.
saveRDS(list(plot=p,theme=theme_get()),file.path(out,'figures','03_tage_shift_relative_to_reference_distributions.rds'))
write_tsv(d,'tables/figure03_source.tsv')
write_tsv(d[,.(n_observations=.N,n_unique_profiles=uniqueN(sample_key),median=median(relative_shift_weeks),q1=quantile(relative_shift_weeks,.25),q3=quantile(relative_shift_weeks,.75)),by=.(model_variant,projection_class)],'tables/figure03_distribution_summary.tsv')

en<-read_tsv(file.path(out,'tables/ovary_elasticnet_features.tsv'))[standardized_coefficient!=0]
pc<-read_tsv(file.path(out,'tables/ovary_pca_backprojected_features.tsv'))
tt<-lapply(c('scaleddiff','yugenediff'),function(v)read_tsv(file.path(out,'tables',paste0('tage_',v,'_features.tsv'))))
universe<-read_tsv(file.path(out,'tables/current_ovary_gene_universe.tsv'))$gene
shared_universe<-Reduce(intersect,list(universe,tt[[1]]$gene,tt[[2]]$gene))
write_tsv(data.table(gene=shared_universe),'tables/feature_overlap_gene_universe.tsv')
details<-list();stats<-list();plots<-list()
for(i in 1:4){
  isen<-i<=2;j<-if(isen)i else i-2
  ov<-if(isen)en else pc;tg<-tt[[j]]
  z<-merge(ov[gene%chin%shared_universe,.(gene,ovary_weight=standardized_coefficient,ovary_raw_scale_weight=coefficient)],tg[,.(gene,tage_weight=coefficient,tage_entrez=entrez)],by='gene')
  cmp<-c('HVG-ElasticNet vs tAge scaled-diff','HVG-ElasticNet vs tAge YuGene','VersionB-PCA all shared vs tAge scaled-diff','VersionB-PCA all shared vs tAge YuGene')[i]
  z[,comparison:=cmp];s<-safe_stats(z$tage_weight,z$ovary_weight);s[,comparison:=cmp]
  stats[[i]]<-s;details[[i]]<-z
  label<-sprintf('n=%s\nPearson=%.2f\nSpearman=%.2f',format(nrow(z),big.mark=','),s$pearson,s$spearman)
  color<-c('#DE9A65','#65AE98')[j]
  plots[[i]]<-ggplot(z,aes(tage_weight,ovary_weight))+geom_hline(yintercept=0,color='#D4DAE2',linewidth=.35)+geom_vline(xintercept=0,color='#D4DAE2',linewidth=.35)+geom_point(color=color,size=if(isen)1.4 else .65,alpha=if(isen).85 else .3)+geom_smooth(method='lm',formula=y~x,se=FALSE,color='#30343B',linewidth=.5)+annotate('text',x=-Inf,y=Inf,label=label,hjust=-.04,vjust=1.1,size=2.5)+labs(title=cmp,x='tAge coefficient',y=if(isen)'Ovary standardized coefficient' else 'Ovary PCA backprojected score')+theme_classic(base_size=9)+theme(plot.title=element_text(size=10),plot.margin=margin(8,8,8,8))
}
stat<-rbindlist(stats);stat[,`:=`(pearson_bh_q=p.adjust(pearson_p,'BH'),spearman_bh_q=p.adjust(spearman_p,'BH'))]
write_tsv(stat,'tables/feature_direction_statistics_all_pca.tsv')
write_tsv(rbindlist(details),'tables/feature_coefficient_concordance_all_pca_source.tsv')
p<-wrap_plots(plots,ncol=2)+plot_annotation(tag_levels='a',caption='Ovary weights use standardized input genes; tAge coefficients retain their native model scale.\nPCA includes all shared input genes. Weight magnitudes are not directly comparable across models.',theme=theme(plot.tag=element_text(face='bold',size=12),plot.caption=element_text(hjust=0,size=8)))
export(p,'feature_coefficient_concordance_all_pca',9,7.5)
write_tsv(data.table(figure=c('10_tage_clock_agreement_readable','03_tage_shift_relative_to_reference_distributions','feature_coefficient_concordance_all_pca'),width_inches=c(12,11.5,9),height_inches=c(8.5,7,7.5),backend='R ggplot2/patchwork',formats='PDF;PNG',png_dpi=300),'provenance/figure_exports.tsv')
message('PASS: three PDF/PNG figure pairs and source/statistics tables exported.')
