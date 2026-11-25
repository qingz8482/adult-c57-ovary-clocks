#!/usr/bin/env Rscript
# Figure contract: explain selection, conditional PC attribution, exact sample
# and gene decompositions, and directional enrichment of the locked clock.
# Quantitative grid: three model panels, a full-width contribution panel, and
# gene/pathway panels. Export 183-mm PDF/SVG plus 600-dpi PNG/TIFF; all sources TSV.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(scales);library(grid)})
arg<-gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),fixed=TRUE)
out<-normalizePath(file.path(dirname(arg),'..'));dd<-file.path(out,'derived_data');panel_dir<-file.path(out,'panels')
obj<-readRDS(file.path(dd,'plot_objects.rds'));list2env(obj,envir=environment())
stopifnot(requireNamespace('ragg',quietly=TRUE))
theme_pub<-theme_classic(base_size=6.5,base_family='Arial')+theme(
 axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3),
 plot.margin=margin(3,4,3,4),legend.title=element_text(size=5.5),legend.text=element_text(size=5),
 plot.title=element_text(size=6.5,face='bold'),plot.subtitle=element_text(size=5.5),
 plot.tag=element_text(size=8,face='bold'),panel.grid=element_blank())
theme_set(theme_pub)
dep_lambda<-lambda
break_ix<-unique(round(seq(1,nrow(cv_dt),length.out=5)))
p_a<-ggplot(cv_dt,aes(log_lambda,mean_mae))+geom_ribbon(aes(ymin=lower,ymax=upper),fill='grey85')+
 geom_line(linewidth=.45)+geom_vline(xintercept=log(dep_lambda),colour='#B2182B',linewidth=.4)+
 geom_vline(xintercept=log(cv_dt$lambda[which.min(cv_dt$mean_mae)]),linetype=2,colour='grey45',linewidth=.3)+
 scale_x_continuous(sec.axis=dup_axis(name='Nonzero PCs in deployment fit',breaks=cv_dt$log_lambda[break_ix],labels=cv_dt$n_nonzero[break_ix]))+
 labs(title='Grouped model selection',subtitle=sprintf('3,000 genes; 10 PC inputs; alpha = %.1f',alpha),
 x='log(lambda)',y='Mean study MAE (weeks)',tag='a')
mark<-path_dt[selected==TRUE, .SD[which.min(abs(lambda-dep_lambda))],by=pc]
p_b<-ggplot(path_dt,aes(log_lambda,coefficient,group=pc))+
 geom_hline(yintercept=0,colour='grey70',linewidth=.25)+
 geom_line(data=path_dt[selected==FALSE],colour='grey78',linewidth=.3)+
 geom_line(data=path_dt[selected==TRUE],aes(colour=coefficient),linewidth=.5)+
 geom_vline(xintercept=log(dep_lambda),linetype=2,linewidth=.3)+
 geom_text(data=mark,aes(label=pc),nudge_x=.7,hjust=0,size=1.8,show.legend=FALSE)+
 scale_colour_gradient2(low='#2166AC',mid='white',high='#B2182B',guide='none')+
 labs(title='Frozen coefficient paths',subtitle=sprintf('Deployed lambda = %.2f; no calibration',dep_lambda),x='log(lambda)',y='PC coefficient (weeks / score)',tag='b')
make_c<-function(d,subtitle){ggplot(d,aes(omitted_pc,delta_mae))+geom_hline(yintercept=0,colour='grey60',linewidth=.3)+
 geom_errorbar(aes(ymin=ci_low,ymax=ci_high),width=.15,linewidth=.4,colour='grey30')+
 geom_point(aes(colour=delta_mae),size=1.8)+scale_colour_gradient2(low='#2166AC',mid='grey50',high='#B2182B',guide='none')+
 labs(title='Conditional PC ablation',subtitle=subtitle,x='Omitted active PC',y='Change in MAE (weeks)',tag='c')}
p_c<-make_c(lopo_summary,'Paired sample bootstrap; fixed PCA space')
p_c_study<-make_c(study_summary,'Equal study weighting; study bootstrap')

d_plot[,age_bin:=factor(age_bin,levels=sort(unique(age_bin)))]
lim_d<-max(abs(d_plot$contribution))
p_d_main<-ggplot(d_plot,aes(age_bin,factor(pc,levels=rev(pcs)),fill=contribution))+geom_tile()+
 scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-lim_d,lim_d),name='Contribution\n(weeks)')+
 guides(fill=guide_colourbar(barheight=unit(14,'mm'),barwidth=unit(2.5,'mm')))+
 scale_x_discrete(labels=NULL,breaks=NULL)+labs(x='136 training profiles in 20 equal-count bins, ordered by chronological age',y=NULL)+
 theme(axis.ticks=element_blank(),axis.line=element_blank())
ann<-unique(d_plot[,.(age_bin,age_week,predicted_age)])
ann<-melt(ann,id.vars='age_bin',variable.name='annotation',value.name='value');ann[,scaled_value:=as.numeric(scale(value)),by=annotation]
fwrite(ann,file.path(dd,'panel_d_age_annotation_tracks.tsv'),sep='\t')
p_d_ann<-ggplot(ann,aes(age_bin,factor(annotation,levels=c('predicted_age','age_week')),fill=scaled_value))+geom_tile()+
 scale_fill_gradient(low='grey95',high='grey25',name='Relative age\n(within-track z)',limits=c(-2,2),oob=scales::squish,breaks=c(-2,0,2))+
 guides(fill=guide_colourbar(barheight=unit(8,'mm'),barwidth=unit(2.5,'mm')))+
 scale_x_discrete(labels=NULL,breaks=NULL)+scale_y_discrete(labels=c('Predicted age','Chronological age'))+
 labs(x=NULL,y=NULL)+theme(axis.ticks=element_blank(),axis.line=element_blank())
p_d<-wrap_elements(full=p_d_ann/p_d_main+plot_layout(heights=c(.6,1)))+labs(tag='d')

e_dt[,gene:=factor(gene,levels=rev(genes_show))];e_dt[,component:=factor(component,levels=c(pcs,'Integrated'))]
lim_e<-max(abs(e_dt$weight))
p_e<-ggplot(e_dt,aes(component,gene,fill=weight))+geom_tile(colour='white',linewidth=.12)+
 scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-lim_e,lim_e),name='Gene weight\n(weeks / ln1p CPM)')+
 guides(fill=guide_colourbar(barheight=unit(18,'mm'),barwidth=unit(3,'mm')))+
 labs(title='Back-projected gene weights',x=NULL,y=NULL,tag='e')+
 theme(axis.text.y=element_text(size=5.6),axis.text.x=element_text(angle=45,hjust=1),axis.line=element_blank(),axis.ticks=element_blank())

slots<-fread(file.path(dd,'gsea_panel_slots.tsv'))
pretty_term<-function(x,width=24){
 x<-gsub('^(HALLMARK|GOBP|REACTOME|WP|CUI)_','',x);x<-tools::toTitleCase(tolower(gsub('_',' ',x)))
 for(pair in list(c('Nk','NK'),c('Il','IL'),c('Tca','TCA'),c('Dna','DNA'),c('E2f','E2F'),c('Cd4','CD4'),c('Cd8','CD8')))x<-gsub(paste0('\\b',pair[1]),pair[2],x)
 paste(strwrap(x,width=width),collapse='\n')
}
make_f<-function(panel_name,compact=FALSE,return_plots=FALSE){
 s<-slots[panel==panel_name][order(display_order)]
 curves<-fread(file.path(dd,paste0('panel_',panel_name,'_enrichment_curves.tsv')))
 stats<-fread(file.path(dd,paste0('panel_',panel_name,'_selected_terms.tsv')))
 pp<-lapply(seq_len(nrow(s)),function(i){
  if(s$status[i]!='selected')return(ggplot()+annotate('text',x=0,y=0,label=paste0('No ',s$category[i],' terms\npass FDR < 0.05'),size=if(compact)1.9 else 2.8)+
   xlim(-1,1)+ylim(-1,1)+labs(title=s$category[i])+theme_void(base_size=6.5)+theme(plot.title=element_text(size=6,face='bold')))
  term<-s$gene_set_name[i];d<-curves[gene_set_name==term];st<-stats[gene_set_name==term]
  colour<-if(st$NES>0)'#B2182B' else '#2166AC'
  xmin<-if(st$NES>0)1 else st$peak_rank;xmax<-if(st$NES>0)st$peak_rank else nrow(d)
  ggplot(d,aes(rank,running_es))+annotate('rect',xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf,fill=colour,alpha=.07)+
   geom_hline(yintercept=0,colour='grey65',linewidth=.25)+geom_line(colour=colour,linewidth=.4)+
   geom_rug(data=d[hit==TRUE],sides='b',linewidth=.13,colour='grey25')+
   geom_vline(xintercept=st$peak_rank,linetype=3,colour=colour,linewidth=.3)+
   scale_x_continuous(breaks=c(0,1500,3000),labels=c('0','1,500','3,000'))+
   labs(title=pretty_term(term,if(compact)22 else 32),subtitle=sprintf('NES %.2f; FDR %.3f\nn=%d; LE=%d',st$NES,st$fdr_within_pc,st$gene_set_size_in_universe,st$n_leading_edge),x='Gene rank',y='ES')+
   theme(plot.title=element_text(size=if(compact)5.3 else 7,face='bold'),plot.subtitle=element_text(size=if(compact)4.8 else 6),
    axis.text=element_text(size=if(compact)4.8 else 6),axis.title=element_text(size=if(compact)5 else 6),plot.margin=margin(4,4,4,4))
 })
 if(return_plots)return(pp)
 wrap_elements(full=wrap_plots(pp,ncol=2))+labs(tag=panel_name)
}
p_f<-make_f('f',TRUE)
# Keep the panel-f reference geometry. All further curve grobs use its exact
# physical inner-panel dimensions, regardless of title length or canvas size.
reference_device<-tempfile(fileext='.pdf')
cairo_pdf(reference_device,width=78/25.4,height=108/25.4)
print(p_f);grid.force()
reference_names<-unique(grid.ls(viewports=TRUE,grobs=FALSE,print=FALSE)$name)
reference_names<-reference_names[grepl('^panel-[0-9]+\\.',reference_names)]
reference_sizes<-rbindlist(lapply(reference_names,function(n){
 seekViewport(n);z<-data.table(width_mm=convertWidth(unit(1,'npc'),'mm',TRUE),
 height_mm=convertHeight(unit(1,'npc'),'mm',TRUE));upViewport(0);z
}))
dev.off();unlink(reference_device)
stopifnot(nrow(reference_sizes)==4L,diff(range(reference_sizes$width_mm))<1e-8,
 diff(range(reference_sizes$height_mm))<1e-8)
curve_width<-reference_sizes$width_mm[1];curve_height<-reference_sizes$height_mm[1]
fwrite(reference_sizes,file.path(dd,'panel_f_reference_curve_dimensions.tsv'),sep='\t')

fixed_curve_grob<-function(p,panel_name,i){
 gt<-ggplotGrob(p);ix<-which(gt$layout$name=='panel');stopifnot(length(ix)==1)
 gt$widths[gt$layout$l[ix]:gt$layout$r[ix]]<-unit(curve_width,'mm')
 gt$heights[gt$layout$t[ix]:gt$layout$b[ix]]<-unit(curve_height,'mm')
 gt$layout$name[ix]<-paste0('curve_panel_',panel_name,'_',i)
 gt
}
make_fixed_f<-function(panel_name){
 pp<-make_f(panel_name,TRUE,return_plots=TRUE)
 is_f4<-panel_name=='f4';w<-if(is_f4)234 else 78;h<-if(is_f4)185 else 54*ceiling(length(pp)/2)
 ss<-slots[panel==panel_name][order(display_order)]
 children<-list()
 if(is_f4)for(g in seq_len(9)){
   col<-(g-1)%%3;row<-(g-1)%/%3
   children[[length(children)+1]]<-textGrob(unique(ss$category)[g],
     x=unit(col*78+39,'mm'),y=unit(h-5-row*60-3,'mm'),
     gp=gpar(fontfamily='Arial',fontface='bold',fontsize=6.5))
 }
 for(i in seq_along(pp)){
   if(is_f4){g<-(i-1)%/%2;col<-g%%3;row<-g%/%3
     x<-col*78+((i-1)%%2+.5)*39;y<-h-5-row*60-6-27
   }else{x<-((i-1)%%2+.5)*39;y<-h-((i-1)%/%2+.5)*54}
   children[[length(children)+1]]<-grobTree(fixed_curve_grob(pp[[i]],panel_name,i),
     vp=viewport(x=unit(x,'mm'),y=unit(y,'mm'),width=unit(39,'mm'),height=unit(54,'mm'),
                 name=paste0('curve_container_',panel_name,'_',i)))
 }
 # Draw the panel tag last so an opaque curve background cannot cover it.
 children[[length(children)+1]]<-textGrob(panel_name,x=unit(1,'mm'),y=unit(h-1,'mm'),just=c('left','top'),
   gp=gpar(fontfamily='Arial',fontface='bold',fontsize=8))
 do.call(grobTree,children)
}
# Resolve Arial metrics on the same Cairo backend used for the final PDFs.
cairo_pdf(reference_device,width=78/25.4,height=108/25.4)
p_f2<-make_fixed_f('f2');p_f3<-make_fixed_f('f3');p_f4<-make_fixed_f('f4');p_f5<-make_fixed_f('f5')
dev.off();unlink(reference_device)
curve_measurements<-list()
draw_plot<-function(p){if(inherits(p,'grob')){grid.newpage();grid.draw(p)}else print(p)}
save_all<-function(p,stem,w,h){
 cairo_pdf(paste0(stem,'.pdf'),width=w/25.4,height=h/25.4,bg='white');draw_plot(p)
 if(basename(stem)%in%c('panel_f','panel_f2','panel_f3','panel_f4','panel_f5')){
   grid.force();vn<-unique(grid.ls(viewports=TRUE,grobs=FALSE,print=FALSE)$name)
   vn<-vn[grepl(if(basename(stem)=='panel_f')'^panel-[0-9]+\\.' else '^curve_panel_',vn)]
   for(n in vn){seekViewport(n);curve_measurements[[length(curve_measurements)+1]]<<-data.table(
     panel=basename(stem),viewport=n,width_mm=convertWidth(unit(1,'npc'),'mm',TRUE),
     height_mm=convertHeight(unit(1,'npc'),'mm',TRUE),canvas_width_mm=w,canvas_height_mm=h);upViewport(0)}
 }
 dev.off()
 ggsave(paste0(stem,'.svg'),p,width=w,height=h,units='mm',device=grDevices::svg,bg='white')
 ragg::agg_png(paste0(stem,'.png'),width=w,height=h,units='mm',res=600,background='white');draw_plot(p);dev.off()
 ragg::agg_tiff(paste0(stem,'.tiff'),width=w,height=h,units='mm',res=600,compression='lzw',background='white');draw_plot(p);dev.off()
}
panels<-list(a=p_a,b=p_b,c=p_c,d=p_d,e=p_e,f=p_f,f2=p_f2,f3=p_f3,f4=p_f4,f5=p_f5,c_study_weighted=p_c_study)
widths<-c(60,60,60,183,105,78,78,78,234,78,75);heights<-c(55,55,55,46,108,108,108,108,185,54,60)
gsea_only<-'--gsea-only'%in%commandArgs(trailingOnly=TRUE)
for(i in seq_along(panels)){
 if(gsea_only&&!names(panels)[i]%in%c('f','f2','f3','f4','f5'))next
 save_all(panels[[i]],file.path(panel_dir,paste0('panel_',names(panels)[i])),widths[i],heights[i])
}
fwrite(rbindlist(curve_measurements),file.path(dd,'gsea_curve_dimensions.tsv'),sep='\t')
if(!gsea_only){
top<-p_a+p_b+p_c+plot_layout(widths=c(1,1,1))
bottom<-p_e+p_f+plot_layout(widths=c(1.35,1))
final_plot<-top/p_d/bottom+plot_layout(heights=c(1,.85,2))
save_all(final_plot,file.path(out,'final/Supplementary_Figure_2_PCA_clock_model_interpretation'),183,220)
}
message('Exported GSEA panels f/f2/f3/f4/f5 with identical curve-area dimensions: ',curve_width,' x ',curve_height,' mm.')
