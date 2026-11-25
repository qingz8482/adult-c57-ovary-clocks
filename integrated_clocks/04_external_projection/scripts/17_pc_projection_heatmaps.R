#!/usr/bin/env Rscript
# Per-dataset frozen-PCA PC-contribution heatmaps (mean contribution per projection group and PC),
# drawn like 07_complete_reference_outputs.R heat_plot(). Tile width and height are fixed to the
# tiles of GSE235021/plots/09_gse235021_pca_pc_contribution_heatmap_bulk_only (11 in wide page); rows follow the
# group order of each dataset's delta-age plot (provenance/complete_group_order_audit.tsv).
# Figures: one per dataset, plus two combined heatmaps (NMN + Cd38: CRA017492 / PRJNA1002222;
# VCD + Mtor: GSE249151 / GSE98497; rice diet + VCD: GSE235021 bulk / GSE249151) with a strip per dataset, and the GSE235021 bulk RNA-seq groups.
# Every figure uses one common 'Weeks' colour scale (fixed limits over all datasets drawn here), so
# figures can be combined with a single colour bar.
# Reads existing tables only; no model fitting or projection.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(grid)})
a<-commandArgs(FALSE);here<-dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE)))
out<-dirname(here)
outdir<-file.path(out,'diagnostics','pc_projection_heatmap');dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
datasets<-c('CRA017492','PRJNA1002222','GSE249151','GSE98497')
figures<-c(
  list(GSE235021_bulk_pca_pc_contribution_heatmap=list(datasets='GSE235021',title='GSE235021 bulk RNA-seq: mean frozen PC contributions')),
  setNames(lapply(datasets,function(ds)list(datasets=ds,title=paste0(ds,': mean frozen PC contributions'))),paste0(datasets,'_pca_pc_contribution_heatmap')),
  list(
    NMN_Cd38_pca_pc_contribution_heatmap=list(datasets=c('CRA017492','PRJNA1002222'),title='NMN (CRA017492) and Cd38 KO (PRJNA1002222): mean frozen PC contributions'),
    VCD_Mtor_pca_pc_contribution_heatmap=list(datasets=c('GSE249151','GSE98497'),title='VCD (GSE249151) and Mtor cKO (GSE98497): mean frozen PC contributions'),
    GSE235021_GSE249151_pca_pc_contribution_heatmap=list(datasets=c('GSE235021','GSE249151'),title='Rice diet (GSE235021 bulk) and VCD (GSE249151): mean frozen PC contributions')
  )
)

# Same theme and heatmap layer as 07_complete_reference_outputs.R
theme_set(theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold'),strip.text=element_text(face='bold'),legend.position='bottom',axis.text.x=element_text(angle=35,hjust=1),plot.margin=margin(12,18,12,12)))
heat_gg<-function(hh,title,lim=NULL){
  g<-ggplot(hh,aes(PC,group_label,fill=contribution))+geom_tile(colour='white')+scale_fill_gradient2(low='#3B6EA8',mid='white',high='#C44E52',midpoint=0,limits=lim)+labs(title=title,x='Frozen selected PC',y='Projection group',fill='Weeks')+theme(axis.text.x=element_text(angle=0),panel.grid=element_blank())
  # Several datasets: one block per dataset (as heat_plot does for several cohorts), shared colour scale
  if('dataset_folder'%in%names(hh)&&uniqueN(hh$dataset_folder)>1)g<-g+facet_grid(dataset_folder~.,scales='free_y',space='free_y')+theme(strip.text.y=element_text(angle=-90,size=7.5),panel.spacing.y=unit(0.25,'lines'))
  g
}
page_w<-11

cc<-fread(file.path(out,'pca_projection_pc_contributions.tsv'))
ann<-unique(fread(file.path(out,'provenance','reference_annotation_dictionary.tsv'))[,.(sample_key,special_group)])
cc<-merge(cc,ann,by='sample_key',all.x=TRUE)
orders<-fread(file.path(out,'provenance','complete_group_order_audit.tsv'))
pcs<-paste0('PC',sort(as.integer(sub('PC','',unique(cc$PC)))))

# GSE235021 bulk RNA-seq groups as in GSE235021/plots/09_gse235021_pca_pc_contribution_heatmap_bulk_only
gse_bulk<-fread(file.path(out,'GSE235021','tables','gse235021_pca_pc_contribution_group_summary.tsv'))[modality=='bulk RNA-seq']
gse_bulk_order<-as.vector(t(outer(paste(c(6,9),'months'),c('Polished rice','Black rice','White/brown rice'),paste,sep=' · ')))
gse_bulk_order<-paste(gse_bulk_order,'Bulk RNA-seq',sep=' · ')
stopifnot(setequal(gse_bulk_order,unique(gse_bulk$display_group)))
group_order<-function(ds)if(ds=='GSE235021')gse_bulk_order else orders[stem==paste0('diagnostics/04_raw_delta_age_by_dataset_treatment_',ds)][order(group_order),plot_group]

group_means<-function(ds){
  if(ds=='GSE235021')return(gse_bulk[,.(dataset_folder='GSE235021',dataset_cohort,modality,group_label=display_group,PC,contribution,sd_contribution,n_samples=n)])
  # Delta-age plot labels: PRJNA1002222 uses its dedicated metadata labels (WT_2M, ...), others simple_group
  d<-cc[dataset_folder==ds];d[,group_label:=if(ds=='PRJNA1002222')special_group else simple_group]
  stopifnot(!anyNA(d$group_label))
  d[,.(contribution=mean(pc_contribution_week),sd_contribution=sd(pc_contribution_week),n_samples=uniqueN(sample_key)),by=.(dataset_folder,dataset_cohort,modality,group_label,PC)]
}

# One colour scale for every figure: limits span all group means drawn in this folder (weeks)
all_means<-rbindlist(lapply(c(datasets,'GSE235021'),group_means))
fill_lim<-c(floor(min(all_means$contribution)),ceiling(max(all_means$contribution)))
message(sprintf('Common Weeks colour scale: %g to %g (midpoint 0)',fill_lim[1],fill_lim[2]))

panel_rows<-function(g)unique(g$layout$t[grepl('^panel',g$layout$name)])
panel_cols<-function(g)unique(g$layout$l[grepl('^panel',g$layout$name)])
fixed<-function(u,skip,conv)sum(vapply(seq_along(u),function(i)if(i%in%skip)0 else conv(u[i],'in',valueOnly=TRUE),numeric(1)))

# Lay out every figure on the device that will draw it: the PDF (cairo_pdf, Arial) and the PNG
# (ragg, wider system font) have different text metrics, as in the reference 07 heatmap.
layout_all<-function(open_dev){
  open_dev();on.exit(dev.off())
  # Reference tile size from GSE235021/plots/09_gse235021_pca_pc_contribution_heatmap_bulk_only
  # (heat_plot: 11 x max(6, .29 * n + 2) in, 6 bulk RNA-seq groups)
  ref<-fread(file.path(out,'GSE235021','tables','gse235021_pca_pc_contribution_group_summary.tsv'))[modality=='bulk RNA-seq']
  ref[,`:=`(group_label=display_group,PC=factor(PC,levels=pcs))]
  n_ref<-uniqueN(ref$group_label);ref_h<-max(6,.29*n_ref+2)
  gref<-ggplotGrob(heat_gg(ref,'GSE235021 bulk RNA-seq: mean frozen PC contributions',fill_lim))
  tile_h<-(ref_h-fixed(gref$heights,panel_rows(gref),convertHeight))/n_ref
  tile_w<-(page_w-fixed(gref$widths,panel_cols(gref),convertWidth))/length(pcs)
  lapply(names(figures),function(stem){
    spec<-figures[[stem]]
    hh<-rbindlist(lapply(spec$datasets,function(ds){
      x<-group_means(ds)
      lv<-group_order(ds)
      stopifnot(length(lv)>0,setequal(lv,unique(x$group_label)))
      x[,row_order:=match(group_label,lv)]
    }))
    stopifnot(!anyDuplicated(unique(hh[,.(dataset_folder,group_label)])$group_label))
    lv_all<-unique(hh[order(match(dataset_folder,spec$datasets),row_order),group_label])
    hh[,`:=`(group_label=factor(group_label,levels=rev(lv_all)),PC=factor(PC,levels=pcs),
             dataset_folder=factor(dataset_folder,levels=spec$datasets))]
    g<-ggplotGrob(heat_gg(hh,spec$title,fill_lim))
    # Panels run top to bottom in dataset order; each gets tile_h per group row
    n_rows<-hh[,.(n=uniqueN(group_label)),by=dataset_folder][order(dataset_folder),n]
    pr<-sort(panel_rows(g));stopifnot(length(pr)==length(n_rows))
    for(k in seq_along(pr))g$heights[pr[k]]<-unit(tile_h*n_rows[k],'in')
    g$widths[panel_cols(g)]<-unit(tile_w*length(pcs),'in')
    # Same page width as the reference: give spare width to the left margin so tiles keep their size
    spare<-page_w-fixed(g$widths,integer(),convertWidth)
    g$widths[1]<-g$widths[1]+unit(max(spare,0),'in')
    list(stem=stem,g=g,h=fixed(g$heights,integer(),convertHeight),w=max(page_w,page_w-spare),
         tile=c(tile_w,tile_h),data=hh[order(dataset_folder,-as.integer(group_label),PC)][,row_order:=NULL])
  })
}
pdf_res<-layout_all(function()cairo_pdf(tempfile(fileext='.pdf'),width=page_w,height=7))
png_res<-layout_all(function()ragg::agg_png(tempfile(fileext='.png'),width=page_w,height=7,units='in',res=300))

for(i in seq_along(figures)){
  rp<-pdf_res[[i]];rg<-png_res[[i]]
  stem<-file.path(outdir,rp$stem)
  ggsave(paste0(stem,'.pdf'),rp$g,width=rp$w,height=rp$h,device=cairo_pdf)
  ggsave(paste0(stem,'.png'),rg$g,width=rg$w,height=rg$h,dpi=300,bg='white',device=ragg::agg_png)
  out_data<-copy(rp$data)[,`:=`(group_label=as.character(group_label),PC=as.character(PC),dataset_folder=as.character(dataset_folder))]
  fwrite(out_data,paste0(stem,'_source.tsv'),sep='\t')
  message(sprintf('%s: %d groups | PDF tile %.3f x %.3f in, page %.2f x %.2f in | PNG tile %.3f x %.3f in, page %.2f x %.2f in',
                  rp$stem,uniqueN(rp$data$group_label),rp$tile[1],rp$tile[2],rp$w,rp$h,rg$tile[1],rg$tile[2],rg$w,rg$h))
}
