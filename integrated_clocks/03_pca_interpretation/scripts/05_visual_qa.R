#!/usr/bin/env Rscript
# R-only contact sheets for checking the report PNGs at their correct aspect ratio.
suppressPackageStartupMessages({library(png);library(grid)})
arg<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
arg<-gsub('~+~',' ',arg,fixed=TRUE)
root<-normalizePath(file.path(dirname(arg),'..'))
paths<-c(sort(list.files(file.path(root,'10_pca_diagnostics'),pattern='png$',full.names=TRUE)),
 sort(list.files(file.path(root,'12_PCA deployment model GSEA'),pattern='png$',full.names=TRUE)))
stopifnot(length(paths)==19)
for(page in seq_len(ceiling(length(paths)/4))) {
 png(file.path(root,'provenance',sprintf('visual_qa_%02d.png',page)),width=2600,height=2200,res=150,type='cairo',bg='white')
 grid.newpage()
 for(i in seq_len(4)) {
  j<-(page-1)*4+i
  if(j>length(paths))next
  im<-readPNG(paths[j])
  pushViewport(viewport(x=if(i%%2==1).25 else .75,y=if(i<=2).75 else .25,width=.49,height=.49))
  grid.text(basename(paths[j]),y=.975,gp=gpar(fontsize=7))
  aspect<-dim(im)[2]/dim(im)[1]
  box_aspect<-(2600*.49)/(2200*.49*.92)
  if(aspect>box_aspect)grid.raster(im,y=.46,width=unit(.98,'npc'),height=unit(.98*2600/(aspect*2200),'npc'),interpolate=TRUE)
  else grid.raster(im,y=.46,height=unit(.90,'npc'),width=unit(.90*(2200*.49)*aspect/(2600*.49),'npc'),interpolate=TRUE)
  popViewport()
 }
 dev.off()
}
cat('Wrote five R-generated contact sheets for 19 figures.\n')
