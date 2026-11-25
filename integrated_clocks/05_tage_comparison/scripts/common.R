suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
setDTthreads(1)
args0 <- commandArgs(FALSE)
script <- normalizePath(gsub('~+~',' ',sub('^--file=','',args0[grepl('^--file=',args0)][1]),fixed=TRUE))
out <- dirname(dirname(script)); project <- dirname(out)
ovary <- dirname(dirname(project))
app <- file.path(project,'lockedmodel_application')
release <- file.path(project,'lockedmodel')  # rebutal4: HVG single-age exclusion + unchanged Version B-PCA
tage <- file.path(ovary,'bulk-pseudobulk_clock/results/external_projection/tAge')
archive <- file.path(ovary,'bulk-pseudobulk_clock/results/external_projection_archive')
models <- c('HVG-ElasticNet','VersionB-PCA/glmnet','tAge scaled-diff','tAge YuGene-diff')
for (d in c('figures','tables','provenance','logs')) dir.create(file.path(out,d),showWarnings=FALSE,recursive=TRUE)
write_tsv <- function(x,p) fwrite(x,file.path(out,p),sep='\t',na='NA')
read_tsv <- function(p) fread(p,na.strings=c('NA',''))
safe_stats <- function(x,y) {
  ok <- is.finite(x)&is.finite(y); x<-x[ok];y<-y[ok]
  if(length(x)<3||sd(x)==0||sd(y)==0) return(data.table(n=length(x),pearson=NA_real_,spearman=NA_real_,pearson_p=NA_real_,spearman_p=NA_real_,direction_agreement=mean(sign(x)==sign(y))))
  a<-cor.test(x,y);b<-suppressWarnings(cor.test(x,y,method='spearman',exact=FALSE))
  data.table(n=length(x),pearson=unname(a$estimate),spearman=unname(b$estimate),pearson_p=a$p.value,spearman_p=b$p.value,direction_agreement=mean(sign(x)==sign(y)))
}
