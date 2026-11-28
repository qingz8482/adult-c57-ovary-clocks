#!/usr/bin/env Rscript
# Input: genes x profiles CSV/TSV in the exact training ln1p(CPM) representation.
a<-commandArgs(FALSE);here<-dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE)))
source(file.path(here,'core.R'));args<-commandArgs(TRUE)
assert(length(args)==3,'Usage: Rscript predict_frozen.R model.rds logcpm.csv output.csv')
model<-readRDS(args[1]);z<-fread(args[2]);g<-z[[1]];assert(!anyDuplicated(g),'Duplicate genes in prediction matrix')
x<-as.matrix(z[,-1,with=FALSE]);storage.mode(x)<-'double';rownames(x)<-g
assert(all(is.na(x)|is.finite(x)&x>=0),'Expect nonnegative ln1p(CPM); NA means unavailable gene')
pred<-predict_model(model,x);fwrite(pred,args[3])
if(any(pred$missing_gene_fraction>.1))warning('More than 10% of model genes missing in at least one profile; inspect transfer suitability.')
