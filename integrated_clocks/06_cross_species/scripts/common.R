suppressPackageStartupMessages({library(data.table);library(glmnet);library(limma);library(Matrix)})
setDTthreads(1)
a<-commandArgs(FALSE);here<-dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE)))
out<-dirname(here);project<-dirname(out);root<-dirname(dirname(project))
dev<-file.path(project,'model_development');release<-file.path(project,'lockedmodel')
# Nested-CV development outputs per locked clock: HVG single-age exclusion (rebutal4) and the unchanged PCA (rebutal3).
dev_by_model<-c('HVG-ElasticNet'=dev,'VersionB-PCA'=file.path(dirname(project),'rebutal3_batcheffects removal_CV_LOO/model_development'))
old<-file.path(root,'bulk-pseudobulk_clock/human_cross_species_stageA_stageB')
gseold<-file.path(root,'bulk-pseudobulk_clock/human_GSE202601_integrated_projection')
source(file.path(release,'scripts/core.R'))
SEED<-20261012L;B<-10000L;set.seed(SEED)
wt<-function(x,p)fwrite(x,file.path(out,p),sep='\t',na='NA')
rt<-function(p)fread(file.path(out,p))
zsafe<-function(x){s<-sd(x);if(!is.finite(s)||s==0)stop('Zero variance or missing score');as.numeric((x-mean(x))/s)}
models<-list('HVG-ElasticNet'=readRDS(file.path(release,'HVG-ElasticNet-exclude-single-age/model.rds')),'VersionB-PCA'=readRDS(file.path(release,'VersionB-PCA-glmnet/model.rds')))
hc3<-function(formula,d,term){
 f<-lm(formula,data=d);X<-model.matrix(f);n<-nrow(X);df<-df.residual(f)
 if(qr(X)$rank<ncol(X)||df<3)return(data.table(n=n,estimate=NA_real_,se=NA_real_,ci_low=NA_real_,ci_high=NA_real_,p=NA_real_,status='not_estimable'))
 h<-hatvalues(f);u<-residuals(f)/pmax(1-h,1e-8);inv<-solve(crossprod(X));V<-inv%*%crossprod(X,X*u^2)%*%inv
 b<-coef(f)[term];se<-sqrt(diag(V))[term];pv<-2*pt(abs(b/se),df,lower.tail=FALSE)
 data.table(n=n,estimate=unname(b),se=unname(se),ci_low=unname(b-qt(.975,df)*se),ci_high=unname(b+qt(.975,df)*se),p=unname(pv),status='estimable',formula=paste(deparse(formula),collapse=' '))
}
wcor<-function(x,y,w){x<-rank(x);y<-rank(y);w<-w/sum(w);x<-x-sum(w*x);y<-y-sum(w*y);sum(w*x*y)/sqrt(sum(w*x*x)*sum(w*y*y))}
readmatrix<-function(p){z<-fread(p);x<-as.matrix(z[,-1,with=FALSE]);storage.mode(x)<-'double';rownames(x)<-z[[1]];x}
