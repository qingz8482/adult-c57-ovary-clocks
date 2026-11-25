source(file.path(dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),fixed=TRUE))),'common.R'))
source(file.path(dirname(script),'permutation_build.R'))
Rcpp::sourceCpp(file.path(dirname(script),'pathway_permutation.cpp'),cacheDir=file.path(out,'cache/Rcpp'))
# An enumerable gene-label null with extensive ties, mixed signs and distinct columns.
toy<-cbind(c(0,0,0,0,0,2,3,4),c(0,0,0,0,-4,-2,1,3),1:8,c(-4,-3,-2,-1,0,1,2,3),8:1,c(0,0,0,0,0,0,0,2),c(1,1,2,2,3,3,4,4),c(-2,-2,-1,-1,1,1,2,2))
for(c in 1:8)toy[,c]<-rank(toy[,c])-4.5
members<-list(c(1L,6L),c(2L,5L,8L),c(1L,2L,3L,4L))
exact<-matrix(NA_real_,3,8)
for(p in 1:3){
  cc<-combn(8,length(members[[p]]))
  for(c in 1:8){values<-colSums(matrix(toy[cc,c],nrow=nrow(cc)));observed<-sum(toy[members[[p]],c]);exact[p,c]<-if(c%%2==1)mean(values>=observed-1e-9)else mean(abs(values)>=abs(observed)-1e-9)}
}
counts<-pathway_permutation(toy,members,20000L,481516L)
empirical<-(counts+1)/20001
error<-max(abs(exact-empirical));stopifnot(error<.02)
write_tsv(data.table(pathway=rep(1:3,8),score_column=rep(1:8,each=3),exact_p=as.vector(exact),monte_carlo_p=as.vector(empirical),absolute_error=as.vector(abs(exact-empirical))),'provenance/permutation_engine_exact_null_check.tsv')
cat('PASS: permutation engine versus exact enumeration; maximum error = ',error,'\n')
