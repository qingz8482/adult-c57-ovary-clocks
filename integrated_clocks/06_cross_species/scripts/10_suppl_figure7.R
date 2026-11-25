# Supplementary Fig. 7: do the frozen mouse ovarian clocks generalize to human ovary?
# Argument order: the test is fair (coverage, human age signal) -> GTEx association and its specificity
# -> independent cohort -> cross-clock consistency -> gene-level direction.
# Reads current-clock result tables only; the human-trained ridge reference is the unchanged cohort-only
# nested-CV analysis. No human refitting or calibration of either mouse clock.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(ggrepel);library(ragg)})
a<-commandArgs(FALSE);here<-dirname(normalizePath(sub('^--file=','',a[grepl('^--file=',a)][1])));out<-dirname(here)
root<-dirname(dirname(dirname(out)))
nested<-file.path(root,'bulk-pseudobulk_clock/human_cross_species_stageA_stageB/4_final_reporting_revision/results/table_nested_cv_predictions_validated.tsv')
src<-file.path(out,'figures/source_data');gsrc<-file.path(src,'GSE202601');sd7<-file.path(src,'Suppl_Figure_7')
dir.create(sd7,showWarnings=FALSE)
rd<-function(p)fread(file.path(out,p))

col_clock<-c('HVG-ElasticNet'='#2878B5','VersionB-PCA'='#D17A22');col_ref<-'#2C8C6B'
lab_clock<-c('HVG-ElasticNet'='ElasticNet (HVG)','VersionB-PCA'='PCA/glmnet')
col_grp<-c('Young'='#3E78B2','Reproductively aged'='#D98C3F')
fmt_p<-function(p)ifelse(p<1e-3,sprintf('%.1e',p),sprintf('%.3f',p))
fisher_ci<-function(r,n){z<-atanh(r);s<-1/sqrt(n-3);tanh(z+c(-1,1)*qnorm(.975)*s)}
hc3<-function(f,d,term){m<-lm(f,data=d);X<-model.matrix(m);u<-residuals(m)/(1-hatvalues(m));inv<-solve(crossprod(X))
  V<-inv%*%crossprod(X,X*u^2)%*%inv;b<-coef(m)[[term]];se<-sqrt(V[term,term]);df<-df.residual(m)
  data.table(estimate=b,ci_low=b-qt(.975,df)*se,ci_high=b+qt(.975,df)*se,p=2*pt(abs(b/se),df,lower.tail=FALSE))}
theme_set(theme_classic(base_size=6.5,base_family='sans')+theme(
  axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3),axis.text=element_text(colour='black'),
  strip.background=element_blank(),strip.text=element_text(face='bold',size=6.5),
  plot.title=element_text(face='bold',size=6.8,margin=margin(0,0,3,0)),plot.tag=element_text(face='bold',size=9),
  legend.key.size=unit(2.8,'mm'),legend.text=element_text(size=5.8),legend.title=element_text(size=6),plot.margin=margin(3,5,3,3)))

# a: design strip
box<-function(x0,x1,lab,fill,col,y0=.3,y1=3.7)list(annotate('rect',xmin=x0,xmax=x1,ymin=y0,ymax=y1,fill=fill,colour=col,linewidth=.4),annotate('text',x=(x0+x1)/2,y=(y0+y1)/2,label=lab,size=1.85,lineheight=.9))
arr<-function(x,xe)annotate('segment',x=x,xend=xe,y=2,yend=2,arrow=arrow(length=unit(1.2,'mm')),linewidth=.3)
pA<-ggplot()+coord_cartesian(xlim=c(0,40),ylim=c(0,4),expand=FALSE)+
  box(.1,7.4,'Frozen mouse ovarian clocks\nElasticNet (HVG), PCA/glmnet\n3,000 input genes each','#E7F1F8','#2878B5')+arr(7.6,8.6)+
  box(8.8,17.2,'Ensembl 116 one-to-one orthologs\nmouse centring, weights, PCA frozen\nno human refitting or calibration','#F1F2F3','#6F7882')+arr(17.4,18.4)+
  box(18.6,24.8,'GTEx V11 ovary\nn = 193, age brackets','#FFF0DF','#D17A22',y0=2.1)+box(18.6,24.8,'GSE202601 ovary\nn = 8, exact ages','#FFF0DF','#D17A22',y1=1.9)+arr(25,26)+
  box(26.2,39.9,'Is the test fair? (b, c)   Gene age directions (d)\nGTEx age association and nulls (e, f)\nIndependent cohort (g)   Clock agreement (h)','white','black')+
  theme_void()+theme(plot.tag=element_text(face='bold',size=9),plot.margin=margin(2,3,2,3))

# b: coverage
cov<-fread(file.path(gsrc,'GTEx_GSE202601_gene_coverage_source_data.tsv'))
cov[,model:=ifelse(model_family=='elasticnet','HVG-ElasticNet','VersionB-PCA')]
cov[,row:=factor(paste0(cohort,' · ',lab_clock[model]),levels=rev(paste0(rep(c('GTEx','GSE202601'),each=2),' · ',lab_clock)))]
covl<-melt(cov,id.vars=c('row','model'),measure.vars=c('human_detected_coverage','absolute_model_weight_coverage'),variable.name='metric',value.name='coverage')
covl[,metric:=factor(metric,labels=c('Input genes measured','Absolute clock weight'))]
pB<-ggplot(covl,aes(coverage,row))+geom_vline(xintercept=.8,linetype=2,colour='#AAAAAA',linewidth=.3)+
  geom_line(aes(group=row),colour='#CCCCCC',linewidth=.5)+geom_point(aes(shape=metric,colour=model),size=1.7)+
  scale_x_continuous(labels=scales::percent_format(accuracy=1),limits=c(.5,1),breaks=c(.5,.75,1))+
  scale_colour_manual(values=col_clock,guide='none')+scale_shape_manual(values=c(16,15),name=NULL)+
  labs(x='Coverage in human data',y=NULL,title='Clock inputs are measured')+
  theme(legend.position='bottom',legend.margin=margin(-4,0,0,0))

# c: human-trained reference (held-out nested CV, GTEx)
ft<-rd('figures/source_data/panel_E_frozen_scores.tsv')[model=='HVG-ElasticNet',.(sample_id,age_bracket,age_order,rin_z,ischemia_z,log_library_z)]
nr<-fread(nested)[specification=='baseline_age_stratified_5fold']
nr<-merge(nr[,.(sample_id,predicted)],ft,by='sample_id');stopifnot(nrow(nr)==193)
nr_rho<-cor.test(nr$predicted,nr$age_order,method='spearman',exact=FALSE)
br<-nr[order(age_order),unique(age_bracket)]
set.seed(7)
pC<-ggplot(nr,aes(age_order,predicted))+
  geom_jitter(width=.15,height=0,colour=col_ref,alpha=.45,size=.6,stroke=0)+
  stat_summary(fun.data=function(x){m<-mean(x);s<-1.96*sd(x)/sqrt(length(x));data.frame(y=m,ymin=m-s,ymax=m+s)},colour='black',size=.25,linewidth=.4)+
  annotate('text',x=.6,y=Inf,label=sprintf('Spearman ρ = %.2f\nn = %d donors, held out',nr_rho$estimate,nrow(nr)),hjust=0,vjust=1.1,size=1.8,lineheight=.9)+
  scale_x_continuous(breaks=seq_along(br),labels=br)+
  labs(x='GTEx age bracket',y='Held-out predicted age (years)',title='GTEx carries learnable age signal',subtitle='Human-trained ridge, nested CV (reference)')+
  theme(axis.text.x=element_text(angle=40,hjust=1),plot.subtitle=element_text(size=5.6,colour='#444444'))

# d: GTEx clock output by bracket
gt<-rd('figures/source_data/panel_E_frozen_scores.tsv');gt[,clock:=factor(lab_clock[model],levels=lab_clock)]
gs<-gt[,.(mean=mean(cross_species_mouse_ovarian_clock_score),se=sd(cross_species_mouse_ovarian_clock_score)/sqrt(.N)),by=.(model,clock,age_order)]
prim<-rd('tables/age_association_sensitivities.tsv')[cohort=='GTEx'&analysis=='primary_adjusted_ordered']
prim[,clock:=factor(lab_clock[model],levels=lab_clock)][,label:=sprintf('slope = %.2f SD/bracket\nP = %s',estimate,fmt_p(p))]
pD<-ggplot(gt,aes(age_order,cross_species_mouse_ovarian_clock_score,colour=model,fill=model))+
  geom_hline(yintercept=0,colour='#DDDDDD',linewidth=.3)+
  geom_jitter(width=.15,height=0,alpha=.3,size=.5,stroke=0)+
  geom_smooth(method='lm',formula=y~x,se=TRUE,linewidth=.55,alpha=.15)+
  geom_errorbar(data=gs,aes(y=mean,ymin=mean-1.96*se,ymax=mean+1.96*se),width=.2,linewidth=.4)+
  geom_point(data=gs,aes(y=mean),size=1.3)+
  geom_text(data=prim,aes(x=.6,y=Inf,label=label),inherit.aes=FALSE,hjust=0,vjust=1.1,size=1.8,lineheight=.9)+
  facet_wrap(~clock,nrow=1)+scale_x_continuous(breaks=seq_along(br),labels=br)+
  scale_y_continuous(expand=expansion(mult=c(.04,.2)))+
  scale_colour_manual(values=col_clock)+scale_fill_manual(values=col_clock)+
  labs(x='GTEx age bracket',y=expression('Clock output ('*italic(z)*', within cohort)'),title='GTEx: only ElasticNet rises with age')+
  theme(legend.position='none',axis.text.x=element_text(angle=40,hjust=1))

# e: specificity of the GTEx slope (all on the technical-adjusted ordered-age slope scale)
ref<-hc3(z~age_order+rin_z+ischemia_z+log_library_z,nr[,z:=(predicted-mean(predicted))/sd(predicted)],'age_order')
nulls<-c(age_label_permutation_sensitivity='Shuffled age labels',mapped_weight_permutation='Clock weights permuted',matched_gene_signature='Matched random genes')
nd<-rbindlist(lapply(names(lab_clock),function(m)rbindlist(lapply(names(nulls),function(n)
  data.table(model=m,null=n,statistic=fread(file.path(out,'null_distributions',sprintf('stageB_%s_%s.tsv',m,n)))$statistic)))))
nb<-rd('tables/stageB_null_benchmarks.tsv')[null_type%in%names(nulls)]
rows<-c('Observed clock',unname(nulls))
nd[,row:=factor(nulls[null],levels=rev(rows))];nb[,row:=factor(nulls[null_type],levels=rev(rows))]
obs<-prim[,.(model,estimate,ci_low,ci_high,row=factor('Observed clock',levels=rev(rows)))]
nb[,plab:=sprintf('P = %s',fmt_p(empirical_one_sided_p))]
for(d in list(nd,nb,obs))d[,clock:=factor(lab_clock[model],levels=lab_clock)]
xr<-c(-.42,max(.68,ref$ci_high+.02))
pE<-ggplot()+
  annotate('rect',xmin=ref$ci_low,xmax=ref$ci_high,ymin=-Inf,ymax=Inf,fill=col_ref,alpha=.12)+
  geom_vline(xintercept=ref$estimate,colour=col_ref,linewidth=.5)+geom_vline(xintercept=0,colour='#BBBBBB',linewidth=.3)+
  geom_violin(data=nd,aes(statistic,row),orientation='y',fill='#D9D9D9',colour=NA,scale='width',width=.75,trim=TRUE)+
  geom_errorbar(data=nb,aes(xmin=null_q025,xmax=null_q975,y=row),orientation='y',width=0,colour='#6B6B6B',linewidth=.45)+
  geom_errorbar(data=obs,aes(xmin=ci_low,xmax=ci_high,y=row,colour=model),orientation='y',width=.25,linewidth=.5)+
  geom_point(data=obs,aes(estimate,row,colour=model),size=2)+
  geom_vline(data=obs,aes(xintercept=estimate,colour=model),linetype=2,linewidth=.3)+
  geom_text(data=nb,aes(x=null_q975+.03,y=row,label=plab),hjust=0,vjust=-.6,size=1.75)+
  annotate('text',x=ref$ci_low-.01,y=-Inf,label='human-trained\nreference',vjust=-.3,hjust=1,size=1.7,colour=col_ref,lineheight=.85)+
  facet_wrap(~clock,nrow=1)+scale_colour_manual(values=col_clock,guide='none')+
  scale_y_discrete(expand=expansion(add=c(.6,1.1)))+coord_cartesian(xlim=xr)+
  labs(x='Technical-adjusted age slope (SD per GTEx bracket)',y=NULL,
       title='ElasticNet slope beats shuffled ages, but not random weights or random genes',
       subtitle='Grey: 10,000 null draws (bar = 95% interval); P = one-sided, observed vs null')+
  theme(plot.subtitle=element_text(size=5.6,colour='#444444'),panel.spacing=unit(4,'mm'))

# f: GSE202601
ge<-fread(file.path(gsrc,'standardized_score_vs_chronological_age_source_data.tsv'))
ge[,model:=ifelse(model_family=='elasticnet','HVG-ElasticNet','VersionB-PCA')][,clock:=factor(lab_clock[model],levels=lab_clock)]
ge[,group:=factor(ifelse(age_group=='young','Young','Reproductively aged'),levels=names(col_grp))]
fp<-rd('tables/frozen_predictions.tsv')[cohort=='GSE202601',.(projection_sample_id=sample_id,model,outside=projection_ood_distance>mouse_ood_q99)]
if(!all(ge$projection_sample_id%in%fp$projection_sample_id))fp<-rd('tables/frozen_predictions.tsv')[cohort=='GSE202601',.(sample_id=donor_id,model,outside=projection_ood_distance>mouse_ood_q99)]
ge<-merge(ge,fp,by=intersect(c('projection_sample_id','sample_id','model'),names(fp)));stopifnot(nrow(ge)==16)
ge[,support:=factor(ifelse(outside,'Outside','Within'),levels=c('Within','Outside'))]
gst<-rd('tables/age_correlations.tsv')[cohort=='GSE202601'];gst[,clock:=factor(lab_clock[model],levels=lab_clock)]
gst<-merge(gst,ge[,.(n_out=sum(outside)),by=model],by='model')
gst[,label:=sprintf('ρ = %.2f, exact P = %.3f\n%d/8 outside mouse support',spearman,spearman_p,n_out)]
pF<-ggplot(ge,aes(age_year,cross_species_mouse_ovarian_clock_score))+
  geom_smooth(method='lm',formula=y~x,se=TRUE,colour='#4A4A4A',fill='#BFC3C7',alpha=.35,linewidth=.4)+
  geom_point(aes(colour=group,shape=support),size=1.6,stroke=.6)+
  geom_text_repel(aes(label=sample_id,colour=group),size=1.7,seed=1,box.padding=.2,min.segment.length=.2,segment.size=.2,segment.colour='#999999',show.legend=FALSE)+
  geom_text(data=gst,aes(x=-Inf,y=Inf,label=label),hjust=-.06,vjust=1.1,size=1.7,lineheight=.9)+
  facet_wrap(~clock,nrow=1)+scale_colour_manual(values=col_grp,name=NULL)+
  scale_shape_manual(values=c(Within=16,Outside=1),name='Mouse support (99%)',drop=FALSE)+
  scale_x_continuous(expand=expansion(mult=c(.08,.12)))+scale_y_continuous(expand=expansion(mult=c(.05,.3)))+
  labs(x='Chronological age (years)',y=expression('Clock output ('*italic(z)*', within cohort)'),title='GSE202601: positive but non-significant (n = 8)')+
  theme(legend.position='right',legend.box='vertical',legend.margin=margin(0,0,0,-4))

# g: cross-clock agreement, mouse vs human
ag<-rd('tables/reference_panels/cross_model_reference_agreement.tsv')
gse_ag<-fread(file.path(gsrc,'GTEx_GSE202601_cross_model_agreement_statistics.tsv'))[cohort=='GSE202601']
ag<-rbind(ag[,.(dataset,n,r=pearson_r)],data.table(dataset='GSE202601',n=gse_ag$n_samples,r=gse_ag$pearson_r))
ag[,label:=c(mouse_deployment_reference_shared_110='Mouse, in-sample',paired_nested_CV_mouse_predictions='Mouse, nested CV (held out)',GTEx_193_donors='Human GTEx',GSE202601='Human GSE202601')[dataset]]
ag[,species:=ifelse(grepl('^Mouse',label),'Mouse','Human')]
ag[,c('lo','hi'):=as.list(fisher_ci(r,n)),by=dataset]
ag[,label:=factor(sprintf('%s (n = %d)',label,n),levels=rev(sprintf('%s (n = %d)',label,n)))]
pG<-ggplot(ag,aes(r,label,colour=species))+geom_vline(xintercept=0,colour='#999999',linewidth=.3)+
  geom_errorbar(aes(xmin=lo,xmax=hi),orientation='y',width=0,linewidth=.5)+geom_point(size=1.8)+
  geom_text(aes(label=sprintf('%.2f',r)),vjust=-.8,size=1.7,show.legend=FALSE)+
  scale_colour_manual(values=c(Mouse='#555555',Human='#B94C4C'),guide='none')+scale_x_continuous(limits=c(-1,1))+
  labs(x='Pearson r, ElasticNet vs PCA/glmnet output (95% CI)',y=NULL,title='Clocks agree in mouse, not in human')

# h: gene-level direction conservation
sa<-rd('tables/stageA_direction_conservation.tsv')
keep<-c(shared_model_genes='Shared by both clocks',hvg_active_genes='ElasticNet nonzero',hvg_input_genes='ElasticNet inputs',pca_input_genes='PCA/glmnet inputs',testable_mouse_universe='All testable genes')
sa<-sa[gene_set%in%names(keep)][order(match(gene_set,names(keep)))]
sa[,lab:=sprintf('%s (%s)',keep[gene_set],format(n_genes,big.mark=',',trim=TRUE))][,lab:=factor(lab,levels=rev(lab))]
sa[,sig:=q<.05]
pH<-ggplot(sa,aes(y=lab))+geom_vline(xintercept=.5,colour='#999999',linewidth=.3,linetype=2)+
  geom_errorbar(aes(xmin=matched_null_q025,xmax=matched_null_q975),orientation='y',width=0,colour='#C8C8C8',linewidth=2.2)+
  geom_point(aes(x=matched_null_mean),shape='|',size=2,colour='#6B6B6B')+
  geom_point(aes(x=sign_concordance,colour=sig),size=1.8)+
  geom_text(aes(x=.64,label=sprintf('q = %s',ifelse(q<.001,sprintf('%.1e',q),sprintf('%.3f',q)))),hjust=1,size=1.7)+
  scale_colour_manual(values=c(`TRUE`='#B94C4C',`FALSE`='#333333'),guide='none')+
  scale_x_continuous(labels=scales::percent_format(accuracy=1),limits=c(.3,.64),breaks=seq(.3,.6,.1))+
  labs(x='Genes with same age direction',y=NULL,title='Gene age directions: near chance',subtitle='Grey: matched-null 95% interval; | null mean')+
  theme(plot.subtitle=element_text(size=5.6,colour='#444444'))

# Rows are assembled independently so long y-axis labels in one panel do not squeeze the others.
row<-function(...,widths)wrap_elements(full=wrap_plots(...,nrow=1,widths=widths))
fig<-wrap_plots(wrap_elements(full=pA+labs(tag='a')),
  row(pB+labs(tag='b'),pC+labs(tag='c'),pH+labs(tag='d'),widths=c(1,1,1.05)),
  row(pD+labs(tag='e'),pE+labs(tag='f'),widths=c(1,1.6)),
  row(pF+labs(tag='g'),pG+labs(tag='h'),widths=c(1.35,1)),
  ncol=1,heights=c(.3,1,1.05,1))+
  plot_annotation(title='Supplementary Fig. 7 | Frozen mouse ovarian clocks show no convergent transfer to human ovary',
    theme=theme(plot.title=element_text(face='bold',size=8.5)))

stem<-file.path(out,'figures/Suppl_Figure_7_human_generalization');w<-200/25.4;h<-225/25.4
grDevices::cairo_pdf(paste0(stem,'.pdf'),width=w,height=h,family='sans');print(fig);invisible(dev.off())
grDevices::svg(paste0(stem,'.svg'),width=w,height=h,family='sans');print(fig);invisible(dev.off())
ragg::agg_png(paste0(stem,'.png'),width=w,height=h,units='in',res=600);print(fig);invisible(dev.off())
ragg::agg_tiff(paste0(stem,'.tiff'),width=w,height=h,units='in',res=600,compression='lzw');print(fig);invisible(dev.off())

unlink(list.files(sd7,full.names=TRUE))
fwrite(cov,file.path(sd7,'panel_b_coverage.tsv'),sep='\t')
fwrite(nr,file.path(sd7,'panel_c_human_ridge_heldout.tsv'),sep='\t')
fwrite(gs,file.path(sd7,'panel_e_GTEx_bracket_means.tsv'),sep='\t')
fwrite(rbind(data.table(row='human_trained_reference',model=NA,ref),obs[,.(row=as.character(row),model,estimate,ci_low,ci_high,p=prim$p)],fill=TRUE),file.path(sd7,'panel_f_slopes.tsv'),sep='\t')
fwrite(nb[,.(model,null_type,observed,null_q025,null_q975,empirical_one_sided_p)],file.path(sd7,'panel_f_null_summary.tsv'),sep='\t')
fwrite(ge,file.path(sd7,'panel_g_GSE202601.tsv'),sep='\t')
fwrite(ag,file.path(sd7,'panel_h_clock_agreement.tsv'),sep='\t')
fwrite(sa,file.path(sd7,'panel_d_gene_direction.tsv'),sep='\t')
message(sprintf('Human-trained reference adjusted slope = %.3f [%.3f, %.3f]; ridge rho = %.3f',ref$estimate,ref$ci_low,ref$ci_high,nr_rho$estimate))
message('Supplementary Fig. 7 written: ',stem)
