#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]))), "common_final.R"))

message("[3/5] Sparse 70-79 bracket and donor-influence audit")
x <- load_inputs(); gtex <- x$gtex; audit <- x$audit; qc <- x$qc
ad0 <- x$de$analysis_data
proj <- fread(pmod("derived_data/stageB/gtex_frozen_zero_shot_scores.tsv"))
gmt <- read_gmt(pbase("external_data/ReactomePathways.gmt"),15L,500L)
symbols <- gtex$age_effects$human_gene_symbol[match(rownames(gtex$logcpm),gtex$age_effects$human_gene_id_stable)]
indices <- lapply(gmt,function(gs)which(symbols%in%gs));indices<-indices[lengths(indices)>=10L&lengths(indices)<=500L]

ad <- data.frame(sample_id=audit$sample_id,age_order=as.integer(audit$age_bracket)-1L,
  age_factor=factor(as.character(audit$age_bracket),levels=age_levels),rin_z=zsafe(audit$rin),
  ischemia_z=zsafe(audit$ischemic_time_minutes),mapping_z=zsafe(audit$mapping_rate),
  exonic_z=zsafe(audit$exonic_read_rate),duplication_z=zsafe(audit$duplicate_mapped_read_rate),
  sequencing_center=factor(audit$sequencing_center))

full_rank <- function(d){q<-qr(d);d[,sort(q$pivot[seq_len(q$rank)]),drop=FALSE]}
partial_r2 <- function(keep,age_factor) {
  Y<-qc$pca$x[keep,1:20,drop=FALSE]
  dd<-droplevels(ad[keep,]);dd$age_sens<-droplevels(factor(age_factor))
  f<-lm(Y~age_sens+rin_z+ischemia_z+mapping_z+exonic_z+duplication_z+sequencing_center,data=dd)
  r<-lm(Y~rin_z+ischemia_z+mapping_z+exonic_z+duplication_z+sequencing_center,data=dd)
  sf<-sum(residuals(f)^2);sr<-sum(residuals(r)^2)
  (sr-sf)/sr
}

clock_slopes <- function(keep,age_num) {
  dd<-droplevels(ad[keep,]);dd$age_sens<-age_num
  out<-list()
  for(fam in c("elasticnet","pca_glmnet")) {
    y<-audit[[paste0("existing_score_",fam)]][keep]
    fit<-lm(y~age_sens+rin_z+ischemia_z+mapping_z+exonic_z+duplication_z+sequencing_center,data=dd)
    co<-summary(fit)$coef["age_sens",]
    out[[fam]]<-data.table(model_family=fam,slope=co[1],se=co[2],p=co[4])
  }
  rbindlist(out)
}

fit_spec <- function(label,keep,age_num,age_factor,type) {
  message("  fitting ",label," (n=",sum(keep),")")
  dd<-droplevels(ad[keep,]);dd$age_sens<-age_num
  design<-full_rank(model.matrix(~age_sens+rin_z+ischemia_z+mapping_z+exonic_z+duplication_z+sequencing_center,dd))
  dge<-DGEList(gtex$counts_filtered[,keep,drop=FALSE]);dge$samples$norm.factors<-gtex$normalization_factors[keep]
  fit<-glmQLFit(dge,design,robust=TRUE);test<-glmQLFTest(fit,coef=match("age_sens",colnames(design)))
  tt<-as.data.table(topTags(test,n=Inf,sort.by="none")$table,keep.rownames="human_gene_id_stable")
  tt[,signed_t:=sign(logFC)*sqrt(pmax(F,0))]
  st<-tt$signed_t[match(rownames(gtex$logcpm),tt$human_gene_id_stable)];st[!is.finite(st)]<-0
  cam<-as.data.table(cameraPR(st,indices,use.ranks=TRUE,sort=FALSE),keep.rownames="pathway")
  clocks<-clock_slopes(keep,age_num)
  summary<-data.table(sensitivity=label,sensitivity_type=type,n_donors=sum(keep),global_sum_squared_t=sum(tt$signed_t^2),
    deg_fdr_lt_0_05=sum(tt$FDR<.05),deg_fdr_lfc_gt_0_25=sum(tt$FDR<.05&abs(tt$logFC)>.25),
    partial_R2=partial_r2(keep,age_factor),pathway_fdr_lt_0_05=sum(cam$FDR<.05))
  summary[,`:=`(elasticnet_slope=clocks[model_family=="elasticnet",slope],elasticnet_p=clocks[model_family=="elasticnet",p],
    pca_glmnet_slope=clocks[model_family=="pca_glmnet",slope],pca_glmnet_p=clocks[model_family=="pca_glmnet",p])]
  tt[,`:=`(sensitivity=label,n_donors=sum(keep))]
  list(summary=summary,genes=tt[,.(sensitivity,n_donors,human_gene_id_stable,logFC,logCPM,F,PValue,FDR,signed_t)],
    pathways=cam[,.(pathway,n_genes=NGenes,direction=Direction,p=PValue,fdr=FDR)][,sensitivity:=label])
}

specs<-list()
add_spec<-function(label,keep,age_num,age_factor,type)specs[[label]]<<-list(keep=keep,age_num=age_num,age_factor=age_factor,type=type)
allkeep<-rep(TRUE,nrow(audit));base_age<-as.integer(audit$age_bracket)-1L;base_fac<-factor(as.character(audit$age_bracket),levels=age_levels)
add_spec("primary_all_bins",allkeep,base_age,base_fac,"primary")
keep70<-as.character(audit$age_bracket)!="70-79"
add_spec("exclude_70_79",keep70,base_age[keep70],droplevels(base_fac[keep70]),"sparse_bin")
merged_fac<-factor(ifelse(as.character(audit$age_bracket)%in%c("60-69","70-79"),"60-79",as.character(audit$age_bracket)),levels=c("20-29","30-39","40-49","50-59","60-79"))
merged_num<-as.integer(merged_fac)-1L
add_spec("merge_60_79",allkeep,merged_num,merged_fac,"sparse_bin")
for(a in age_levels){k<-as.character(audit$age_bracket)!=a;add_spec(paste0("leave_bin_",gsub("-","_",a)),k,base_age[k],droplevels(base_fac[k]),"leave_one_age_bin")}
old_ids<-audit[as.character(age_bracket)=="70-79",sample_id]
for(id in old_ids){k<-audit$sample_id!=id;add_spec(paste0("leave_oldest_donor_",id),k,base_age[k],droplevels(base_fac[k]),"leave_one_70_79_donor")}

## Avoid duplicate refit: leave_bin_70_79 is the same data/model as exclude_70_79.
results<-list();gene_tabs<-list();path_tabs<-list()
for(nm in names(specs)) {
  if(nm=="leave_bin_70_79")next
  s<-specs[[nm]];z<-fit_spec(nm,s$keep,s$age_num,s$age_factor,s$type)
  results[[nm]]<-z$summary;gene_tabs[[nm]]<-z$genes;path_tabs[[nm]]<-z$pathways
}
dup<-copy(results[["exclude_70_79"]]);dup[,`:=`(sensitivity="leave_bin_70_79",sensitivity_type="leave_one_age_bin")]
results[["leave_bin_70_79"]]<-dup
gd<-copy(gene_tabs[["exclude_70_79"]]);gd[,sensitivity:="leave_bin_70_79"];gene_tabs[["leave_bin_70_79"]]<-gd
pd<-copy(path_tabs[["exclude_70_79"]]);pd[,sensitivity:="leave_bin_70_79"];path_tabs[["leave_bin_70_79"]]<-pd

summary<-rbindlist(results,fill=TRUE);genes<-rbindlist(gene_tabs,fill=TRUE);paths<-rbindlist(path_tabs,fill=TRUE)
summary[,specification_order:=match(sensitivity,names(specs))];setorder(summary,specification_order);summary[,specification_order:=NULL]
write_tsv(summary,pfinal("results/table_sparse_oldest_sensitivity_summary.tsv"))
write_tsv(genes,pfinal("results/table_sparse_oldest_sensitivity_gene_results.tsv"))
write_tsv(paths,pfinal("results/table_sparse_oldest_sensitivity_pathway_results.tsv"))

## Donor influence and Cook's distance for both clock scores and expression PCs.
message("  computing donor influence and Cook's distances")
cook_rows<-list();loo_rows<-list()
for(fam in c("elasticnet","pca_glmnet")) {
  y<-audit[[paste0("existing_score_",fam)]]
  fit<-lm(y~age_order+rin_z+ischemia_z+mapping_z+exonic_z+duplication_z+sequencing_center,data=ad)
  cd<-cooks.distance(fit);base_slope<-coef(fit)["age_order"]
  cook_rows[[fam]]<-data.table(sample_id=audit$sample_id,donor_id=audit$donor_id,age_bracket=audit$age_bracket,
    model=fam,cooks_distance=cd,cook_flag_gt_4_over_n=cd>4/nrow(audit))
  loo<-vapply(seq_len(nrow(audit)),function(i){coef(update(fit,subset=-i))["age_order"]},numeric(1))
  loo_rows[[fam]]<-data.table(sample_id=audit$sample_id,donor_id=audit$donor_id,age_bracket=audit$age_bracket,
    metric=paste0(fam,"_age_slope"),full_estimate=base_slope,leave_one_out_estimate=loo,influence_change=loo-base_slope)
}
pc_cooks<-sapply(1:20,function(j)cooks.distance(lm(qc$pca$x[,j]~age_order+rin_z+ischemia_z+mapping_z+exonic_z+duplication_z+sequencing_center,data=ad)))
cook_pc<-data.table(sample_id=audit$sample_id,donor_id=audit$donor_id,age_bracket=audit$age_bracket,model="expression_PC1_20",
  cooks_distance=rowMeans(pc_cooks),cooks_distance_max=rowMaxs(pc_cooks),cook_flag_gt_4_over_n=rowMaxs(pc_cooks)>4/nrow(audit))
cooks<-rbindlist(c(cook_rows,list(cook_pc)),fill=TRUE)
write_tsv(cooks,pfinal("results/table_donor_cooks_distance.tsv"))

## Leave-one-donor partial R2 influence (multivariate, inexpensive).
full_r2<-partial_r2(allkeep,base_fac)
loo_r2<-vapply(seq_len(nrow(audit)),function(i){k<-seq_len(nrow(audit))!=i;partial_r2(k,droplevels(base_fac[k]))},numeric(1))
loo_rows[["partial_R2"]]<-data.table(sample_id=audit$sample_id,donor_id=audit$donor_id,age_bracket=audit$age_bracket,
  metric="expression_PC_partial_R2",full_estimate=full_r2,leave_one_out_estimate=loo_r2,influence_change=loo_r2-full_r2)
influence<-rbindlist(loo_rows,fill=TRUE)
write_tsv(influence,pfinal("results/table_donor_leave_one_out_influence.tsv"))

## Sensitivity figures ----------------------------------------------------------
metric_cols<-c("global_sum_squared_t","deg_fdr_lt_0_05","partial_R2","pathway_fdr_lt_0_05","elasticnet_slope","pca_glmnet_slope")
sm<-melt(summary,id.vars=c("sensitivity","sensitivity_type","n_donors"),measure.vars=metric_cols,variable.name="metric",value.name="value")
sm[,primary_value:=value[sensitivity=="primary_all_bins"],by=metric]
sm[,relative_to_primary:=fifelse(abs(primary_value)>1e-12,value/primary_value,value-primary_value)]
p_sens<-ggplot(sm,aes(sensitivity,relative_to_primary,colour=sensitivity_type))+geom_hline(yintercept=1,linetype=2,colour="#888888")+
  geom_point(size=1.25)+facet_wrap(~metric,scales="free_y",ncol=2)+coord_flip()+
  scale_colour_brewer(palette="Dark2")+labs(x=NULL,y="Ratio to primary (slope ratio where applicable)",colour=NULL,title="Sparse-oldest-bin robustness")+
  theme(legend.position="top",axis.text.y=element_text(size=5.2))
save_pub_r(p_sens,pfinal("figures/12_sparse_oldest_sensitivity_metrics"),183,180)
write_tsv(sm,pfinal("figures/source_data/sparse_oldest_sensitivity_metrics.tsv"))

old<-summary[sensitivity_type=="leave_one_70_79_donor"]
ol<-melt(old,id.vars=c("sensitivity","n_donors"),measure.vars=c("elasticnet_slope","pca_glmnet_slope"),variable.name="model",value.name="slope")
p_old<-ggplot(ol,aes(sensitivity,slope,colour=model,group=model))+geom_hline(data=melt(summary[sensitivity=="primary_all_bins"],measure.vars=c("elasticnet_slope","pca_glmnet_slope"),variable.name="model",value.name="slope"),aes(yintercept=slope,colour=model),linetype=2)+
  geom_point(size=1.5)+coord_flip()+scale_colour_manual(values=c(elasticnet_slope=unname(palette_contract["signal_blue"]),pca_glmnet_slope=unname(palette_contract["accent_red"])))+
  labs(x=NULL,y="Adjusted score slope per age bracket",colour=NULL,title="Leave-one-70-79-donor clock slopes")+theme(axis.text.y=element_text(size=5.4),legend.position="top")
save_pub_r(p_old,pfinal("figures/11_leave_one_oldest_donor_clock_slopes"),145,85)

p_cook<-ggplot(cooks,aes(reorder(sample_id,cooks_distance),cooks_distance,colour=age_bracket))+geom_point(size=.8,alpha=.8)+facet_wrap(~model,scales="free_y",ncol=1)+
  scale_colour_manual(values=palette_age)+labs(x="Donors ordered within model",y="Cook's distance",colour="Age",title="Donor influence audit")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="top")
save_pub_r(p_cook,pfinal("figures/10_donor_cooks_distance"),145,125)

saveRDS(list(summary=summary,genes=genes,paths=paths,cooks=cooks,influence=influence),pfinal("derived_data/sparse_oldest_sensitivity_objects.rds"))
message("Sparse-oldest sensitivity complete: ",nrow(summary)," specifications")
