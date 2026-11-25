theme_set(theme_context())

qc<-readRDS(file.path(old,'3_gtex_age_signal_adequacy/derived_data/gtex_age_signal_qc_objects.rds'));audit<-copy(qc$audit)
for(ff in c('elasticnet','pca_glmnet')){d<-cp[cohort=='GTEx'&model_family==ff];ii<-match(audit$sample_id,d$sample_id);stopifnot(!anyNA(ii));audit[[paste0('existing_score_',ff)]]<-d$cross_species_mouse_ovarian_clock_score[ii];audit[[paste0('projection_ood_distance_',ff)]]<-d$projection_ood_distance[ii]}
ad <- data.frame(sample_id=audit$sample_id,age_order=as.integer(audit$age_bracket)-1L,
  age_factor=factor(as.character(audit$age_bracket),levels=age_levels),rin_z=zsafe(audit$rin),
  ischemia_z=zsafe(audit$ischemic_time_minutes),mapping_z=zsafe(audit$mapping_rate),
  exonic_z=zsafe(audit$exonic_read_rate),duplication_z=zsafe(audit$duplicate_mapped_read_rate),
  sequencing_center=factor(audit$sequencing_center))

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


summary<-fread(file.path(old,'4_final_reporting_revision/results/table_sparse_oldest_sensitivity_summary.tsv'))
for(nm in names(specs)){ss<-specs[[nm]];cs<-clock_slopes(ss$keep,ss$age_num);summary[sensitivity==nm,`:=`(elasticnet_slope=cs[model_family=='elasticnet',slope],elasticnet_p=cs[model_family=='elasticnet',p],pca_glmnet_slope=cs[model_family=='pca_glmnet',slope],pca_glmnet_p=cs[model_family=='pca_glmnet',p])]}
write_tsv(summary,file.path(ref_tables,'table_sparse_oldest_sensitivity_summary.tsv'))
cook_rows<-list()
for(ff in c('elasticnet','pca_glmnet')){y<-audit[[paste0('existing_score_',ff)]];fit<-lm(y~age_order+rin_z+ischemia_z+mapping_z+exonic_z+duplication_z+sequencing_center,data=ad);cd<-cooks.distance(fit);cook_rows[[ff]]<-data.table(sample_id=audit$sample_id,donor_id=audit$donor_id,age_bracket=audit$age_bracket,model=ff,cooks_distance=cd,cook_flag_gt_4_over_n=cd>4/nrow(audit))}
cook_pc<-fread(file.path(old,'4_final_reporting_revision/results/table_donor_cooks_distance.tsv'))[model=='expression_PC1_20'];cooks<-rbindlist(c(cook_rows,list(cook_pc)),fill=TRUE)
write_tsv(cooks,file.path(ref_tables,'table_donor_cooks_distance.tsv'))
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

