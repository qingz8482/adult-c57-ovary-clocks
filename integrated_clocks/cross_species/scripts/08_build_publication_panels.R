#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1] |> sub("^--file=", "", x = _))), "utilities.R"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(ragg) })

pal <- c(elasticnet="#2878B5", pca_glmnet="#D17A22", neutral="#6F7882", negative="#B94C4C", positive="#2C8C6B")
theme_pub <- function(base_size=7.5) theme_classic(base_size=base_size, base_family="sans") +
  theme(axis.line=element_line(linewidth=.35),axis.ticks=element_line(linewidth=.35),
        plot.title=element_text(face="bold",size=base_size+1),strip.text=element_text(face="bold"),
        legend.key.height=unit(3,"mm"),legend.key.width=unit(4,"mm"),plot.margin=margin(4,4,4,4))
theme_set(theme_pub())

save_one <- function(p, stem, width_mm=86, height_mm=70, transparent=FALSE) {
  w<-width_mm/25.4; h<-height_mm/25.4; bg<-if(transparent)"transparent" else "white"
  grDevices::svg(paste0(stem,".svg"),width=w,height=h,bg=bg,family="sans"); print(p); dev.off()
  grDevices::cairo_pdf(paste0(stem,".pdf"),width=w,height=h,family="sans",bg=bg); print(p); dev.off()
  ragg::agg_png(paste0(stem,".png"),width=w,height=h,units="in",res=600,background=bg); print(p); dev.off()
  ragg::agg_tiff(paste0(stem,".tiff"),width=w,height=h,units="in",res=600,background=bg,compression="lzw"); print(p); dev.off()
}

# Figure contract: the frozen mouse clocks do not jointly support a conserved
# human ovarian ageing axis; the apparent ElasticNet trend fails matched-null
# and OOD safeguards, while PCA is directionally inconsistent.
panel_dir<-path_module("1_publication/panels")

pA <- ggplot() + xlim(0,10)+ylim(0,6) +
  annotate("rect",xmin=.2,xmax=2.3,ymin=3.7,ymax=5.4,fill="#E7F1F8",colour=pal["elasticnet"],linewidth=.5) +
  annotate("text",x=1.25,y=4.75,label="Frozen mouse clocks\nElasticNet + PCA",fontface="bold",size=2.4) +
  annotate("segment",x=2.4,xend=3.5,y=4.55,yend=4.55,arrow=arrow(length=unit(2,"mm")),linewidth=.4) +
  annotate("rect",xmin=3.6,xmax=5.6,ymin=3.7,ymax=5.4,fill="#F1F2F3",colour=pal["neutral"],linewidth=.5) +
  annotate("text",x=4.6,y=4.75,label="Ensembl 116\nhigh-confidence 1:1",size=2.35) +
  annotate("segment",x=5.7,xend=6.8,y=4.55,yend=4.55,arrow=arrow(length=unit(2,"mm")),linewidth=.4) +
  annotate("rect",xmin=6.9,xmax=9.8,ymin=3.7,ymax=5.4,fill="#FFF0DF",colour=pal["pca_glmnet"],linewidth=.5) +
  annotate("text",x=8.35,y=4.75,label="GTEx V11 ovary\n193 donors; bracketed age",size=2.35) +
  annotate("segment",x=4.6,xend=4.6,y=3.6,yend=2.6,arrow=arrow(length=unit(2,"mm")),linewidth=.4) +
  annotate("rect",xmin=2,xmax=7.2,ymin=.8,ymax=2.5,fill="white",colour="#444444") +
  annotate("text",x=4.6,y=1.65,label="Stage A: direction conservation\nStage B: frozen mouse-clock projection into GTEx",fontface="bold",size=2.2) +
  annotate("text",x=8.2,y=1.65,label="HPA (n=7)\nStage A exploratory only",size=2.2,colour=pal["neutral"]) +
  annotate("text",x=5,y=.25,label="No refitting • no human calibration • mouse-week output audit-only",size=2.15) +
  theme_void() + ggtitle("Analysis design")

cov<-fread(path_module("results/tables/table_orthology_coverage.tsv"))
cov[,gene_set:=factor(gene_set,levels=rev(gene_set))]
covlong<-melt(cov,id.vars="gene_set",measure.vars=c("raw_gene_coverage","human_detected_coverage","absolute_model_weight_coverage"),variable.name="metric",value.name="coverage")
covlong<-covlong[is.finite(coverage)]
pB<-ggplot(covlong,aes(coverage,gene_set,shape=metric,colour=metric))+geom_point(size=2)+
  geom_vline(xintercept=c(.6,.8),linetype=c(3,2),colour="#AAAAAA")+
  scale_x_continuous(labels=scales::percent_format(),limits=c(0,1))+
  scale_colour_manual(values=c("#777777","#4E9B88","#222222"))+
  labs(x="Coverage",y=NULL,colour=NULL,shape=NULL,title="Orthology and detectability")+theme(legend.position="bottom")

master<-fread(path_module("derived_data/stageA/mouse_human_gene_effect_master.tsv"))
c1<-master[elasticnet_hvg_2000==TRUE & is.finite(mouse_age_beta)&is.finite(human_age_beta_per_decade),.(mouse_age_beta,human_age_beta_per_decade,model_family="ElasticNet")]
c2<-master[pca_versionB_2000==TRUE & is.finite(mouse_age_beta)&is.finite(human_age_beta_per_decade),.(mouse_age_beta,human_age_beta_per_decade,model_family="PCA/glmnet")]
csrc<-rbind(c1,c2)
pC<-ggplot(csrc,aes(mouse_age_beta,human_age_beta_per_decade,colour=model_family))+geom_hline(yintercept=0,colour="#BBBBBB",linewidth=.3)+geom_vline(xintercept=0,colour="#BBBBBB",linewidth=.3)+geom_point(alpha=.28,size=.7)+facet_wrap(~model_family,scales="free")+
  scale_colour_manual(values=c("ElasticNet"=pal["elasticnet"],"PCA/glmnet"=pal["pca_glmnet"]))+
  labs(x="Mouse expression–age slope (logCPM per week)",y="GTEx expression–age slope (logCPM per decade)",title="Observed age-direction conservation")+theme(legend.position="none")

null_panel<-fread(path_module("results/tables/table_random_signature_benchmarks.tsv"))
if(!"signature" %in% names(null_panel)) { null_panel<-null_panel[stage=="Stage A"]; setnames(null_panel,"test","signature") }
enr<-merge(fread(path_module("results/tables/table_stageA_enrichment.tsv")),null_panel,by="signature")
enr[,signature:=factor(signature,levels=rev(signature))]
pD<-ggplot(enr,aes(observed_weighted_signed_t,signature,colour=grepl("elasticnet",as.character(signature))))+geom_vline(xintercept=0,colour="#999999")+geom_segment(aes(x=null_mean,xend=null_q95,yend=signature),colour="#B5B5B5",linewidth=1)+geom_point(size=2)+
  scale_colour_manual(values=c("TRUE"=pal["elasticnet"],"FALSE"=pal["pca_glmnet"]),guide="none")+
  labs(x="Weighted signed human age t",y=NULL,title="Model-direction enrichment vs matched null")

proj<-fread(path_module("derived_data/stageB/gtex_frozen_zero_shot_scores.tsv"))
proj[,model_label:=factor(model_family,levels=c("elasticnet","pca_glmnet"),labels=c("ElasticNet","PCA/glmnet"))]
proj[,age_order:=frank(age_midpoint_sensitivity,ties.method="dense")]
age_labels<-proj[order(age_order),unique(age_bracket)]
summ<-proj[,.(mean=mean(cross_species_mouse_ovarian_clock_score),se=sd(cross_species_mouse_ovarian_clock_score)/sqrt(.N),n=.N),by=.(model_label,age_order,age_bracket)]
pE<-ggplot(proj,aes(age_order,cross_species_mouse_ovarian_clock_score,colour=model_label,fill=model_label))+
  geom_smooth(method="lm",formula=y~x,se=TRUE,level=.95,linewidth=.65,alpha=.10)+
  geom_jitter(width=.12,height=0,alpha=.28,size=.65)+
  geom_errorbar(data=summ,aes(ymin=mean-1.96*se,ymax=mean+1.96*se,y=mean),width=.12,linewidth=.45)+
  geom_point(data=summ,aes(y=mean),size=1.8)+facet_wrap(~model_label,ncol=1)+
  scale_x_continuous(breaks=seq_along(age_labels),labels=age_labels)+
  scale_colour_manual(values=c("ElasticNet"=unname(pal["elasticnet"]),"PCA/glmnet"=unname(pal["pca_glmnet"])))+
  scale_fill_manual(values=c("ElasticNet"=unname(pal["elasticnet"]),"PCA/glmnet"=unname(pal["pca_glmnet"])))+
  labs(x="GTEx public age bracket",y=expression("Standardized predicted age ("*italic(z)[i]*")"),title="Standardized predicted age by bracket")+theme(legend.position="none",axis.text.x=element_text(angle=35,hjust=1))

wide<-dcast(proj,donor_id+age_bracket~model_family,value.var=c("cross_species_mouse_ovarian_clock_score","projection_ood_distance"))
cc<-cor(wide$cross_species_mouse_ovarian_clock_score_elasticnet,wide$cross_species_mouse_ovarian_clock_score_pca_glmnet)
pF<-ggplot(wide,aes(cross_species_mouse_ovarian_clock_score_elasticnet,cross_species_mouse_ovarian_clock_score_pca_glmnet,colour=projection_ood_distance_elasticnet))+geom_hline(yintercept=0,colour="#CCCCCC")+geom_vline(xintercept=0,colour="#CCCCCC")+geom_point(alpha=.65,size=1)+scale_colour_viridis_c(option="C",name="ElasticNet\nOOD distance")+
  annotate("text",x=-Inf,y=Inf,label=sprintf("r = %.02f",cc),hjust=-.1,vjust=1.2,size=2.4)+
  labs(x="ElasticNet standardized\npredicted age (zᵢ)",y="PCA/glmnet standardized\npredicted age (zᵢ)",title="Cross-model agreement and OOD")

plots<-list(A=pA,B=pB,C=pC,D=pD,E=pE,F=pF)
sizes<-list(A=c(183,70),B=c(90,70),C=c(183,75),D=c(90,70),E=c(120,105),F=c(90,75))
for(nm in names(plots)) save_one(plots[[nm]],file.path(panel_dir,paste0("panel_",nm)),sizes[[nm]][1],sizes[[nm]][2])

tagged <- lapply(list(A=pA,B=pB,C=pC,D=pD,E=pE,F=pF), function(z) z + theme(plot.tag=element_text(face="bold",size=9)))
tagged$A <- tagged$A + labs(tag="A")
tagged$B <- tagged$B + labs(tag="B")
tagged$C <- tagged$C + labs(tag="C")
tagged$D <- tagged$D + labs(tag="D")
tagged$E <- tagged$E + labs(tag="E")
tagged$F <- tagged$F + labs(tag="F")
main <- (tagged$A / (tagged$B | tagged$D) / tagged$C / (tagged$E | tagged$F)) +
  plot_layout(heights=c(.75,.9,1,1.3))
assembled<-path_module("1_publication/assembled/human_cross_species_main/human_cross_species_main")
save_one(main,assembled,183,230)

write_tsv(cov,path_module("1_publication/source_data/panel_B_orthology_coverage.tsv"))
write_tsv(csrc,path_module("1_publication/source_data/panel_C_age_direction.tsv"))
write_tsv(enr,path_module("1_publication/source_data/panel_D_enrichment.tsv"))
write_tsv(proj,path_module("1_publication/source_data/panel_E_frozen_scores.tsv"))
write_tsv(wide,path_module("1_publication/source_data/panel_F_model_ood.tsv"))
manifest<-rbindlist(lapply(names(plots),function(nm)data.table(panel=nm,question=c("design","coverage","age-direction","enrichment","score trend","model/OOD")[match(nm,names(plots))],source_data=if(nm=="A")"schematic" else normalizePath(path_module("1_publication/source_data",paste0("panel_",nm,c(B="_orthology_coverage.tsv",C="_age_direction.tsv",D="_enrichment.tsv",E="_frozen_scores.tsv",F="_model_ood.tsv")[[nm]]))),pdf=normalizePath(file.path(panel_dir,paste0("panel_",nm,".pdf"))))))
write_tsv(manifest,path_module("1_publication/figure_manifest.tsv"))
writeLines(c("Figure. Cross-species conservation and strict zero-shot transfer audit.","A, Analysis design. B, Official one-to-one orthology and human detectability coverage. C, Gene-wise mouse expression-age slopes (logCPM per week) versus GTEx expression-age slopes (logCPM per decade). These are adjusted regression slopes, not clock coefficients or correlations. D, Signed model-direction effects relative to expression-, detection-, and variance-matched nulls. E, GTEx-standardized mouse-calibrated predicted age across public GTEx age brackets. Small points denote donors, large points and error bars show bracket means and 95% confidence intervals, and the model-coloured line and soft ribbon show a descriptive linear fit and its 95% confidence interval across ordered age brackets. Values are shown in standard-deviation units and are not predicted human age, human delta age or age acceleration. F, Cross-model agreement and ElasticNet OOD distance. All projections use frozen mouse centers, scales, coefficients, PCA rotations, and mouse calibration; no human refitting or calibration was performed."),path_module("1_publication/captions/human_cross_species_main_caption.txt"))
message("Publication panels and assembled figure exported in R.")
