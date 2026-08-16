#!/usr/bin/env Rscript
module_dir <- dirname(normalizePath(commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))][1] |> sub("^--file=","",x=_)))
source(file.path(module_dir,"scripts","utilities.R"))

checks <- list()
add <- function(name, pass, detail) checks[[length(checks)+1L]] <<- data.table(check=name,status=if(pass)"PASS" else "FAIL",detail=as.character(detail))

arts <- readRDS(path_module("derived_data/mouse_models/authoritative_artifacts.rds"))
add("exactly_two_authoritative_artifacts",length(arts)==2L,paste(names(arts),collapse=","))
add("unique_2000_model_genes",all(vapply(arts,function(a)length(a$selected_genes)==2000L&&!anyDuplicated(a$selected_genes),logical(1))),"2,000 unique genes per artifact")
add("centers_scales_aligned",all(vapply(arts,function(a)identical(names(a$gene_center),a$selected_genes)&&identical(names(a$gene_scale),a$selected_genes),logical(1))),"exact name/order match")
add("pca_pc1_to_pc7",identical(arts$pca_glmnet$model$selected_pc_names,paste0("PC",1:7)),paste(arts$pca_glmnet$model$selected_pc_names,collapse=","))
repro<-fread(path_module("results/qc/mouse_artifact_reproduction.tsv"))
add("artifact_reproduction_1e8",all(repro$status=="PASS" & repro$max_abs_predict_vs_formula<=1e-8 & (is.na(repro$max_abs_pc_vs_gene_space)|repro$max_abs_pc_vs_gene_space<=1e-8)),max(c(repro$max_abs_predict_vs_formula,repro$max_abs_pc_vs_gene_space),na.rm=TRUE))

g<-readRDS(path_module("derived_data/gtex/gtex_ovary_prepared.rds"))
add("gtex_raw_counts_primary",is.matrix(g$counts_filtered)&&all(g$counts_filtered>=0),paste(dim(g$counts_filtered),collapse="x"))
add("gtex_unique_donors",!anyDuplicated(g$metadata$donor_id),nrow(g$metadata))
add("public_bracket_ages_only",all(grepl("^[0-9]+-[0-9]+$",as.character(g$metadata$age_bracket))),paste(unique(g$metadata$age_bracket),collapse=","))
h<-readRDS(path_module("derived_data/hpa/hpa_ovary_prepared.rds"))
add("hpa_independent_sample_data",nrow(h$metadata)==7L&&!any(grepl("consensus|gtex",names(h),ignore.case=TRUE)),"seven official HPA samples; no consensus/GTEx track")
orth<-fread(path_module("derived_data/orthology/mouse_human_orthology_primary_one2one.tsv"))
add("primary_mapping_one2one",all(orth$orthology_type=="ortholog_one2one"&orth$orthology_confidence==1),nrow(orth))

proj<-fread(path_module("derived_data/stageB/gtex_frozen_zero_shot_scores.tsv"))
forbidden<-c("human_predicted_age","human_delta_age","age_acceleration")
add("no_forbidden_human_outputs",!any(forbidden%in%names(proj)),paste(intersect(forbidden,names(proj)),collapse=","))
add("primary_endpoint_named_correctly","cross_species_mouse_ovarian_clock_score"%in%names(proj),"present")
add("audit_mouse_week_labeled","mouse_week_equivalent_output_audit_only"%in%names(proj),"present")

src<-list.files(path_module("1_publication/source_data"),pattern="[.]tsv$",full.names=TRUE)
fig<-list.files(path_module("1_publication/assembled/human_cross_species_main"),pattern="human_cross_species_main[.](pdf|svg|png|tiff)$",full.names=TRUE)
add("figure_source_data",length(src)>=5L&&all(file.info(src)$size>0),length(src))
add("assembled_exports",length(fig)==4L&&all(file.info(fig)$size>0),paste(basename(fig),collapse=","))
grade<-fread(path_module("results/tables/table_claim_grade.tsv"))
add("combined_claim_grade",grade[scope=="combined_stageA_stageB",grade]=="INSUFFICIENT_OR_NEGATIVE",grade[scope=="combined_stageA_stageB",grade])
add("no_single_cell_inputs",!any(grepl("single.cell|snRNA|scRNA|GSE202601",fread(path_module("data_manifest/download_manifest.tsv"))$label,ignore.case=TRUE)),"bulk GTEx/HPA only")

out<-rbindlist(checks)
write_tsv(out,path_module("results/qc/validation_summary.tsv"))
if(any(out$status=="FAIL")){print(out);stop("Output validation failed.")}
message("All ",nrow(out)," validation checks PASS.")
