from pathlib import Path
import shutil,csv,hashlib
root=Path(__file__).resolve().parents[4]
out=root/'rebutal/rebutal4/generalization_tohuman'
old=root/'bulk-pseudobulk_clock/human_cross_species_stageA_stageB'
gse=root/'bulk-pseudobulk_clock/human_GSE202601_integrated_projection'
archive=out/'provenance/pre_reference_style_20261012';archive.mkdir(exist_ok=True)
if False:  # rebutal4 starts from a clean folder; there are no earlier figures to archive
 shutil.move(out/'figures',archive/'figures')
 shutil.copy2(out/'scripts/04_figures.R',archive/'04_figures.R')
 shutil.copytree(out/'reports',archive/'reports',dirs_exist_ok=True)
 shutil.copy2(out/'README.md',archive/'README.md')
for d in ['figures/source_data','inputs/reference_adapter/cross_species','inputs/reference_adapter/gse/projection','scripts/reference_renderers','provenance/reference_sources','tables/reference_panels']:(out/d).mkdir(parents=True,exist_ok=True)
changed={'pca_by_elasticnet_ood_distance','pca_by_pca_glmnet_ood_distance','expression_pca_qc_panel','gtex_age_signal_adequacy_integrated','10_donor_cooks_distance','11_leave_one_oldest_donor_clock_slopes','12_sparse_oldest_sensitivity_metrics'}
rows=[]
for ref in [old,gse]:
 for f in sorted(ref.rglob('*.pdf')):
  rel=f.relative_to(ref);stem=f.stem
  dependent=(ref==gse or rel.parts[0] in ['1_publication','2_revision_defensibility'] or stem in changed)
  mode='recomputed_current_models_original_plot_code' if dependent else 'unchanged_cohort_analysis_exact_copy'
  rows.append(dict(reference_module=ref.name,reference_path=str(rel),output_stem=stem,status=mode,reference_sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
  if not dependent:
   for ext in ['.pdf','.png','.svg','.tiff']:
    p=f.with_suffix(ext)
    if p.exists():shutil.copy2(p,out/'figures'/p.name)
assert len(rows)==len({r['output_stem'] for r in rows})
with (out/'provenance/reference_panel_manifest.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
for part in ['3_gtex_age_signal_adequacy','4_final_reporting_revision']:
 shutil.copytree(old/part/'figures/source_data',out/'figures/source_data'/part,dirs_exist_ok=True)
refs=[old/'scripts/08_build_publication_panels.R',old/'2_revision_defensibility/scripts/01_run_revision.R',old/'3_gtex_age_signal_adequacy/scripts/01_metadata_expression_qc.R',old/'3_gtex_age_signal_adequacy/scripts/04_figures_reports.R',old/'4_final_reporting_revision/scripts/03_sparse_oldest_sensitivity.R']+list((gse/'scripts').glob('plot_*.R'))
for i,f in enumerate(refs):shutil.copy2(f,out/'provenance/reference_sources'/f'{i:02d}_{f.name}')
print(len(rows),'reference figure stems;',sum(r['status'].startswith('recomputed') for r in rows),'require current-model rendering')
scripts=out/'scripts/reference_renderers'
s=(old/'scripts/08_build_publication_panels.R').read_text()
s=s[s.index('pal <-'):].replace('hvg_2000','hvg_3000').replace('versionB_2000','versionB_3000')
s=s.replace('slopes (logCPM per week)', 'slopes (ln1pCPM per week)')
s=s.replace('Mouse expression–age slope (logCPM per week)','Mouse expression–age slope (ln1pCPM per week)').replace('GTEx-standardized mouse-calibrated predicted age','GTEx-standardized frozen mouse-clock output').replace('and mouse calibration','with no calibration')
(scripts/'01_publication.R').write_text(s)
s=(gse/'scripts/plot_predicted_vs_chronological_age.R').read_text();s=s[s.index('x <- fread(projection_path)'):]
head="projection_path<-file.path(gcompat,'projection/GSE202601_two_integrated_clock_projections.tsv');plot_dir<-figdir;source_dir<-file.path(figdir,'source_data/GSE202601');dir.create(source_dir,recursive=TRUE,showWarnings=FALSE);set.seed(202601)\n"
(scripts/'02_gse_age.R').write_text(head+s)
s=(gse/'scripts/plot_gene_coverage.R').read_text();s=s[s.index('x <- fread(coverage_path)'):s.index('expected <- data.table(')].replace('4000L','6000L').replace('hvg_2000','hvg_3000').replace('versionB_2000','versionB_3000')
head="coverage_path<-file.path(gcompat,'projection/GSE202601_model_gene_coverage.tsv');gtex_coverage_path<-file.path(figdir,'source_data/panel_B_orthology_coverage.tsv');plot_dir<-figdir;source_dir<-file.path(figdir,'source_data/GSE202601')\n"
extra="""
z<-copy(coverage_long[cohort=='GSE202601']);z[,row_label:=factor(ifelse(model_family=='elasticnet','ElasticNet/HVG (3,000 genes)','PCA/glmnet (3,000 genes)'),levels=c('PCA/glmnet (3,000 genes)','ElasticNet/HVG (3,000 genes)'))]
z[,metric:=factor(as.character(metric),levels=levels(coverage_long$metric),labels=c('High-confidence 1:1 orthologue','Detected in GSE202601','Absolute model-weight coverage'))]
p2<-p %+% z;p2$layers<-p2$layers[-2];p2<-p2+labs(title='GSE202601 orthology and detectability')+scale_colour_manual(values=setNames(unname(palette),levels(z$metric)))+scale_shape_manual(values=setNames(unname(shapes),levels(z$metric)))
save_bundle(p2,file.path(figdir,'03_GSE202601_gene_coverage'),90,80)
write_tsv(z,file.path(source_dir,'GSE202601_gene_coverage_source_data.tsv'))
"""
(scripts/'03_gse_coverage.R').write_text(head+s+extra)
s=(gse/'scripts/plot_cross_model_agreement.R').read_text();s=s[s.index('gtex <- fread(gtex_path)'):s.index('qa <- data.table(')]
head="gtex_path<-file.path(figdir,'source_data/panel_F_model_ood.tsv');gse_path<-file.path(gcompat,'projection/GSE202601_two_integrated_clock_projections.tsv');plot_dir<-figdir;source_dir<-file.path(figdir,'source_data/GSE202601')\n"
(scripts/'04_gse_agreement.R').write_text(head+s)
s=(old/'2_revision_defensibility/scripts/01_run_revision.R').read_text();s=s[s.index('## Figures '):s.index('## Reports ')]
s=s.replace('elasticnet_hvg_2000','elasticnet_hvg_3000').replace('pca_versionB_2000','pca_versionB_3000')
s=s.replace('coord_cartesian(ylim=c(-.075,.01))','coord_cartesian(ylim=range(c(-.075,.01,pcdat$spearman_correlation*1.25)))').replace('coord_cartesian(ylim=c(-.7,.15))','coord_cartesian(ylim=range(c(-.7,.15,pfdat$value*1.2)))')
s=s.replace('PCA/glmnet PC1-PC7 contributions','PCA/glmnet PC1-PC10 contributions').replace('Calibrated age-slope contribution','Frozen output age-slope contribution')
(scripts/'05_supplements.R').write_text("source(file.path(here,'reference_supplement_data.R'),local=TRUE)\nset.seed(SEED)\n"+s)
# Preserve actual script source and modifications for transparent reproduction.
print('Generated original-code renderers for publication, GSE202601 and supplements.')
# Use the exact cohort palette and themes in each source module.
a=(old/'3_gtex_age_signal_adequacy/scripts/common.R').read_text()

base="""
qc<-readRDS(file.path(old,'3_gtex_age_signal_adequacy/derived_data/gtex_age_signal_qc_objects.rds'));audit<-copy(qc$audit)
for(ff in c('elasticnet','pca_glmnet')){d<-cp[cohort=='GTEx'&model_family==ff];ii<-match(audit$sample_id,d$sample_id);stopifnot(!anyNA(ii));audit[[paste0('existing_score_',ff)]]<-d$cross_species_mouse_ovarian_clock_score[ii];audit[[paste0('projection_ood_distance_',ff)]]<-d$projection_ood_distance[ii]}
"""
theme=a[a.index('palette_age <-'):a.index('theme_set(theme_nature_contract())')+len('theme_set(theme_nature_contract())')]
s=(old/'3_gtex_age_signal_adequacy/scripts/01_metadata_expression_qc.R').read_text();s=s[s.index('pdat <- merge(pcs'):s.index('## Library/detected-gene')]
s='\n'.join(l for l in s.splitlines() if not l.startswith(('save_pub_r(p_age,','save_pub_r(p_rin,','save_pub_r(p_isch,','save_pub_r(p_center,','save_pub_r(p_det,')))
feature="""
gtex<-readRDS(file.path(old,'derived_data/gtex/gtex_ovary_prepared.rds'));logcpm<-gtex$logcpm
ordered_primary<-readRDS(file.path(old,'3_gtex_age_signal_adequacy/derived_data/gtex_age_de_objects.rds'))$ordered_primary$table
stat_primary<-ordered_primary$signed_quasi_likelihood_t[match(rownames(logcpm),ordered_primary$human_gene_id_stable)]
symbols<-gtex$age_effects$human_gene_symbol[match(rownames(logcpm),gtex$age_effects$human_gene_id_stable)]
stagea<-copy(master0);stagea[,deployment_shared_artifact:=hvg_selected&pca_selected]
"""
f=(old/'3_gtex_age_signal_adequacy/scripts/03_multivariate_pathway_prediction.R').read_text();f=f[f.index('universe <- data.table'):f.index('## Direct histology/pathology')].replace('2000','3000').replace('artifact_shared_297',f'artifact_shared_current').replace('sample(pools[[s]], 1L)','pools[[s]][sample.int(length(pools[[s]]), 1L)]')
f=f.replace('paudit("results/', 'file.path(ref_tables,"')
fig=(old/'3_gtex_age_signal_adequacy/scripts/04_figures_reports.R').read_text();reads=fig[fig.index('age_counts <-'):fig.index('technical_failure <-')]
reads=reads.replace('fread(paudit("results/table_clock_feature_set_age_information.tsv"))','copy(feature_summary)').replace('fread(paudit("results/table_clock_feature_set_matched_random_nulls.tsv"))','rbindlist(feature_null)')
plot=fig[fig.index('## Integrated quantitative-grid'):fig.index('## Reports ')]
(scripts/'06_adequacy.R').write_text(theme+'\n'+base+"pcs<-data.table(sample_id=audit$sample_id,qc$pca$x);pve<-qc$pve\np_string<-function(x)ifelse(x<.001,format(x,scientific=TRUE,digits=2),sprintf('%.3f',x))\n"+s+feature+f+reads+plot)
sp=(old/'4_final_reporting_revision/scripts/03_sparse_oldest_sensitivity.R').read_text()
ad=sp[sp.index('ad <- data.frame'):sp.index('full_rank <-')]
clock=sp[sp.index('clock_slopes <-'):sp.index('fit_spec <-')]
spec=sp[sp.index('specs<-list()'):sp.index('## Avoid duplicate')]
recompute="""
summary<-fread(file.path(old,'4_final_reporting_revision/results/table_sparse_oldest_sensitivity_summary.tsv'))
for(nm in names(specs)){ss<-specs[[nm]];cs<-clock_slopes(ss$keep,ss$age_num);summary[sensitivity==nm,`:=`(elasticnet_slope=cs[model_family=='elasticnet',slope],elasticnet_p=cs[model_family=='elasticnet',p],pca_glmnet_slope=cs[model_family=='pca_glmnet',slope],pca_glmnet_p=cs[model_family=='pca_glmnet',p])]}
write_tsv(summary,file.path(ref_tables,'table_sparse_oldest_sensitivity_summary.tsv'))
cook_rows<-list()
for(ff in c('elasticnet','pca_glmnet')){y<-audit[[paste0('existing_score_',ff)]];fit<-lm(y~age_order+rin_z+ischemia_z+mapping_z+exonic_z+duplication_z+sequencing_center,data=ad);cd<-cooks.distance(fit);cook_rows[[ff]]<-data.table(sample_id=audit$sample_id,donor_id=audit$donor_id,age_bracket=audit$age_bracket,model=ff,cooks_distance=cd,cook_flag_gt_4_over_n=cd>4/nrow(audit))}
cook_pc<-fread(file.path(old,'4_final_reporting_revision/results/table_donor_cooks_distance.tsv'))[model=='expression_PC1_20'];cooks<-rbindlist(c(cook_rows,list(cook_pc)),fill=TRUE)
write_tsv(cooks,file.path(ref_tables,'table_donor_cooks_distance.tsv'))
"""
finaltheme=(old/'4_final_reporting_revision/scripts/common_final.R').read_text()

plot=sp[sp.index('## Sensitivity figures'):sp.index('saveRDS(list(summary')]
(scripts/'07_sparse_oldest.R').write_text("theme_set(theme_context())\n"+base+ad+clock+spec+recompute+plot)
