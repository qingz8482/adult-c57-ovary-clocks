from pathlib import Path
import csv, hashlib, json
out=Path(__file__).resolve().parents[1]
root=out.parents[2]
old=root/'bulk-pseudobulk_clock/human_cross_species_stageA_stageB'
def read(rel):
    with (out/rel).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def f(x,d=3):
    try:return f'{float(x):.{d}g}'
    except (ValueError,TypeError):return str(x)
def table(rows,cols):
    return '| '+' | '.join(c[0] for c in cols)+' |\n| '+' | '.join('---' for _ in cols)+' |\n'+'\n'.join('| '+' | '.join(f(r[c[1]]) if c[1] in r else '' for c in cols)+' |' for r in rows)
a=read('tables/age_association_sensitivities.tsv');primary=[r for r in a if r['analysis']=='primary_adjusted_ordered']
c=read('tables/age_correlations.tsv');n=read('tables/stageB_null_benchmarks.tsv');coverage=read('tables/coverage_summary.tsv');m=read('provenance/model_manifest.tsv');v=read('provenance/validation_checks.tsv');hpa=read('tables/HPA_exploratory_direction.tsv');st=read('tables/stageA_direction_conservation.tsv');agree=read('tables/cross_model_agreement.tsv')
res='''# Current frozen mouse clocks: human generalization

Completed 2026-10-12 (Asia/Singapore). The rebutal4 locked release was tested: HVG-ElasticNet trained without single-age studies (110 mouse profiles; new) and VersionB-PCA (136 profiles; unchanged from rebutal3, whose results it reproduces exactly; see provenance/rebutal3_pca_reuse_check.tsv). Both use 3,000 input genes. The 201 human donors generated 402 frozen projections.

## Interpretation

The two clocks **disagree, so there is no convergent evidence of human age generalization**. HVG-ElasticNet shows a positive GTEx age association that persists after technical adjustment (adjusted slope 0.19 SD per age bracket, HC3 P = 8×10⁻⁵; Freedman–Lane permutation P = 0.0002; all sensitivities P < 0.001), but it is not specific to the clock's gene weights: it does not exceed the random-weight (P = 0.067) or matched-gene (P = 0.19) null signatures. VersionB-PCA does not show a positive GTEx age association. Both GSE202601 correlations are positive (HVG ρ = 0.67, exact P = 0.083; PCA ρ = 0.48), but neither exact donor-level test is significant in this eight-donor cohort. This is a transfer analysis, not validation of human chronological-age prediction. There is nevertheless partial gene-level conservation: the 414 testable genes shared by the two clocks show 56.0% direction agreement (matched-null BH q = 0.0098). Whole-clock transfer and gene-level conservation are distinct endpoints.

## GTEx primary analysis

193 unique ovarian donors; ordered public age bracket plus RIN, ischemic time and log library size. Estimates are within-cohort score SD per ordered bracket, with HC3 standard errors and t-based 95% intervals. BH correction is across the two prespecified model tests.

'''+table(primary,[('Clock','model'),('Slope','estimate'),('95% CI low','ci_low'),('95% CI high','ci_high'),('P','p'),('BH q','q_primary_two_models')])+'''

## Age correlations

GTEx uses public age brackets (midpoints preserve rank); GSE202601 uses exact ages. GSE202601 Spearman P values enumerate all 40,320 donor permutations. Pearson confidence intervals and leave-one-donor-out sensitivity are provided separately. No nuclei are treated as independent replicates.

'''+table(c,[('Cohort','cohort'),('Clock','model'),('Donors','n'),('Pearson r','pearson'),('Spearman rho','spearman'),('Spearman P','spearman_p')])+'''

## Feature and model-weight coverage

Primary projection uses every measured raw-input gene, including zero counts. The full prefilter library total is computed before orthology/model subsetting. Genes unavailable through unique one-to-one mappings receive the frozen mouse mean. An additional sensitivity excludes genes failing the historical human count filter while keeping the same full-library denominator.

'''+table(coverage,[('Cohort','cohort'),('Clock','model'),('Measured genes','n_measured'),('Input genes','n_input'),('Gene fraction','raw_coverage'),('Absolute-weight fraction','absolute_weight_coverage')])+'''

## Random-signature and permutation benchmarks

Each null has 10,000 fixed-seed draws. The observed and randomized statistics use the same technical-adjusted ordered-age slope of a within-cohort standardized score. The matched-gene benchmark forms actual donor-level random signatures, preserving inter-gene covariance. One-sided P tests enrichment in the prespecified positive-age direction; an adverse observed direction cannot establish successful transfer. Two-sided values are available in the table file. Residual-permutation inference is approximate and assumes exchangeable residuals.

'''+table(n,[('Clock','model'),('Null','null_type'),('Observed slope','observed'),('One-sided P','empirical_one_sided_p'),('BH q (two clocks per null)','q_within_null_across_two_models')])+'''

## Stage A direction conservation

Mouse gene-age effects were recomputed from the 136-profile development reference with source-batch adjustment, over the union of the two clocks' eligible genes measured in all 136 profiles (the clocks were filtered on different training sets). Human GTEx voom/limma effects are reused from the previous cohort preparation because these gene-level outcomes are independent of the mouse clock release. The matched background is the current mouse-eligible, one-to-one, human-testable universe. Null strata use expression, detection and variance, never human-age effects. Gene-set rows are exploratory, overlap heavily, and are not independent replications. Saturated top-N sets can duplicate the active-gene set; BH correction includes all reported rows and is consequently conservative.

'''+table([r for r in st if r['gene_set'] in ['testable_mouse_universe','shared_model_genes','hvg_input_genes','hvg_active_genes','pca_input_genes']],[('Gene set','gene_set'),('Genes','n_genes'),('Signed-t rho','spearman_signed_t'),('Sign concordance','sign_concordance'),('Matched-null P','empirical_p'),('BH q across sets','q')])+'''

Weighted older/younger signatures, top-absolute-weight subsets, standardized gene effects and current-feature human age information are in separate tables. Top-N sets saturate at the number of available active genes; their actual n is reported. HPA is confined to exploratory Stage A, with all-seven and adult-only analyses:

'''+table(hpa,[('Subset','subset'),('Clock','model'),('Donors','n_donors'),('Genes','n_genes'),('Mouse-HPA rho','mouse_hpa_spearman'),('Sign concordance','sign_concordance')])+'''

## Domain and agreement diagnostics

Low-dimensional Mahalanobis distances and nearest-reference distances use mouse-fitted axes only. HVG diagnostic PCA is fitted to the mouse reference solely for QC; it does not change its clock. VersionB-PCA uses its frozen ten axes. The orthogonal expression residual measures signal outside those axes. All 193 GTEx donors fall within the mouse 99th-percentile low-dimensional threshold for both clocks. GSE202601 has 8/8 within that threshold for HVG but only 2/8 for PCA. All human projections exceed the mouse 99th-percentile orthogonal-residual threshold. Being within a low-dimensional mouse threshold therefore does not establish full-expression similarity or reliable human age prediction. Missing-mean imputation can also pull projections toward the reference.

'''+table(agree,[('Cohort','cohort'),('Profiles','n'),('Pearson r','pearson'),('Spearman rho','spearman')])+'''

Mouse reference agreement is in-sample and labelled as such. The reference-style supplement separately reports current paired nested-CV recipe agreement. Human agreement is between two projections of the same samples, not an independent accuracy endpoint. Domain metrics, mouse 95th/99th-percentile support, Cook's distance, leave-one-age-bracket, leave-one-sequencing-batch and GSE202601 leave-one-donor results are exported. PCA contributions retain PC1, PC2, PC5 and PC8 as the active predictors; all ten axes are available for audit.

## Scope and limitations

- Frozen mouse outputs are in mouse weeks and are **audit-only**, not predicted human ages. No human MAE, human delta age or identity line is used. Human ages never enter the frozen prediction calculation.
- Human scores are standardized only after frozen prediction for within-cohort display. Input genes are never human-recentered or human-rescaled; PCA, coefficients and features remain frozen. Scores cannot be compared in magnitude between cohorts.
- Current model normalization is natural-log `ln(1 + CPM)` from full available prefilter counts. The historical TMM logCPM projection normalization is incompatible with these artifacts and was not reused for projection. Historical human expression filtering is a separate sensitivity, not a reason to impute a measured zero.
- GTEx chronological age is bracketed; technical and postmortem characteristics, cellular composition and donor selection may confound associations. HPA is small and includes nonadult samples in the all-seven view. GSE202601 has only eight donors with young/aged sampling and is vulnerable to individual influence.
- Sparse model weights and small matched strata limit null diversity. Matched genes are sampled per feature slot with replacement; repeated draws add their coefficients. Eligibility and matching do not use age.
- Dataset-wide GTEx adequacy and nested human-predictor audits from the previous analysis retain their original cohort-only figures unchanged; they are not rerun or relabelled as validation of the current clocks. All current-feature adequacy comparisons are recomputed. The supplementary mouse CV comparison uses paired nested-CV predictions from the current development recipes; these are not predictions from the final frozen artifacts.
- Model-level results do not establish that ageing biology is universally conserved or unconserved across species.

## Reproduction and validation

Run `bash scripts/run_all.sh` from this directory. Inputs and original models are read-only. All new results are written here. Source SHA-256 checksums, the release manifest, R session information, figure source data and numerical validation are under `provenance/` and `figures/source_data/`.

'''+f"Validation: {sum(r['pass']=='TRUE' for r in v)}/{len(v)} checks passed. Checks include source integrity, exact formula/PC replay, missing-value behavior, donor-independent frozen prediction, full-library normalization, null P-value replay, exact n=8 Spearman enumeration and figure output presence.\n"
(out/'reports/results_and_interpretation.md').write_text(res)
methods='''# Methods and differences from the previous human analyses

1. **Frozen release:** verify the locked release checksums; load the HVG-ElasticNet (single-age exclusion) and VersionB-PCA/glmnet artifacts from the rebutal4 release `../lockedmodel`. No fitting or recalibration of either clock occurs. Both representations have 3,000 genes; HVG was fitted to 110 mouse profiles (79 nonzero genes) and PCA to 136 (ten inputs, four active PCs).
2. **Inputs:** GTEx V11 ovary raw gene counts, the previous 193-donor deduplication/metadata, the existing 42,568-nucleus GSE202601 raw donor pseudobulk sums (eight donors), Ensembl release 116 high-confidence one-to-one orthology, and exploratory HPA gene TPM. The GSE202601 aggregate input is reused and its donor totals/count integerness verified; the original Seurat aggregation is not repeated.
3. **Normalization:** sum full prefilter gene counts per donor; compute ln(1 + 10^6 counts / total). GTEx is keyed by stable human Ensembl ID; GSE202601 by unique human gene symbol. Ambiguous symbols are excluded rather than selected by row order. Missing mapped features are NA before the frozen prediction function imputes mouse centers. Measured zeros remain observed zeros. Historical filtering is evaluated as a sensitivity with the same denominator.
4. **Prediction:** apply exact artifact gene order, frozen mouse centers/scales, frozen PCA rotation where applicable, and glmnet lambda. Independently reconstruct using raw coefficients/intercept and standardized coefficients/intercept. Sum PC contributions plus intercept to verify PCA predictions. The release has no calibration. Cohort z scores are presentation endpoints calculated after prediction.
5. **Stage A:** current mouse expression-age limma models adjust for source batch; use the model-eligible background. Compare standardized effects, signed t ranks and sign concordance to existing human GTEx voom effects. Test expression/detection/variance-matched nulls in 10 x 5 x 10 quantile strata, with replacement. All-gene conservation sets are unweighted to avoid assigning arbitrary weight to zero-coefficient HVGs. Separate directional compatibility tests use absolute raw coefficients for positive/negative signatures. HPA uses all seven and age >=18 subsets; leave-one-donor sign stability is recomputed.
6. **Age association:** GTEx primary is score ~ ordered age bracket + standardized RIN + ischemic time + log10 library size. HC3 covariance is computed explicitly and inference uses residual-df t statistics. Across-model BH q values are reported for the two primary tests. Sensitivities include historical RIN-only, unadjusted, midpoint, robust MASS regression, oldest-bracket exclusion, 2.5/97.5% winsorization, filter sensitivity, and leave-one-bracket/batch. GSE202601 uses exact 8! Spearman permutations, Pearson intervals, HC3 per-decade association, unpaired Welch young/aged comparison, and leave-one-donor influence. Secondary tests are descriptive unless their multiplicity family is explicitly given.
7. **Stage B nulls:** (a) permute frozen standardized weights among measured standardized features; (b) replace each active gene with a measured one-to-one gene matched by human ln1pCPM expression, detection and variance, preserving its raw coefficient and actual donor covariance; (c) Freedman-Lane permutation of score residuals from the technical-covariate model; (d) historical age-label permutation as sensitivity. Each statistic uses the same adjusted age slope and cohort score normalization, 10,000 draws, seed 20261012. One-sided positive-direction and centered two-sided empirical P values include the +1 correction. BH adjustment is across the two models within each null family. These nulls diagnose specificity; they do not fit a human clock.
8. **QC:** mouse-derived Mahalanobis thresholds, nearest-reference distances, fractions of observed standardized genes outside mouse ranges, standardized extremes, and orthogonal residual RMS complement gene/weight coverage. Mouse thresholds are descriptive reference percentiles, not validated human acceptability cutoffs. Cross-model agreement includes paired donor bootstrap intervals (2,000 draws) and current mouse in-sample reference agreement.
9. **Figures:** R ggplot2/patchwork, Cairo PDF and ragg PNG. Points are donor-level; boxplots show median/IQR and 1.5-IQR whiskers. Scatter ribbons are descriptive ordinary least-squares 95% mean-fit intervals, not HC3 primary intervals. Sensitivity/errorbar plots explicitly use computed HC3 intervals except the labelled MASS robust sensitivity. Null shading/segments are empirical quantiles. Source data accompany plots.

## Corrections required for applying the new models

The earlier Stage B scripts mixed a technical-adjusted observed slope with age-only null slopes in one benchmark, and a later matched-gene approximation used a different scale from the observed score slope. This reproduction computes the same donor-level statistic for both sides and preserves gene covariance for matched signatures. Historical numerical null P values therefore must not be transplanted. The older 2,000-gene models, 117-profile mouse reference, calibration terms and active PC identities are also replaced by the current artifacts. No interpretation from the old frozen-clock results is assumed.
'''
(out/'reports/methods.md').write_text(methods)
(out/'README.md').write_text('''# Human generalization: current frozen ovarian clocks

Completed 2026-10-12. HVG-ElasticNet (single-age exclusion) and VersionB-PCA from the rebutal4 locked release (3,000 genes each), applied to **193 GTEx + 8 GSE202601 donors**. HPA is exploratory Stage A only.

**Finding:** no convergent human-age transfer. HVG has a positive GTEx age trend robust to technical adjustment (P = 8×10⁻⁵) but not specific to its gene weights (random-weight null P = 0.067); PCA lacks a positive GTEx trend. Both GSE202601 exact Spearman tests are nonsignificant (n = 8).

- [Results and interpretation](reports/results_and_interpretation.md)
- [Methods and changes from the previous analyses](reports/methods.md)
- [Main figure](figures/human_cross_species_main.pdf) and [PNG preview](figures/human_cross_species_main.png)
- [All 402 frozen predictions](tables/frozen_predictions.tsv)
- [Age correlations](tables/age_correlations.tsv), [adjusted analyses](tables/age_association_sensitivities.tsv), [null benchmarks](tables/stageB_null_benchmarks.tsv)
- [Gene coverage](tables/coverage_summary.tsv), [Stage A](tables/stageA_direction_conservation.tsv), [HPA](tables/HPA_exploratory_direction.tsv)
- [Validation checks](provenance/validation_checks.tsv), [scope audit](provenance/reference_scope.tsv), [source checksums](provenance/source_checksums.tsv)

Run `bash scripts/run_all.sh`. The existing R environment is used; all new outputs stay in this folder. **Mouse-week outputs are audit-only, not predicted human ages.** Cohort-standardized scores are not comparable across cohorts.
''')
scope=[
('original Stage A orthology, direction, signed signatures','recomputed','tables/stageA_direction_conservation.tsv'),
('original Stage B GTEx frozen projection','recomputed using current frozen preprocessing','tables/frozen_predictions.tsv'),
('random signature and permutation benchmarks','recomputed with matched observed/null statistics','tables/stageB_null_benchmarks.tsv'),
('revision standardized effects/top-weight/shared sets','recomputed using the 136-profile mouse reference, union eligible-gene background','tables/stageA_gene_effect_master.tsv'),
('revision age/OOD/influence/support sensitivities','recomputed','tables/age_association_sensitivities.tsv'),
('HPA all-seven and adult-only direction','recomputed, exploratory only','tables/HPA_exploratory_direction.tsv'),
('PCA contribution diagnostics','recomputed for current PC1/PC2/PC5/PC8','tables/pca_scores_contributions.tsv'),
('GSE202601 projection, age plots, coverage, agreement','recomputed from existing raw donor sums','tables/frozen_predictions.tsv'),
('current-feature human age-information matched nulls','recomputed','tables/current_features_human_age_information.tsv'),
('3_gtex_age_signal_adequacy cohort-wide/pathway audits','historical context: cohort unchanged; not a current-clock result',str(old/'3_gtex_age_signal_adequacy')),
('4_final_reporting_revision cohort-wide nested prediction/centroid/Hardy audit','historical context: no human model training in this reproduction',str(old/'4_final_reporting_revision')),
('previous held-out clock agreement','replaced with paired current nested-CV recipe predictions; final-artifact reference remains in-sample','tables/cross_model_agreement.tsv')]
with (out/'provenance/reference_scope.tsv').open('w') as fobj:
    wr=csv.writer(fobj,delimiter='\t');wr.writerow(['reference_analysis','status','output_or_context']);wr.writerows(scope)
context=[]
for rel in ['3_gtex_age_signal_adequacy/reports','4_final_reporting_revision/reports']:
    for p in sorted((old/rel).glob('*.md')):context.append((str(p),hashlib.sha256(p.read_bytes()).hexdigest(),'historical cohort-only context; not recomputed'))
with (out/'provenance/historical_context_manifest.tsv').open('w') as fobj:
    wr=csv.writer(fobj,delimiter='\t');wr.writerow(['path','sha256','scope']);wr.writerows(context)
# Reference figure guide is generated from the one-to-one inventory.
manifest=read('provenance/reference_panel_manifest.tsv')
intro="""# Reference-matched figure guide

All 58 reference figure stems are retained in `figures/`, with PDF, SVG, PNG and TIFF. The original plotting code defines dimensions, panel arrangement, colours, shapes, facets and statistical displays. Twenty-six model-dependent figures use the current frozen HVG-ElasticNet and VersionB-PCA results; 32 cohort-only figures are exact copies. Values, labels for 3,000 genes / ten PCs, and data-dependent axis limits necessarily reflect the new models. The original input directories remain unchanged.

- **human_cross_species_main** and **panel_A–F**: original schematic, orthology coverage, mouse–human slope scatter, weighted directional enrichment, GTEx score-by-age and cross-model/OOD arrangement.
- **supplement_***: original six-panel summary, four null families, bracket omission, age/OOD, mouse/human agreement, HPA and gene/PC contributions.
- **GSE202601**: original raw mouse-output and standardized-score age panels, gene-coverage panels and pooled GTEx/GSE agreement. Short donor labels, young/aged colours and shapes, exact Spearman tests, and descriptive OLS ribbons are retained. Raw mouse output is audit-only; it is not a human-age estimate.
- **Adequacy / final reporting**: cohort-wide expression, pathways and human nested-prediction panels are preserved exactly. Current clock feature sets, OOD colours, Cook's distances and sparse-oldest clock slopes are updated.

GTEx n=193; GSE202601 n=8. Points represent donors. Scores are standardized within model/cohort after frozen prediction. No human fitting or calibration of either mouse clock occurs. Null distributions use 10,000 draws; the four-family reference panel additionally includes independent weight/sign permutation. Its observed and randomized statistics both use the current adjusted slope. As in the reference, a degenerate zero-variance randomized signature has score zero (and slope zero); sparse HVG weights can produce such draws. HPA weight-direction panels include active mapped genes only. Current nested mouse CV predictions describe the development recipe, not held-out evaluation of the final artifact. Contribution slopes are reconstructed in the current ln1pCPM units and checked to sum to the frozen output slope.

The sparse-oldest sensitivity panels intentionally retain the reference regression adjustment (RIN, ischemia, mapping, exonic and duplication rates, sequencing center), with ordinary-lm inference. This is a secondary analysis; the primary HC3 regression in the main results retains RIN, ischemia and log library size. Feature-set adequacy uses the original cohort-wide count-aware age statistics and current gene membership. Thus plot form matches the reference while model-dependent quantities are recomputed.

`figures/source_data/`, `tables/reference_panels/` and `inputs/reference_adapter/` provide the plotted values. The SHA-256 inventory below identifies exact-copy versus updated figures. The former figure design is archived under `provenance/pre_reference_style_20261012/`.

## One-to-one figure inventory

| Figure stem | Reference module/path | Treatment |
| --- | --- | --- |
"""
intro+='\n'.join('| '+r['output_stem']+' | '+r['reference_module']+'/'+r['reference_path']+' | '+r['status']+' |' for r in manifest)+'\n'
(out/'reports/figure_legends.md').write_text(intro)
p=out/'README.md'
p.write_text(p.read_text()+"\n## Reference-matched plots\n\nAll **58 corresponding figure stems** are in [figures](figures/): **26 updated for the current clocks**, **32 exact copies of unchanged cohort-only analyses**. Original layouts, dimensions, colours and panel structure are retained. Each has PDF, SVG, PNG and TIFF.\n\n- [Figure guide and one-to-one inventory](reports/figure_legends.md)\n- [GSE202601 standardized-score panels](figures/02_standardized_clock_score_vs_chronological_age.pdf)\n- [Six-panel supplement](figures/supplement_complete_stageA_stageB_A_to_F.pdf)\n- [Reference matching validation](provenance/reference_figure_validation.tsv)\n\nTo rebuild plots only: `bash scripts/run_reference_figures.sh`.\n")

# Machine-readable complete output inventory; hash after report generation.
with (out/'provenance/output_inventory.tsv').open('w') as fobj:
    wr=csv.writer(fobj,delimiter='\t');wr.writerow(['path','bytes','sha256'])
    for p in sorted(out.rglob('*')):
        if p.is_file() and p.name!='output_inventory.tsv':wr.writerow([str(p.relative_to(out)),p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest()])
print('PASS: report, methods, scope audit, README and output inventory written.')
