# Figure assembly

Most panels are drawn by the analysis steps under `integrated_clocks/` and `metacell_clocks/`. The scripts here assemble figures from those saved outputs without refitting any model. Each checks that its inputs come from the locked release `hvg_exclusion_vbpca_20261012`. The complete panel-to-script map is in [`../docs/figure_code_map.tsv`](../docs/figure_code_map.tsv).

| Folder | Figure | Entry point |
|---|---|---|
| `figure1_supplementary_figure1/` | Fig. 1, Supplementary Fig. 1 | `scripts/make_figures.R` (main panels), `make_sfig1_v2.R` (final Supplementary Fig. 1 layout), `make_fig1d_tall.R`. Gene–age concordance, overlap and diagnostics come from the other `make_*.R` scripts. `enrich_elasticnet_nonzero_genes.R` and `plot_ora_significant.R` run the EN non-zero-gene ORA; `enrich_elasticnet_nonzero_genes.R` needs clusterProfiler and org.Mm.eg.db. |
| `supplementary_figure2/` | Supplementary Fig. 2 | `scripts/make_sfig2.R`, re-plotted from `integrated_clocks/03_pca_interpretation/supple2_PCA` outputs |
| `supplementary_figure3/` | Supplementary Fig. 3, plus the contrast-effect and agreement panels of Fig. 3 | `scripts/make_sfig3.R`, from `integrated_clocks/05_tage_comparison` outputs |
| `figure4e/` | Fig. 4e | `make_fig4e_update.R`: internal out-of-fold metacell-CV performance from the v2 metacell clocks |

`legend_notes.md` in each folder records the panel order and the numbers quoted in the legend at the time the figure was generated.
