#!/usr/bin/env Rscript
# Side-by-side plot of the BH-significant ORA terms from enrich_elasticnet_nonzero_genes.R:
# younger (coefficient < 0) on the left, older (coefficient > 0) on the right, on a shared x scale.
# Reads pathway_enrichment/ora_significant_BH0.05.tsv; nothing is re-tested.
# Run with: Rscript <this>
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

a <- commandArgs(FALSE)
script <- gsub('~+~', ' ', sub('^--file=', '', a[grepl('^--file=', a)][1]), fixed = TRUE)
fig_dir <- dirname(dirname(normalizePath(script)))
out_dir <- file.path(fig_dir, 'diagnostics/09_selected_gene_diagnostics/elasticnet_hvg_coefficients/pathway_enrichment')
q_cut <- 0.05

save_plot <- function(p, stem, width, height) {
  ggsave(paste0(stem, '.png'), p, width = width, height = height, dpi = 300, bg = 'white', limitsize = FALSE)
  ggsave(paste0(stem, '.pdf'), p, width = width, height = height, bg = 'white', limitsize = FALSE, device = cairo_pdf)
}

sig <- fread(file.path(out_dir, 'ora_significant_BH0.05.tsv'))
stopifnot(nrow(sig) > 0, all(sig$p.adjust < q_cut))
mapping <- fread(file.path(out_dir, 'gene_id_mapping_summary.tsv'))
n_univ <- mapping[set == 'universe', n_entrez]
db_lab <- c(GO_BP = 'BP', GO_CC = 'CC', GO_MF = 'MF', KEGG = 'KEGG', Reactome = 'Reactome', Hallmark = 'Hallmark')
sig[, `:=`(neg_log_q = -log10(p.adjust), db = factor(db_lab[database], levels = db_lab))]
# Grouped by database, then most significant first.
setorder(sig, direction, db, p.adjust, Description)
sig[, term := sprintf('%s (%s)', Description, db)]
n_dir <- sig[, .N, by = direction][, setNames(N, direction)]
for (d in c('older', 'younger')) if (is.na(n_dir[d])) n_dir[[d]] <- 0L
x_max <- max(sig$neg_log_q) * 1.22
dir_col <- c(older = '#F04E37', younger = '#4DB6CC')

ora_panel <- function(d, title) {
  dt <- sig[direction == d]
  dt[, term := factor(term, levels = rev(term))]
  pad <- max(n_dir) - nrow(dt)  # empty rows at the bottom keep both panels top-aligned
  ggplot(dt, aes(neg_log_q, term)) +
    geom_col(width = 0.72, fill = dir_col[[d]], colour = 'black', linewidth = 0.25) +
    geom_text(aes(label = GeneRatio), hjust = 0, nudge_x = 0.08, size = 2.1, colour = 'black') +
    geom_vline(xintercept = -log10(q_cut), linetype = 'dashed', linewidth = 0.3, colour = 'grey40') +
    scale_x_continuous(limits = c(0, x_max), expand = expansion(mult = c(0, 0.02))) +
    scale_y_discrete(expand = expansion(add = c(0.6 + pad, 0.6))) +
    labs(title = title, x = expression(-log[10]~'BH-adjusted p'), y = NULL) +
    theme_classic(base_size = 8) +
    theme(axis.text = element_text(colour = 'black'), axis.text.y = element_text(size = 7),
          axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.35),
          plot.title = element_text(size = 8.5, face = 'bold', hjust = 0.5, colour = dir_col[[d]]),
          plot.title.position = 'plot')
}
mapped <- setNames(mapping$n_entrez, mapping$set)
p <- ora_panel('younger', sprintf('Toward younger predicted age (%d terms)', n_dir[['younger']])) +
  ora_panel('older', sprintf('Toward older predicted age (%d terms)', n_dir[['older']])) +
  plot_annotation(
    title = 'Pathway over-representation of HVG–EN nonzero genes',
    subtitle = sprintf(paste0('Terms with BH q < %.2f (hypergeometric ORA; GO, KEGG, Reactome, Hallmark tested). Gene lists: ',
                              '%d younger / %d older Entrez-mapped genes;\nuniverse = %d HVG candidates. ',
                              'Labels = genes in term / genes in list; dashed line = q %.2f. BP/CC/MF = GO ontology.'),
                       q_cut, mapped[['younger']], mapped[['older']], n_univ, q_cut),
    theme = theme(plot.title = element_text(size = 10, face = 'bold'), plot.subtitle = element_text(size = 6.8)))
fwrite(sig[, !c('db', 'term')], file.path(out_dir, 'ora_significant_older_vs_younger_plot_source.tsv'), sep = '\t')
save_plot(p, file.path(out_dir, 'ora_significant_older_vs_younger'), 8.6, 1.3 + 0.17 * max(n_dir))
cat('done\n')
