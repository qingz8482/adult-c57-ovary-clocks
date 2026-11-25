#!/usr/bin/env Rscript
# Pathway over-representation of all nonzero coefficient genes of the locked HVG-EN clock
# (HVG-ElasticNet-exclude-single-age), split by sign:
#   older   coefficient > 0  (higher expression -> higher predicted age)
#   younger coefficient < 0  (higher expression -> lower predicted age)
# Universe = the 3000 HVG candidate genes the model could select from (coefficients.csv), not the genome.
# Databases: GO BP/MF/CC (org.Mm.eg.db), KEGG (mmu, online), Reactome and Hallmark (msigdbr, mouse-native MSigDB).
# Run with: Rscript <this>  (scRNA_r4 lacks clusterProfiler/org.Mm.eg.db)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(clusterProfiler); library(org.Mm.eg.db); library(msigdbr)
})

a <- commandArgs(FALSE)
script <- gsub('~+~', ' ', sub('^--file=', '', a[grepl('^--file=', a)][1]), fixed = TRUE)
fig_dir <- dirname(dirname(normalizePath(script)))
lock <- file.path(dirname(fig_dir), 'lockedmodel')
en_dir <- file.path(fig_dir, 'diagnostics/09_selected_gene_diagnostics/elasticnet_hvg_coefficients')
out_dir <- file.path(en_dir, 'pathway_enrichment')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

min_gs <- 5; max_gs <- 500; q_cut <- 0.05

save_plot <- function(p, stem, width, height) {
  ggsave(paste0(stem, '.png'), p, width = width, height = height, dpi = 300, bg = 'white', limitsize = FALSE)
  ggsave(paste0(stem, '.pdf'), p, width = width, height = height, bg = 'white', limitsize = FALSE, device = cairo_pdf)
}

# ---------------------------------------------------------------- gene lists
en_coef <- fread(file.path(lock, 'HVG-ElasticNet-exclude-single-age/coefficients.csv'))
stopifnot(nrow(en_coef) == 3000, !anyDuplicated(en_coef$gene))
id_map <- suppressMessages(bitr(en_coef$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Mm.eg.db))
id_map <- as.data.table(id_map)[!duplicated(SYMBOL)]
en_coef <- merge(en_coef, id_map, by.x = 'gene', by.y = 'SYMBOL', all.x = TRUE, sort = FALSE)

nz <- en_coef[coefficient != 0][order(-abs(coefficient), gene)]
nz[, `:=`(rank = .I, direction = fifelse(coefficient > 0, 'older', 'younger'))]
setcolorder(nz, c('rank', 'gene', 'ENTREZID', 'direction', 'coefficient', 'standardized_coefficient'))
fwrite(nz, file.path(out_dir, 'all_nonzero_elasticnet_hvg_coefficients.tsv'), sep = '\t')

universe <- en_coef[!is.na(ENTREZID), ENTREZID]
lists <- list(older = nz[direction == 'older' & !is.na(ENTREZID), ENTREZID],
              younger = nz[direction == 'younger' & !is.na(ENTREZID), ENTREZID])
mapping <- data.table(set = c('universe', 'older', 'younger'),
                      n_genes = c(nrow(en_coef), nz[direction == 'older', .N], nz[direction == 'younger', .N]),
                      n_entrez = c(length(universe), lengths(lists)))
print(mapping)
fwrite(mapping, file.path(out_dir, 'gene_id_mapping_summary.tsv'), sep = '\t')

# ---------------------------------------------------------------- gene sets
msig_t2g <- function(...) {
  m <- as.data.table(msigdbr(species = 'Mus musculus', db_species = 'MM', ...))
  unique(m[, .(gs_name, ENTREZID = as.character(ncbi_gene))])
}
t2g <- list(Reactome = msig_t2g(collection = 'M2', subcollection = 'CP:REACTOME'),
            Hallmark = msig_t2g(collection = 'MH'))

run_ora <- function(genes, db) {
  args <- list(gene = genes, universe = universe, minGSSize = min_gs, maxGSSize = max_gs,
               pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = 'BH')
  res <- switch(db,
    GO_BP = do.call(enrichGO, c(args, OrgDb = org.Mm.eg.db, ont = 'BP', readable = TRUE)),
    GO_MF = do.call(enrichGO, c(args, OrgDb = org.Mm.eg.db, ont = 'MF', readable = TRUE)),
    GO_CC = do.call(enrichGO, c(args, OrgDb = org.Mm.eg.db, ont = 'CC', readable = TRUE)),
    KEGG  = tryCatch(setReadable(do.call(enrichKEGG, c(args, organism = 'mmu')), org.Mm.eg.db, 'ENTREZID'),
                     error = function(e) { message('KEGG skipped: ', conditionMessage(e)); NULL }),
    do.call(enricher, c(args, TERM2GENE = list(t2g[[db]]))))
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
  if (db %in% names(t2g)) res <- setReadable(res, org.Mm.eg.db, 'ENTREZID')
  as.data.table(as.data.frame(res))[, database := db]
}

dbs <- c('GO_BP', 'GO_MF', 'GO_CC', 'KEGG', 'Reactome', 'Hallmark')
ora <- rbindlist(lapply(names(lists), function(d) rbindlist(lapply(dbs, function(db)
  run_ora(lists[[d]], db)), fill = TRUE)[, direction := d]), fill = TRUE)
ora[, c('k', 'n') := tstrsplit(GeneRatio, '/', type.convert = TRUE)]
ora[, c('K', 'N') := tstrsplit(BgRatio, '/', type.convert = TRUE)]
ora[, fold_enrichment := (k / n) / (K / N)]
setorder(ora, direction, database, pvalue)
keep <- intersect(c('direction', 'database', 'ID', 'Description', 'GeneRatio', 'BgRatio', 'fold_enrichment',
                    'pvalue', 'p.adjust', 'qvalue', 'Count', 'geneID'), names(ora))
ora <- ora[, ..keep]
fwrite(ora, file.path(out_dir, 'ora_all_terms.tsv'), sep = '\t')
fwrite(ora[p.adjust < q_cut], file.path(out_dir, 'ora_significant_BH0.05.tsv'), sep = '\t')
summ <- ora[, .(tested_terms = .N, BH_lt_0.05 = sum(p.adjust < q_cut), nominal_p_lt_0.01 = sum(pvalue < 0.01),
                top_term = Description[1], top_p = pvalue[1], top_BH = p.adjust[1]), by = .(direction, database)]
print(summ)
fwrite(summ, file.path(out_dir, 'ora_summary.tsv'), sep = '\t')

# ---------------------------------------------------------------- plot: top terms per direction x database
top <- ora[Count >= 2][order(pvalue), head(.SD, 5), by = .(direction, database)]
top[, label := sprintf('%s [%s]', substr(Description, 1, 60), database)]
top[, key := factor(paste(direction, label, sep = '||'), levels = paste(direction, label, sep = '||')[order(direction, -pvalue)])]
top[, sig := fifelse(p.adjust < q_cut, 'BH q < 0.05', 'BH q ≥ 0.05')]
p <- ggplot(top, aes(-log10(pvalue), key)) +
  geom_segment(aes(x = 0, xend = -log10(pvalue), yend = key), colour = 'grey70') +
  geom_point(aes(size = Count, fill = sig), shape = 21, colour = 'black') +
  facet_grid(direction ~ ., scales = 'free_y', space = 'free_y',
             labeller = labeller(direction = c(older = sprintf('Older (coef > 0, n = %d)', length(lists$older)),
                                               younger = sprintf('Younger (coef < 0, n = %d)', length(lists$younger))))) +
  scale_fill_manual(values = c('BH q < 0.05' = '#B2182B', 'BH q ≥ 0.05' = 'grey80'), name = NULL) +
  scale_y_discrete(labels = function(x) sub('^.*\\|\\|', '', x)) +
  scale_size_continuous(range = c(1.5, 4.5), breaks = c(2, 4, 6, 8)) +
  labs(x = expression(-log[10]~'nominal p (hypergeometric)'), y = NULL,
       title = 'HVG-EN nonzero genes: pathway over-representation',
       subtitle = sprintf('Top 5 terms (≥ 2 genes) per database; universe = %d HVG candidates', length(universe))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = 'bold'), plot.title.position = 'plot', strip.text = element_text(face = 'bold'),
        legend.position = 'bottom', panel.grid.minor = element_blank())
save_plot(p, file.path(out_dir, 'ora_top_terms_older_vs_younger'), 9, 0.5 + 0.17 * nrow(top) + 1.5)
writeLines(capture.output(sessionInfo()), file.path(out_dir, 'sessionInfo.txt'))
