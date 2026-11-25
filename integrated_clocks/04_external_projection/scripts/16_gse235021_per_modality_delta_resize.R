#!/usr/bin/env Rscript
# Redraw the GSE235021 bulk and pseudobulk raw delta-age plots with the panel layout and
# page size of the mixed-modality relative plot (script 15). Values and Welch tests are
# read unchanged from the script 07 outputs.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
setDTthreads(1)
args <- commandArgs(FALSE)
script <- gsub('~+~', ' ', sub('^--file=', '', args[grepl('^--file=', args)][1]), fixed=TRUE)
out <- dirname(dirname(normalizePath(script)))
rd <- function(f) fread(file.path(out, f))
models <- c('HVG-ElasticNet', 'VersionB-PCA/glmnet')
model_labels <- c('HVG-ElasticNet'='Elasticnet', 'VersionB-PCA/glmnet'='PCA/glmnet')
groups <- c('6m polished', '6m black', '6m brown', '9m polished', '9m black', '9m brown')
cohorts <- c(bulk_rnaseq='bulk RNA-seq', scrna_pseudobulk='scRNA-derived pseudobulk')
src <- rd('provenance/complete_plot_source_data.tsv')
tests <- rd('GSE235021/02_gse235021_raw_delta_by_group_ttest.tsv')
format_p <- function(p) ifelse(p<1e-4, '****', ifelse(p<.001, '***', ifelse(p<.01, '**', ifelse(p<.05, '*', sprintf('p = %.3f', p)))))
cols <- c(diet='#E69F00', 'reference control'='grey50')

for (co in names(cohorts)) {
  stem <- paste0('GSE235021/plots/GSE235021__', co, '_all_groups_delta_age')
  d <- unique(src[figure_stem==stem & value_type=='delta_age_week', .(sample_key, model, plot_group, value)])
  stopifnot(!anyDuplicated(d[, .(sample_key, model)]), setequal(d$plot_group, groups),
            setequal(d$model, models))
  d[, `:=`(model=factor(model, levels=models), plot_group=factor(plot_group, levels=groups),
           plot_class=factor(ifelse(grepl('polished$', plot_group), 'reference control', 'diet'),
                             levels=c('diet', 'reference control')))]
  setorder(d, model, plot_group, sample_key)
  set.seed(310)
  d[, point_x:=as.integer(plot_group)+runif(.N, -.12, .12)]

  tt <- tests[modality==cohorts[[co]], .(model, Reference, Target, welch_p)]
  stopifnot(nrow(tt)==8L)
  tt[, `:=`(model=factor(model, levels=models), annotation=format_p(welch_p),
            xmin=match(Reference, groups), xmax=match(Target, groups),
            tier=ifelse(grepl('black$', Target), 1L, 2L))]
  rng <- d[, .(hi=max(value), span=diff(range(c(0, value)))), by=model]
  tt <- merge(tt, rng, by='model')
  tt[, `:=`(y=hi+span*(.12+.16*(tier-1)), tip=hi+span*(.09+.16*(tier-1)))]
  tt[, label_y:=y+span*.015]
  d[, model_label:=factor(model_labels[as.character(model)], levels=model_labels)]
  tt[, model_label:=factor(model_labels[as.character(model)], levels=model_labels)]

  # Circles for bulk, triangles for pseudobulk, as in the mixed-modality plots.
  pt_shape <- if (co=='scrna_pseudobulk') 17L else 16L
  g <- ggplot(d, aes(as.integer(plot_group), value, colour=plot_class)) +
    geom_hline(yintercept=0, linetype='dashed', colour='grey45') +
    geom_boxplot(aes(group=plot_group), width=.55, outlier.shape=NA, linewidth=.6) +
    geom_point(aes(x=point_x), shape=pt_shape, size=2.2, alpha=.9) +
    geom_segment(data=tt, aes(x=xmin, xend=xmax, y=y, yend=y), inherit.aes=FALSE, linewidth=.45) +
    geom_segment(data=tt, aes(x=xmin, xend=xmin, y=tip, yend=y), inherit.aes=FALSE, linewidth=.45) +
    geom_segment(data=tt, aes(x=xmax, xend=xmax, y=tip, yend=y), inherit.aes=FALSE, linewidth=.45) +
    geom_text(data=tt, aes(x=ifelse(grepl('^\\*', annotation), (xmin+xmax)/2, xmin),
                           y=label_y, label=annotation,
                           hjust=ifelse(grepl('^\\*', annotation), .5, 0)),
              inherit.aes=FALSE, size=4.2, vjust=0) +
    facet_wrap(~model_label, nrow=1, scales='free_y') +
    scale_colour_manual(values=cols, name='Projection class') +
    scale_x_continuous(breaks=seq_along(groups), labels=groups, expand=expansion(add=.5)) +
    scale_y_continuous(expand=expansion(mult=c(.05, .06))) +
    labs(x=NULL, y='Raw delta age (weeks)') +
    theme_bw(base_size=13) +
    theme(panel.grid.minor=element_blank(), strip.text=element_text(face='bold', size=13),
          strip.background=element_rect(fill='grey95'), axis.title.y=element_text(face='bold'),
          axis.text.x=element_text(angle=30, hjust=1, size=12, colour='black'),
          axis.text.y=element_text(face='bold', size=12, colour='black'),
          legend.title=element_text(face='bold'), legend.position='right',
          legend.key.height=unit(1.1, 'lines')) +
    guides(colour=guide_legend(order=1, override.aes=list(shape=pt_shape)))
  # Same page size as script 15.
  width <- 11; height <- 4.6
  ggsave(file.path(out, paste0(stem, '.pdf')), g, width=width, height=height, device=cairo_pdf)
  ggsave(file.path(out, paste0(stem, '.png')), g, width=width, height=height, dpi=300, bg='white')
  print(tt[, .(model, Reference, Target, welch_p, annotation)])
}
message('PASS: GSE235021 per-modality raw delta plots redrawn at 11 x 4.6 in.')
