#!/usr/bin/env Rscript
# Overlay both modalities while preserving their original within-modality tests.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
setDTthreads(1)
args <- commandArgs(FALSE)
script <- gsub('~+~', ' ', sub('^--file=', '', args[grepl('^--file=', args)][1]), fixed=TRUE)
out <- dirname(dirname(normalizePath(script)))
rd <- function(f) fread(file.path(out, f))
wr <- function(x, f) fwrite(x, file.path(out, f), sep='\t', na='NA')
stem <- 'GSE235021/plots/GSE235021__mixed_modality_all_groups_delta_age'
tab <- 'GSE235021/tables/mixed_modality_delta_age_'
models <- c('HVG-ElasticNet', 'VersionB-PCA/glmnet')
groups <- c('6m polished', '6m black', '6m brown', '9m polished', '9m black', '9m brown')
modalities <- c('bulk RNA-seq', 'scRNA-derived pseudobulk')
protected <- c('external_projection_predictions.tsv', 'tables/external_projection_contrast_statistics.tsv',
               'tables/external_projection_relative_to_controls.tsv', 'inputs/model_gene_expression.rds')
before <- tools::md5sum(file.path(out, protected))
source_stems <- paste0('GSE235021/plots/GSE235021__', c('bulk_rnaseq', 'scrna_pseudobulk'), '_all_groups_delta_age')
d <- rd('provenance/complete_plot_source_data.tsv')[figure_stem %in% source_stems]
stopifnot(nrow(d)==70L, uniqueN(d$sample_key)==35L, !anyDuplicated(d[,.(sample_key, model)]),
          all(d$value_type=='delta_age_week'), setequal(d$plot_group, groups))
meta <- unique(rd('external_projection_predictions.tsv')[dataset_folder=='GSE235021',.(sample_key, modality)])
d <- merge(d, meta, by='sample_key')
d[, `:=`(model=factor(model, levels=models), plot_group=factor(plot_group, levels=groups),
         modality=factor(modality, levels=modalities))]
stopifnot(!anyNA(d$modality), !anyNA(d$model), !anyNA(d$plot_group))
d[, `:=`(modality_label=ifelse(modality==modalities[1], 'Bulk', 'Pseudobulk'),
         point_shape=ifelse(modality==modalities[1], 16L, 17L),
         x=as.integer(plot_group)+ifelse(modality==modalities[1], -.17, .17))]
setorder(d, model, plot_group, modality, sample_key)
set.seed(310)
d[, point_x:=x+runif(.N, -.055, .055)]
tt <- rd('GSE235021/02_gse235021_raw_delta_by_group_ttest.tsv')
stopifnot(nrow(tt)==16L, !anyDuplicated(tt[,.(analysis_contrast_id, model)]))
tt[, `:=`(model=factor(model, levels=models), modality=factor(modality, levels=modalities))]
tt[, `:=`(modality_label=ifelse(modality==modalities[1], 'Bulk', 'PB'),
          xmin=match(Reference, groups)+ifelse(modality==modalities[1], -.17, .17),
          xmax=match(Target, groups)+ifelse(modality==modalities[1], -.17, .17))]
format_p <- function(p) ifelse(p<1e-4, '****', ifelse(p<.001, '***', ifelse(p<.01, '**', ifelse(p<.05, '*', sprintf('p = %.3f', p)))))
stopifnot(identical(tt$annotation, format_p(tt$welch_p)))
replay <- rbindlist(lapply(seq_len(nrow(tt)), function(i) {
  s <- tt[i]
  ref_keys <- strsplit(s$reference_sample_keys, ';', fixed=TRUE)[[1]]
  target_keys <- strsplit(s$target_sample_keys, ';', fixed=TRUE)[[1]]
  a <- d[model==s$model & sample_key %in% ref_keys]
  b <- d[model==s$model & sample_key %in% target_keys]
  stopifnot(nrow(a)==s$n_control, nrow(b)==s$n_experimental,
            all(a$modality==s$modality), all(b$modality==s$modality),
            all(a$plot_group==s$Reference), all(b$plot_group==s$Target),
            all(abs(a$x-s$xmin)<1e-12), all(abs(b$x-s$xmax)<1e-12))
  z <- t.test(b$value, a$value)
  data.table(analysis_contrast_id=s$analysis_contrast_id, model=s$model,
             p_error=abs(z$p.value-s$welch_p), t_error=abs(unname(z$statistic)-s$t_statistic),
             df_error=abs(unname(z$parameter)-s$df), membership_verified=TRUE)
}))
stopifnot(max(replay$p_error)<1e-10, max(replay$t_error)<1e-9, max(replay$df_error)<1e-9)
# Two comparisons per age, each with separate bulk and pseudobulk annotation tiers.
tt[, tier:=ifelse(grepl('black$', Target), 1L, 3L)+ifelse(modality==modalities[1], 0L, 1L)]
hi <- max(d$value); span <- diff(range(c(0, d$value)))
tt[, `:=`(y=hi+span*(.10+.14*(tier-1)), tip=hi+span*(.078+.14*(tier-1)),
          label_y=hi+span*(.115+.14*(tier-1)), label=paste0(modality_label, ': ', annotation))]
counts <- unique(d[,.(sample_key, plot_group, modality)])[,.(n=.N),by=.(plot_group, modality)]
ns <- dcast(counts, plot_group~modality, value.var='n')
ns <- ns[match(groups, as.character(plot_group))]
xlabels <- paste0(groups, '\n(n = ', ns[[modalities[1]]], '/', ns[[modalities[2]]], ')')
g <- ggplot(d, aes(x, value, colour=plot_class)) +
  geom_hline(yintercept=0, linetype='dashed', colour='grey45') +
  geom_boxplot(aes(group=interaction(plot_group, modality)), width=.28, outlier.shape=NA, linewidth=.45) +
  geom_point(aes(x=point_x, shape=modality), size=2.1, alpha=.9) +
  facet_wrap(~model, nrow=1) +
  scale_colour_manual(values=c(diet='#E69F00', reference_control='#2A9D8F')) +
  scale_shape_manual(values=setNames(c(16,17), modalities), labels=c('Bulk RNA-seq (Bulk)', 'scRNA pseudobulk (PB)')) +
  scale_x_continuous(breaks=seq_along(groups), labels=xlabels, limits=c(.4,6.6), expand=expansion(mult=0)) +
  scale_y_continuous(expand=expansion(mult=c(.06,.09))) +
  geom_segment(data=tt, aes(x=xmin,xend=xmax,y=y,yend=y), inherit.aes=FALSE, linewidth=.35) +
  geom_segment(data=tt, aes(x=xmin,xend=xmin,y=tip,yend=y), inherit.aes=FALSE, linewidth=.35) +
  geom_segment(data=tt, aes(x=xmax,xend=xmax,y=tip,yend=y), inherit.aes=FALSE, linewidth=.35) +
  geom_text(data=tt, aes(x=(xmin+xmax)/2,y=label_y,label=label), inherit.aes=FALSE,
            fontface='bold', size=3.0, vjust=0) +
  labs(title='Raw delta age by projection group: GSE235021',
       subtitle='Mixed-modality display | Two frozen clocks | Circles: bulk; triangles: pseudobulk',
       x='Projection group', y='Raw delta age (weeks)', colour='Projection class', shape='Modality',
       caption=paste0('Boxes: median, IQR, whiskers to 1.5 × IQR, separately by modality. Points: sample profiles; n = bulk/pseudobulk.\n',
                      'Brackets: within-modality two-sided Welch t-tests, unadjusted p; Bulk = bulk RNA-seq, PB = scRNA pseudobulk.\n',
                      '* p < 0.05; ** < 0.01; *** < 0.001; **** < 0.0001; otherwise numerical p. Delta age = predicted age − chronological age.')) +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(), plot.title=element_text(face='bold'),
        strip.text=element_text(face='bold', size=10), legend.position='bottom', legend.box='vertical',
        axis.text.x=element_text(angle=35,hjust=1,size=10), plot.margin=margin(12,18,12,12),
        plot.caption=element_text(hjust=0,size=8)) +
  guides(colour=guide_legend(order=1,override.aes=list(shape=16)),
         shape=guide_legend(order=2,override.aes=list(colour='grey25',alpha=1,size=2.5)))
width <- 10.2; height <- 7.2
ggsave(file.path(out,paste0(stem,'.pdf')),g,width=width,height=height,device=cairo_pdf)
ggsave(file.path(out,paste0(stem,'.png')),g,width=width,height=height,dpi=300,bg='white')
wr(d,paste0(tab,'source.tsv')); wr(tt,paste0(tab,'annotations.tsv'))
wr(counts,paste0(tab,'sample_counts.tsv')); wr(replay,paste0(tab,'validation.tsv'))
wr(data.table(stem=stem,width_inches=width,height_inches=height,dpi=300,n_profiles=35,n_points=70,n_tests=16),
   paste0(tab,'export_manifest.tsv'))
stopifnot(identical(before,tools::md5sum(file.path(out,protected))))
writeLines(c('# GSE235021 mixed-modality delta age', '',
  'The two panels show the same frozen-clock projections as the separate bulk and pseudobulk figures, together on shared axes.',
  'Circles denote bulk RNA-seq (18 profiles); triangles denote scRNA-derived pseudobulk (17 profiles). Boxes are side by side within each projection group, bulk on the left and pseudobulk on the right. Colour denotes projection class.',
  'Groups remain in this order: 6m polished, 6m black, 6m brown, 9m polished, 9m black, 9m brown. Axis sample counts are bulk/pseudobulk; the 9m polished group has 3/2 profiles and the other groups each have 3/3.',
  'Boxes show the median, interquartile range and whiskers extending to observations within 1.5 IQR. All sample profiles are shown. Delta age is predicted age minus chronological age, in weeks.',
  'Each bracket retains its original within-modality two-sided, unpaired Welch t-test; Bulk and PB identify the modality. Stars denote raw P <0.05, <0.01, <0.001 and <0.0001; nonsignificant comparisons show numeric P. Exact P, t, df, confidence intervals and existing adjusted q-values are in the annotation table. No pooled-modality or between-modality test is introduced.',
  'Numerical validation: 35 unique profiles, 105 sample-by-clock values, 24 independently replayed tests, exact modality-specific bracket membership and unchanged protected numerical inputs.'),
  file.path(out,paste0(tab,'caption.md')))
message('PASS: mixed-modality GSE235021 plot; 35 profiles, 70 points, 16 existing within-modality tests; frozen numerical inputs unchanged.')
