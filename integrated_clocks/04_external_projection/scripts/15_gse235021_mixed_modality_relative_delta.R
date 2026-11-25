#!/usr/bin/env Rscript
# Mixed-modality relative delta age for GSE235021: one box per group pooling bulk and
# pseudobulk, after centring each modality on its own age-matched polished-rice controls.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
setDTthreads(1)
args <- commandArgs(FALSE)
script <- gsub('~+~', ' ', sub('^--file=', '', args[grepl('^--file=', args)][1]), fixed=TRUE)
out <- dirname(dirname(normalizePath(script)))
rd <- function(f) fread(file.path(out, f))
wr <- function(x, f) fwrite(x, file.path(out, f), sep='\t', na='NA')
stem <- 'GSE235021/plots/GSE235021__mixed_modality_relative_delta_age'
tab <- 'GSE235021/tables/mixed_modality_relative_delta_age_'
models <- c('HVG-ElasticNet', 'VersionB-PCA/glmnet')
model_labels <- c('HVG-ElasticNet'='Elasticnet', 'VersionB-PCA/glmnet'='PCA/glmnet')
groups <- c('6m polished', '6m black', '6m brown', '9m polished', '9m black', '9m brown')
modalities <- c('bulk RNA-seq', 'scRNA-derived pseudobulk')

r <- rd('GSE235021/gse235021_relative_to_controls.tsv')
# Reference rows repeat once per contrast; keep one row per sample and clock.
d <- unique(r[, .(sample_key, model, modality, plot_group=new_group, plot_class=new_class,
                  delta_age_week, control_mean_delta, relative_delta_age_week)])
stopifnot(!anyDuplicated(d[, .(sample_key, model)]), uniqueN(d$sample_key)==35L, nrow(d)==70L,
          setequal(d$plot_group, groups), setequal(d$modality, modalities))
# Recompute the centring to confirm it is per modality x age against polished controls.
d[, age:=sub(' .*', '', plot_group)]
cm <- d[grepl('polished$', plot_group), .(control_mean_check=mean(delta_age_week)), by=.(model, modality, age)]
d <- merge(d, cm, by=c('model', 'modality', 'age'))
stopifnot(max(abs(d$control_mean_check-d$control_mean_delta))<1e-10,
          max(abs(d$delta_age_week-d$control_mean_delta-d$relative_delta_age_week))<1e-10)
d[, `:=`(model=factor(model, levels=models), plot_group=factor(plot_group, levels=groups),
         modality=factor(modality, levels=modalities),
         plot_class=factor(ifelse(grepl('polished$', plot_group), 'reference control', 'diet'),
                           levels=c('diet', 'reference control')),
         value=relative_delta_age_week)]
setorder(d, model, plot_group, modality, sample_key)
set.seed(310)
d[, point_x:=as.integer(plot_group)+runif(.N, -.12, .12)]

# Pooled-modality Welch t-tests on relative delta age, diet vs age-matched polished.
cmp <- data.table(Reference=rep(c('6m polished', '9m polished'), each=2),
                  Target=c('6m black', '6m brown', '9m black', '9m brown'))
tt <- rbindlist(lapply(models, function(m) rbindlist(lapply(seq_len(nrow(cmp)), function(i) {
  a <- d[model==m & plot_group==cmp$Reference[i], value]
  b <- d[model==m & plot_group==cmp$Target[i], value]
  z <- t.test(b, a)
  data.table(model=m, Reference=cmp$Reference[i], Target=cmp$Target[i],
             n_Reference=length(a), n_Target=length(b), mean_Reference=mean(a), mean_Target=mean(b),
             effect_weeks=mean(b)-mean(a), ci95_low=z$conf.int[1], ci95_high=z$conf.int[2],
             t_statistic=unname(z$statistic), df=unname(z$parameter), welch_p=z$p.value,
             test='two-sided unpaired Welch t-test on relative delta age, bulk and pseudobulk pooled')
}))))
format_p <- function(p) ifelse(p<1e-4, '****', ifelse(p<.001, '***', ifelse(p<.01, '**', ifelse(p<.05, '*', sprintf('p = %.3f', p)))))
tt[, `:=`(model=factor(model, levels=models), annotation=format_p(welch_p),
          xmin=match(Reference, groups), xmax=match(Target, groups),
          tier=ifelse(grepl('black$', Target), 1L, 2L))]
rng <- d[, .(hi=max(value), span=diff(range(c(0, value)))), by=model]
tt <- merge(tt, rng, by='model')
tt[, `:=`(y=hi+span*(.12+.16*(tier-1)), tip=hi+span*(.09+.16*(tier-1)))]
tt[, label_y:=y+span*.015]
d[, model_label:=factor(model_labels[as.character(model)], levels=model_labels)]
tt[, model_label:=factor(model_labels[as.character(model)], levels=model_labels)]

cols <- c(diet='#E69F00', 'reference control'='grey50')
g <- ggplot(d, aes(as.integer(plot_group), value, colour=plot_class)) +
  geom_hline(yintercept=0, linetype='dashed', colour='grey45') +
  geom_boxplot(aes(group=plot_group), width=.55, outlier.shape=NA, linewidth=.6) +
  geom_point(aes(x=point_x, shape=modality), size=2.2, alpha=.9) +
  geom_segment(data=tt, aes(x=xmin, xend=xmax, y=y, yend=y), inherit.aes=FALSE, linewidth=.45) +
  geom_segment(data=tt, aes(x=xmin, xend=xmin, y=tip, yend=y), inherit.aes=FALSE, linewidth=.45) +
  geom_segment(data=tt, aes(x=xmax, xend=xmax, y=tip, yend=y), inherit.aes=FALSE, linewidth=.45) +
  geom_text(data=tt, aes(x=ifelse(annotation=='*' | grepl('^\\*', annotation), (xmin+xmax)/2, xmin),
                         y=label_y, label=annotation,
                         hjust=ifelse(grepl('^\\*', annotation), .5, 0)),
            inherit.aes=FALSE, size=4.2, vjust=0) +
  facet_wrap(~model_label, nrow=1, scales='free_y') +
  scale_colour_manual(values=cols, name='Projection class') +
  scale_shape_manual(values=setNames(c(16, 17), modalities),
                     labels=c('Bulk-seq', 'ScRNA-seq derived\npseudobulk'), name='Modality') +
  scale_x_continuous(breaks=seq_along(groups), labels=groups, expand=expansion(add=.5)) +
  scale_y_continuous(expand=expansion(mult=c(.05, .06))) +
  labs(x=NULL, y='Relative delta age (weeks)') +
  theme_bw(base_size=13) +
  theme(panel.grid.minor=element_blank(), strip.text=element_text(face='bold', size=13),
        strip.background=element_rect(fill='grey95'), axis.title.y=element_text(face='bold'),
        axis.text.x=element_text(angle=30, hjust=1, size=12, colour='black'),
        axis.text.y=element_text(face='bold', size=12, colour='black'),
        legend.title=element_text(face='bold'), legend.position='right',
        legend.key.height=unit(1.1, 'lines')) +
  guides(colour=guide_legend(order=1, override.aes=list(shape=16)),
         shape=guide_legend(order=2, override.aes=list(colour='black', alpha=1, size=2.5)))
width <- 11; height <- 4.6
ggsave(file.path(out, paste0(stem, '.pdf')), g, width=width, height=height, device=cairo_pdf)
ggsave(file.path(out, paste0(stem, '.png')), g, width=width, height=height, dpi=300, bg='white')
wr(d[, .(sample_key, model, modality, plot_group, plot_class, delta_age_week, control_mean_delta,
         relative_delta_age_week)], paste0(tab, 'source.tsv'))
wr(tt[, !c('hi', 'span', 'y', 'tip', 'label_y', 'tier', 'model_label')], paste0(tab, 'tests.tsv'))
print(tt[, .(model, Reference, Target, n_Reference, n_Target, effect_weeks, welch_p, annotation)])
message('PASS: mixed-modality relative delta plot; 35 profiles, 70 points, 8 pooled Welch tests.')
