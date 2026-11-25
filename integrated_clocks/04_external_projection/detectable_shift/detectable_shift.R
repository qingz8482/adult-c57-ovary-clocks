#!/usr/bin/env Rscript
# Minimum detectable delta-age shift for the locked-clock matched-control contrasts (reviewer 2, point 7).
# Spread: SD of delta age within each distinct control group of the prespecified, estimated, non-placement
# contrasts. Detectable shift: two-sided two-sample t-test, alpha = 0.05, power = 0.8, equal n per group,
# using each clock's median control SD. Reads saved application outputs only.
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(FALSE); here <- dirname(normalizePath(sub('^--file=', '', a[grepl('^--file=', a)][1])))
app <- dirname(here)
cs <- fread(file.path(app, 'external_projection_contrast_statistics.tsv'))[status == 'estimated' & placement_only == FALSE]
pr <- fread(file.path(app, 'external_projection_predictions_all.tsv'))
g <- unique(cs[, .(model, reference_sample_keys, n_control)])
g[, control_sd := sapply(seq_len(.N), function(i) sd(pr[model == g$model[i] & sample_key %in% strsplit(reference_sample_keys[i], ';')[[1]], delta_age_week]))]
spread <- g[, .(n_control_groups = .N, median_n_control = median(n_control), median_control_sd = median(control_sd, na.rm = TRUE),
                q25_sd = quantile(control_sd, .25, na.rm = TRUE), q75_sd = quantile(control_sd, .75, na.rm = TRUE)), by = model]
mde <- spread[, .(n_per_group = c(3, 4, 5, 10)), by = .(model, median_control_sd)]
mde[, detectable_shift_weeks := mapply(function(n, s) power.t.test(n = n, sd = s, sig.level = 0.05, power = 0.8)$delta, n_per_group, median_control_sd)]
fwrite(g[, .(model, n_control, control_sd)], file.path(here, 'control_group_delta_sd.tsv'), sep = '\t')
fwrite(spread, file.path(here, 'control_spread_summary.tsv'), sep = '\t')
fwrite(mde, file.path(here, 'detectable_shift.tsv'), sep = '\t')
print(spread); print(mde)
