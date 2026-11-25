# Side-by-side figure: locked-clock raw delta age (07) and compare_tAge reference-centred tAge shift (compare_tAge 03).
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
here<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
out<-dirname(here);project<-dirname(out)
a<-readRDS(file.path(out,'plots/02_delta_age_by_projection_class.rds'))
tb<-readRDS(file.path(project,'compare_tAge/figures/03_tage_shift_relative_to_reference_distributions.rds'))
# Panel b inherits its theme from the compare_tAge session; panel a carries a complete theme itself.
theme_set(tb$theme);b<-tb$plot
p<-(a|b)+plot_annotation(tag_levels='a')&theme(plot.tag=element_text(face='bold',size=14))
stem<-file.path(out,'diagnostics/03_raw_delta_age_and_tage_shift_combined')
ggsave(paste0(stem,'.pdf'),p,width=23,height=7,device=cairo_pdf,limitsize=FALSE)
ggsave(paste0(stem,'.png'),p,width=23,height=7,dpi=300,bg='white',limitsize=FALSE)
message('Rendered ',stem)
