## helper function for ggplot
library(ggplot2)
tilt_xlab <- theme(axis.text.x=element_text(angle=45, hjust=1))
tilt_ylab <- theme(axis.text.y=element_text(angle=90, hjust=1))

legend_top <- theme(legend.position="top")

legend_top_ver <- theme(legend.position="top", legend.direction = "vertical")

legend_off <- theme(legend.position="none")

facet_label_horizontal <- theme(strip.text.y = element_text(angle = 0))
## correlation matrix
## http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

## add n-observations to boxplot
boxplot_n_obs <- stat_summary(fun.data = function(x) c(y = median(x)*1.05, label = length(x)),
                              geom = "text", fun.y = median) 


save_plot_mul <- function(basename, formats, ...) {
  library(cowplot)
  for (format in formats) {
    save_plot(paste0(basename, ".", format), ...)
  }
}
