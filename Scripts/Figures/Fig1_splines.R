#'---
#' title: Create figure 1c
#' author: Ziga Avsec
#' wb:
#'  input: []
#'---
library(knitr)
opts_chunk$set(echo=FALSE, cache=F, results = 'hide', messages = FALSE)
options(width=140)
#+ 
library(motifp)
library(cowplot)
plt_dir <- "./data/plots/Concise/Paper/"

## get the output dir

dir.create(plt_dir, showWarnings = FALSE)
#'
#'
#' ## Plot
gam_obj <- get_gam_splines(seq(0, 100, 0.1), n_bases = 10, spline_order = 3)

x_ext <- seq(-50, 150, 0.1)
gam_obj_ext <- predict_gam_splines(gam_obj$mgcv_smoothcon_obj,
                                   x=x_ext, col_names=paste0("spline", 1:10))
dtm_ext<- data.table(x = x_ext, gam_obj_ext) %>% melt(id.vars ="x")
dtm_ext[, i := str_replace(variable, "spline", "")]
dtm_ext[, bi := paste0("b[", i, "](x)")]


dtm <- data.table(x = gam_obj$x, gam_obj$X) %>% melt(id.vars ="x")
dtm[, i := str_replace(variable, "spline", "")]
dtm[, bi := paste0("b[", i, "](x)")]



## get the peak annotations

w <- c(1, 1, 1, 1.5, 1.1, 0.9, 0.9, 1, 1, 1)
dtm[, value_w := value * w[as.integer(i)]]
dtm_xy <- data.table(x = gam_obj$x, y = (gam_obj$X %*% w)[,1])
dtm_max <- dtm[, .SD[which.max(value)], by = variable]

dtm_knots <- data.table(i = 1:length(gam_obj$knots), knot_x= gam_obj$knots)


plt <- ggplot(dtm, aes(x = x, y = value_w, group=bi, color=bi)) +
  geom_line()+ 
  geom_hline(yintercept = 1, lty="dashed", alpha = 0.4) + 
  ## geom_text(aes(x = x + 2.99, y = value + 0.03, label = bi), dtm_max, parse=TRUE) +
  geom_line(aes(x, y, group=NULL, color=NULL), dtm_xy) +
  geom_point(aes(x = knot_x, y = 0, group=NULL, color=NULL), dtm_knots) +
  xlim(-5, 105) + 
  geom_text(aes(x = knot_x, y = -0.04, label=i, group=NULL, color=NULL), dtm_knots) + 
  annotate("text", x = 0, y = -0.04, size=4, label="Knot id   ", hjust="outward") +
  annotate("text", x = 6, y = 0.78, size=4.5, label="Spline basis functions") +
  annotate("text", x = 2, y = 1.15, size=4.5, label="f[S](x) == sum(w[i] * b[i](x), i=1, B)",
           parse=TRUE) +
  annotate("text", x = 100, y = 1.03, size=3.5,
           hjust="inward",
           label="w = [1, 1, 1, 1.5, 1.1, 0.9, 0.9, 1, 1, 1]") +
  annotate("text", x = 5, y = 1.27, size=4.5, label="Spline transformation") +
  theme(legend.position="none") +
  ylab(NULL) +
  ggtitle("P-splines")
plt

save_plot_mul(file.path(plt_dir, "fig1c"), c("pdf", "png"), plt,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 1,
          base_height=6,
          base_aspect_ratio = 1.3
          )

#+
a=1
