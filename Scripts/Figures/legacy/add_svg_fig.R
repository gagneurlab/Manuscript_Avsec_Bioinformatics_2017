## TODO - add overall performance?

## plt_global <- plot_grid(plot_grid(plt1,
##                                   plot_grid(perf_boxplot, loss_plt,
##                                             ncol=1,
##                                             labels = c("b", "c"),
##                                             rel_heights=c(1.2,1)),
##                                   labels = "a",
##                                   rel_widths = c(1, 0.7)
##                                   ),
##                         lgd,
##                         ncol=1,
##                         rel_heights = c(1,0.01))
## plt_global
a <- 1


## install.packages("grImport2")
## devtools::install_github("sjp/grConvert")

## gg_read_svg <- function(path) {
##   library(grImport2)
##   library(grConvert)
##   ## create tmp dir
##   tmp_file <- tempfile(fileext=".svg")
##   convertPicture(path, tmp_file)
##   fig_raw <- readPicture(tmp_file)
##   fig_grob <- gTree(children=gList(pictureGrob(fig_raw)))
##   fig <- grid.arrange(fig_grob)
##   return(fig)
## }

## fig2a <- gg_read_svg("~/droak/work/concise/paper/fig2a.svg")
## fig2b <- gg_read_svg("~/droak/work/concise/paper/fig2b-plain.svg")
## TODO - why is there an error?
## plot_grid(fig2a, fig2b)


## plot_grid(grid.arrange(fig1a_grob))
## TODO - assemble it in a meaningful way 
## TODO - inject pdf, svg plots using ggplot & cowplot?
