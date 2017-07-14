#'---
#' title: Create figure 2
#' author: Žiga Avsec
#' wb:
#'   input: ["data/eclip/processed/predictions/UPF1/kmer_glmnet-no_position.csv",
#'            "data/eclip/processed/predictions/AARS/DeepNN_2.csv",
#'             # ... all rbps, all methods (RBP_ALL, RBP_LIST)
#'           "data/eclip/processed/protein_peak_overlaps.rds",
#'           "data/eclip/processed/peak_center-gene_mapping.rds",
#'           "data/annotation/gencode.v25.annotation.gtf.rds",
#'           "data/eclip/processed/design_matrix/train/AARS_extended.csv"
#'           ]
#'---
##' ## Goal
##'
##'
##' ---------------------------------------------------------------
##' 
##+ set_workdir, echo=F
library(ggrepel)
library(ggsignif)
library(knitr)
library(rmarkdown)
library(cowplot)
library(gridExtra)
opts_knit$set(root.dir = getwd())
opts_chunk$set(echo=FALSE, cache=F, results = 'hide', messages = FALSE)
options(width=140)
## 
plt_dir <- "data/plots/RBP/Eclip/"
dir.create(plt_dir, showWarnings = FALSE, recursive = TRUE)
source("Scripts/Figures/config.R")

## ---------------------------------
## color scheme
cols <- c(col_glmnet, col_branchpointer, col_concise_shallow, col_concise_deep)

rbps <- c("UPF1", "PUM2", "DDX3X", "NKRF", "TARDBP", "SUGP2")
## ---------------------------------
dtb <- get_dtb_eclip_subset() # dtb_auc, dtb_auprc

q_auc <- ggplot(dtb[metric == "auc"], aes(x = Method, y = value)) +
  geom_boxplot(aes(color=Method)) +
  geom_jitter(aes(color=Method), alpha=0.1, size=0.5) +
  scale_color_manual(labels=c("glmnet", "glmnet w/ rel. distances",
                              "DNN", "DNN w/ rel. distances"),
                     values=cols) + 
  geom_signif(comparisons = list(c("DeepNN", "DeepNN_pos_spline-t"),
                                 c("kmer-glmnet_pos", "DeepNN"),
                                 c("kmer-glmnet_pos", "DeepNN_pos_spline-t")
                                 ),
              ## map_signif_level = TRUE,
              ## hjust = 0.9,
              map_signif_level = TRUE,
              tip_length=0.01,
              step_increase = c(0, 0, 0.03)
) +
  facet_grid(.~rbp) +
  ylab("Area under ROC curve") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) +
  ## guides(color=guide_legend(nrow=2)) +
  ## theme(legend.position = "bottom") + 
  theme(legend.justification = c(1, 0), legend.position = c(1, 0), legend.direction = "horizontal") 
q_auc

q_auprc <- (q_auc + ylab("Area under Precision-Recall curve")) %+% dtb[metric == "auprc"]

q_auprc + theme(panel.grid.major.y=element_line("grey50"))

plt_auc <- plot_grid(q_auc, labels="")
plt_auprc <- plot_grid(q_auprc, labels="a")

lgd <- get_legend(q_auc)
plt_boxplot <- plot_grid(plot_grid(q_auc + legend_off, q_auprc + legend_off,
                                   labels=c("a", "b"), align="v"),
                       lgd, ncol=1, rel_heights=c(1, 0.07))

## --------------------------------------------
## add the torkler plot
source_all("Scripts/RBP/Eclip/plot_functions")
library(scales)
tp <- top_genes_tork_plot(rbps,
                          bp_bucket_size = 1000,
                          gene_bucket_size = 100, saturate_at = 10,
                          cylab = "Gene index",
                          line_alpha = 0.2)

scientific_10 <- function(x, digits=3) {
  sci_format <- scales::scientific(x, digits=digits)
  ## remove +
  sci_format <- gsub("\\+", "", sci_format)
  ifelse(grepl("e00", sci_format),
         parse(text=gsub("e00", "", sci_format)),
         parse(text=gsub("e", " %*% 10^", sci_format))
         )
}

tp <- tp + scale_y_continuous(labels = scales::comma) +  scale_x_continuous(breaks = c(0, 50000), labels = scales::comma)
## tp <- tp + scale_y_continuous(labels = scientific_10) +  scale_x_continuous(breaks = c(0, 50000), labels = scientific_10)
tp

# scientific y-axis
plt_ab <- plot_grid(tp, q_auprc, ncol=1, labels = c("a", "b"), align="v", rel_heights = c(0.5, 1))

## --------------------------------------------
##' ## iDeep performance plot
dtm <- get_metric_comparison(1e-4)
df1 <- dtm[task != "branch-point" & metric == "auprc"]

## statistics for the paper
df1[task == "eClip" & subtask =="UPF1"]
##      task subtask metric       gam      relu    no_pos  p.gam__relu p.gam__no_pos signif.gam__relu signif.gam__no_pos
##    <char>  <char> <char>     <num>     <num>     <num>        <num>         <num>           <lgcl>             <lgcl>
## 1:  eClip    UPF1    auc 0.9444418 0.9392438 0.7992773 1.440169e-44  4.830856e-67             TRUE               TRUE

## > df1[task == "eClip", gam - no_pos] %>% summary
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## -0.002254  0.027951  0.051840  0.055303  0.071242  0.239969 

## df1[task == "eClip", wilcox.test(gam, no_pos, paired = TRUE)]
## V = 6306, p-value < 2.2e-16
## > dtm[task == "iClip" & metric == "auc", wilcox.test(gam, no_pos, paired = TRUE)]

## 	Wilcoxon signed rank test

## data:  gam and no_pos
## V = 488, p-value = 2.328e-08
## alternative hypothesis: true location shift is not equal to 0

## > dtm[task == "iClip" & metric == "auc", gam- no_pos] %>% summary
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## -0.004433  0.007605  0.017141  0.021911  0.029563  0.090431   
## --------------------------------------------
## shared RBP's: hnRNPC, Nsun2, PUM2, QKI, SRSF1, TAF15, TIA1, U2AF2

## dtm[task == "eCLIP" & metric == "auprc", .(subtask, signif.gam__no_pos)] %>% print(112)
## dtm[task == "eCLIP" & metric == "auprc", .(subtask, signif.gam__no_pos)][, table(signif.gam__o_pos)]
## CLIP not signif: ESWR1, PUM2, SRSF1, TAF15

## 104 / 112 - PUM2, SRSF1, TAF15 are all significant
## 10/11 for unique CLIP
## In total: 114 / 123
in_top <- function(x, n = 20) {
  return(rank(- x) <= n)
}
df1[, add_label := (subtask %in% rbps) |  !signif.gam__no_pos|
        (task == "eClip" & in_top(gam- no_pos, 10)) |
        (task == "iClip" & in_top(gam- no_pos, 10))
  , by = task] # TOP 20 from 
df1[, from_before := subtask %in% rbps & task != "iClip"]

df1[, task_name := ifelse(task == "iClip", "iClip (iDeep)", task)]

## TODO - split into two pieces
## TODO - show just one color
pl_scatter <- function(df1) {
  qpl1 <- ggplot(df1,
                 aes(x = no_pos,
                     y = gam,
                     color = !signif.gam__no_pos))+#factor(sign(gam - no_pos) * signif.gam__no_pos, levels = c(1, 0, -1)))) +
    geom_abline(alpha=0.2) +
    geom_point(size=2, alpha=0.7) +
    geom_text_repel(data = df1[add_label==TRUE], aes(label=subtask),
                    color="black",
                    box.padding = unit(0.4, "lines"),
                    point.padding = unit(0.4, "lines"),
                    alpha=0.5) +
    facet_wrap(~task_name, scales="free") + 
    scale_color_manual(values=c("black", "grey", "black",
                                mypal[3], "black", mypal[1], "black", "red"), # mypal[7] - alternatively for grey
                       ## limits = c(3, 4, 5),
                       guide = guide_legend("Wilcoxon test"),
                       drop=FALSE,
                       labels = c(bquote(list(p[adj] < 10^-4)),
                                  "Not significant",
                                  bquote(list(p[adj] < 0.0001, bar(x) > bar(y))),
                                  "", "")) + 
    xlim(c(0.29, .97)) + 
    ylim(c(0.29, .97)) + 
    ## xlim(c(0.68, 1)) +
    ## ylim(c(0.68, 1)) +
    ylab("auPR with spline transformation") +
    xlab("auPR default model (rel. distances not used)") + 
    legend_top_ver +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0)) 
  qpl1  
}
pl_eclip <- pl_scatter(df1[task =="eCLIP"]) + legend_off + xlab(NULL) + ylab(NULL) + white_facets
pl_iclip <- pl_scatter(df1[task =="CLIP"]) +
  xlim(c(.5, .97)) + 
  ylim(c(.5, .97)) + 
  xlab(NULL) +
  ylab(NULL)+ white_facets

xl <- ggdraw()  + cowplot::draw_label("auPR default model (rel. distances not used)")
yl <- ggdraw()  + cowplot::draw_label("auPR with rel. distances (spline tr.)", angle=90)

plt_test_perf <- plot_grid(plot_grid(yl, pl_eclip, pl_iclip,
                                     labels= c("c", "", "d"),
                                     rel_widths = c(0.05, 1, 1), nrow=1),
                           xl,
                           rel_heights = c(1, 0.05),
                           ncol=1)

## plot_grid(qpl2_auprc, qpl1_auprc)

fig2 <- plot_grid(plot_grid(tp, q_auprc, ncol=1, labels = c("a", "b"),
                            align="v", rel_heights = c(0.5, 1)),
                  plt_test_perf, rel_heights=c(1.3, 1), ncol = 1)
fig2

overleaf_plt_dir <- "~/projects-work/spline_trans/plots"

## eps will be made from inkscape, re-arranging the labels
save_plot_mul(file.path(overleaf_plt_dir, "fig2"), c("pdf", "png"), fig2,
          ncol = 1,
          nrow = 3,
          base_aspect_ratio = 2.8,
          base_height=3.5
          # each individual subplot should have an aspect ratio of 1.3
          )

## --------------------------------------------
## supplementary figures
##
## 2nd

dir_list <- list.files("data/eclip/processed/design_matrix/train", full.names = T) %>%
  grep("extended.csv$", ., value=TRUE)
dt_all <- lapply(dir_list, fread) %>% rbindlist

rbp_subset <- df1[task == "eCLIP" & add_label][order(-from_before,
                                                     -signif.gam__no_pos,
                                                     no_pos)][, subtask]
rbp_subset <- c(rbps, rbp_subset[!rbp_subset %in% rbps])
features <- c("gene_start", "TSS", "start_codon",
              "exon_intron", "intron_exon",
              "stop_codon", "polya", "gene_end")

dtm_all <- melt(dt_all, id.vars=c("rbp", "binding_site"), measure.vars=features)
dtm_all_pl <- dtm_all[rbp %in% rbp_subset]
dtm_all_pl[, rbp := fct_relevel(rbp, rbp_subset)]
dtm_all_pl[, variable := fct_relevel(variable, features)]

fig2_sup <- ggplot(dtm_all_pl, aes(x = - sign(value) * log10(abs(value)), color = binding_site)) +
  geom_density() +
  facet_grid(rbp~variable, scales="free_y") + 
  xlab("Signed log10 distance to landmark") + 
  theme(strip.text.y = element_text(angle=0)) +
  geom_vline(xintercept = 0, lty="dashed", alpha = 0.1)

#+ fig.height=15, fig.width=10
fig2_sup

save_plot_mul(file.path(overleaf_plt_dir, "fig2_sup2"), c("png", "eps", "pdf"), fig2_sup,
          ncol = 1,
          nrow = 1,
          base_aspect_ratio = 0.8,
          base_height=15
          )


## --------------------------------------------
## 1st
dt_genes <- readRDS("data/eclip/processed/peak_center-gene_mapping.rds")
dt_t <- dt_genes[, c(t.test(TSS_distance[binding_site == TRUE]/width[binding_site == TRUE],
                            TSS_distance[binding_site == FALSE]/width[binding_site == FALSE]
                            )[c("statistic", "p.value")], "N"= .N), by = rbp] %>%
  .[order(-abs(statistic),-N)]

dt_genes_wp <- merge(dt_genes, dt_t, by="rbp")
dt_genes_wp[, rbp := fct_reorder(rbp, -statistic)]
fig2_sup1 <- ggplot(dt_genes_wp ,aes(x = TSS_distance/ width, color = binding_site)) +
  geom_density() +
  facet_wrap(~rbp, ncol = 8) +
  xlab("Relative gene position") + 
  ## ggtitle("Poly-a distance") +
  legend_top + 
  tilt_xlab

save_plot_mul(file.path(overleaf_plt_dir, "fig2_sup1"), c("png", "eps", "pdf"), fig2_sup1,
          ncol = 1,
          nrow = 1,
          base_aspect_ratio = 0.8,
          base_height=15
          )
