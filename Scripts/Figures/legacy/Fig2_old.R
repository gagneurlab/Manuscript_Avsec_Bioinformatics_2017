#'---
#' title: Create figure 2
#' author: Žiga Avsec
#' wb:
#'   input: ["data/encode/eclip/processed/predictions/UPF1/kmer_glmnet-no_position.csv"]
#'---
##' ## Goal
##'
##' 
##' 
##' ---------------------------------------------------------------
##' 
##+ set_workdir, echo=F
library(ggrepel)
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
source("Scripts/RBP/Eclip/plot_functions/pr_roc.R")
source("Scripts/RBP/Eclip/config.R")
source("Scripts/Concise/Paper/func/fig2.R")

## ---------------------------------
## color scheme
library(cowplot)
library(ggsci)
mypal = pal_locuszoom()(7)
col_concise_shallow <- mypal[4]
col_concise_deep <- mypal[5]
col_relu_shallow <- mypal[3]
col_relu_deep <- mypal[1]
col_branchpointer <- mypal[2]
col_glmnet <- mypal[7]

cols <- c(col_glmnet, col_relu_shallow, mypal[2], col_relu_deep)
rbps <- c("UPF1", "PUM2", "DDX3X", "NKRF", "TARDBP", "SUGP2")
## ---------------------------------
attach(get_dtb()) # dtb_auc, dtb_auprc

## remove relu 
dtb_auc <- dtb_auc[!grepl("relu", as.character(Method))]
dtb_auprc <- dtb_auprc[!grepl("relu", as.character(Method))]

q_auc <- ggplot(dtb_auc, aes(x = Method, y = auc)) +
  geom_boxplot(aes(color=Method)) +
  geom_jitter(aes(color=Method), alpha=0.1, size=0.5) +
  scale_color_manual(values=cols) + 
  geom_signif(comparisons = list(c("DeepNN", "DeepNN_pos_spline-t"),
                                 c("kmer-glmnet_pos", "DeepNN_pos_spline-t")
                                 ),
              ## map_signif_level = TRUE,
              hjust = 0.9,
              step_increase = 0.07) +
  facet_grid(.~rbp) +
  ylab("Area under ROC curve") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) +
  ## guides(color=guide_legend(nrow=2)) +
  theme(legend.position = "bottom")
q_auc

q_auprc <- (q_auc + aes(y = auprc) + ylab("Area under Precision-Recall curve")) %+% dtb_auprc

plt_auc <- plot_grid(q_auc, labels="")
plt_auprc <- plot_grid(q_auprc, labels="a")

lgd <- get_legend(q_auc)
plt_boxplot <- plot_grid(plot_grid(q_auc + legend_off, q_auprc + legend_off,
                                   labels=c("a", "b"), align="v"),
                       lgd, ncol=1, rel_heights=c(1, 0.07))
#+fig.width=12, fig.height=6
plt_boxplot

save_plot(file.path(plt_dir, "roc_pr_boxplot.pdf"), plt_boxplot,
          ncol = 2,
          nrow = 1,
          base_aspect_ratio = 1.4,
          base_height=5,
          # each individual subplot should have an aspect ratio of 1.3
          )

save_plot(file.path(plt_dir, "roc_pr_boxplot.png"), plt_boxplot,
          ncol = 2,
          nrow = 1,
          base_aspect_ratio = 1.4,
          base_height=5,
          )


save_plot(file.path(plt_dir, "pr_boxplot.pdf"), plt_auprc,
          ncol = 1,
          nrow = 1,
          base_aspect_ratio = 2,
          base_height=5,
          # each individual subplot should have an aspect ratio of 1.3
          )

save_plot(file.path(plt_dir, "pr_boxplot.png"), plt_auprc,
          ncol = 1,
          nrow = 1,
          base_aspect_ratio = 2,
          base_height=5,
          )

save_plot(file.path(plt_dir, "roc_boxplot.pdf"), plt_auc,
          ncol = 1,
          nrow = 1,
          base_aspect_ratio = 2,
          base_height=5,
          # each individual subplot should have an aspect ratio of 1.3
          )

save_plot(file.path(plt_dir, "roc_boxplot.png"), plt_auc,
          ncol = 1,
          nrow = 1,
          base_aspect_ratio = 2,
          base_height=5,
          )

## --------------------------------------------
## add the torkler plot
source_all("Scripts/RBP/Eclip/plot_functions")
tp <- top_genes_tork_plot(rbps,
                          bp_bucket_size = 1000,
                          gene_bucket_size = 100, saturate_at = 10,
                          cylab = "Gene index",
                          line_alpha = 0.2)

plt_ab <- plot_grid(tp, q_auprc, ncol=1, labels = c("a", "b"), align="v", rel_heights = c(0.5, 1))

## --------------------------------------------
## position mapping
dt_genes <- readRDS("data/encode/eclip/processed/peak_center-gene_mapping.rds")

dtp <- dt_genes[rbp %in% rbps]
dtp[, rbp := fct_relevel(rbp, rbps)]
melt(dtp, measure.vars = c("TSS_distance", "polya_distance"))

ggplot(dtp, aes(x = log10(abs(TSS_distance)), color = binding_site)) +
  geom_density() +
  facet_grid(~rbp)
  
#' 
##' ## iDeep performance plot

dti <- fread("data/encode/eclip/processed/predictions/iClip-iDeep_boostrap_Amin.csv")
dti[, V1 := NULL]
dti %>% setnames("method", "Method")
dti[, Method := str_replace(Method, "_scaler", "")]
dti[, Method := str_replace(Method, "_position", "_pos")]
dti[,i := 1:.N , by = .(Method, rbp)]


dti_auprc<- fread("data/encode/eclip/processed/predictions/iClip-iDeep_boostrap_Amin_auprc.csv")
dti_auprc[, V1 := NULL]
dti_auprc %>% setnames("method", "Method")
dti_auprc[, Method := str_replace(Method, "_scaler", "")]
dti_auprc[, Method := str_replace(Method, "_position", "_pos")]
dti_auprc[,i := 1:.N , by = .(Method, rbp)]

library(ggsignif)

## q_auc <- ggplot(dti, aes(x = Method, y = auc)) +
##   geom_boxplot(aes(color=Method)) +
##   geom_jitter(aes(color=Method), alpha=0.1, size=0.5) +
##   scale_color_manual(values=cols) + 
##   geom_signif(comparisons = list(c("iDeep", "iDeep_pos_gam"),
##                                  c("iDeep_pos_relu", "iDeep_pos_gam")
##                                  ),
##               ## map_signif_level = TRUE,
##               hjust = 0.9,
##               step_increase = 0.07) +
##   facet_wrap(~rbp, nrow=2) +
##   ylab("Area under ROC curve") + 
##   theme(axis.title.x=element_blank(),
##         axis.text.x=element_blank(),
##         axis.ticks.x = element_blank()) +
##   ## guides(color=guide_legend(nrow=2)) +
##   theme(legend.position = "bottom")
## q_auc

## do the scatterplot

dti_x <- data.table::dcast(dti, rbp  + i ~ Method, value.var="auc")

dti_iDeep_vs_iDeep_pos_gam <- dti_x[, .(pval = wilcox.test(iDeep, iDeep_pos_gam)$p.value,
                                        iDeep = mean(iDeep),
                                        iDeep_std = std(iDeep),
                                        iDeep_min = min(iDeep),
                                        iDeep_max = max(iDeep),
                                        iDeep_pos_gam = mean(iDeep_pos_gam),
                                        iDeep_pos_gam_std = std(iDeep_pos_gam),
                                        iDeep_pos_gam_min = min(iDeep_pos_gam),
                                        iDeep_pos_gam_max = max(iDeep_pos_gam)
                                        ) , by = rbp]

dti_iDeep_pos_relu_vs_iDeep_pos_gam <- dti_x[, .(pval = wilcox.test(iDeep_pos_relu, iDeep_pos_gam)$p.value,
                                        iDeep_pos_relu = mean(iDeep_pos_relu),
                                        iDeep_pos_relu_std = std(iDeep_pos_relu),
                                        iDeep_pos_relu_min = min(iDeep_pos_relu),
                                        iDeep_pos_relu_max = max(iDeep_pos_relu),
                                        iDeep_pos_gam = mean(iDeep_pos_gam),
                                        iDeep_pos_gam_std = std(iDeep_pos_gam),
                                        iDeep_pos_gam_min = min(iDeep_pos_gam),
                                        iDeep_pos_gam_max = max(iDeep_pos_gam)
                                        ) , by = rbp]
## auprc
dti_auprc_x <- data.table::dcast(dti_auprc, rbp  + i ~ Method, value.var="auprc")

dti_auprc_iDeep_vs_iDeep_pos_gam <- dti_auprc_x[, .(pval = wilcox.test(iDeep, iDeep_pos_gam)$p.value,
                                                    iDeep = mean(iDeep),
                                                    iDeep_std = std(iDeep),
                                                    iDeep_min = min(iDeep),
                                                    iDeep_max = max(iDeep),
                                                    iDeep_pos_gam = mean(iDeep_pos_gam),
                                                    iDeep_pos_gam_std = std(iDeep_pos_gam),
                                                    iDeep_pos_gam_min = min(iDeep_pos_gam),
                                                    iDeep_pos_gam_max = max(iDeep_pos_gam)
                                                    ) , by = rbp]

dti_auprc_iDeep_pos_relu_vs_iDeep_pos_gam <- dti_auprc_x[, .(pval = wilcox.test(iDeep_pos_relu, iDeep_pos_gam)$p.value,
                                                             iDeep_pos_relu = mean(iDeep_pos_relu),
                                                             iDeep_pos_relu_std = std(iDeep_pos_relu),
                                                             iDeep_pos_relu_min = min(iDeep_pos_relu),
                                                             iDeep_pos_relu_max = max(iDeep_pos_relu),
                                                             iDeep_pos_gam = mean(iDeep_pos_gam),
                                                             iDeep_pos_gam_std = std(iDeep_pos_gam),
                                                             iDeep_pos_gam_min = min(iDeep_pos_gam),
                                                             iDeep_pos_gam_max = max(iDeep_pos_gam)
                                                             ) , by = rbp]

qpl1 <- ggplot(dti_iDeep_vs_iDeep_pos_gam, aes(x = iDeep_pos_gam,
                                       xmin = iDeep_pos_gam - iDeep_pos_gam_std,
                                       xmax = iDeep_pos_gam + iDeep_pos_gam_std,
                                       y = iDeep,
                                       ymin = iDeep - iDeep_std,
                                       ymax = iDeep + iDeep_std,
                                       color = as.factor(sign(iDeep_pos_gam - iDeep) * (pval < 0.05 / 31)))) +
  geom_text_repel(aes(label=rbp), color="black", alpha=0.4) +
  geom_point(size=2) +
  ## geom_errorbar(size=0.5) +
  ## geom_errorbarh(size=0.5) +
  scale_color_manual(values=c(mypal[3], mypal[7], mypal[1]),
                     guide = guide_legend("Wilcoxon test"),
                     labels = c(bquote(list(p < 0.0016, bar(x) < bar(y))),
                                "NS",
                                bquote(list(p < 0.0016, bar(x) > bar(y))))) + 
  geom_abline() +
  xlim(c(0.68, 1)) +
  ylim(c(0.68, 1)) +
  xlab("AUC iDeep_pos_gam") +
  ylab("AUC iDeep") + 
  legend_top_ver +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0))


qpl2 <- ggplot(dti_iDeep_pos_relu_vs_iDeep_pos_gam, aes(x = iDeep_pos_gam,
                                       xmin = iDeep_pos_gam - iDeep_pos_gam_std,
                                       xmax = iDeep_pos_gam + iDeep_pos_gam_std,
                                       y = iDeep_pos_relu,
                                       ymin = iDeep_pos_relu - iDeep_pos_relu_std,
                                       ymax = iDeep_pos_relu + iDeep_pos_relu_std,
                                       color = as.factor(sign(iDeep_pos_gam - iDeep_pos_relu) * (pval < 0.05 / 31)))) +
  geom_text_repel(aes(label=rbp), color="black", alpha=0.4) +
  geom_point(size=2) +
  ## geom_errorbar(size=0.5) +
  ## geom_errorbarh(size=0.5) +
  scale_color_manual(values=c(mypal[3], mypal[7], mypal[1]),
                     guide = guide_legend("Wilcoxon test"),
                     labels = c(bquote(list(p < 0.0016, bar(x) < bar(y))),
                                "NS",
                                bquote(list(p < 0.0016, bar(x) > bar(y))))) + 
  geom_abline() +
  xlim(c(0.68, 1)) +
  ylim(c(0.68, 1)) +
  xlab("AUC iDeep_pos_gam") +
  ylab("AUC iDeep_pos_relu") + 
 legend_top_ver +
    theme(legend.justification = c(1, 0), legend.position = c(1, 0))

## plot_grid(qpl1 + legend_off, qpl2)

## TODO - save the plot

qpl1_auprc <- qpl1 %+% dti_auprc_iDeep_vs_iDeep_pos_gam + legend_off +
  xlab("auPRC iDeep_pos_gam") +
  ylab("auPRC iDeep")  +
  xlim(c(0.5, 1)) + 
  ylim(c(0.5, 1))

qpl2_auprc <- qpl2 %+% dti_auprc_iDeep_pos_relu_vs_iDeep_pos_gam  +
  xlab("auPRC iDeep_pos_gam") +
  ylab("auPRC iDeep_pos_relu")  +
  xlim(c(0.5, 1))  + 
  ylim(c(0.5, 1))

plt_ab <- plot_grid(tp, q_auprc, ncol=1, labels = c("a", "b"), align="v", rel_heights = c(0.5, 1))

fig2 <- plot_grid(plt_ab, 
                  plot_grid(qpl2_auprc, qpl1_auprc), rel_heights=c(1.5, 1), ncol = 1)

## TODO add the scatter-plot from eClip
fig2

save_plot(file.path(plt_dir, "fig2.pdf"), fig2,
          ncol = 1,
          nrow = 3,
          base_aspect_ratio = 2.8,
          base_height=4,
          # each individual subplot should have an aspect ratio of 1.3
          )

save_plot(file.path(plt_dir, "fig2.png"), plt_ab,
          ncol = 1,
          nrow = 3,
          base_aspect_ratio = 2.2,
          base_height=4,
          )
