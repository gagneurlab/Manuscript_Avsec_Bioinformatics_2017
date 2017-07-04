#'---
#' title: Analyze mercer HC features
#' author: Žiga Avsec
#' wb:
#'  input: ["data/Concise/Splice_branchpoints/processed/branchpointer/train/branchpoint_df_HCN.csv"]
#'---
#'
#' ## References
#'
#' - Dean et. al.  http://web.mit.edu/vdean/www/MLCB2016_paper.pdf
#' - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4315302/
#' - ppt tract https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4315302/figure/F4/
#'
#' ## Conclusions
#'
#' - dataset fits our expectations
#' - Positions have a very strong effect on the response
#'
#' ---------------
library(cowplot)
library(seqLogo)
library(gglogo)

#+ cache=TRUE
dt <- fread("data/Concise/Splice_branchpoints/processed/branchpointer/train/branchpoint_df_HCN.csv")
dt[, V1:= NULL]
dt[, is_hc := as.integer(set == "HC")]
#+

#' 
#' ## Distance dependencies

dt_dist1 <- dt[, .(is_hc= mean(is_hc), N = .N), by = dist.1]
pl_d1 <- qplot(dist.1, is_hc ,data = dt_dist1[dist.1 < 10000], color = N, alpha = I(0.3)) +
  ylab("Fraction of positive branchpoints") + 
  xlab("Donor distance") +
  ggtitle("Donor distance ")

dt_dist2 <- dt[, .(is_hc= mean(is_hc), N = .N), by = dist.2]
pl_d2 <- qplot(-dist.2, is_hc ,data = dt_dist2, color = N, size = N) +
  ylab("Fraction of positive branchpoints") + 
  xlab("Acceptor distance") +
  ggtitle("Acceptor distance ")

dt_ppts <- dt[, .(is_hc= mean(is_hc), N = .N), by = ppt_start]
pl_ppts <- qplot(ppt_start, is_hc ,data = dt_ppts, color = N, size = N) +
  ylab("Fraction of positive branchpoints") + 
  xlab("PPT start") +
  ggtitle("PPT start ") 
  ## ylim(c(0, 0.1))

dt_pptl <- dt[, .(is_hc= mean(is_hc), N = .N), by = ppt_run_length]
pl_pptl <- qplot(ppt_run_length, is_hc ,data = dt_pptl, color = N, size = N) +
  ylab("Fraction of positive branchpoints") + 
  xlab("PPT length") +
  ggtitle("PPT length ") +
  ylim(c(0, 0.1))

#+ fig.width =10, fig.height=8
plot_grid(pl_d1, pl_d2,
          pl_ppts, pl_pptl, ncol=2
          )
plot_cannon <- function(i =1) {
  id <- paste0("canon_hit",i)
  dt_canon_hit<- dt[, .(is_hc= mean(is_hc), N = .N), by = id]
  setnames(dt_canon_hit, paste0("canon_hit",i), "canon_hit")
  qplot(canon_hit, is_hc, data = dt_canon_hit, color = N, size = N) +
    ylab("Fraction of positive branchpoints") + 
    xlab(paste("Cannonical AG distance", i)) +
    ggtitle(paste("Cannonical AG distance", i))  
}

#'
#' ## Cannonical AG distances
#'
#' Left unfiltered, right filtered

#+ fig.width = 10, fig.height = 15
plot_grid(
  plot_cannon(1),
  plot_cannon(1) + ylim(c(0, .15)) + xlim(c(0, 75)),
  plot_cannon(2),
  plot_cannon(2) + ylim(c(0, .15)) + xlim(c(0, 150)),
  plot_cannon(3),
  plot_cannon(3) + ylim(c(0, .10)) + xlim(c(0, 150)),
  plot_cannon(4),
  plot_cannon(4) + ylim(c(0, .10)) + xlim(c(0, 150)),
  plot_cannon(5),
  plot_cannon(5) + ylim(c(0, .10)) + xlim(c(0, 150)),
  ncol = 2
)

#'
#' ## Correlation between cannonical AG distances and acceptor distance
#+ fig.width =10, fig.height = 6
par(mfrow=c(2,3))
for (i in 1:5) {
  ch <- paste0("canon_hit", i)
  with(dt[set == "HC"], heatscatter(x = dist.2,
                                    y = get(ch),
                                    ylab = ch,
                                    main = ch))
}
par(mfrow=c(1,1))

#'
#' ## Sequence logos

seq_cols <- grep("^seq_", colnames(dt), val=T)
dt[, seq := do.call(paste0, c(.SD[, seq_cols, with = F]))]

cm <- dt[set == "HC", seq] %>% consensusMatrix(as.prob=TRUE)
cm_neg <- dt[set == "NEG", seq] %>% consensusMatrix(as.prob=TRUE)

#'
#' ### HC sites
#+ fig.width = 4, fig.height = 4
seqLogo(cm, ic.scale = FALSE)
seqLogo(cm, ic.scale = TRUE)
#'
#' ### NEG sites
#+ fig.width = 4, fig.height = 4
seqLogo(cm_neg, ic.scale = FALSE)
seqLogo(cm_neg, ic.scale = TRUE)
