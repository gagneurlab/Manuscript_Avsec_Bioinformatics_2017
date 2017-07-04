## make a scatterplot with logistic regression decision boundary
scatterplot_w_logit_boundary <- function(y_var, x1_var, x2_var, dt, with_interaction= TRUE,
                                         xlim = c(NA, NA),
                                         ylim = c(NA, NA), ...) {
  decision_func <- function(coef) {
    dummy <- function(x) {
      if(length(coef)==3) coef[4] <- 0
      return((- coef[1] - coef[2] * x)/(coef[3] + coef[4] * x))
    }
    return(dummy)
  }

  if(with_interaction == TRUE) {
    form = paste0(y_var, "~ ", x1_var, "*", x2_var)
  } else {
    form = paste0(y_var, "~ ", x1_var, "+", x2_var)
  }
  
  coef <- glm(as.formula(form), data = dt, family = binomial)$coef

  qplot(get(x1_var), get(x2_var), data =  dt,
        xlab = x1_var,
        ylab = x2_var,
        main = form,
        xlim = as.numeric(xlim),
        ylim = as.numeric(ylim)
        ) + aes_string(...)+
    scale_colour_manual(values = c("black", "red")) + stat_function(fun= decision_func(coef))
}

##' Produce a "Torkler" plot
##'
##' @param ... see arguments to torkler_plot_preproc

torkler_plot <- function(line_alpha=1,
                         use_theme=theme_bw(), ...) {
  summary_table <- torkler_plot_table(...)
  ## make a nice ggplot
  ggplot(summary_table, aes(position_bucket, gene_bucket)) +
    geom_raster(aes(fill = signal), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    xlab("Position wrt gene start [bp]") +
    ylab("Genes sorted by length") +
    use_theme +
    geom_line(aes(gene_end_bucket, gene_bucket), alpha=line_alpha)+
    geom_line(aes(gene_start_bucket, gene_bucket), alpha=line_alpha)
    ## scale_x_discrete(expand = c(0, 0)) +
    ## scale_y_discrete(expand = c(0, 0)) 
}


##' Produce torkler plot summary table
##' 
##' @param genes_gr (GRanges) Gene annotation
##' @param signal_gr (GRanges)Signal GRanges
##' @param signal_colname (character) Column name denoting the signal column. If NULL, a column of ones is used.
##' @param agg (function) Aggregation function.
##' @param fill (numeric) Signal value for bins without any data.
##' @param gene_bucket_size (integer) How many genes should we merge (useful for anti-aliasing)
##' @param bp_bucket_size 
##' @param n_upstream By how many bases to extend the ranges upstream of gene boundaries.
##' @param n_downstream By how many bases to extend the ranges downstream of gene boundaries.
##' @param position_bucket_size (integer) Number of base-pairs to be considered as one bucket.
##' @author Žiga Avsec <avsec@in.tum.de>
torkler_plot_table <- function(genes_gr, signal_gr, signal_colname = NULL,
                         agg = length, saturate_at = Inf, gene_length_quantile = .8, 
                         gene_bucket_size = 100, bp_bucket_size = 1000,
                         n_upstream = 5000, n_downstream = 5000) {

  ## consider only say 80% of genes
  width_limit <- quantile(width(genes_gr), probs = gene_length_quantile)
  genes_gr <- genes_gr[width(genes_gr) < width_limit]

  ## create genes table
  genes_gr2 <- copy(genes_gr)
  mcols(genes_gr2) <- NULL
  genes_gr2$gene_id = 1:length(genes_gr2)
  dt_genes <- genes_gr2 %>% as.data.frame %>% as.data.table
  dt_genes[, gene_length := width]
  dt_genes <- dt_genes[order(gene_length)][, c("gene_id", "gene_length")]
  ## create gene bucket
  dt_genes[, gene_bucket := cut_floor(1:.N,breaks = seq(0, .N + gene_bucket_size + 1, by = gene_bucket_size))]
  dt_genes[, gene_length_bucket := as.numeric(median(gene_length)), by = gene_bucket]
  dt_genes[, gene_start_bucket := 0]
  dt_genes[, gene_end_bucket := as.numeric(median(gene_length)), by = gene_bucket]

  ## create signal table
  dt <- map_ranges_to_genes(genes_gr, signal_gr,
                            n_upstream = n_upstream, n_downstream = n_downstream
                            )

  ## choose the signal column
  if (is.null(signal_colname)) {
    dt[, signal := 1]
  } else {
    stopifnot(length(signal_colname) == 1)
    dt[, (signal) := get(signal_colname)]
  }

  range_min <- -n_upstream
  range_max <- dt[, max(width) ]

  ## generate buckets
  ## 
  ## 1. position bucket wrt TSS
  dt[, position_bucket := cut_floor(position_wrt_TSS,
                                    breaks = seq(range_min, range_max + bp_bucket_size + 1, by = bp_bucket_size))]
  dt <- merge(dt, dt_genes, by = c("gene_id", "gene_length"))

  ## summarize the signal wrt buckets
  summary_table <- dt[, .(signal = agg(signal)), by = .(gene_bucket, gene_length_bucket,
                                                        gene_start_bucket, gene_end_bucket,
                                                        position_bucket
                                                        )]
  summary_table[, signal := pmin(signal, saturate_at)] #saturate signal

  ## make a nice ggplot
  return(summary_table)
}




##' TODO - Make it more general than only for genes
##' 
##' Map genome-wide signal to units like genes.
##' 
##' @param genes_gr (GRanges) Gene annotation
##' @param signal_gr (GRanges) Signal GRanges
##' @param n_upstream By how many bases to extend the ranges upstream of gene boundaries.
##' @param n_downstream By how many bases to extend the ranges downstream of gene boundaries.
##' @author Žiga Avsec <avsec@in.tum.de>
##'
##' @return Data.table with features like:
##' seqnames, start, end, width, strand: range of extended gene
##' gene_id, gene_start, gene_end, gene_length: gene related features (from genes_gr)
##' position: position centers of signal_gr
##' ...: mcols from signal_gr
##' position_wrt_*, gene_start_wrt_polya ..: same as position, gene_start, gene_end,
##'                                          but with respect to TSS or poly-a site (gene_start/gene_end)
map_ranges_to_genes <- function(genes_gr, signal_gr,
                                n_upstream = 5000, n_downstream = 5000) {

  ## restrict signal_gr only to the middle
  signal_gr <- resize(signal_gr, width = 1, fix = "center")
  
  genes_gr <- copy(genes_gr)
  mcols(genes_gr) <- DataFrame(gene_id= 1:length(genes_gr),
                               gene_start = start(genes_gr),
                               gene_end = end(genes_gr),
                               gene_length = width(genes_gr))
  ggr <- extend_range(genes_gr, n_upstream = n_upstream, n_downstream = n_downstream)

  ## merge genes and signal into one data.table
  ol <- findOverlaps(ggr, signal_gr)
  message("Fraction of data unmapped to genes: ", signif( 1 - uniqueN(subjectHits(ol)) / length(signal_gr), 3))
  ggr_ol <- ggr[queryHits(ol)]
  signal_ol <- signal_gr[subjectHits(ol)]
  ggr_ol$position <- start(signal_ol)
  mcols(ggr_ol) <- DataFrame(mcols(ggr_ol), mcols(signal_ol))

  dt <- ggr_ol %>% as.data.frame %>% as.data.table

  ## <colname>_x is the position with respect to the TSS
  dt[, position_wrt_TSS := ifelse(strand == "+", position - gene_start, gene_end - position)]
  dt[, gene_start_wrt_TSS := 0]
  dt[, gene_end_wrt_TSS := ifelse(strand == "+", gene_end - gene_start, gene_end - gene_end)]

  dt[, position_wrt_polya := ifelse(strand == "+", position - gene_end, gene_start - position)]
  dt[, gene_start_wrt_polya := 0]
  dt[, gene_end_wrt_polya := ifelse(strand == "+", gene_end - gene_end, gene_start - gene_end)]
  return(dt)
}
