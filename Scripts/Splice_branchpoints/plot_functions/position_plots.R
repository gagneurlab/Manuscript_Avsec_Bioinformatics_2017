#'---
#' title: Position plots
#'---

#  input: ["data/Concise/Splice_branchpoints/processed/branchpointer/train/branchpoint_df_HCN.csv",
#          "data/Concise/Splice_branchpoints/interpret/position/concise_shallow.csv"]

DIST_FEATURES <- c("dist1" = "Donor distance",
                   "dist2" = "Acceptor distance",
                   "ppt_start" = "PPT start",
                   "ppt_run_length" = "PPT length", 
                   "canon_hit1" = "Cannonical AG distance 1",
                   "canon_hit2" = "Cannonical AG distance 2",
                   "canon_hit3" = "Cannonical AG distance 3",
                   "canon_hit4" = "Cannonical AG distance 4",
                   "canon_hit5" = "Cannonical AG distance 5")


tidy_position_measured <- function() {
  dt <- fread("data/Concise/Splice_branchpoints/processed/branchpointer/train/branchpoint_df_HCN.csv")
  dt[, V1:= NULL]
  dt[, is_hc := as.integer(set == "HC")]
  setnames(dt, "dist.1", "dist1")
  setnames(dt, "dist.2", "dist2")
  response <- "is_hc"

  dtm <- melt(dt, id.vars = response, measure.vars = names(DIST_FEATURES))
  dtm <- dtm[, .(p = mean(is_hc), N = .N), by = .(variable, value)]
  ## dtm[, variable_name := DIST_FEATURES[variable]]
  return(dtm)
}

tidy_position <- function() {
  dtm <- tidy_position_measured()
  dti <- tidy_position_inferred()

  setnames(dtm, "value", "position")

  setnames(dti, "x", "position")
  setnames(dti, "y", "p")
  setnames(dti, "feature", "variable")
  dtm[, method := "measured"]

  dt <- rbindlist(list(dtm, dti), fill=TRUE)
  dt[, variable_name := DIST_FEATURES[variable]]
  dt[, primary := variable %in% c("dist2", "ppt_start")] #, "canon_hit1", "canon_hit2")]
  return(dt)
}

tidy_position_inferred <- function() {
  methods <- c("concise_shallow",
               "shallow_relu_tid=4250")
  hash <- c("concise_shallow" = "shallow_gam",
            "shallow_relu_tid=4250" = "shallow_relu")


  dtpos <- lapply(methods, function(method) {
    path <- paste0("data/Concise/Splice_branchpoints/interpret/position/", method, ".csv")
    dt <- fread(path)
    if (method %in% names(hash)) {
      dt[, method := hash[method]]    
    } else {
      dt[, method := method]
    }
    dt[, V1 := NULL]
    dt[, y := sigmoid(y)]
    return(dt)
  }) %>% rbindlist
  return(dtpos)
}


get_position_plots <- function(show_title=TRUE, aes_str=aes_string(alpha="N")) {
  dt <- fread("data/Concise/Splice_branchpoints/processed/branchpointer/train/branchpoint_df_HCN.csv")
  dt[, V1:= NULL]
  dt[, is_hc := as.integer(set == "HC")]


  my_ggtitle <- function(label) {
    if (isTRUE(show_title)) {
      return(ggtitle(label))
    } else {
      return(NULL)
    }
  }

  fn_marginal <- function(x) {
    ## input = p
    return(x)
    ## return(- log(x) / log(1-x))
  }
  ## 2663
  ## methods <- c("shallow_gam_tid=2663",
  ##              "shallow_relu_tid=4250")
  ## hash <- c("shallow_gam_tid=2663" = "shallow_gam",
  ##           "shallow_relu_tid=4250" = "shallow_relu")
  methods <- c("concise_shallow",
               "shallow_relu_tid=4250")
  hash <- c("concise_shallow" = "shallow_gam",
            "shallow_relu_tid=4250" = "shallow_relu")


  dtpos <- lapply(methods, function(method) {
    path <- paste0("data/Concise/Splice_branchpoints/interpret/position/", method, ".csv")
    dt <- fread(path)
    if (method %in% names(hash)) {
      dt[, method := hash[method]]    
    } else {
      dt[, method := method]
    }
    dt[, V1 := NULL]
    dt[, y := sigmoid(y)]
    return(dt)
  }) %>% rbindlist

  #' 
  #' ## Distance dependencies
  #'
  #' ### Left = computed fraction, right = inferred by concise model
  #'
  #' 
  #' - on the left, points color and size match to the number of points with that
  #' particular x position (large blue circle means high numbers)
  #' - y axis = faction of branchpoints at position x
  #'
  #' Example code for the left part:
  #+ eval = FALSE, echo =TRUE
  dt_dist2 <- dt[, .(is_hc= fn_marginal(mean(is_hc)), N = .N), by = dist.2]
  pl_d2 <- qplot(-dist.2, is_hc ,data = dt_dist2) + aes_str
  pl_d2
  #+
  no_legend <- theme(legend.position="none")
  yl <- ylab("BP fraction") 
  ym <- ylab("Inferred effect")

  ## acceptor
  dt_dist2 <- dt[, .(is_hc= fn_marginal(mean(is_hc)), N = .N), by = dist.2]
  pl_d2 <- qplot(-dist.2, is_hc ,data = dt_dist2) + aes_str + 
    yl + 
    xlab("Acceptor distance") +
    my_ggtitle("Acceptor distance ") +
    no_legend

  m_pl_d2 <- qplot(-x, y, data = dtpos[feature == "dist2"], color = method, geom = 'line') +
    ym + 
    xlab("Acceptor distance") +
    my_ggtitle("Acceptor distance ") +
    no_legend

  ## donor
  dt_dist1 <- dt[, .(is_hc= fn_marginal(mean(is_hc)), N = .N), by = dist.1]
  pl_d1 <- qplot(dist.1, is_hc ,data = dt_dist1[dist.1 < 10000]) +
    aes_str + 
    yl + 
    xlab("Donor distance") +
    my_ggtitle("Donor distance ") +
    no_legend

  m_pl_d1 <- qplot(x, y, data = dtpos[feature == "dist1"], color = method, geom = 'line') +
    ym + 
    xlab("Donor distance") +
    my_ggtitle("Donor distance ")  +
    no_legend

  ## ----

  dt_ppts <- dt[, .(is_hc= fn_marginal(mean(is_hc)), N = .N), by = ppt_start]
  pl_ppts <- qplot(ppt_start, is_hc ,data = dt_ppts) +
    aes_str + 
    yl + 
    xlab("PPT start") +
    my_ggtitle("PPT start ") +
    no_legend
  ## ylim(c(0, 0.1))

  m_pl_ppts <- qplot(x, y, data = dtpos[feature == "ppt_start"], color = method, geom = 'line') +
    ym + 
    xlab("PPT start") +
    my_ggtitle("PPT start ")  +
    no_legend


  dt_pptl <- dt[, .(is_hc= fn_marginal(mean(is_hc)), N = .N), by = ppt_run_length]
  pl_pptl <- qplot(ppt_run_length, is_hc ,data = dt_pptl) + aes_str + 
    yl + 
    xlab("PPT length") +
    my_ggtitle("PPT length ") +
    ylim(c(0, 0.1)) +
    no_legend

  m_pl_pptl <- qplot(x, y, data = dtpos[feature == "ppt_run_length"], color = method, geom = 'line') +
    ym + 
    xlab("PPT length") +
    my_ggtitle("PPT length ") +
    no_legend

  plot_cannon <- function(i =1) {
    if (i == 1) {
      xl <- xlim(c(0, 75))
    } else {
      xl <- xlim(c(0, 150))
    }
    id <- paste0("canon_hit",i)
    dt_canon_hit<- dt[, .(is_hc= fn_marginal(mean(is_hc)), N = .N), by = id]
    setnames(dt_canon_hit, paste0("canon_hit",i), "canon_hit")
    qplot(canon_hit, is_hc, data = dt_canon_hit) +
      aes_str + 
      yl + 
      xlab(paste("Cannonical AG distance", i)) +
      my_ggtitle(paste("Cannonical AG distance", i)) +
      xl + 
      no_legend
  }

  m_plot_cannon <- function(i = 1) {
    if (i == 1) {
      xl <- xlim(c(0, 75))
    } else {
      xl <- xlim(c(0, 150))
    }
    id <- paste0("canon_hit",i)
    qplot(x, y, data = dtpos[feature == id], color = method, geom = 'line') +
      ym + 
      xlab(paste("Cannonical AG distance", i)) +
      my_ggtitle(paste("Cannonical AG distance", i)) +
      xl  +
      no_legend
  }

  plot_list <- list(pl_d2=pl_d2,
       m_pl_d2=m_pl_d2,
       pl_d1=pl_d1,
       m_pl_d1=m_pl_d1,
       pl_ppts=pl_ppts,
       m_pl_ppts=m_pl_ppts,
       pl_pptl=pl_pptl,
       m_pl_pptl=m_pl_pptl,
       pl_c1=plot_cannon(1) + ylim(c(0, .15)),
       m_pl_c1=m_plot_cannon(1),
       pl_c2=plot_cannon(2) + ylim(c(0, .15)),
       m_pl_c2=m_plot_cannon(2),
       pl_c3=plot_cannon(3) + ylim(c(0, .10)),
       m_pl_c3=m_plot_cannon(3),
       pl_c4=plot_cannon(4) + ylim(c(0, .10)),
       m_pl_c4=m_plot_cannon(4),
       pl_c5=plot_cannon(5) + ylim(c(0, .10)),
       m_pl_c5=m_plot_cannon(5)
       )

  return(plot_list)
      
}

