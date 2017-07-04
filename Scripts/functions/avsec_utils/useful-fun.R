## Other many useful functions

##' Get basename without extension
basename_core <- function(path) {
  tools::file_path_sans_ext(basename(path))
}

##' Test whether x is either:
##' - NULL
##' - NA
##' - NaN
##' - vector of length 0
##' - ""
## http://stackoverflow.com/questions/19655579/a-function-that-returns-true-on-na-null-nan-in-r
is.blank <- function(x, false.triggers=FALSE){
  if(is.function(x)) return(FALSE) # Some of the tests below trigger
  # warnings when used on functions
  return(
    is.null(x) ||                # Actually this line is unnecessary since
      length(x) == 0 ||            # length(NULL) = 0, but I like to be clear
      all(is.na(x)) ||
      all(x=="") ||
      (false.triggers && all(!x))
  )
}

##' is the vector unique
##' @param vec
##' @return TRUE or FALSE
is_unique <- function(vec) {
  !any(duplicated(vec))
}

#' Fill NA values assuming the values increase or decrease by some value
#'
#' @examples
#' x <- c(NA, 7, NA, 9)
#' y <- c(NA, 9, NA, 7)
#' z <- c(NA, 7, NA, 11)
#' all.equal(fill_incr_na(x), 6:9)
#' all.equal(fill_incr_na(y, -1), 10:7)
#' all.equal(fill_incr_na(z, 2), seq(5,11,2))
fill_incr_na <- function(x, incr_by=1) {
  if (incr_by >= 0) {
    ref <- as.integer(seq(1, length(x) * incr_by, by = incr_by))
  } else {
    ref <- as.integer(seq(length(x)*abs(incr_by), 1, by = incr_by))
  }
  first_non_na <- which(!is.na(x))[1] # position of first non-na
  diff <- x[first_non_na] - ref[first_non_na] 
  recovered <- ref + diff
  ## sanity_check
  stopifnot(mean(recovered == x, na.rm=TRUE) == 1)
  return(recovered)
}

##' Replacement to which.max - it returns all the maximal values
##' 
##' @return a vector of indicies where the value is maximal
which_max <- function(vec) {
  which(vec == max(vec, na.rm = TRUE))
}

##' Replacement to which.min - it returns all the minimal values
##' 
##' @return a vector of indicies where the value is maximal
which_min <- function(vec) {
  which(vec == min(vec, na.rm = TRUE))
}



## list the duplicated values
duplicated_values <- function(vec) {
  unique(vec[duplicated(vec)])
}

## list the duplicated values - boolean vector
duplicated_values_bool <- function(vec) {
  vec %in% vec[duplicated(vec)]
}

## Rectified linear function
relu <- function(x) {
  pmax(x,0)
}

## useful function - transform list names into factor levels and append them to the dataframe
df_list2dt <- function(df_list,colname="which",to.numeric=FALSE,ordered.factor=FALSE){
  ## we have to copy the table in order to avoid trouble
  df_list=lapply(df_list,function(x) as.data.table(copy(x)))

  ##' we don't have any list names
  if(is.null(names(df_list))) {
    names(df_list) <- 1:length(df_list)
  }
  
  lapply(df_list,function(df) if(ncol(df)==1 && "V1" %in% names(df)) setnames(df,"V1","vec"))
  
  lapply(df_list,function(df)
    if(any(names(df)=="dummy")) warning("'dummy' variable is present in the list!")
    )
  df_list=sapply(names(df_list),function(dummy) df_list[[dummy]][,dummy:=dummy],simplify=FALSE)
  df=rbindlist(df_list)

  if(ordered.factor) df[,dummy:=factor(dummy,levels = names(df_list))]
  else     df[,dummy:=factor(dummy)]

  if(to.numeric) {
    dummy <- df[,dummy]
    if(!any(is.na(dummy))){
      ## check if we get an error after conversion
      tt <- tryCatch(res <- as.numeric(as.character(dummy)),warning=function(w) w)
      ## if there is no error - convert it to numeric
      if(!is(tt,"warning")) {
        df[,dummy:=res]
        res <- as.numeric(as.character(dummy))
      }}
  }
  setnames(df,"dummy",colname)
  df
}    

## do df_list2dt recursively 
df_list2dt_rec <- function(df_list,colnames=c("which"),to.numeric=FALSE,ordered.factor=FALSE){
  ldepth=listdepth(df_list)
  if(length(colnames)!=ldepth) stop(paste0(
    "Colnames have to have length",ldepth,"yours is",length(colnames)))
  if(ldepth==0) return("listdepth of df_list is 0")
  if(ldepth==1) return(df_list2dt(df_list,colnames[1],to.numeric,ordered.factor))
  if(ldepth>1) return( df_list2dt(lapply(df_list,df_list2dt_rec,colnames[-1],
                                         to.numeric,ordered.factor) ,colnames[1],to.numeric,ordered.factor) )
}

## returns the list depth
listdepth <- function(somelist){
  listdepth_rec<- function(somelist){
    if(all(class(somelist)!="list")) return(1)
    else return(1+listdepth_rec(somelist[[1]]))
  }
  if(all(class(somelist)!="list")) return(0)
  else return(listdepth_rec(somelist)-1)
}

##' Cut the vector into classes similar to base::cut,
##' however cut_floor returns the interval's lower boundary (numeric)
##' 0.5 -> [0, 1) ; 0 -> [0, 1)
##' @example
##' .bincode(c(-1, 0,.5, 1, 100), c(-1,0,1), include.lowest = TRUE, right = FALSE)
##' ## [1]  1  2  2 NA
cut_floor <- function(x, breaks) {
  breaks[.bincode(x, breaks = breaks, include.lowest = TRUE, right = FALSE)]
}

##' Transform the data.table with a variable of type "list of atomic vectors" into a long data.table, with one row per element (extend the list of atomic vectors)
##' 
##' @param dt data.table
##' @param variables variable names from dt to be expanded (variables should be of class list of atomic vectors)
##' @param id_name how to name the unique_col_id
##' @param position_name how to name the position column
##' @param Class Character vector of the corresponding classes. Can be also of length one.
##' @param Separator Character vector of separators. Can be also of length one.
##' @return data.table with one additional column "pos" (position) and changed `get(variable)` column
dt_expand_vector_list <- function(dt, variables, id_name = "row_id", position_name = "pos", Class = "numeric", separator = ";") {
  dt <- copy(dt)

  stopifnot(all(variables %in% names(dt)))

  Class <- rep(Class,length.out = length(variables))
  separator <- rep(separator, length.out = length(variables))
  
  ## repeat
  for (i in 1:length(variables)) {
    variable <- variables[i]
    if (is.character(dt[[variable]]) | is.factor(dt[[variable]])) {
      message("variable '", variable, "' is a character or factor. Using function char_vec2vec_list to convert it to a list of vectors.")
      dt[, (variable):= char_vec2vec_list(get(variable), Class = Class[i], separator = separator[i])]
    }
  }

  ## if more variables, they all have to have the same length
  if (length(variables) > 1) {
    X <- sapply(variables, function(variable) dt[[variable]] %>% sapply(length))
    if (!all_identical(as.data.frame(X))) {
      stop("Not all the length are the same")
    }
  }
  
  dt[, unique_col_id := 1:.N]
  ## multi-column support
  dt_expanded <- dt[, c(list(pos = 1:length(get(variable[1])[[1]])), lapply(variables, function(x) get(x)[[1]])), by = unique_col_id]
  setnames(dt_expanded, c("unique_col_id", position_name, variables))
  dtm <- merge(dt[, -variables, with = F], dt_expanded, by = "unique_col_id")
  setnames(dtm, "unique_col_id", id_name)

  ## for (i in 1:length(variables)) {
  ##   message("i = ", i)
  ##   if (Class[i] == "character") {
  ##     if (!all(nchar(as.character(dtm[[variables[i]]])) == 1)) {
  ##       warning("Not all character variables in variable: ", variable, " have final length 1")
  ##     }
  ##   }
  ## }
  return(dtm)
}

##' Make a vector of valid names:
make_names <- function(name) {
  ## replace space or -
  regex <- "-| |\\."
  gsub(regex,"_", name)
}

##' test if the elements of a list are identical
##'
##' @param elem_list List of elements to be compared with identical
##' @param ... additional arguments passed to identical
##' @return Boolean
all_identical <- function(elem_list, ...) {
  ## http://stackoverflow.com/a/18813590
  if (!is.list(elem_list)) stop("elem_list has to be a list-like object")
  if (length(elem_list) == 1) stop("elem_list has to have length > 1")
  
  ## compare all elements to the first element
  for (i in 2:length(elem_list)) {
    identical_elem <- identical(elem_list[[1]], elem_list[[i]], ...)
    
    if (identical_elem == FALSE) return(FALSE)
  }

  return(TRUE)
}

## transform a data.table of p-values to a data.table of fdr's; keep the first column as it is
p_to_fdr <- function(res_path_dis=res_path_ibd){
  result=foreach(i = 1:length(res_path_dis), pvalues=res_path_dis) %do% {
    if(i==1) pvalues else p.adjust(pvalues, 'BH')
  }
  names(result)=names(res_path_dis)
  return(as.data.table(result))
}

## convert the confusion matrix to a data.table
confusionMatrix2dt <- function(cm){
  data.table(TP=cm[1,1],FN=cm[1,2],FNR=cm[1,3],
             FP=cm[2,1],TN=cm[2,2],FPR=cm[2,3])
}

## suggested by Leo; mclapply wrapper - nice printing, 
mcadply <- function(X, FUN2,...) {
  ## Runs multicore lapply with progress indicator and transformation to
  ## data.table output
  ##
  ## Arguments (same as lapply):
  ## X:   Vector
  ## FUN: Function to apply to each value of X
  ##
  ## Output: data.table stack of each mclapply return value
  ##         Note FUN is transformed to a data.frame return if necessary
  ##
  ## Progress bar code based on http://stackoverflow.com/a/10993589

  require(multicore)
  require(plyr)

  local({
    f <- fifo(tempfile(), open="w+b", blocking=T)
    if (inherits(fork(), "masterProcess")) {
      ## Child
      progress <- 0
      print.progress <- 0
      while (progress < 1 && !isIncomplete(f)) {
        msg <- readBin(f, "double")
        progress <- progress + as.numeric(msg)
        ## Print every 1%
        if(progress >= print.progress + 0.01) {
          cat(sprintf("Progress: %.0f%%\n", progress * 100))
          print.progress <- floor(progress * 100) / 100
        }
      }
      exit()
    }

    newFun <- function(...) {
      ret = FUN2(...)
      writeBin(1 / length(X), f)
      return(ret)
    }
    ##result = lapply(X = X, FUN = newFun,...)
    result <- mclapply(X = X, FUN = newFun,...,mc.preschedule=FALSE)
    close(f)
    cat("Done\n")
    return(result)
  })
}
############################################
## other useful functions
############################################

## [force attach] works similar like attach, but deletes all the variables in the list x first and then attaches them
attach.all <- function (x, overwrite = NA, name = "attach.all")  {
  rem <- names(x) %in% ls(.GlobalEnv)
  if (!any(rem)) overwrite <- FALSE
  rem <- names(x)[rem]
  if (is.na(overwrite)) {
    question <- paste("The following objects in .GlobalEnv will mask\nobjects in the attached database:\n", paste(rem, collapse = ", "), "\nRemove these objects from .GlobalEnv?", sep = "")
    if (interactive()) {
      if (.Platform$OS.type == "windows")  overwrite <- "YES" == winDialog(type = "yesno",  question)
      else overwrite <- 1 == menu(c("YES", "NO"), graphics = FALSE, title = question)
    }
    else overwrite <- FALSE
  }
  if (overwrite) remove(list = rem, envir = .GlobalEnv)
  attach(x, name = name)
}

##' Hash a vector values using the hash_list
##' 
##' @param vec Vector of values we would like to hash
##' @param hash_list list representing the hash function key:value
##' @return vector of the same length as `vec`, but with hashed values
##' 
hash_vector <- function(vec, hash_list) {
  if (!all(as.character(unique(vec)) %in% names(hash_list))) {
    stop("Some values in the vector were not present in the hash_list")
  }
  return(unlist(hash_list[as.character(vec)]))
}


## takes a factor, changes the level order upside down and returns the releveled factor
reverse_level <- function(f){
  if(is.character(f)) {
    f=as.factor(f)
  }
  nlevel=length(levels(f))
  f=factor(f,levels(f)[nlevel:1])
  return(f)
}

##' Put factor levels to the end
##'
##' @param f Factor
##' @param levels Levels to be put back
##' @return re-leveled factor f
##' @author Žiga Avsec
factor_put_back <- function(f, levels) {
  stopifnot(all(levels %in% levels(f)))

  f <- factor(f, levels = c(levels(f)[!levels(f) %in% levels], levels))
  return(f)
}

refactor <- function(f) {
  if(class(levels(f))=="character") return(as.factor(as.character(f)))
  if(class(f)=="character") return(f)
  stop("You should use this function only for factors with character levels or character vectors")

}

to_numeric <- function(f){
  if(class(f)=="numeric") return(f)
  if(class(f) %in% c("integer","character","logical")) return(as.numeric(f))
  if(class(f)=="factor") return( as.numeric(f) )
  stop("class not recognized")
}

to_integer <- function(f){
  if(class(f)=="integer") return(f)
  if(class(f) %in% c("numeric","character","logical")) return(as.integer(f))
  if(class(f)=="factor") return( as.integer(f) )
  stop("class not recognized")
}

as.factor_keep_order <- function(f){
  return(factor(f,unique(f)))
}


##' Threshold integers
threshold_int<- function(x, n = 2) {
  ifelse(x >= n, paste0(">=", n), as.character(x))
}
threshold_int_fctr <- function(x, n = 2) {
  factor(threshold_int(x, n = n), levels = threshold_int(min(x):n, n = n))
}

##' Cut by quantiles but with nice names
cut_quantile <- function(x, n_quantiles) {
  p <- quantile(x, probs = seq(0, 1, length = n_quantiles + 1))
  cut(x, breaks = p, labels = names(p)[-1], include.lowest = TRUE)
}

##' Convert a vector of characters into an integer vector with
##' subsequent occurrence number
##' @param char_vec Character vector or factor
##' 
##' @examples
##' c("a", "a", "b", "c", "a") %>% char_id
##' ## [1] 1 1 2 3 1
char_id <- function(char_vec) {
  v <- as.factor_keep_order(char_vec)
  return(as.integer(v))
}

## list of vectors to a character vector (elements separated by ;)
## [[1]] c(1,2) -> c("1;2")
vec_list2char_vec <- function(vec_list, separator = ";") {
  char_vec <- sapply(vec_list, paste, collapse = separator)
  return(char_vec)
}
## inverse operation of vec_list2char_vec
## c("1;2") -> [[1]] c(1,2)
char_vec2vec_list <- function(char_vec, Class = "numeric", separator = ";") {
  if (is.factor(char_vec)) {
    warning("char_vec is a factor. Converting to a character vector.")
    char_vec <- as.character(char_vec)
  }
  
  vec_list <- strsplit(char_vec, split = separator)
  vec_list <- lapply(vec_list, as, Class = Class)
  return(vec_list)
}

## create folds for cross validation    
cv_folds = function(n, folds = 10, seed = NULL, test_and_train = FALSE){
  if (!is.null(seed)) set.seed(seed)
  res <- split(sample(1:n), rep(1:folds, length = n))
  if (isTRUE(test_and_train)) {
    res <- lapply(res, function(x) list(test = (1:n)[x], train = (1:n)[-x]))
  }
  return(res)
}

## takes a list of data.tables, finds common column names and does rbindlist on them
rbindlist_match <- function(dtlist){
  common_names <- Reduce(intersect,lapply(dtlist,names))
  rbindlist(lapply(dtlist,function(dt) dt[,common_names,with=FALSE]))
}

## pretty print in knitr - uncomment if you want to switch it on
## from
## http://cran.rstudio.com/web/packages/knitr/vignettes/knit_print.html
## knit_print.data.frame <- function(x, ...) {
##     require(knitr)
##     res = paste(c("", "", kable(x)), collapse = "\n")
##     asis_output(res)
## }

## lload - lazy load
#' don't load the object if all of the following conditions hold:
#' - the file was already loaded once with lload
#' - the file hasn't changed since the last call of lload on this file
#' - all of the variables from the file are still present in the global enviroment
#'
#' verbose - print the message if the file was actually loaded or not
#' Useful for fun2source in iterargs
#'
#' Take care: if the loaded file was changed, and you want to reload it, then use `load`
lload <- function(file,verbose=FALSE,force=FALSE){
  
  if(!file.exists(file)) stop(paste("No such file:", file))

  ## create the new enviroment
  if(!exists(".lload.env")) local(.lload.env <- new.env(),envir=.GlobalEnv)

  ## create the loaded file if it doesn't exist
  if(!exists("loaded",envir=.lload.env)) assign("loaded",list(),envir=.lload.env)

  ## pass the file to the enviroment
  assign("file",file,envir=.lload.env)
  assign("verbose",verbose,envir=.lload.env)

  assign("force",force,envir=.lload.env)
  
  local({
    load_the_file=TRUE
    
    if(file %in% names(loaded)){
      if(loaded[[file]]$mtime == file.info(file)$mtime &
                                                 all(loaded[[file]]$vars %in% ls(envir=.GlobalEnv))) {
        if(verbose) message("File already loaded.")
        load_the_file=FALSE
      }}
    
    if(force) load_the_file=TRUE
    if(load_the_file) {
      if(verbose) message("Actually loading the file with load().")
      ## if not, load the file
      vars <- load(file=file,envir=.GlobalEnv)

      loaded[[file]]$mtime <- file.info(file)$mtime
      loaded[[file]]$vars <- vars
    }
  }
 ,envir=.lload.env)
  return(
    invisible(local(loaded[[file]]$vars,envir=.lload.env))
  )
}

## Remove all the load history of lload function. Files calles with `lload` will be loaded from the file again using `load`.
clear.lload <- function(){
  rm(list=ls(envir=.lload.env),envir=.lload.env)
}


## get the repository file
## get_dirname=TRUE - get the full file path
## get_dirname=TRUE - get only directory name of the file

## in case the file name can not be found, the function will display a warning and return NULL
thisFile <- function(get_dirname=FALSE) {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    ## Rscript
    thisfile <- normalizePath(sub(needle, "", cmdArgs[match]))
  } else {
    ## 'source'd via R console
    if(!is.null(sys.frames()[[1]]$ofile)) {
      thisfile <- normalizePath(sys.frames()[[1]]$ofile)
    } else {
      warning("Can't get the file name, returning NULL")
      return(NULL)
    }
  }

  if(get_dirname) return(dirname(thisfile))
  return(thisfile)
}
############################################
## knitr  Rmd -> pdf without using render()
############################################

## rmd.convert("2-state-of-the-art.Rmd","pdf")
rmd.convert <- function(fname, output=c('latex', 'word', 'html', "pdf"),rm_md_log=TRUE){
  ## Thanks to Robert Musk for helpful additions to make this run better on Windows

  require(knitr)
  require(tools)
  
  thedir <- file_path_as_absolute(dirname(fname))
  thefile <- (basename(fname)) 
  
  create_latex <- function(f){
    knit(f, 'tmp-outputfile.md'); 
    newname <- paste0(file_path_sans_ext(f), ".tex")
    mess <- paste('pandoc -f markdown -t latex -s -o', shQuote(newname), 
                  "tmp-outputfile.md")
    system(mess)
    cat("The Latex file is", file.path(thedir, newname), 
        "\nIf transporting do not forget to include the folder", file.path(thedir, "figure"), "\n")
    mess <- paste('rm tmp-outputfile.md')
    system(mess)
  }

  create_word <- function(f){
    knit(f, 'tmp-outputfile.md');
    newname <- paste0(file_path_sans_ext(f),".docx")
    mess <- paste('pandoc -f markdown -t docx -o', shQuote(newname), "tmp-outputfile.md")
    system(mess)
    cat("The Word (docx) file is", file.path(thedir, newname), "\n")
    mess <- paste('rm tmp-outputfile.md')
    system(mess)
  }
  
  create_html <- function(f){
    knit2html(f)
    cat("The main HTML file is", file.path(thedir, paste0(file_path_sans_ext(f), ".html")), 
        "\nIf transporting do not forget to include the folder", file.path(thedir, "figure"), "\n")
  }

  create_pdf <- function(f){
    knit(f, 'tmp-outputfile.md');
    newname <- paste0(file_path_sans_ext(f),".pdf")
    mess <- paste('pandoc -f markdown -o', shQuote(newname), "tmp-outputfile.md")
    system(mess)
    cat("The PDF file is", file.path(thedir, newname), "\n")
    mess <- paste('rm tmp-outputfile.md')
    system(mess)
  }

  origdir <- getwd()  
  tryCatch({
    setwd(thedir) ## put us next to the original Rmarkdown file
    out <- match.arg(output)
    switch(out,
           latex=create_latex(thefile),
           html=create_html(thefile),
           pdf=create_pdf(thefile),
           word=create_word(thefile)
           )}, finally=setwd(origdir))

}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


## this function will remove all the cached figures in order to avoid duplicated plots (if you rename a chunk)
rm_cached_figures <- function(fig.path=knitr::opts_chunk$get("fig.path")){
  if(is.null(fig.path)) fig.path <- "./spin-ex_files/figure-html"
  if(file.exists(fig.path)) invisible(sapply(dir(fig.path, full.names = TRUE), file.remove))
}


multiline_kable <- function(dt,nc_per_line=10,...){
  df <- as.data.frame(dt)
  nc <- ncol(df)

  nrow <- if(nc %% nc_per_line ==0) nc/nc_per_line else nc %/% nc_per_line +1
  set.seed(1)
  for(i in 1:(nrow-1)) print(knitr::kable(df[,(1+(i-1)*nc_per_line):((i)*nc_per_line)],...))
  print(knitr::kable(df[,(1+(nrow-1)*nc_per_line):nc],...))

}


##' Show only figures from R output and no code or printed output
knitr_only_plots <- knitr::opts_chunk$set(echo=FALSE,results="hide", cache=TRUE, warning = FALSE)


save_ggplot <- function(plot,path,width=NA,height=NA){
  if(is.na(width)) width <- 7
  if(is.na(height)) height <- 7
  
  pdf(path,width=width,height=height)
  print(plot)
  dev.off()
}


## use with care - I don't know if this works correctly
## py_pred - probability p(y=1)
## y - true class
roc.curve <- function(py_pred,y){
  thresh_vec <- seq(1,0,length.out = 100)
  lapply(thresh_vec,function(s){
    Ps <- ifelse(py_pred>s,1,0)
    FPR=sum(Ps==1 & y==0)/sum(y==0)
    TPR=sum(Ps==1 & y==1)/sum(y==1)
    data.table(TPR,FPR,s)
  }) %>% rbindlist
}

## make (print) a nice scatterplot of a data.frame using the ggpairs::GGally function
my_scatterplot <- function(dt,axisLabels="none"){
  library(GGally)
  scatterplot <- ggpairs(dt,axisLabels=axisLabels)+
    theme(legend.position = "none", 
          panel.grid.major = element_blank(),
          ## axis.ticks = element_blank(), 
          panel.border = element_rect(
            ## linetype = "dashed",
            colour = "black",
            fill = NA)
          )
  
  ## print(scatterplot,left=0.4,bottom=0.2,spacing = 0.05)
  print(scatterplot)
}

## set column order in data.table - specify the first columns:
setcolorder_first <- function(x,startorder){
  colnames <- names(x)
  if(!all(startorder %in% colnames)) stop("some columns don't exist in the initial data.table")

  neworder <- c(startorder,colnames[!(colnames %in% startorder)])
  setcolorder(x,neworder)
}

setcolorder_last <- function(x, endorder){
  colnames <- names(x)
  if(!all(endorder %in% colnames)) stop("some columns don't exist in the initial data.table")

  neworder <- c(colnames[!(colnames %in% endorder)], endorder)
  setcolorder(x,neworder)
}

pretty_print_list <- function(mylist) {
  max_nchar <- max(nchar(names(mylist)))

  mylist_char <- sapply(names(mylist), function(elem_name) {
    elem <- mylist[[elem_name]]
    if (all(is.numeric(elem))) {
      elem_str <- sapply(elem, function(x) sprintf("%.3g", x))
    } else {
      elem_str <- elem
    }
    ph <- paste0("%-", max_nchar + 1, "s")
    elem_str <- paste(elem_str, collapse = ", ")
    paste0(sprintf(ph, elem_name), ": ", elem_str)
  })

  message(paste(mylist_char, collapse = "\n"))
  
}


## to be used for timing your code, start_timer(); stop_timer()
## run this for all the samples
start_timer <- function(){
  message(paste("Timer start: ",Sys.time()))
  ## create a new environment if it doesn't exist
  if(!exists(".timer.env")) {
    local(.timer.env <- new.env(),envir=.GlobalEnv)
    .timer.env$t_list <- list()
  }
  ## append current time to t_list
  t0 <- Sys.time()
  .timer.env$t0 <- t0
  .timer.env$t_list[[length(.timer.env$t_list) + 1]] <- t0
  return(t0)
}
stop_timer <- function(){
  ## time environment has to exist
  if(!exists(".timer.env")) {
    warning("You first have to run start_timer()!")
    return(invisible(NULL))
  }

  ## t_list has have some elements
  if(length(.timer.env$t_list) < 1 ) {
    warning("Too few start_timer() instances. Using the last instance of start_timer.")
    t0 <- .timer.env$t0
  } else {
    ## get the latest time entry
    t0 <- .timer.env$t_list[[length(.timer.env$t_list)]]
    ## remove the latest time entry
    .timer.env$t_list[[length(.timer.env$t_list)]] <- NULL
  }
  ## compute the time difference with now
  t1 <- Sys.time()
  td <- difftime(t1,t0,units="secs")
  ## get the number of days 
  tdays <- floor(difftime(t1,t0,units="days"))
  tdformat <- format(.POSIXct(td,tz="GMT"), "%H:%M:%S")
  ## append day counts
  tdformat <- paste0(as.character(tdays),"-",tdformat)
  message(paste("\nEnd:   ",Sys.time()))
  message(paste("\nExecution time:",tdformat))
  return(tdformat)
}

reset_timer <- function(){
  rm(".timer.env", envir = .GlobalEnv)
}

## find a common substring always starting at the beginning
## example: words <- c("abc123","abc432")
## return should be: abc
find_common_substring_beginning <- function(words) {
  ## took from:
  ## 
  ## http://stackoverflow.com/questions/26285010/r-find-largest-common-substring-starting-at-the-beginning
  ## Roland solution
  ##
  
  if(length(words)==1) return(words)
  ##extract substrings from length 1 to length of shortest word
  subs <- sapply(seq_len(min(nchar(words))), 
                 function(x, words) substring(words, 1, x), 
                 words=words)
  #max length for which substrings are equal
  neqal <- max(cumsum(apply(subs, 2, function(x) length(unique(x)) == 1L)))
  #return substring
  substring(words[1], 1, neqal)
}

## Chris' function print_log
print_log <- function(..., print_stack = TRUE){
  msg <- list(...)
  
  # check if print_stack is given
  if(!is.logical(print_stack)){
    msg <- append(msg, print_stack, 0)
    print_stack = TRUE
  }
  
  # combine all input
  if(length(msg) > 0){
    log_msg <- paste(sapply(msg, function(x) eval(x)), collapse = "")
  } else {
    log_msg <- "No message was given!"
  }
  
  # create stack trace
  stack_trace <- ""
  if(print_stack){
    calls <- sys.calls()
    call_names <- sapply(calls, function(x) x[[1]])
    call_names <- call_names[1:(length(call_names) - 1)]
    stack_trace <- paste(paste(call_names, collapse = "() -> "), "()", sep="")
  }
  
  # create message
  log_msg <- gsub("[\n]", "\n\t", log_msg)
  log_msg <- paste(date(), ": ", stack_trace, "\n\tlog: ", log_msg, "\n", sep="")
  
  # print message
  message(log_msg)
}

## print data.table in a nice format, so that you can use in the warning messages
print_and_capture <- function(x) {
  paste(capture.output(print(x)),collapse = "\n")
}

## does a vector have only one value?
has_single_value <- function(vec) {
  length(unique(vec))==1
}

## does a vector have only multiple values?
has_multiple_values <- function(vec) {
  !has_single_value(vec)
}


## get the columns that have multiple values
mult_value_col <- function(dt){
  names(dt)[sapply(dt,has_multiple_values)]
}

single_value_col <- function(dt) {
  mvc <- mult_value_col(dt)
  return(setdiff(names(dt), mvc))
}

## get the duplicated "by_variables"
## variables that are not unique:
subset_dt_by_duplicated_variables <- function(dt,by_variables){
  dt_unique <- dt[,by_variables,with=FALSE]
  dt_unique <- dt_unique[duplicated(dt_unique)]
  setkey(dt,"sample")
  dt[dt_unique]
}

## get the columns that cause the values of by_variables to not be unique
## i.e. for each "by_variables", look at which columns have multiple values
show_mult_value_col <- function(dt,by_variables){
  dt_duplicated <- subset_dt_by_duplicated_variables(dt,by_variables)

  ## which variables are causing this?
  show_variables <- dt_duplicated[,.(vars=mult_value_col(.SD)),by=c(by_variables)][,unique(vars)]


  ## display only this table
  dt_duplicated <- unique(dt_duplicated[,by_variables,with=FALSE])
  setkey(dt)
  dt[dt_duplicated][,c(by_variables,show_variables),with=FALSE]
}

## concatenate together cells with multiple elements in data.table
## use case:
## dt[,variable:=concat(variable)]
concat <- function(variable){
  sapply(variable,function(x) paste(x,collapse= ' '))
}

##
solve_multiple_value_col <- function(dt,by_variables){
  ## helper function for solve_multiple_value_col
  solve_only_one_non_NA <- function(vec){
    if(has_single_value(vec)) return(unique(vec))

    ## if only 1 is not NA, return this value
    values <- unique(vec)
    if(sum(!is.na(values))==1) return(values[!is.na(values)])

    ## in case you're not able to solve it, return the vector
    return(vec)
  }

  
  problematic_dt <- subset_dt_by_duplicated_variables(dt,by_variables)
  problematic_samples <- problematic_dt[,by_variables,with=FALSE] %>% unique

  problematic_dt_solution <- problematic_dt[,lapply(.SD,solve_only_one_non_NA),by=c(by_variables)]

  ## integrate the solution into data.table
  setkeyv(dt,by_variables)
  dt_remove_problems <- dt[!problematic_dt]
  dt_solved <- rbindlist(list(dt_remove_problems,problematic_dt_solution), use.names=TRUE)

  
  if(nrow(subset_dt_by_duplicated_variables(dt_solved,by_variables))==0) {
    message("\nProblem with sample_data solved.\n")
    return(dt_solved)
  }
  
  return(NULL)
}

## list functionf to debug
## http://stackoverflow.com/questions/6950602/listing-functions-with-debug-flag-set-in-r
ls.deb  <- function(items = search ()){
  .ls.deb <-  function (i){
    f <- ls (i)
    f <- mget (f, as.environment (i), mode = "function",

               ## return a function that is not debugged
               ifnotfound = list (function (x) function () NULL)
               )

    if (length (f) == 0)
      return (NULL)

    f <- f [sapply (f, isdebugged)]
    f <- names (f)

    ## now check whether the debugged function is masked by a not debugged one
    masked <- !sapply (f, function (f) isdebugged (get (f)))

    ## generate pretty output format:
    ## "package::function"  and "(package::function)" for masked debugged functions
    if (length (f) > 0) {
      if (grepl ('^package:', i)) {
        i <- gsub ('^package:', '', i)
        f <- paste (i, f, sep = "::")
      }

      f [masked] <- paste ("(", f [masked], ")", sep = "")

      f
    } else {
      NULL
    }
  }


  functions <- lapply (items, .ls.deb)
  unlist (functions)
}

##' merge a list of data.tables by merge_by vector of variables
##'
##' @param dt_list A list of data.tables of data.frames.
##' @param merge_by Vector of variables to merge by.
##' @return Merged data.table.
mergelist <- function(dt_list,merge_by, all=FALSE){
  Reduce(function(...) merge(...,by=merge_by, all=all),dt_list)
}

## split the 
split_dt <- function(dt, N){
  split(dt, sample(1:N, nrow(dt), replace=T))
}

##' autosize all columns in an excel `file` (for all sheets)
xls_autosizecol <- function(file) {
  library(xlsx)
  wb <- loadWorkbook(file)
  sheets <- getSheets(wb)
  # set widths to 20
  for(csheet in sheets) autoSizeColumn(csheet, colIndex = 1:csheet$getLastRowNum() )
  saveWorkbook(wb,file)
}


##' Get the most common level from a factor or character vector
most_freq <- function(f) {
  names(which.max(table(f)))
}

##' Handle NA values
##' @param NA_handling
##' - "none" - do nothing
##' - "-1"   - replace numeric NA's with -1 and for factor NA's add another level
##' - "median" - replace numeric NA's with median factor NA's with the most common class
##' @param exclude_col which column should be ignored
handle_NA <- function(dt, NA_handling = "none", exclude_col = "y"){
  dt <- copy(dt)
  message("NA_handling = ", NA_handling)
  if(NA_handling == "none") return(dt)
  
  if(NA_handling  %in% c("-1", "median")) {
    for(col in names(dt)){
      ## skip if the column is in exclude_co
      if(col %in% exclude_col) next()
      
      ## numeric + integer
      if(all(class(dt[[col]]) %in%  c("integer","numeric"))) {
        ## -1
        if(NA_handling == "-1"){
          dt[, (col) := ifelse(is.na(get(col)),
                               -1,      #impute -1
                               get(col))]
        }
        ## median
        if(NA_handling == "median"){
          dt[, (col) := ifelse(is.na(get(col)),
                               median(get(col), na.rm = TRUE), #impute median
                               get(col))]
        }
      }
      ## factors
      if(all(class(dt[[col]]) %in%  c("factor"))) {
        if(NA_handling == "-1"){
          dt[, (col) := addNA(get(col), ifany = TRUE)] #add NA as another level
        }
        ## median
        if(NA_handling == "median"){
          x <- dt[,get(col)]
          x[is.na(x)] <- most_freq(x)
          dt[, (col) := as.factor(x)]
        }
        
      }
    }
    
    return(dt)
  }

}


##' Read multiple data.tables and stack them into one bigger data.table
##' @param path_vec Character vector of file paths to a table.
##' @param colname Name of the column containing the file information
##' @param fold_names How should we name each element (rows in colname)
##' @param ... Additional parameters to be passed to fread
fread_multiple <- function(path_vec, colname = "fold", fold_names = NULL, ...) {
  library(data.table)
  ## input check
  if (is.null(fold_names)) fold_names <- path_vec
  stopifnot(length(fold_names)==length(path_vec))
  
  dt_list <- sapply(1:length(path_vec) , function(i) {
    path <- path_vec[i]
    fold_name <- fold_names[i]
    dt <- fread(path, ...)
    dt[,(colname) := fold_name]
    return(dt)
  }, simplify = FALSE)
  
  dt <- dt_list %>% rbindlist
  return(dt)
}

##' Sort the unique vector values decreasingly by their table() values
sort_by_prevalance <- function(f) {
  names(sort(table(f),decreasing=TRUE))
}
