## global config functions

## prettier print
options(datatable.prettyprint.char = 25)

## not working
options(datatable.print.class = TRUE)

options(repos=structure(c(CRAN="http://ftp5.gwdg.de/pub/misc/cran/")))
options(browser = "google-chrome")

## for emacs
rmarkdown_open_chrome <- function(input, ...){
    ## substitute .R with .html
    html_input <- gsub("\\.R$","\\.html",input)
    browseURL(html_input, browser="google-chrome")
}

rmarkdown_open_firefox <- function(input, ...){
    ## substitute .R with .html
    html_input <- gsub("\\.R$","\\.html",input)
    browseURL(html_input)
}

render_pdf <- function(f){
    rmarkdown::render(f,output_format="pdf_document")
}

## get the bioc functions
install_bioc <- function(...){
  source("https://bioconductor.org/biocLite.R")
  biocLite(...)
}
