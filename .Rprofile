## my Rprofile. I will put a link of it also into src/r.

message("Loading packages and configuration from src/r/config.R...")
suppressMessages(source("./Scripts/config.R")) # load your own functions
message("Done")

## set the appropriate path
options(repos=structure(c(CRAN="http://ftp5.gwdg.de/pub/misc/cran/")))
options(browser = "google-chrome")


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

install_bioc <- function(...){
  source("http://bioconductor.org/biocLite.R")
  biocLite(...)
}
