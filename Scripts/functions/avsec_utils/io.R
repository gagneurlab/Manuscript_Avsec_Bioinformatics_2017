## useful read, write functions

##' Read json file in R
##'
##' @param file Path to .json
##' @param fast_and_innacurate If TRUE, use `rjson::fromJSON`, else use `jsonlite::fromJSON(txt = file, simplifyDataFrame = FALSE)`
##' @return List
##' @author Žiga Avsec
read_json <- function(file, fast_and_innacurate = FALSE) {
  if (isTRUE(fast_and_innacurate)) {
    json_data <- rjson::fromJSON(file = file)
  } else {
    json_data <- jsonlite::fromJSON(txt = file, simplifyDataFrame = FALSE)
  }

  return(json_data)
}

#' @param obj List to write
#' @param path File path to write
write_json <- function(obj, path) {
  write(jsonlite::toJSON(obj, pretty = TRUE, auto_unbox = TRUE), path)
}
