eclip_load_modelling_data <- function(rbp_name) {
  info <<- jsonlite::fromJSON(paste0("data/encode/eclip/processed/design_matrix/meta_info/", rbp_name, ".json"))
  stopifnot(info$rbp == rbp_name)
  ## train the model
  dt_train <<- fread(info$path_train)
  dt_valid <<- fread(info$path_valid)
  dt_train_valid <<- rbind(dt_train, dt_valid)
  dt_test <<- fread(info$path_test)
  y_train <<- dt_train[[info$response]]
  y_valid <<- dt_valid[[info$response]]
  y_train_valid <<- dt_train_valid[[info$response]]
  y_test <<- dt_test[[info$response]]
}

eclip_modelling_get_output_file <- function(script_name, rbp_name,
                                            position_type="no_position") {
  output_dir <- file.path("data/encode/eclip/processed/design_matrix/predictive_models/", tools::file_path_sans_ext(script_name))
  output_file <- file.path(output_dir, paste0(rbp_name, "-", position_type, ".rds"))
  return(output_file)
}

eclip_modelling_get_output_file_csv <- function(script_name, rbp_name,
                                                subinfo="no_position") {
  output_dir <- file.path("data/encode/eclip/processed/predictions/", rbp_name)
  dir.create(output_dir, showWarnings = FALSE)
  output_file <- file.path(output_dir,
			   paste0(tools::file_path_sans_ext(script_name), "-", subinfo, ".csv"))
  return(output_file)
}
