require(rlang)
require(purrr)
require(assertthat)
require(readr)

#' Read data file
#'
#' Read an RDS data file.
#'
#' @param file_name File name
#'
#' @return Content of the RDS file
sh_read_file <- function(file_name) {
  file_path <- file.path(CONFIG$data_path, file_name)
  assertthat::assert_that(file.exists(file_path))
  readr::read_rds(file_path)
}


#' Read selected data sets
#'
#' @param sets A character vector with set names, e.g. c("data", "de")
#'
#' @return A list of data
sh_read_data <- function(sets, with_progress = TRUE) {
  step <- 1 / length(sets)
  purrr::map(sets, function(s) {
    fname <- CONFIG$file_names[[s]]
    res <- sh_read_file(fname)
    if (with_progress)
      incProgress(step)
    return(res)
  }) |> 
    rlang::set_names(sets)
}
