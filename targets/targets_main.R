targets_main <- function() {

  read_data <- list(
    tar_target(metadata, make_metadata("info/samples.txt")),
    tar_target(set, read_spectronaut_data(SP_FILE, metadata))
  )



  c(
    read_data
  )

}
