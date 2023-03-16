write_rds_name <- function(obj) {
  path <- file.path("shiny", "data")
  if (!dir.exists(path)) dir.create(path)
  obj_name <- deparse(substitute(obj))
  file_name <- file.path(path, str_glue("{obj_name}.rds"))
  write_rds(obj, file_name, compress = "xz")
}

save_data_for_shiny <- function(dataset, de, terms, fterms, gse) {
  
  gns <- dataset$info |> 
    select(id, description = protein_descriptions) |> 
    mutate(description = str_remove(description, ";.*$"))
  
  de <- de |> 
    rename(p_value = PValue, fdr = FDR, log_fc = logFC, log_exp = AveExpr)
  data <- dataset$dat
  metadata <- dataset$metadata
  
  features <- de |>
    select(id, name = gene_symbol) |>
    distinct() |> 
    left_join(gns, by = "id")
  
  write_rds_name(data)
  write_rds_name(metadata)
  write_rds_name(features)
  write_rds_name(de)
  write_rds_name(fterms)
  write_rds_name(gse)
  
  return(1)
}
