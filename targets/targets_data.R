targets_data <- function() {
  annotations <- list(
    tar_target(mart, useEnsembl(biomart = "ensembl", dataset = ENSEMBL_DATASET, version = ENSEMBL_VERSION)),
    tar_target(bm_genes, biomart_fetch_genes(mart) |> filter(chr %in% CHROMOSOMES)),
    tar_target(uni_gene,  download_uniprot_mapping(UNIPROT_MAPPING_FILE)),
    tar_target(terms, get_functional_terms(bm_genes, covid$id_prot_gene, species = SPECIES)),
    tar_target(gterms, prepare_terms_gsea(terms, all_ids)),
    tar_target(fterms, prepare_terms_fenr(terms, all_ids)),
    tar_target(contaminants, read_contaminants(CONTAMINANT_FILE))
  )
  
  read_data <- list(
    tar_target(metadata_covid, read_covid_metadata(COVID_METADATA_FILE, COVID_BAD_SAMPLES)),
    tar_target(covid, read_spectronaut_long_data(SPECTRONAUT_FILE, metadata_covid, uni_gene, contaminants))
  )
  
  selections <- list(
    tar_target(detected_proteins, select_detected_proteins(covid)),
    tar_target(genes, detected_proteins$gene_symbol |> unique()),
    tar_target(all_genes, covid$id_prot_gene$gene_symbol |> unique()),
    tar_target(all_ids, covid$info$id),
    tar_target(participants, collect_participant_stats(covid$metadata))
  )
  
  save_shiny <- list(
    tar_target(save_shiny_terms, write_rds(fterms, "shiny/terms.rds", compress = "xz")),
    tar_target(save_shiny_genes, write_rds(bm_genes, "shiny/genes.rds", compress = "xz")),
    tar_target(save_shiny_da_day29, write_rds(da_day29, "shiny/da_day29.rds", compress = "xz")),
    tar_target(save_shiny_data, write_rds(covid, "shiny/data.rds", compress = "xz"))
  )
  
  
  c(
    annotations,
    read_data,
    selections,
    save_shiny
  )
}
