targets_data <- function() {
  annotations <- list(
    tar_target(mart, useEnsembl(biomart = "ensembl", dataset = ENSEMBL_DATASET, version = ENSEMBL_VERSION)),
    tar_target(bm_genes, biomart_fetch_genes(mart)),
    tar_target(go_terms,  bm_fetch_go(mart, all_genes)),
    tar_target(re_terms, fetch_reactome(mart, all_genes)),
    tar_target(kg_terms, get_kegg(species = KEGG_SPECIES, bm_genes = bm_genes)),
    tar_target(uni_gene,  download_uniprot_mapping(UNIPROT_MAPPING_FILE)),
    tar_target(all_terms, list(go = go_terms, re = re_terms, kg = kg_terms))
  )
  
  read_data <- list(
    tar_target(set, read_spectronaut_long_data(SPECTRONAUT_FILE, METADATA_FILE, uni_gene)),
    tar_target(metadata, set$metadata)
  )
  
  selections <- list(
    tar_target(detected_proteins, select_detected_proteins(set)),
    tar_target(genes, detected_proteins$gene_name |> unique()),
    tar_target(all_genes, set$id_prot_gene$gene_name |> unique()),
    tar_target(participants, collect_participant_stats(set$metadata))
  )
  
  save_shiny <- list(
    tar_target(save_shiny_terms, write_rds(all_terms, "shiny/terms.rds", compress = "xz")),
    tar_target(save_shiny_genes, write_rds(bm_genes, "shiny/genes.rds", compress = "xz")),
    tar_target(save_shiny_da_day29, write_rds(da_day29, "shiny/da_day29.rds", compress = "xz")),
    tar_target(save_shiny_data, write_rds(set, "shiny/data.rds", compress = "xz"))
  )
  
  
  c(
    annotations,
    read_data,
    selections,
    save_shiny
  )
}
