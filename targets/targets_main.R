targets_main <- function() {
  
  biomart <- list(
    tar_target(mart, useEnsembl(biomart = "ensembl", dataset = ENSEMBL_DATASET, version = ENSEMBL_VERSION)),
    #tar_target(bm_genes, biomart_fetch_genes(mart)),
    tar_target(go_terms,  bm_fetch_go(mart, genes)),
    tar_target(re_terms, fetch_reactome(mart, genes))
    #tar_target(kg_terms, get_kegg(species = KEGG_SPECIES, bm_genes = bm_genes)),
    #tar_target(all_terms, list(go = go_terms, gs = gs_terms, re = re_terms, kg = kg_terms))
  )

  read_data <- list(
    tar_target(metadata, make_metadata("info/samples.txt")),
    tar_target(set, read_spectronaut_data(SP_FILE, metadata, min_pep = MIN_PEPTIDES))
  )
  
  selections <- list(
    tar_target(detected_proteins, select_detected_proteins(set)),
    tar_target(genes, detected_proteins$gene_name %>% unique())
  )

  overview <- list(
    tar_target(fig_unique_pep, plot_unique_pep(set)),
    tar_target(fig_detection, plot_detection(set)),
    tar_target(fig_sample_detection, plot_sample_detection(set)),
    tar_target(png_sample_dist, plot_sample_distributions(set, x_breaks = c(-4, 0, 4), text_size = 7) %>% gs("sample_dist", 10, 10)),
    tar_target(png_kernels, plot_kernel_comparison(set) %>% gs("kernels", 6, 3)),
    tar_target(fig_clustering, plot_clustering(set)),
    tar_target(fig_cormat, plot_distance_matrix(set)),
    tar_target(fig_pca, plot_pca(set)),
    tar_target(fig_umap, plot_umap(set, n_neighbours = 20, min_dist = 0.1)),
    tar_target(png_big_heat, plot_big_heatmap(set) %>% gs("big_heatmap", 20, 49))
  )


  c(
    biomart,
    read_data,
    selections,
    overview
  )

}
