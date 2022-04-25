targets_main <- function() {
  
  biomart <- list(
    tar_target(mart, useEnsembl(biomart = "ensembl", dataset = ENSEMBL_DATASET, version = ENSEMBL_VERSION)),
    tar_target(bm_genes, biomart_fetch_genes(mart)),
    tar_target(go_terms,  bm_fetch_go(mart, all_genes)),
    tar_target(re_terms, fetch_reactome(mart, all_genes)),
    tar_target(kg_terms, get_kegg(species = KEGG_SPECIES, bm_genes = bm_genes))
    #tar_target(all_terms, list(go = go_terms, gs = gs_terms, re = re_terms, kg = kg_terms))
  )

  read_data <- list(
    tar_target(set, read_spectronaut_long_data(SPECTRONAUT_FILE, METADATA_FILE, UNIPROT_MAPPING_FILE)),
    tar_target(metadata, set$metadata)
  )
  
  selections <- list(
    tar_target(detected_proteins, select_detected_proteins(set)),
    tar_target(genes, detected_proteins$gene_name %>% unique()),
    tar_target(all_genes, set$id_prot_gene$gene_name %>% unique())
  )

  overview <- list(
    tar_target(fig_detection, plot_detection(set)),
    tar_target(fig_sample_detection, plot_sample_detection(set)),
    tar_target(png_sample_dist, plot_sample_distributions(set, x_lim = c(0, 12), x_breaks = c(0, 6, 12), text_size = 7) %>% gs("sample_dist", 10, 10)),
    tar_target(png_kernels, plot_kernel_comparison(set) %>% gs("kernels", 6, 3)),
    tar_target(fig_clustering, plot_clustering(set, colour_var = "batch")),
    tar_target(fig_cormat, plot_distance_matrix(set)),
    tar_target(fig_pca, plot_pca(set, colour_var = "batch", shape_var = "day")),
    tar_target(fig_umap, plot_umap(set, n_neighbours = 20, min_dist = 0.1, colour_var = "batch", shape_var = "day")),
    tar_target(png_big_heat, plot_big_heatmap(set) %>% gs("big_heatmap", 20, 49))
  )
  
  differential_abundance <- list(
    tar_target(da, limma_de_f(set, "~ treatment + day + batch")),
    tar_target(fig_volcano, plot_volcano(da, fdr_limit = FDR_LIMIT)),
    tar_target(fig_pdist, plot_pdist(da))
  )

  profiles <- list(
    tar_target(prot_umap, make_prot_umap(set, n_neighbours = 30, min_dist = 0.1)),
    tar_target(time_prof, make_time_profiles(set, treat = "drug")),
    tar_target(time_prof_cl6, cluster_time_profiles(time_prof, n_clust = 6)),
    tar_target(fig_profiles_cl6, plot_cluster_profiles(time_prof_cl6)),
    tar_target(fig_prof_heatmap, plot_prof_heatmap(time_prof)),
    tar_target(fig_gradients_eigen, plot_gradients(time_prof))
  )
  
  for_report <- list(
    tar_target(n_all_proteins, nrow(set$info)),
    tar_target(n_hit_proteins, set$dat$id %>% unique() %>% length()),
    tar_target(n_full_detection, detection_samples(set)$n[1]),
    tar_target(n_samples, nrow(set$metadata))
  )

  c(
    biomart,
    read_data,
    selections,
    overview,
    for_report,
    differential_abundance,
    profiles
  )

}
