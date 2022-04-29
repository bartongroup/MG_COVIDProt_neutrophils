targets_main <- function() {
  
  biomart <- list(
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
    tar_target(genes, detected_proteins$gene_name %>% unique()),
    tar_target(all_genes, set$id_prot_gene$gene_name %>% unique())
  )

  overview <- list(
    tar_target(prot_umap, make_prot_umap(set, n_neighbours = 30, min_dist = 0.1)),
    tar_target(fig_detection, plot_detection(set)),
    tar_target(fig_sample_detection, plot_sample_detection(set)),
    tar_target(png_sample_dist, plot_sample_distributions(set, x_lim = c(0, 12), x_breaks = c(0, 6, 12), text_size = 7) %>% gs("sample_dist", 10, 10)),
    tar_target(png_kernels, plot_kernel_comparison(set) %>% gs("kernels", 6, 3)),
    tar_target(fig_clustering, plot_clustering(set, colour_var = "batch")),
    tar_target(fig_cormat, plot_distance_matrix(set)),
    tar_target(fig_pca, plot_pca(set, colour_var = "batch", shape_var = "day")),
    tar_target(fig_umap, plot_umap(set, n_neighbours = 20, min_dist = 0.1, colour_var = "batch", shape_var = "day")),
    tar_target(fig_participants, plot_participants(set$metadata)),
    tar_target(fig_counts, plot_meta_numbers(set$metadata))
    #tar_target(png_big_heat, plot_big_heatmap(set) %>% gs("big_heatmap", 20, 49))
  )
  
  differential_abundance <- list(
    tar_target(base_filter, "completion & batch %in% c(3, 4, 5)"),
    tar_target(da_full, limma_de_f(set, "~ treatment + day + batch + age_group + sex", filt = base_filter)),
    tar_target(da_batch, limma_de(set, contrasts = c("x5-x4", "x5-x3", "x4-x3"), group_var = "batch")),
    tar_target(da_day1, limma_de_f(set, "~ treatment + batch + age_group + sex", filt = paste(base_filter, "& day == 1"))),
    tar_target(da_day29, limma_de_f(set, "~ treatment + batch + age_group + sex", filt = paste(base_filter, "& day == 29"))),
    tar_target(dal_day29, de_list(da_day29, "contrast", split_up_down = TRUE)),
    
    tar_target(fig_ma_full, plot_ma(da_full) + ylim(-10, 10)),
    tar_target(fig_ma_batch, plot_ma(da_batch) + ylim(-10, 10)),
    tar_target(fig_ma_day1, plot_ma(da_day1) + ylim(-10, 10)),
    tar_target(fig_ma_day29, plot_ma(da_day29) + ylim(-10, 10)),
    
    tar_target(lograt, logfc_days(set)),
    tar_target(dd_drug, limma_de_ratio(lograt, filt = "treatment == 'drug'")),
    tar_target(dd_placebo, limma_de_ratio(lograt, filt = "treatment == 'placebo'")),
    tar_target(dd_all, bind_rows(dd_drug %>% mutate(contrast = "drug"), dd_placebo %>% mutate(contrast = "placebo"))),
    tar_target(fig_volcano_dd, plot_volcano(dd_all)),

    tar_target(fig_batch_examples, plot_protein(set, pids = BATCH_EXAMPLES))
  )
  
  
  fgsea <- list(
    tar_target(gse_dd, fgsea_all_terms(add_genes(dd_all, set$info), all_terms, prefix = "dd")),
    tar_target(gse_batch5, fgsea_all_terms(add_genes(da_full %>% filter(contrast == "batch5"), set$info), all_terms, prefix = "batch5"))
  )

  profiles <- list(
    tar_target(time_prof, make_time_profiles(set, treat = "drug")),
    tar_target(time_prof_cl6, cluster_time_profiles(time_prof, n_clust = 6)),
    tar_target(fig_profiles_cl6, plot_cluster_profiles(time_prof_cl6)),
    tar_target(fig_prof_heatmap, plot_prof_heatmap(time_prof)),
    tar_target(fig_gradients_eigen, plot_gradients(time_prof))
  )
  
  for_report <- list(
    tar_target(spectronaut_columns, get_file_columns(SPECTRONAUT_FILE)),
    tar_target(max_qvalue, max(set$qc$qvalue, na.rm = TRUE)),
    tar_target(n_all_proteins, nrow(set$info)),
    tar_target(n_hit_proteins, set$dat$id %>% unique() %>% length()),
    tar_target(n_full_detection, detection_samples(set)$n[1]),
    tar_target(n_samples, nrow(set$metadata))
  )
  
  save_shiny <- list(
    tar_target(save_shiny_terms, write_rds(all_terms, "shiny/terms.rds", compress = "xz")),
    tar_target(save_shiny_genes, write_rds(bm_genes, "shiny/genes.rds", compress = "xz")),
    tar_target(save_shiny_da_day29, write_rds(da_day29 %>% add_genes(set$info), "shiny/da_day29.rds", compress = "xz")),
    tar_target(save_shiny_data, write_rds(set, "shiny/data.rds", compress = "xz"))
  )

  c(
    biomart,
    read_data,
    selections,
    overview,
    for_report,
    differential_abundance,
    fgsea,
    save_shiny,
    profiles
  )

}
