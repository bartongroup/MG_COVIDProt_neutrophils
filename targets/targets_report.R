targets_report <- function() {
  
  overview_figures <- list(
    tar_target(fig_participat_stats, plot_paricipant_stats(participants)),
    
    tar_target(prot_umap, make_prot_umap(set, n_neighbours = 30, min_dist = 0.1)),
    tar_target(fig_detection, plot_detection(set)),
    tar_target(fig_sample_detection, plot_sample_detection(set)),
    tar_target(png_sample_dist, plot_sample_distributions(set, x_lim = c(0, 12), x_breaks = c(0, 6, 12), text_size = 7) |> gs("sample_dist", 10, 10)),
    tar_target(png_kernels, plot_kernel_comparison(set) |> gs("kernels", 6, 3)),
    tar_target(fig_clustering, plot_clustering(set, colour_var = "batch")),
    tar_target(fig_clustering_circular, plot_clustering_circular(set, colour_var = "batch")),
    tar_target(fig_cormat, plot_distance_matrix(set)),
    tar_target(fig_pca, plot_pca(set, colour_var = "batch", shape_var = "day")),
    tar_target(fig_umap, plot_umap(set, n_neighbours = 20, min_dist = 0.1, colour_var = "batch", shape_var = "day")),
    tar_target(fig_participant_1_29, plot_participant_1_29(set$metadata)),
    tar_target(fig_counts_batch, plot_meta_numbers(set$metadata, fill_var = "batch")),
    tar_target(fig_batch_heatmap, plot_batch_heatmap(set$metadata)),
    tar_target(fig_counts_run, plot_meta_numbers(set$metadata, fill_var = "run_index")),
    tar_target(fig_batch_age, plot_batch_dependence(set$metadata, "age_group", BASE_FILTER)),
    tar_target(fig_batch_sex, plot_batch_dependence(set$metadata, "sex", BASE_FILTER))
  )

  
  da_figures <- list(
    tar_target(fig_dpp1, plot_protein(set, 2126)),
    tar_target(fig_ddx58, plot_protein(set, 694)),
    
    tar_target(fig_ma_full, plot_ma(da_full, fdr_limit = FDR_LIMIT) + ylim(-10, 10)),
    tar_target(fig_pdist_full, plot_pdist(da_full)),
    
    tar_target(fig_ma_batch, plot_ma(da_batch, fdr_limit = FDR_LIMIT) + ylim(-10, 10)),
    tar_target(fig_ma_day1, plot_ma(da_day1, fdr_limit = FDR_LIMIT) + ylim(-10, 10)),
    tar_target(fig_ma_day29, plot_ma(da_day29, fdr_limit = FDR_LIMIT) + ylim(-10, 10)),
    tar_target(fig_ma_day29_b4, plot_ma(da_day29_b4, fdr_limit = FDR_LIMIT) + ylim(-10, 10)),
    tar_target(fig_ma_day29_b5, plot_ma(da_day29_b5, fdr_limit = FDR_LIMIT) + ylim(-10, 10)),
    
    tar_target(fig_volcano_dl, plot_volcano(dl_all, fdr_limit = FDR_LIMIT)),
    tar_target(fig_volcano_dl_drug_vs_placebo, plot_volcano(dl_drug_vs_placebo, fdr_limit = FDR_LIMIT)),
    
    tar_target(fig_batch_examples, plot_protein(set, pids = BATCH_EXAMPLES, colour_var = "run_index", shape_var = "batch")),
    tar_target(fig_day29_best, plot_protein(set, pids = dal_day29_treatment_best_ids, ncol = 5)),
    tar_target(fig_day29_best_part1132, plot_protein(set, pids = dal_day29_treatment_best_ids, ncol = 5, sample_sel = c("D5_01-35", "D5_29-18"))),
    tar_target(fig_dl_drug_vs_placebo, plot_lograt_protein(lograt, pids = dal_dl_drug_vs_placebo$treatmentdrug, ncol = 5))
  )
  
  upsets <- list(
    tar_target(ups_day29, list(batch4 = dal_day29_b4$treatmentdrug, batch5 = dal_day29_b5$treatmentdrug, both = dal_day29$treatmentdrug))
  )
  
  numbers_and_tables <- list(
    tar_target(spectronaut_nrow, get_file_nrow(SPECTRONAUT_FILE)),
    tar_target(spectronaut_head, get_file_head(SPECTRONAUT_FILE, n = 100)),
    tar_target(spectronaut_columns, get_file_columns(SPECTRONAUT_FILE)),
    tar_target(max_qvalue, max(set$qc$qvalue, na.rm = TRUE)),
    tar_target(n_all_proteins, nrow(set$info)),
    tar_target(n_hit_proteins, set$dat$id |> unique() |> length()),
    tar_target(n_full_detection, detection_samples(set)$n[1]),
    tar_target(n_samples, nrow(set$metadata)),
    
    tar_target(dal_day29_treatment_best_ids, Reduce(intersect, ups_day29)),
    tar_target(dal_day29_treatment_best, set$info |> filter(id %in% dal_day29_treatment_best_ids)),
    
    tar_target(dal_dl_drug_vs_placebo_proteins, set$info |> filter(id %in% dal_dl_drug_vs_placebo$treatmentdrug)),
    
    tar_target(sav_csv, save_data_csv(set, "tab/full_data.csv"))
  )
  
  
  c(
    da_figures,
    overview_figures,
    upsets,
    numbers_and_tables
  )
}
