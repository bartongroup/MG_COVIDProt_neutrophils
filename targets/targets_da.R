targets_differential_abundance <- function() {
  
  direct_da <- list(
    tar_target(da_full, limma_de_f(set, "~ treatment + day + batch + age_group + sex", filt = BASE_FILTER)),
    tar_target(da_batch, limma_de(set, contrasts = c("x5-x4", "x5-x3", "x4-x3"), group_var = "batch", filt = BASE_FILTER)),
    tar_target(da_day1, limma_de_f(set, "~ treatment + batch + age_group + sex", filt = paste(BASE_FILTER, "& day == 1"))),
    tar_target(da_day29, limma_de_f(set, "~ treatment + batch + age_group + sex", filt = paste(BASE_FILTER, "& day == 29"))),
    tar_target(da_day29_b4, limma_de_f(set, "~ treatment + age_group + sex", filt = paste(BASE_FILTER, "& batch == 4 & day == 29"))),
    tar_target(da_day29_b5, limma_de_f(set, "~ treatment + age_group + sex", filt = paste(BASE_FILTER, "& batch == 5 & day == 29")))
  )
  
  lograt_da <- list(
    tar_target(lograt, logfc_days(set)),
    tar_target(dl_drug, limma_de_ratio(lograt, filt = "treatment == 'drug'")),
    tar_target(dl_placebo, limma_de_ratio(lograt, filt = "treatment == 'placebo'")),
    tar_target(dl_all, bind_rows(dl_drug |> mutate(contrast = "drug"), dl_placebo |> mutate(contrast = "placebo"))),
    tar_target(dl_drug_vs_placebo, limma_de_f(lograt, formula = "~ treatment + age_group + sex", what = "logFC", names = "participant_id", filt = "completion"))
  )
  
  da_lists <- list(
    tar_target(dal_day29, de_list(da_day29, "contrast", split_up_down = FALSE, fdr_limit = FDR_LIMIT)),
    tar_target(dal_day29_b4, de_list(da_day29_b4, "contrast", split_up_down = FALSE, fdr_limit = FDR_LIMIT)),
    tar_target(dal_day29_b5, de_list(da_day29_b5, "contrast", split_up_down = FALSE, fdr_limit = FDR_LIMIT)),
    tar_target(dal_dl_drug_vs_placebo, de_list(dl_drug_vs_placebo, "contrast", split_up_down = FALSE, fdr_limit = FDR_LIMIT))
  )
  
  fgsea <- list(
    tar_target(gse_da_day29, fgsea_all_terms(add_genes(da_day29, set$info), all_terms, prefix = "d29")),
    tar_target(gse_dl, fgsea_all_terms(add_genes(dl_all, set$info), all_terms, prefix = "dd")),
    tar_target(gse_batch5, fgsea_all_terms(add_genes(da_full |> filter(contrast == "batch5"), set$info), all_terms, prefix = "batch5")),
    
    tar_target(gse_da_day29_sig, get_significant_fgsea(gse_da_day29, da_day29, fdr_limit = 0.05) |> filter(!str_detect(contrast, "batch"))),
    tar_target(png_gse_examples, plot_volcano_enrichment(gse_da_day29_sig, da_day29, all_terms, ncol = 3) |> gs("gse_examples", 12, 25))
  )
  
  
  c(
    direct_da,
    lograt_da,
    da_lists,
    fgsea
  )
  
  
  
}
