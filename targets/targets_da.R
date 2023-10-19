targets_differential_abundance <- function() {
  
  direct_da <- list(
    tar_target(da_full, limma_de_f(covid, "~ treatment + day + batch + age_group + sex", filt = BASE_FILTER)),
    tar_target(da_batch, limma_de(covid, contrasts = c("x5-x4", "x5-x3", "x4-x3"), group_var = "batch", filt = BASE_FILTER)),
    tar_target(da_day1, limma_de_f(covid, "~ treatment + batch + age_group + sex", filt = paste(BASE_FILTER, "& day == 1"))),
    tar_target(da_day29, limma_de_f(covid, "~ treatment + batch + age_group + sex", filt = paste(BASE_FILTER, "& day == 29"))),
    tar_target(da_day29_b4, limma_de_f(covid, "~ treatment + age_group + sex", filt = paste(BASE_FILTER, "& batch == 4 & day == 29"))),
    tar_target(da_day29_b5, limma_de_f(covid, "~ treatment + age_group + sex", filt = paste(BASE_FILTER, "& batch == 5 & day == 29")))
  )
  
  lograt_da <- list(
    tar_target(lograt, logfc_days(covid)),
    tar_target(dl_drug, limma_de_ratio(lograt, what = "logFC_quant", filt = "treatment == 'Brensocatib'")),
    tar_target(dl_placebo, limma_de_ratio(lograt, what = "logFC_quant", filt = "treatment == 'Placebo'")),
    tar_target(dl_all, bind_rows(dl_drug |> mutate(contrast = "Brensocatib"), dl_placebo |> mutate(contrast = "Placebo"))),
    tar_target(dl_drug_vs_placebo, limma_de_f(lograt, formula = "~ treatment + age_group + sex", what = "logFC_quant",
                                              names = "participant_id", filt = "completion", convert_log = FALSE))
  )
  
  da_lists <- list(
    tar_target(da_list_day29, de_list(da_day29, "contrast", split_up_down = FALSE, fdr_limit = FDR_LIMIT)),
    tar_target(da_list_day29_b4, de_list(da_day29_b4, "contrast", split_up_down = FALSE, fdr_limit = FDR_LIMIT)),
    tar_target(da_list_day29_b5, de_list(da_day29_b5, "contrast", split_up_down = FALSE, fdr_limit = FDR_LIMIT)),
    tar_target(dl_list_drug_vs_placebo, de_list(dl_drug_vs_placebo, "contrast", split_up_down = FALSE, fdr_limit = FDR_LIMIT))
  )
  
  fgsea <- list(
    tar_target(gse_da_day29, fgsea_all_terms(da_day29, gterms)),
    tar_target(gse_dl, fgsea_all_terms(dl_all, gterms)),
    tar_target(gse_dl_drug_vs_placebo, fgsea_all_terms(dl_drug_vs_placebo, gterms)),
    tar_target(gse_batch5, fgsea_all_terms(add_genes(da_full |> filter(contrast == "batch5"), covid$info), gterms)),
    
    tar_target(gse_da_day29_sig, get_significant_fgsea(gse_da_day29, da_day29, fdr_limit = 0.05) |> filter(!str_detect(contrast, "batch"))),
    tar_target(png_gse_examples, plot_volcano_enrichment(gse_da_day29_sig, da_day29, gterms, ncol = 3) |> gs("gse_examples", 12, 25))
  )
  
  
  c(
    direct_da,
    lograt_da,
    da_lists,
    fgsea
  )
  
  
  
}
