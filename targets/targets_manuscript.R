targets_manuscript <- function() {
  
  selections <- tar_plan(
    ids_significant = dl_list_drug_vs_placebo$treatmentBrensocatib,
    genes_to_show = c("AZU1", "CTSG", "DEFA3", "TUBB", "ELANE", "PRTN3"),
    ids_to_show = covid$info |> filter(gene_names %in% genes_to_show) |> pull(id),
    additional_genes = c("CTSC", "C5AR1", "CXCR2"),
    additional_ids = covid$info |> filter(gene_names %in% additional_genes) |> pull(id),
  )
  
  calculations <- tar_plan(
    enr_go = functional_enrichment(as.character(all_ids), as.character(ids_significant), fterms$go, id2gene)
  )
  
  figures <- tar_plan(
    pdf_volcano_brensocatib = mn_plot_volcano(dl_drug_vs_placebo, ctr = "treatmentBrensocatib") |> gp("volcano_drug_vs_placebo_29-1", 4, 4),
    pdf_significant_heatmap = mn_plot_diff_heatmap(covid, id_sel = ids_significant) |> gp("heatmap_significant", 3, 4),
    pdf_sel_prot_timelines = mn_plot_protein_means(covid, ids_to_show) |> gp("timelines_selection", 5, 5),
    pdf_add_prot_timelines = mn_plot_protein_means(covid, additional_ids) |> gp("timelines_supplment", 5, 3),
    pdf_lograt_sig = mn_plot_lograt_protein(lograt, ids_significant) |> gp("lograt_significant", 5, 5),
    pdf_go_terms = mn_plot_enrichment(enr_go) |> gp("go_enrichment_significant", 6.5, 3.8)
  )
  
  c(
    selections,
    calculations,
    figures
  )
  
}