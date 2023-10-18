targets_manuscript <- function() {
  
  selections <- tar_plan(
    genes_selection = c("AZU1", "CTSG", "DEFA3", "TUBB", "ELANE", "PRTN3"),
    supplement_genes = c("CTSC", "C5AR1", "CXCR2"),
    ids_selection = covid$info |> filter(gene_names %in% genes_selection) |> pull(id),
    ids_supplement = covid$info |> filter(gene_names %in% supplement_genes) |> pull(id),
    ids_significant = dl_list_drug_vs_placebo$treatmentBrensocatib
  )
  
  calculations <- tar_plan(
    enr_go = functional_enrichment(as.character(all_ids), as.character(ids_significant), fterms$go, id2gene)
  )
  
  figures <- tar_plan(
    pdf_volcano_brensocatib = mn_plot_volcano(dl_drug_vs_placebo, ctr = "treatmentBrensocatib") |> gp("volcano_drug_vs_placebo_29-1", 4, 4),
    pdf_heatmap_sig = mn_plot_diff_heatmap(covid, id_sel = ids_significant) |> gp("heatmap_significant", 3, 4),
    pdf_timelines_sel = mn_plot_protein_means(covid, ids_selection) |> gp("timelines_selection", 5, 5),
    pdf_timelines_add = mn_plot_protein_means(covid, ids_supplement) |> gp("timelines_supplment", 5, 3),
    pdf_lograt_sig = mn_plot_lograt_protein(lograt, ids_significant) |> gp("lograt_significant", 6, 9),
    pdf_lograt_add = mn_plot_lograt_protein(lograt, ids_selection) |> gp("lograt_selection", 4.5, 4.5),
    pdf_lograt_sel = mn_plot_lograt_protein(lograt, ids_supplement) |> gp("lograt_supplement", 4.5, 2.8),
    pdf_go_terms = mn_plot_enrichment(enr_go) |> gp("go_enrichment_significant", 6.5, 3.8)
  )
  
  c(
    selections,
    calculations,
    figures
  )
  
}