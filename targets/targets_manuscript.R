targets_manuscript <- function() {
  
  selections <- tar_plan(
    genes_selection = c("AZU1", "CTSG", "DEFA3", "TUBB", "ELANE", "PRTN3"),
    supplement_genes = c("CTSC", "C5AR1", "CXCR2"),
    ids_selection = covid$info |> filter(gene_names %in% genes_selection) |> pull(id),
    ids_supplement = covid$info |> filter(gene_names %in% supplement_genes) |> pull(id),
    ids_significant = dl_list_drug_vs_placebo$treatmentBrensocatib
  )
  
  calculations <- tar_plan(
    enr_go = functional_enrichment(as.character(all_ids), as.character(ids_significant), fterms$go, id2gene) |> add_column(ont = "go"),
    enr_re = functional_enrichment(as.character(all_ids), as.character(ids_significant), fterms$re, id2gene) |> add_column(ont = "re"),
    enr_kg = functional_enrichment(as.character(all_ids), as.character(ids_significant), fterms$kg, id2gene) |> add_column(ont = "kg"),
    enr = bind_rows(enr_go, enr_re, enr_kg) |> filter(p_adjust < 0.05) |> unite(term_name, c(ont, term_name), sep = ": ")
  )
  
  figures <- tar_plan(
    pdf_volcano_brensocatib = mn_plot_volcano(dl_drug_vs_placebo, ctr = "treatmentBrensocatib") |> gp("volcano_drug_vs_placebo_29-1", 4, 4),
    pdf_heatmap_sig = mn_plot_diff_heatmap(covid, id_sel = ids_significant) |> gp("heatmap_significant", 2.8, 4),
    pdf_timelines_sel = mn_plot_protein_means(covid, ids_selection) |> gp("timelines_selection", 5, 5),
    pdf_timelines_add = mn_plot_protein_means(covid, ids_supplement) |> gp("timelines_supplement", 5, 3),
    pdf_lograt_sig = mn_plot_lograt_protein(lograt, ids_significant) |> gp("lograt_significant", 6, 9),
    pdf_lograt_add = mn_plot_lograt_protein(lograt, ids_selection) |> gp("lograt_selection", 4.5, 4.5),
    pdf_lograt_sel = mn_plot_lograt_protein(lograt, ids_supplement) |> gp("lograt_supplement", 4.5, 2.8),
    pdf_go_terms = mn_plot_enrichment(enr) |> gp("enrichment_significant", 6.5, 3.8)
  )
  
  c(
    selections,
    calculations,
    figures
  )
  
}