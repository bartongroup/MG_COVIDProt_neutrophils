targets_manuscript <- function() {
  
  figures = tar_plan(
    fig_volcano_brensocatib = mn_plot_volcano(dl_drug_vs_placebo, ctr = "treatmentBrensocatib") |> gp("volcano", 4, 4)
  )
  
  c(
    figures
  )
  
}