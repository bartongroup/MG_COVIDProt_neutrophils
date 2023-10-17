th <- ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank()
  )

gp <- function(gg, name, width, height) {
  if(!dir.exists("manuscript_fig")) dir.create("manuscript_fig")
  ggsave(filename = file.path("manuscript_fig", paste0(name, ".pdf")), plot = gg, device = "pdf",
         width = width, height = height)
}


mn_plot_volcano <- function(de, ctr, fdr_limit = 0.01) {
  d <- de |> 
    filter(contrast == ctr) |> 
    mutate(logp = -log10(PValue))
  d_nsig <- d |> 
    filter(FDR >= fdr_limit)
  d_sig <- d |> 
    filter(FDR < fdr_limit)
  
  ggplot() +
    th +
    geom_point(aes(x = logFC, y = logp), data = d_nsig, size = 0.3, colour = "grey70") +
    geom_point(aes(x = logFC, y = logp), data = d_sig, size = 0.4, colour = "black") +
    geom_vline(xintercept = 0, linewidth = 0.1, alpha = 0.5) +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)))
}