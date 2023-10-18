th <- ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    strip.background = element_rect(fill = "grey90")
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


mn_plot_diff_heatmap <- function(set, what = "abu_norm", max_fc = 2, id_sel = NULL, text.size = 8) {
  
  meta <- set$metadata
  
  d <- set$dat |> 
    mutate(val = get(what)) |>
    filter(id %in% id_sel) |> 
    left_join(set$metadata, by = "sample") |> 
    group_by(id, treatment, day) |> 
    summarise(val = mean(val)) |> 
    pivot_wider(id_cols = c(id, day), names_from = treatment, values_from = val) |> 
    mutate(log_fc = (Brensocatib - Placebo) / log10(2)) |> 
    left_join(set$info, by = "id")

  gene_ord <- d |> 
    filter(day == 29) |> 
    arrange(-log_fc) |> 
    pull(gene_names)
  
  mx <- max(abs(d$log_fc))
  
  d |>  
    ggplot(aes(x = day, y = factor(gene_names, levels = gene_ord), fill = log_fc)) +
    th +
    theme(
      legend.position = "right"
    ) +
    geom_tile() +
    scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-mx, mx),
                         name = expression(log[2]~FC)) +
    labs(x = "Day", y = NULL)
}


mn_plot_protein_means <- function(set, pids, what = "abu_norm", min_n = 3, ncol = 3,
                                  exclude_days = 4, dodge_width = 1) {
  info <- set$info |> 
    filter(id %in% pids) |> 
    mutate(prot = gene_names)
  
  d <- set$dat |> 
    filter(id %in% pids) |> 
    left_join(set$metadata, by = "sample") |> 
    left_join(info, by = "id") |> 
    mutate(val = get(what)) |> 
    filter(!(day %in% exclude_days)) |> 
    select(prot, day, treatment, val) |> 
    mutate(day = as.numeric(as.character(day))) |> 
    group_by(prot, treatment, day) |> 
    summarise(
      m = mean(val),
      s = sd(val),
      n = n(),
      se = s / sqrt(n),
      t_crit = qt(0.975, df = n - 1),
      ci = se * t_crit
    ) |> 
    filter(n >= min_n)
  days <- unique(d$day) |> sort()
  pd <- position_dodge(width = dodge_width)
  d |> 
    ggplot(aes(x = day, y = m, ymin = m - ci, ymax = m + ci, colour = treatment, fill = treatment)) +
    th +
    theme(legend.position = "top") +
    geom_line(alpha = 0.5, position = pd) +
    geom_errorbar(position = pd, width = 1) +
    geom_point(position = pd, shape = 21, colour = "grey30") +
    facet_wrap(~prot, labeller = label_wrap_gen(), scales = "free_y", ncol = ncol) +
    scale_colour_manual(values = okabe_ito_palette, name = "Treatment") +
    scale_fill_manual(values = okabe_ito_palette, name = "Treatment") +
    scale_x_continuous(breaks = days) +
    labs(x = "Day", y = "Normalised abundance")
}


mn_plot_lograt_protein <- function(df, pids, ncol = NULL, what = "logFC", filt = "completion") {
  info <- df$info |> 
    filter(id %in% pids) |> 
    mutate(prot = gene_names)
  d <- df$dat |> 
    filter(id %in% pids) |> 
    left_join(df$metadata, by = "participant_id") |> 
    filter(!bad & !!rlang::parse_expr(filt)) |> 
    mutate(val = get(what)) |> 
    left_join(info, by = "id") |> 
    mutate(x = as_factor(treatment), xi = as.integer(x))
  dm <- d |> 
    group_by(id, prot, xi) |> 
    summarise(M = mean(val))
  gene_ord <- dm |> 
    pivot_wider(id_cols = prot, names_from = xi, values_from = M) |> 
    mutate(diff = `2` - `1`) |> 
    arrange(diff) |> 
    pull(prot)
  d <- d |> 
    mutate(prot = factor(prot, levels = gene_ord))
  dm <- dm |> 
    mutate(prot = factor(prot, levels = gene_ord))
  
  rm(df)
  
  ggplot(d, aes(x = x, y = val)) +
    th +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_hline(yintercept = 0, colour = "grey50") +
    ggbeeswarm::geom_quasirandom(width = 0.2) +
    scale_colour_manual(values = okabe_ito_palette) +
    geom_segment(data = dm, aes(x = xi - 0.3, y = M, xend = xi + 0.3, yend = M), linewidth = 1, colour = "brown") +
    facet_wrap(~ prot, labeller = label_wrap_gen(), ncol = ncol) +
    labs(x = NULL, y = expression(log[2]~I[29]/I[1]))
}


mn_plot_enrichment <- function(enr, fdr_limit = 0.01) {
  enr |> 
    filter(p_adjust < fdr_limit) |> 
    separate_longer_delim(ids, delim = ",") |> 
    mutate(ids = str_remove(ids, " ")) |> 
    group_by(term_id) |> 
    mutate(n_term =  n()) |> 
    ungroup() |> 
    group_by(ids) |> 
    mutate(n_ids =  n()) |> 
    ungroup() |> 
  ggplot(aes(x = fct_reorder(ids, -n_ids), y = fct_reorder(term_name, n_term), size = odds_ratio, fill = log10(p_value))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    geom_point(shape = 21) +
    scale_fill_viridis_c(option = "cividis") +
    labs(x = NULL, y = NULL, size = "Odds ratio", fill = "log P value")
}