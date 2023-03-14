okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


sh_enrichment <- function(genes_all, genes_sel, term_data, gene2name = NULL,
                                 min_count = 3, sig_limit = 0.05) {
  gene2term <- term_data$gene2term
  term_info <- term_data$terms
  
  # select only terms represented in our gene set
  gene2term <- gene2term |>
    filter(gene_symbol %in% genes_all)
  
  # all terms present in the selection
  terms <- gene2term |>
    filter(gene_symbol %in% genes_sel) |>
    pull(term_id) |>
    unique()
  
  # number of selected genes
  Nsel <- length(genes_sel)
  # size of the universe
  Nuni <- length(genes_all)
  
  # empty line for missing terms
  na_term <- term_info |> slice(1) |> mutate_all(~NA)
  
  res <- map_dfr(terms, function(term) {
    info <- term_info |>
      filter(term_id == term)
    # returns NAs if no term found
    if (nrow(info) == 0) info <- na_term |> mutate(term_id = term)
    
    # all genes with the term
    tgenes <- gene2term |>
      filter(term_id == term) |>
      pull(gene_symbol)
    # genes from selection with the term
    tgenes_sel <- intersect(tgenes, genes_sel)
    
    nuni <- length(tgenes)
    nsel <- length(tgenes_sel)
    
    expected <- nuni * Nsel / Nuni
    fish <- matrix(c(nsel, nuni - nsel, Nsel - nsel, Nuni + nsel - Nsel - nuni), nrow = 2)
    ft <- fisher.test(fish, alternative = "greater")
    p <- as.numeric(ft$p.value)
    
    if (!is.null(gene2name)) tgenes_sel <- gene2name[tgenes_sel] |> unname()
    
    bind_cols(
      info,
      tibble(
        tot = nuni,
        sel = nsel,
        expect = expected,
        enrich = nsel / expected,
        ids = paste(tgenes_sel, collapse = ","),
        P = p
      )
    )
  }) |>
    mutate(P = p.adjust(P, method = "BH")) |>
    filter(sel >= min_count & P <= sig_limit) |>
    arrange(desc(enrich)) |>
    mutate(enrich = round(enrich, 1), expect = round(expect, 2))
  
  res
}



sh_plot_protein <- function(set, pids, what = "abu_norm", colour_var = "batch") {
  info <- set$info |> 
    filter(id %in% pids) |> 
    mutate(protein_descriptions = str_remove(protein_descriptions, ";.+$")) |> 
    mutate(prot = glue::glue("{gene_names}: {protein_descriptions}"))
  d <- set$dat |>
    mutate(val = get(what)) |> 
    filter(id %in% pids) |>
    left_join(set$metadata, by = "sample") |> 
    rename(colvar = !!colour_var) |> 
    arrange(treatment, day) |> 
    unite(x, c(treatment, day)) |> 
    mutate(x = as_factor(x), xi = as.integer(x)) |> 
    left_join(info, by = "id") |> 
    select(id, prot, x, xi, val, colvar)
  dm <- d |> 
    group_by(id, prot, xi) |> 
    summarise(M = mean(val))
  rm(set)
  ggplot(d, aes(x = x, y = val, colour = colvar)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "bottom"
    ) +
    scale_colour_manual(values = okabe_ito_palette, name = colour_var) +
    ggbeeswarm::geom_quasirandom(width = 0.2, size = 1, alpha = 0.8) +
    geom_segment(data = dm, aes(x = xi - 0.3, y = M, xend = xi + 0.3, yend = M), linewidth = 1, colour = "brown") +
    labs(x = NULL, y = what)
}


sh_plot_volma <- function(res, point_size, point_alpha) {
  ggplot(res, aes(x = x, y = y, colour = sel)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    geom_hline(yintercept = 0, colour = "darkgreen") +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_colour_manual(values = c("grey70", "black"))
}

sh_plot_ma <- function(res,  point_size = 0.5, point_alpha = 0.5) {
  res |> 
    sh_plot_volma(point_size, point_alpha) +
    geom_hline(yintercept = 0, linewidth = 0.1, alpha = 0.5) +
    labs(x = expression(log[10]~Intensity), y = expression(log[2]~FC))
}

sh_plot_volcano <- function(res, point_size = 0.5, point_alpha = 0.5) {
  res |> 
    sh_plot_volma(point_size, point_alpha) +
    geom_vline(xintercept = 0, linewidth = 0.1, alpha = 0.5) +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)))
}
