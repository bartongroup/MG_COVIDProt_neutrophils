fgsea_run <- function(term_data, res, min.size = 3) {
  res <- res |>
    filter(!is.na(value) & !is.na(feature_id))
  ranks <-  set_names(res$value, res$id)
  fgsea::fgsea(pathways = term_data$term2feature, stats = ranks, nproc = 6, minSize = min.size, eps = 0) |>
    as_tibble() |> 
    rename(term_id = pathway) |> 
    mutate(term_name = term_data$term2name[term_id]) |> 
    arrange(NES) |>
    select(term_id, term_name, p_value = pval, fdr = padj, nes = NES, size, leading_edge = leadingEdge)
}

fgsea_groups <- function(d, term_data, feature_var, value_var, group_var) {
  d |>
    mutate(value = get(value_var), feature_id = get(feature_var)) |>
    group_split(!!sym(group_var)) |>
    map_dfr(function(w) {
      fgsea_run(term_data, w) |>
        mutate(!!group_var := dplyr::first(w[[group_var]]))
    })
}

fgsea_all_terms <- function(d, terms, feature_var = "id", value_var = "logFC", group_var = "contrast") {
  ontologies <- names(terms)
  map(ontologies, function(ont) {
    cat(str_glue("  Computing fgsea for {ont}\n\n"))
    fgsea_groups(d, terms[[ont]], feature_var, value_var, group_var)
  }) |>
    set_names(ontologies)
}

plot_fgsea_enrichment <- function(term_id, res, trms, valvar = "logFC") {
  lst <- trms$term2feature[[term_id]]
  rnks <- set_names(res[[valvar]], res$gene_id)
  fgsea::plotEnrichment(lst, rnks)
}

split_genes_fgsea <- function(se, fg, groupvar = "contrast") {
  fg |> 
    filter(padj < 0.05) |> 
    group_split(term, !!sym(groupvar)) |> 
    map_dfr(function(w) {
      term <- as.character(w$term)
      gr <- as.character(w[[groupvar]])
      genes <- w$leading_edge[[1]]
      se |> 
        filter(gene_id %in% genes & !!sym(groupvar) == gr) |> 
        mutate(term_id = term, .before = "gene_symbol")
    })
}


gsea_de <- function(gse, res, fdr_limit = 0.01) {
  contrasts <- unique(res$contrast) |> 
    as.character()
  ontologies <- names(gse)
  
  map_dfr(ontologies, function(ont) {
    map_dfr(contrasts, function(ctr) {
      res_sel <- res |> 
        filter(contrast == ctr)
      gse_sel <-  gse[[ont]] |> 
        filter(contrast == ctr)
      sig_genes <- res_sel |> 
        filter(sig) |> 
        pull(gene_id)
      
      gse_sel |> 
        filter(fdr < fdr_limit) |> 
        unnest(leading_edge) |> 
        filter(leading_edge %in% sig_genes) |> 
        rename(gene_id = leading_edge) |> 
        left_join(res_sel, by = "gene_symbol") |> 
        select(term_id, term_name, nes, fdr, gene_id, gene_symbol, logFC, logCPM, PValue, FDR) |> 
        add_column(ontology = ont, .before = 1) |> 
        add_column(contrast = ctr)
    })
  })
}


get_terms_str <- function(gso, query, fdr_limit = 0.05) {
  gso |> 
    filter(str_detect(term_name, query) & padj < fdr_limit) |>
    pull(term_id)
}


get_significant_fgsea <- function(gse, da, fdr_limit = 0.05, n_top = 10) {
  map_dfr(names(gse), function(ont) {
    gse[[ont]] |> 
      filter(fdr < fdr_limit) |> 
      unnest(leading_edge) |> 
      mutate(leading_edge = as.integer(leading_edge)) |> 
      left_join(da, by = c("contrast", "leading_edge" = "id")) |> 
      group_by(term_id, term_name, p_value, fdr, nes, size, contrast) |> 
      arrange(-sign(nes) * logFC) |>  
      summarise(top_genes = gene_symbol |> head(n_top) |> str_c(collapse = ", ")) |>
      ungroup() |> 
      add_column(ontology = ont, .before = 1) |> 
      mutate(across(c(p_value, fdr, nes), ~signif(.x, 3))) |> 
      select(ontology, term_id, term_name, fdr, nes, size, contrast, top_genes)
    #mutate(ontology = recode(ontology, go = "GO", re = "Reactome", kg = "KEGG"))
  }) |> 
    arrange(nes)
}