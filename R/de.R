limma_de <- function(set, formula = "~ 0 + group", what = "abu_norm", info_cols = NULL) {
  X <- dat2mat(set$dat, what)

  meta <- set$metadata
  conditions <- unique(meta$group)

  design <- model.matrix(as.formula(formula), meta)
  colnames(design) <- conditions

  ctrs <- expand_grid(x = as_factor(conditions), y = as_factor(conditions)) %>%
    filter(as.integer(x) < as.integer(y)) %>%
    unite(contrast, c(y, x), sep = "-") %>%
    pull(contrast)
  contrast_mat <- makeContrasts(contrasts = ctrs, levels = design)

  fit <- lmFit(X, design) %>%
    contrasts.fit(contrasts = contrast_mat) %>%
    eBayes()


  res <- map_dfr(ctrs, function(cf) {
    topTable(fit, coef = cf, number = 1e6, sort.by = "none") %>%
      as_tibble(rownames = "id") %>%
      mutate(across(c(id), as.integer)) %>% 
      mutate(contrast = cf) %>%
      rename(FDR = adj.P.Val, PValue = P.Value) %>%
      select(-c(t, B))
  }) %>%
    drop_na() %>%
    mutate(contrast = factor(contrast, levels = ctrs))

  if (!is.null(info_cols)) {
    info <- set$info %>%
      select(c("id", all_of(info_cols)))
    res <- res %>%
      left_join(info, by = "id")
  }
  res
}

# DE with formula
limma_de_f <- function(set, formula, what = "abu_norm", info_cols = NULL) {
  tab <- dat2mat(set$dat, what)
  meta <- set$metadata %>% 
    drop_na()
  tab <- tab[, meta$sample]
  
  design_mat <- model.matrix(as.formula(formula), data = set$metadata)
  coefs <- colnames(design_mat)[-1]
  
  fit <- tab %>%
    lmFit(design_mat) %>%
    eBayes()
  
  res <- map_dfr(coefs, function(cf) {
    topTable(fit, coef = cf, number = 1e6, sort.by = "none") %>%
      as_tibble(rownames = "id") %>%
      mutate(id = as.integer(id)) %>% 
      mutate(contrast = cf) %>%
      rename(FDR = adj.P.Val, PValue = P.Value) %>%
      select(-c(t, B))
  }) %>% 
    drop_na() %>% 
    mutate(contrast = factor(contrast, levels = coefs))
  
  if (!is.null(info_cols)) {
    info <- set$info %>%
      select(c("id", all_of(info_cols)))
    res <- res %>%
      left_join(info, by = "id")
  }
  res
}






de_list <- function(res, group_var, fdr = "FDR", logfc = "logFC", logfc_limit = 0, fdr_limit = 0.05, name = NULL, split_up_down = FALSE) {
  d <- res %>%
    filter(!!sym(fdr) < fdr_limit & abs(!!sym(logfc)) >= logfc_limit) %>%
    mutate(group = !!sym(group_var))
  if (split_up_down) {
    d <- d %>%
      mutate(direction = if_else(!!sym(logfc) > 0, "up", "down")) %>%
      unite(group, c(group, direction), sep = ":")
  }
  d <- d %>%
    select(group, id) %>%
    group_by(group)
  kname <- ifelse(is.null(name), "", paste0(name, ":"))
  ks <- paste0(kname, group_keys(d)[[1]])
  d %>%
    distinct() %>%
    group_map(~pull(.x, id)) %>%
    set_names(ks)
}


pull_proteins <- function(des) {
  des %>%
    select(protein, gene_name) %>%
    distinct() %>%
    arrange(gene_name) %>%
    filter(!is.na(gene_name))
}

make_de_genes <- function(de, fdr_limit = 0.05, logfc_limit = 1) {
  list(
    up = pull_proteins(de %>% filter(FDR < fdr_limit & logFC >= logfc_limit)),
    down = pull_proteins(de %>% filter(FDR < fdr_limit & logFC <= -logfc_limit))
  )
}

make_de_table <- function(de, info) {
  de %>%
    #filter(FDR < fdr_limit & abs(logFC) >= logfc_limit) %>%
    select(id, multi, logFC, PValue, FDR) %>%
    left_join(info, by = "id") %>%
    select(id, multi, logFC, FDR, proteins, protein_names, gene_name, aa = amino_acid, in_pep = position_in_peptide, in_prot = position, charge) %>%
    mutate(across(c(logFC, FDR), ~signif(.x, 3)))
}
