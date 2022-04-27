# DE for selected contrasts; if not specified, all pairs of contrasts 
limma_de <- function(set, contrasts = NULL, group_var = "treatment", what = "abu_norm", filt = "TRUE") {
  meta <- set$metadata %>% 
    filter(!bad & !!rlang::parse_expr(filt)) %>% 
    mutate(group = get(group_var)) %>% 
    filter(!is.na(group)) %>% 
    droplevels()
  groups <- unique(as.character(meta$group)) %>% 
    janitor::make_clean_names()
  design_mat <- model.matrix(~ 0 + group, data = meta)
  colnames(design_mat) <- groups
  
  tab <- dat2mat(set$dat, what)[, meta$sample] / log10(2)
  
  if (is.null(contrasts)) {
    contrasts <- expand_grid(x = as_factor(groups), y = as_factor(groups)) %>%
      filter(as.integer(x) < as.integer(y)) %>%
      unite(contrast, c(y, x), sep = "-") %>%
      pull(contrast)
  }
  contrast_mat <- limma::makeContrasts(contrasts = contrasts, levels = design_mat)

  fit <- tab %>%
    limma::lmFit(design_mat) %>%
    limma::contrasts.fit(contrasts = contrast_mat) %>%
    limma::eBayes()
  
  map_dfr(contrasts, function(ctr) {
    limma::topTable(fit, coef = ctr, number = 1e6, sort.by = "none") %>%
      as_tibble(rownames = "id") %>%
      mutate(across(c(id), as.integer)) %>% 
      rename(FDR = adj.P.Val, PValue = P.Value) %>%
      select(-c(t, B)) %>% 
      add_column(contrast = ctr)
  }) %>% 
    drop_na() %>% 
    mutate(contrast = factor(contrast, levels = contrasts))
}

# DE with formula
limma_de_f <- function(set, formula, what = "abu_norm", filt = "TRUE") {
  tab <- dat2mat(set$dat, what) / log10(2)
  meta <- set$metadata %>% 
    filter(!bad & !!rlang::parse_expr(filt)) %>% 
    droplevels()
  tab <- tab[, meta$sample]
  
  design_mat <- model.matrix(as.formula(formula), data = meta)
  coefs <- colnames(design_mat)[-1]
  
  fit <- tab %>%
    limma::lmFit(design_mat) %>%
    limma::eBayes()
  
  res <- map_dfr(coefs, function(cf) {
    limma::topTable(fit, coef = cf, number = 1e6, sort.by = "none") %>%
      as_tibble(rownames = "id") %>%
      mutate(id = as.integer(id)) %>% 
      mutate(contrast = cf) %>%
      rename(FDR = adj.P.Val, PValue = P.Value) %>%
      select(-c(t, B))
  }) %>% 
    drop_na() %>% 
    mutate(contrast = factor(contrast, levels = coefs))
  
  res
}

# one-sample limma against zero
limma_de_ratio <- function(df, what = "logFC", id_var = "participant_id", filt = "TRUE") {
  meta <- df$metadata %>% 
    filter(!bad & !!rlang::parse_expr(filt)) %>% 
    droplevels()
  
  tab <- dat2mat(df$dat, what = what, names = id_var)
  tab <- tab[, as.character(meta[[id_var]])]
  
  design_mat <- cbind(Intercept = rep(1, ncol(tab)))
  fit <- tab %>%
    limma::lmFit(design_mat) %>%
    limma::eBayes()
  
  res <- limma::topTable(fit, number = 1e6, sort.by = "none") %>%
      as_tibble(rownames = "id") %>%
      mutate(id = as.integer(id)) %>% 
      rename(FDR = adj.P.Val, PValue = P.Value) %>%
      select(-c(t, B)) %>% 
      add_column(contrast = "ratio") %>% 
      drop_na()
  
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
