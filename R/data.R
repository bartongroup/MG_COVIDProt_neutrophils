read_spectronaut_data <- function(file, meta, min_pep) {
  d <- read_csv(file, show_col_types = FALSE) %>% 
    rename(id = 1) %>% 
    mutate(id = row_number())  # prefer starting from 1
  info <- d %>% 
    select(id, starts_with("PG")) %>% 
    rename_with(~str_remove(.x, "PG.")) %>% 
    janitor::clean_names()
  # columns are complicated, need to simplify naming
  cols <- tibble(name = colnames(d)) %>% 
    filter(name != "id" & !(str_detect(name, "^PG"))) %>%
    mutate(new = name) %>% 
    mutate(new = str_replace(new, "Razor", "PG.Razor")) %>% 
    mutate(new = str_replace(new, ".PG.", ":"))
  
  d_ren <- d %>% 
    select(c(id, all_of(cols$name))) %>% 
    set_names(c("id", cols$new))
  
  # Coverage is a string, cannot pivot it with numerical values
  cover <- d_ren %>% 
    select(c(id, contains("Coverage"))) %>% 
    pivot_longer(-id, names_to = c("raw_sample", "quantity"), names_pattern = "(.*):(.*)") %>% 
    pivot_wider(id_cols = c(id, raw_sample), names_from = quantity, values_from = value) %>% 
    rename(coverage = Coverage) %>% 
    mutate(coverage = str_remove_all(coverage, "%")) %>% 
    mutate(coverage = str_remove(coverage, ";.+$")) %>% # with multiple hits select the first
    mutate(coverage = as.numeric(coverage) / 100)

  full <- d_ren %>% 
    select(-contains("Coverage")) %>% 
    pivot_longer(-id, names_to = c("raw_sample", "quantity"), names_pattern = "(.*):(.*)") %>% 
    pivot_wider(id_cols = c(id, raw_sample), names_from = quantity, values_from = value) %>% 
    rename(abu = Quantity, abu_norm = Quantity.Normalised, unique_pep = Razor.Unique.Peptides, q_value = Qvalue, single_hit = IsSingleHit) %>% 
    left_join(cover, by = c("id", "raw_sample")) %>% 
    mutate(single_hit = as.logical(single_hit)) %>% 
    mutate(abu = log10(abu), abu_norm = log10(abu_norm)) %>% 
    inner_join(select(meta, sample, raw_sample), by = "raw_sample") %>% 
    select(-raw_sample) %>% 
    mutate(sample = as_factor(sample)) %>% 
    relocate(sample, .after = id)
  
  
  qc <- full %>% 
    select(id, sample, unique_pep, q_value, coverage, single_hit)
  dat <- full %>% 
    filter(unique_pep >= min_pep) %>% 
    select(-c(unique_pep, q_value, coverage, single_hit)) %>% 
    drop_na()
  id_prot_gene <- info %>% 
    select(id, gene_name = genes, uniprot = protein_accessions) %>% 
    separate_rows(c(gene_name, uniprot), sep = ";") %>% 
    mutate(uniprot = str_remove(uniprot, "-\\d+$")) %>%   # redundant isoforms
    distinct()
  
  list(
    info = info,
    id_prot_gene = id_prot_gene,
    qc = qc,
    dat = dat,
    metadata = meta
  )
}


select_detected_proteins <- function(set) {
  set$dat %>% 
    select(id) %>% 
    distinct() %>% 
    left_join(set$id_prot_gene, by = "id")
}


dat2mat <- function(dat, what="abu_norm") {
  dat %>% 
    pivot_wider(id_cols = id, names_from = sample, values_from = !!what) %>% 
    column_to_rownames("id") %>% 
    as.matrix()
}


detection_samples <- function(set) {
  set$dat %>% 
    group_by(id) %>% 
    summarise(n_det = n()) %>% 
    group_by(n_det) %>% 
    tally() %>% 
    arrange(desc(n_det))
}



make_prot_umap <- function(set, what = "abu_norm", seed = 4654,
                           n_neighbours = 15, min_dist = 0.01, n_dim = 10) {
  d <- set$dat %>% 
    mutate(val = get(what)) %>% 
    group_by(id) %>% 
    mutate(fc = val - mean(val)) %>% 
    ungroup()
  
  tab <- dat2mat(d, "fc")
  tab <- tab[apply(tab, 1, function(v) sum(is.na(v)) == 0), ]
  
  set.seed(seed)
  u <- uwot::umap(tab, n_neighbors = n_neighbours, min_dist = min_dist, n_components = n_dim)
  colnames(u) <- paste0("X", 1:n_dim)
  u %>% 
    as_tibble(rownames = "id") %>% 
    mutate(id = as.integer(id))
}

cluster_dim <- function(d, n_clust) { 
  X <- d %>% 
    select(starts_with("X")) %>%
    as.matrix()
  km <- kmeans(X, centers = n_clust)
  d %>% 
    mutate(clust = km$cluster %>% as_factor())
}


