read_spectronaut_data <- function(file, meta, min_pep) {
  d <- read_csv(file, show_col_types = FALSE) %>% 
    rename(id = 1) %>% 
    mutate(id = row_number())  # prefer starting from 1
  info <- d %>% 
    select(id, starts_with("PG")) %>% 
    rename_with(~str_remove(.x, "PG."))
  # columns are complicated, need to simplify naming
  cols <- tibble(name = colnames(d)) %>% 
    filter(name != "id" & !(str_detect(name, "^PG")) & !str_detect(name, "Coverage")) %>%
    mutate(new = name) %>% 
    mutate(new = str_replace(new, "Razor", "PG.Razor")) %>% 
    mutate(new = str_replace(new, ".PG.", ":"))
  
  full <- d %>% 
    select(c(id, all_of(cols$name))) %>% 
    set_names(c("id", cols$new)) %>% 
    pivot_longer(-id, names_to = c("raw_sample", "quantity"), names_pattern = "(.*):(.*)") %>% 
    pivot_wider(id_cols = c(id, raw_sample), names_from = quantity, values_from = value) %>% 
    rename(abu = Quantity, abu_norm = Quantity.Normalised, unique_pep = Razor.Unique.Peptides, q_value = Qvalue, single_hit = IsSingleHit) %>% 
    mutate(single_hit = as.logical(single_hit)) %>% 
    mutate(abu = log10(abu), abu_norm = log10(abu_norm)) %>% 
    left_join(select(meta, sample, raw_sample), by = "raw_sample") %>% 
    select(-raw_sample) %>% 
    mutate(sample = as_factor(sample)) %>% 
    relocate(sample, .after = id)
  
  qc <- full %>% 
    select(id, sample, unique_pep, q_value, single_hit)
  dat <- full %>% 
    filter(unique_pep >= min_pep) %>% 
    select(-c(unique_pep, q_value, single_hit)) %>% 
    drop_na()
  id_prot_gene <- info %>% 
    select(id, gene_name = Genes, uniprot = ProteinAccessions) %>% 
    separate_rows(c(gene_name, uniprot), sep = ";") %>% 
    mutate(uniprot = str_remove(uniprot, "-\\d+$")) %>%   # redundand isoforms
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

