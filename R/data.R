read_spectronaut_data <- function(file, meta) {
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
  
  dat <- d %>% 
    select(c(id, all_of(cols$name))) %>% 
    set_names(c("id", cols$new)) %>% 
    pivot_longer(-id, names_to = c("raw_sample", "quantity"), names_pattern = "(.*):(.*)") %>% 
    pivot_wider(id_cols = c(id, raw_sample), names_from = quantity, values_from = value) %>% 
    rename(abun = Quantity, abun_norm = Quantity.Normalised, unique_pep = Razor.Unique.Peptides, q_value = Qvalue, single_hit = IsSingleHit) %>% 
    mutate(single_hit = as.logical(single_hit)) %>% 
    drop_na() %>% 
    mutate(abun_norm = log10(abun_norm)) %>% 
    left_join(select(meta, sample, raw_sample), by = "raw_sample") %>% 
    select(-raw_sample) %>% 
    mutate(sample = as_factor(sample)) %>% 
    relocate(sample, .after = id)
  
  qc <- dat %>% 
    select(id, sample, unique_pep, q_value, single_hit)
  dat <- dat %>% 
    select(-c(unique_pep, q_value, single_hit))
  
  list(
    info = info,
    qc = qc,
    dat = dat,
    metadata = meta
  )
}


dat2mat <- function(dat, value="abun_norm") {
  dat %>% 
    pivot_wider(id_cols = id, names_from = sample, values_from = !!value) %>% 
    drop_na() %>% 
    column_to_rownames("id") %>% 
    as.matrix()
}

