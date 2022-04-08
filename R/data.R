read_spectronaut_data <- function(file) {
  d <- read_csv(file, show_col_types = FALSE) %>% 
    rename(id = 1)
  info <- d %>% 
    select(id, starts_with("PG")) %>% 
    rename_with(~str_remove(.x, "PG."))
  # columns are complicated, need to simplfy naming
  cols <- tibble(name = colnames(d)) %>% 
    filter(name != "id" & !(str_detect(name, "^PG")) & !str_detect(name, "Coverage")) %>%
    mutate(new = name) %>% 
    mutate(new = str_replace(new, "Razor", "PG.Razor")) %>% 
    mutate(new = str_replace(new, ".PG.", ":"))
  
  dat <- d %>% 
    select(c(id, all_of(cols$name))) %>% 
    set_names(c("id", cols$new)) %>% 
    pivot_longer(-id, names_to = c("sample", "quantity"), names_pattern = "(.*):(.*)") %>% 
    pivot_wider(id_cols = c(id, sample), names_from = quantity, values_from = value) %>% 
    rename(abun = Quantity, abun_norm = Quantity.Normalised, unique_pep = Razor.Unique.Peptides, q_value = Qvalue, single_hit = IsSingleHit) %>% 
    mutate(sample = as_factor(sample), single_hit = as.logical(single_hit)) %>% 
    drop_na() %>% 
    mutate(abun_norm = log10(abun_norm))
  
  list(
    info = info,
    dat = dat
  )
}
