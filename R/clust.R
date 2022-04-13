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

make_time_profiles <- function(set, what = "abu_norm") {
  set$dat %>%
    mutate(val = get(what) / log10(2)) %>% 
    left_join(set$metadata, by = "sample") %>%
    filter(group == "drug" & !is.na(day)) %>%
    group_by(id, day) %>%
    summarise(m = mean(val)) %>%
    ungroup() %>%
    group_by(id) %>%
    arrange(day) %>%
    mutate(dif = m - m[1]) %>%
    ungroup() %>% 
    pivot_wider(id_cols = id, names_from = day, values_from = dif)
}

cluster_time_profiles <- function(tpr, n_clust = 5) {
  tpr <- drop_na(tpr)
  X <- tpr %>% 
    select(-id) %>% 
    as.matrix()
  km <- kmeans(X, centers = n_clust)
  tpr %>% 
    mutate(clust = km$cluster %>% as_factor()) %>% 
    pivot_longer(-c(id, clust), names_to = "day", values_to = "log_fc") %>% 
    mutate(day = as.numeric(day))
}


plot_cluster_profiles <- function(tpc) {
  tpc %>% 
    ggplot(aes(x = day, y = log_fc, group = id)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_line(alpha = 0.1) +
    geom_hline(yintercept = 0, colour = "brown") +
    facet_wrap(~ clust) +
    labs(x = "Day", y = expression(log[2]~FC))
}