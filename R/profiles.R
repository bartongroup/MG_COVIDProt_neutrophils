make_time_profiles <- function(set, what = "abu_norm") {
  set$dat %>%
    mutate(val = get(what) / log10(2)) %>% 
    left_join(set$metadata, by = "sample") %>%
    filter(treatment == "drug" & !is.na(day)) %>%
    droplevels() %>% 
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
    geom_point(size = 0.3, alpha = 0.5) +
    geom_hline(yintercept = 0, colour = "brown") +
    facet_wrap(~ clust) +
    labs(x = "Day", y = expression(log[2]~FC))
}


plot_prof_heatmap <- function(tpr, max_fc = 2) {
  tpr %>% 
    drop_na() %>% 
    column_to_rownames("id") %>% 
    ggheatmap(order.col = FALSE, dendro.line.size = 0.1, max.fc = max_fc)
}



# Gradient-gradient plot with pricipal componend axes
plot_gradients <- function(tpr) {
  days <- colnames(tpr)[-1]
  n <- length(days)
  tb <- tibble(d1 = days[1:(n - 1)], d2 = days[2:n]) %>% 
    rowwise() %>% 
    mutate(pair = list(c(d1, d2)))
  xy <- map2_dfr(tb$pair[1:(n - 2)], tb$pair[2:(n - 1)], function(p1, p2) {
    tpr %>% 
      transmute(
        r1 = get(p1[2]) - get(p1[1]),
        r2 = get(p2[2]) - get(p2[1])
      ) %>% 
      add_column(rat = str_c(p1, p2, sep = "-", collapse = ":")) %>% 
      drop_na()
  })

  # Calculate principal components' slopes and intercepts
  eig <- map_dfr(unique(xy$rat), function(r) {
    d <- xy %>% filter(rat == r) %>% select(r1, r2)
    eigen_vec <- prcomp(d)$rotation
    slopes <- eigen_vec[2, ] / eigen_vec[1, ]
    tibble(
      intercpt = mean(d$r2),
      slp1 = slopes[1],
      slp2 = slopes[2]
    ) %>% 
      add_column(rat = r)
  })
  
  ggplot() + 
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point(data = xy, aes(x = r1, y = r2), size = 0.3) +
    stat_ellipse(data = xy, aes(x = r1, y = r2), colour = "red", size = 1, level = 0.99) +
    geom_abline(data = eig, aes(slope = slp1, intercept = intercpt), colour = "red") +
    geom_abline(data = eig, aes(slope = slp2, intercept = intercpt), colour = "red", alpha = 0.4) +
    facet_wrap(~ rat, nrow = 1) +
    coord_fixed()
}