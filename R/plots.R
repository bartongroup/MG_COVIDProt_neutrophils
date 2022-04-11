okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

gs <- function(gg, name, width, height) {
  ggsave(filename = file.path("fig", paste0(name, ".png")), plot = gg, device = "png",
         width = width, height = height, dpi = 300)
}

plot_sample_distributions <- function(set, bins = 100, ncol = 15, value = "abun_norm", x_breaks = c(-5, 0, 5)) {
  d <- set$dat %>% 
    mutate(val = get(value)) %>% 
    #left_join(set$metadata, by = "sample") %>% 
    select(sample, val, group = 1) 
    #mutate(sample = factor(sample, levels = set$metadata$sample))
  ggplot(d, aes(x = val, y = after_stat(density), fill = group)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.background = element_blank()
    ) +
    geom_histogram(bins = bins) +
    geom_hline(yintercept = 0) +
    facet_wrap(~sample, scales = "free_y", ncol = ncol) +
    scale_fill_manual(values = okabe_ito_palette) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x = value, y = "density")
}

plot_clustering <- function(set, text_size = 10, value = "abun_norm", dist.method = "euclidean", clust.method = "complete") {
  tab <- dat2mat(set$dat, value)
  
  dendr <- t(tab) %>% 
    dist(method = dist.method) %>% 
    hclust(method = clust.method) %>%
    dendsort::dendsort() %>% 
    ggdendro::dendro_data()
  
  seg <- ggdendro::segment(dendr)
  meta <- set$metadata %>% mutate(sample = as.character(sample))
  labs <- left_join(dendr$labels %>% mutate(label = as.character(label)), meta, by = c("label" = "sample")) %>% 
    mutate(colour = okabe_ito_palette[as_factor(group)])
  theme.d <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = text_size, colour = labs$colour),
    axis.line.y = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_line(size = 0.5),
    axis.ticks.y = ggplot2::element_blank()
  )
  g <- ggplot() +
    theme.d +
    coord_flip() +
    geom_segment(data = seg, aes_(x = ~x, y = ~y, xend = ~xend, yend = ~yend)) +
    scale_x_continuous(breaks = seq_along(labs$label), labels = labs$label) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(seg$y) * 1.03)) +
    scale_colour_manual(values = okabe_ito_palette) +
    labs(x = NULL, y = "Distance", title = value)
  # This dramatically reduces file size in serialization used by targets
  cowplot::plot_grid(g)
}


plot_distance_matrix <- function(set, value = "abun_norm", text_size = 10) {
  tab <- dat2mat(set$dat, value)
  
  d <- cor(tab, use = "complete.obs") %>% 
    as_tibble(rownames = "sample") %>%
    pivot_longer(-sample) %>% 
    mutate(sample = factor(sample, levels = set$metadata$sample)) %>% 
    mutate(name = factor(name, levels = set$metadata$sample))
  
  g <- ggplot(d, aes(x = sample, y = name)) +
    geom_tile(aes(fill = value)) +
    scale_fill_viridis(option = "cividis") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = text_size),
      axis.text.y = element_text(size = text_size)
    ) +
    labs(x = NULL, y = NULL, fill = "Correlation")
  # This dramatically reduces file size in serialization used by targets
  cowplot::plot_grid(g)
}


plotXY <- function(dat, point_size = 1) {
  ggplot(dat, aes(x = x, y = y, colour = group, shape = day)) +
    theme_bw() +
    geom_point(size = point_size) +
    scale_color_manual(values = okabe_ito_palette) +
    scale_shape_manual(values = c(15:18, 0:9)) +
    labs(x = NULL, y = NULL) +
    theme(
      panel.grid = element_blank()
    ) 
    #geom_text_repel(aes(label = sample), size = 3, colour = "black")
}

pca2xy <- function(pc, meta) {
  pc$x %>% 
    as_tibble(rownames = "sample") %>% 
    select(x = PC1, y = PC2, sample) %>% 
    left_join(meta, by = "sample")
}

umap2xy <- function(um, meta) {
  colnames(um) <- c("x", "y")
  um %>% 
    as_tibble() %>% 
    mutate(sample = meta$sample) %>% 
    left_join(meta, by = "sample")
}

plot_pca <- function(set, point_size = 2, value = "abun_norm") {
  tab <- dat2mat(set$dat, value)
  
  # remove rows with zero variance
  tab <- tab[apply(tab, 1, function(v) sum(is.na(v)) == 0), ]
  tab <- tab[apply(tab, 1, sd) > 0, ]
  
  pca <- prcomp(t(tab), scale. = TRUE, center = TRUE)
  var.perc <- 100 * (pca$sdev)^2 / sum((pca$sdev)^2)
  pca1 <- sprintf("PCA1 (%5.1f%%)", var.perc[1])
  pca2 <- sprintf("PCA2 (%5.1f%%)", var.perc[2])
  pca2xy(pca, set$metadata) %>% 
    plotXY(point_size = point_size) + labs(x = pca1, y = pca2, title = value)
}


plot_umap <- function(set, value = "abun_norm", point_size = 2, seed = 1, n_neighbours = 15, min_dist = 0.01) {
  tab <- dat2mat(set$dat, value)
  
  set.seed(seed)
  tab <- tab[apply(tab, 1, function(v) sum(is.na(v)) == 0), ]
  uwot::umap(t(tab), n_neighbors = n_neighbours, min_dist = min_dist) %>% 
    umap2xy(set$metadata) %>% 
    plotXY(point_size = point_size)
}


