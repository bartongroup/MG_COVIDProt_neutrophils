okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

gs <- function(gg, name, width, height) {
  ggsave(filename = file.path("fig", paste0(name, ".png")), plot = gg, device = "png",
         width = width, height = height, dpi = 300)
}

plot_volma <- function(res, p, fdr, group, fdr_limit, point_size, point_alpha) {
  r <- res %>%
    mutate(
      p = get(p),
      sig = get(fdr) < fdr_limit,
      group = get(group)
    ) %>% 
    select(x, y, sig, group)
  rm(res)  # Minimise environment for serialisation
  ggplot(r, aes(x = x, y = y, colour = sig)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_colour_manual(values = c("grey70", "black")) +
    facet_grid(. ~ group) 
}

plot_ma <- function(res, a = "AveExpr", fc = "logFC", p = "PValue", fdr = "FDR", group = "contrast",
                    fdr_limit = 0.05, point_size = 0.5, point_alpha = 0.5) {
  res %>% 
    mutate(
      x = get(a),
      y = get(fc),
    ) %>% 
    plot_volma(p, fdr, group, fdr_limit, point_size, point_alpha) +
    geom_hline(yintercept = 0, size = 0.1, alpha = 0.5) +
    labs(x = expression(log[10]~Intensity), y = expression(log[2]~FC))
}

plot_volcano <- function(res, fc = "logFC", p = "PValue", fdr = "FDR", group = "contrast",
                         fdr_limit = 0.05, point_size = 0.5, point_alpha = 0.5) {
  res %>% 
    mutate(
      x = get(fc),
      y = -log10(get(p)),
    ) %>% 
    plot_volma(p, fdr, group, fdr_limit, point_size, point_alpha) +
    geom_vline(xintercept = 0, size = 0.1, alpha = 0.5) +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)))
}

plot_pdist <- function(res, p = "PValue", group = "contrast", n_bins = 50) {
  brks <- seq(0, 1, length.out = n_bins)
  r <- res %>% 
    mutate(p = get(p), grp = get(group)) %>% 
    select(p, grp)
  rm(res)
  ggplot(r, aes(x = p, y = after_stat(density))) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_histogram(breaks = brks) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    facet_wrap(~grp)
}


make_hist <- function(x, bins) {
  h <- hist(x, breaks = bins, plot = FALSE)
  tibble(
    x = h$mids,
    y = h$density
  )
}


plot_sample_distributions <- function(set, bins = 100, ncol = 15, what = "abu_norm",
                                      x_breaks = c(-5, 0, 5), text_size = 10, x_lim = c(-5, 5)) {
  d <- set$dat %>% 
    mutate(val = get(what)) %>% 
    left_join(set$metadata, by = "sample") %>% 
    select(sample, val, group) %>% 
    nest(data = val) %>% 
    mutate(hist = map(data, ~make_hist(.x$val, bins = bins))) %>% 
    select(-data) %>% 
    unnest(hist)
  w <- d$x[2] - d$x[1]
  
  rm(set)
  
  ggplot(d, aes(x = x, y = y, fill = group)) +
    theme_bw() +
    theme(
      text = element_text(size = text_size),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.background = element_blank()
    ) +
    geom_col(width = w) +
    geom_hline(yintercept = 0, colour = "brown") +
    geom_vline(xintercept = 0, colour = "brown") +
    facet_wrap(~sample, scales = "free_y", ncol = ncol) +
    scale_fill_manual(values = okabe_ito_palette) +
    scale_x_continuous(breaks = x_breaks, limits = x_lim) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x = what, y = "Density")
}

plot_kernels <- function(set, what = "abu_norm", xlab = "Normalised abundance") {
  d <- set$dat %>% 
    mutate(val = get(what)) %>% 
    select(sample, val) %>% 
    nest(data = val) %>% 
    mutate(
      krn = map(data, ~density(.x$val, bw = "SJ") %>% tidy())
    ) %>% 
    select(-c(data)) %>% 
    unnest(krn)
  rm(set)
  ggplot(d, aes(x = x, y = y, group = sample)) +
    theme_bw() +
    theme(
      panel.grid = element_blank()
    ) +
    geom_line(alpha = 0.1) +
    labs(x = xlab, y = "Density") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)))
}

plot_kernel_comparison <- function(set) {
  g1 <- plot_kernels(set, "abu", "Abundance")
  g2 <- plot_kernels(set, "abu_norm", "Normalised abundance")
  cowplot::plot_grid(g1, g2, nrow = 1)
}

plot_clustering <- function(set, text_size = 10, what = "abu_norm", dist.method = "euclidean", clust.method = "complete") {
  tab <- dat2mat(set$dat, what)
  
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
  rm(set, tab, dendr)
  ggplot() +
    theme.d +
    coord_flip() +
    geom_segment(data = seg, aes_(x = ~x, y = ~y, xend = ~xend, yend = ~yend)) +
    scale_x_continuous(breaks = seq_along(labs$label), labels = labs$label) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(seg$y) * 1.03)) +
    scale_colour_manual(values = okabe_ito_palette) +
    labs(x = NULL, y = "Distance", title = what)
}


plot_distance_matrix <- function(set, what = "abu_norm", text_size = 10) {
  tab <- dat2mat(set$dat, what)
  
  d <- cor(tab, use = "complete.obs") %>% 
    as_tibble(rownames = "sample") %>%
    pivot_longer(-sample) %>% 
    mutate(sample = factor(sample, levels = set$metadata$sample)) %>% 
    mutate(name = factor(name, levels = set$metadata$sample))
  
  rm(set, tab)
  ggplot(d, aes(x = sample, y = name)) +
    geom_tile(aes(fill = value)) +
    scale_fill_viridis(option = "cividis") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = text_size),
      axis.text.y = element_text(size = text_size)
    ) +
    labs(x = NULL, y = NULL, fill = "Correlation")
}


plot_xy <- function(dat, point_size = 1) {
  ggplot(dat, aes(x = x, y = y, colour = group, shape = day)) +
    theme_bw() +
    geom_point(size = point_size) +
    scale_color_manual(values = okabe_ito_palette) +
    scale_shape_manual(values = c(15:18, 0, 1, 2, 5, 6), na.value = 4) +
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

plot_pca <- function(set, point_size = 2, what = "abu_norm") {
  tab <- dat2mat(set$dat, what)
  
  # remove rows with zero variance
  tab <- tab[apply(tab, 1, function(v) sum(is.na(v)) == 0), ]
  tab <- tab[apply(tab, 1, sd) > 0, ]
  
  pca <- prcomp(t(tab), scale. = TRUE, center = TRUE)
  var.perc <- 100 * (pca$sdev)^2 / sum((pca$sdev)^2)
  pca1 <- sprintf("PCA1 (%5.1f%%)", var.perc[1])
  pca2 <- sprintf("PCA2 (%5.1f%%)", var.perc[2])
  pca2xy(pca, set$metadata) %>% 
    plot_xy(point_size = point_size) + labs(x = pca1, y = pca2, title = what)
}


plot_umap <- function(set, what = "abu_norm", point_size = 2, seed = 1,
                      n_neighbours = 15, min_dist = 0.01) {
  tab <- dat2mat(set$dat, what)
  
  set.seed(seed)
  tab <- tab[apply(tab, 1, function(v) sum(is.na(v)) == 0), ]
  uwot::umap(t(tab), n_neighbors = n_neighbours, min_dist = min_dist) %>% 
    umap2xy(set$metadata) %>% 
    plot_xy(point_size = point_size)
}

plot_k_umap <- function(ku, point_size = 2) {
  ku %>% 
    plot_xy(point_size = point_size)
}

plot_detection <- function(set) {
  d <- set$dat %>% 
    group_by(id) %>% 
    summarise(n_detect = n()) %>% 
    group_by(n_detect) %>% 
    tally() %>% 
    mutate(cum = cumsum(n))
  rm(set)
  ggplot(d, aes(x = max(cum) - cum + 1, y = n_detect)) +
    theme_bw() +
    theme(
    ) +
    geom_step(direction = "vh") +
    labs(x = "Proteins", y = "Detected in that many samples") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA))
}

plot_sample_detection <- function(set) {
  d <- set$dat %>% 
    group_by(sample) %>% 
    tally()
  rm(set)
  ggplot(d, aes(x = n, y = fct_reorder(sample, n))) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_point() +
    labs(x = "Number of detected proteins", y = "Sample")
}

plot_unique_pep <- function(set) {
  d <- set$qc %>% 
    group_by(unique_pep) %>%
    tally()
  rm(set)
  ggplot(d, aes(x = unique_pep, y = n, fill = unique_pep < MIN_PEPTIDES)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    geom_point(shape = 21) +
    scale_fill_manual(values = c("black", "white")) +
    scale_x_continuous(trans = "log1p", breaks = c(0, 1, 2, 5, 10, 50, 100, 200, 300)) +
    scale_y_continuous(trans = "log10", labels = scales::comma) +
    labs(x = "Number of unique peptides", y = "Count")
}


plot_big_heatmap <- function(set, what = "abu_norm", min_n = 100, max_fc = 2,
                             id_sel = NULL, sample_sel = NULL) {
  d <- set$dat %>% 
    mutate(val = get(what))
  if (!is.null(sample_sel))
    d <- d %>% filter(sample %in% sample_sel)
  if (!is.null(id_sel))
    d <- d %>% filter(id %in% id_sel)
  
  d <- d %>% 
    group_by(id) %>% 
    mutate(
      n = n(),
      fc = val - mean(val)
    ) %>% 
    filter(n > min_n)
  tab <- dat2mat(d, "fc")
  ggheatmap(tab, with.x.text = FALSE, with.y.text = FALSE, max.fc = max_fc, legend.name = "logFC")
}


plot_protein <- function(set, pid, what = "abu_norm") {
  d <- set$dat %>%
    mutate(val = get(what)) %>% 
    filter(id == pid) %>%
    left_join(set$metadata, by = "sample") %>% 
    arrange(group, day) %>% 
    unite(x, c(group, day)) %>% 
    mutate(x = as_factor(x), xi = as.integer(x))
  dm <- d %>% 
    group_by(xi) %>% 
    summarise(M = mean(val))
  rm(set)
  ggplot(d, aes(x = x, y = val)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    geom_quasirandom(width = 0.2, size = 1, alpha = 0.8) +
    geom_segment(data = dm, aes(x = xi - 0.3, y = M, xend = xi + 0.3, yend = M), size = 1, colour = "brown") +
    labs(x = NULL, y = what)
}