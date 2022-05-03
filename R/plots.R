okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "grey80", "grey30", "black")

gs <- function(gg, name, width, height) {
  ggsave(filename = file.path("fig", paste0(name, ".png")), plot = gg, device = "png",
         width = width, height = height, dpi = 300)
}

plot_volma <- function(res, group, point_size, point_alpha) {
  r <- res %>%
    mutate(
      group = get(group)
    ) %>% 
    select(x, y, sel, group)
  r_sel <- r %>% filter(sel)
  r_nsel <- r %>% filter(!sel)

  rm(res, r)  # Minimise environment for serialisation
    
  g <- ggplot() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    geom_point(data = r_nsel, aes(x = x, y = y), colour = "grey70",
               size = point_size, alpha = point_alpha) +
    facet_wrap( ~ group) 
  
  if (nrow(r_sel) > 0) {
    g <- g + geom_point(data = r_sel, aes(x = x, y = y), colour = "black",
                        size = point_size, alpha = point_alpha)
  }
  g
}

plot_ma <- function(res, a = "AveExpr", fc = "logFC", p = "PValue", fdr = "FDR", group = "contrast",
                    fdr_limit = 0.05, point_size = 0.5, point_alpha = 0.5) {
  res %>% 
    mutate(
      x = get(a),
      y = get(fc),
      sel = get(fdr) < fdr_limit
    ) %>% 
    plot_volma(group, point_size, point_alpha) +
    geom_hline(yintercept = 0, size = 0.1, alpha = 0.5) +
    labs(x = expression(log[10]~Intensity), y = expression(log[2]~FC))
}

plot_volcano <- function(res, fc = "logFC", p = "PValue", fdr = "FDR", group = "contrast",
                         fdr_limit = 0.05, point_size = 0.5, point_alpha = 0.5) {
  res %>% 
    mutate(
      x = get(fc),
      y = -log10(get(p)),
      sel = get(fdr) < fdr_limit
    ) %>% 
    plot_volma(group, point_size, point_alpha) +
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


plot_up_down <- function(res, fc = "logFC", fdr = "FDR", group = "contrast", fdr_limit = 0.05) {
  r <- res %>%
    mutate(
      direction = if_else(get(fc) > 0, "up", "down"),
      sig = get(fdr) < fdr_limit,
      group = get(group)
    ) %>% 
    filter(sig) %>% 
    group_by(group, direction, .drop = FALSE) %>% 
    tally() %>%
    mutate(x = if_else(direction == "down", -n, n)) %>% 
    mutate(adj = -1.3 * sign(x) / 2 + 0.5)
  rm(res)
  ggplot(r, aes(x = x, y = group, fill = direction)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    geom_vline(xintercept = 0, colour = "grey50") +
    geom_col() +
    scale_fill_manual(values = okabe_ito_palette) +
    geom_text(aes(label = n, hjust = adj)) +
    scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
    labs(x = "Count", y = "Coefficient")
}


make_hist <- function(x, bins) {
  h <- hist(x, breaks = bins, plot = FALSE)
  tibble(
    x = h$mids,
    y = h$density
  )
}


plot_sample_distributions <- function(set, bins = 100, ncol = 15, what = "abu_norm", colour_var = "treatment",
                                      x_breaks = c(-5, 0, 5), text_size = 10, x_lim = c(-5, 5)) {
  d <- set$dat %>% 
    left_join(set$metadata, by = "sample") %>% 
    mutate(val = get(what), colvar = get(colour_var)) %>% 
    select(sample, val, colvar) %>% 
    nest(data = val) %>% 
    mutate(hist = map(data, ~make_hist(.x$val, bins = bins))) %>% 
    select(-data) %>% 
    unnest(hist)
  w <- d$x[2] - d$x[1]
  
  rm(set)
  
  ggplot(d, aes(x = x, y = y, fill = colvar)) +
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
    #geom_vline(xintercept = 0, colour = "brown") +
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

plot_clustering <- function(set, text_size = 10, what = "abu_norm", dist.method = "euclidean",
                            clust.method = "complete", colour_var = "treatment") {
  tab <- dat2mat(set$dat, what)
  
  dendr <- t(tab) %>% 
    dist(method = dist.method) %>% 
    hclust(method = clust.method) %>%
    dendsort::dendsort() %>% 
    ggdendro::dendro_data()
  
  seg <- ggdendro::segment(dendr)
  meta <- set$metadata %>%
    mutate(
      colvar = get(colour_var),
      sample = as.character(sample)
    )
  labs <- left_join(dendr$labels %>% mutate(label = as.character(label)), meta, by = c("label" = "sample")) %>% 
    mutate(colour = okabe_ito_palette[as_factor(colvar)])
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
    labs(x = NULL, y = "Distance")
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
    viridis::scale_fill_viridis(option = "cividis") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = text_size),
      axis.text.y = element_text(size = text_size)
    ) +
    labs(x = NULL, y = NULL, fill = "Correlation")
}


plot_xy <- function(dat, colour_var, shape_var, point_size = 1) {
  dat %>% 
    rename(colvar = !!colour_var, shapevar = !!shape_var) %>% 
  ggplot(aes(x = x, y = y, colour = colvar, shape = shapevar)) +
    theme_bw() +
    geom_point(size = point_size) +
    scale_color_manual(values = okabe_ito_palette, name = colour_var) +
    scale_shape_manual(values = c(15:18, 0, 1, 2, 5, 6), na.value = 4, name = shape_var) +
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

plot_pca <- function(set, point_size = 2, what = "abu_norm", colour_var = "treatment", shape_var = "day") {
  tab <- dat2mat(set$dat, what)
  
  # remove rows with zero variance
  tab <- tab[apply(tab, 1, function(v) sum(is.na(v)) == 0), ]
  tab <- tab[apply(tab, 1, sd) > 0, ]
  
  pca <- prcomp(t(tab), scale. = TRUE, center = TRUE)
  var.perc <- 100 * (pca$sdev)^2 / sum((pca$sdev)^2)
  pca1 <- sprintf("PCA1 (%5.1f%%)", var.perc[1])
  pca2 <- sprintf("PCA2 (%5.1f%%)", var.perc[2])
  pca2xy(pca, set$metadata) %>% 
    plot_xy(colour_var, shape_var, point_size) +
    labs(x = pca1, y = pca2)
}


plot_umap <- function(set, what = "abu_norm", point_size = 2, seed = 1,
                      n_neighbours = 15, min_dist = 0.01,
                      id_sel = NULL, sample_sel = NULL,
                      colour_var = "treatment", shape_var = "day") {
  tab <- dat2mat(set$dat, what)
  meta <- set$metadata
  
  if (!is.null(id_sel))
    tab <- tab[as.character(id_sel), ]
  if (!is.null(sample_sel)) {
    tab <- tab[, sample_sel]
    meta <- meta %>% filter(sample %in% sample_sel)
  }
  
  set.seed(seed)
  tab <- tab[apply(tab, 1, function(v) sum(is.na(v)) == 0), ]
  uwot::umap(t(tab), n_neighbors = n_neighbours, min_dist = min_dist) %>% 
    umap2xy(meta) %>% 
    plot_xy(colour_var, shape_var, point_size)
}

plot_k_umap <- function(ku, point_size = 2) {
  ku %>% 
    plot_xy(point_size = point_size)
}

plot_detection <- function(set) {
  cumcurve <- function(d, what) {
    d %>%
      group_by(get(what)) %>% 
      tally() %>%
      arrange(desc(n)) %>% 
      mutate(x = row_number())
  }
  plotcurve <- function(d) {
    ggplot(d, aes(x = x, y = n)) +
      theme_bw() +
      geom_step(direction = "vh") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA))
  }
  
  dp <- cumcurve(set$dat, "id")
  ds <- cumcurve(set$dat, "sample")
  
  rm(set)
  
  g1 <- plotcurve(dp) + labs(x = "Proteins", y = "Detected in that many samples")
  g2 <- plotcurve(ds) + labs(x = "Samples", y = "Contains that many proteins")

  cowplot::plot_grid(g1, g2, align = "h")
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

plot_unique_pep <- function(set, min_pep) {
  d <- set$qc %>% 
    group_by(unique_pep) %>%
    tally() %>% 
    mutate(uni = unique_pep < min_pep)
  rm(set)
  ggplot(d, aes(x = unique_pep, y = n, fill = uni)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    geom_point(shape = 21) +
    scale_fill_manual(values = c("black", "white")) +
    scale_x_continuous(trans = "log1p", breaks = c(0, 1, 2, 5, 10, 50, 100, 200, 300)) +
    scale_y_continuous(trans = "log10", labels = scales::comma) +
    labs(x = "Number of unique peptides", y = "Count of protein-samples")
}


plot_big_heatmap <- function(set, what = "abu_norm", min_n = 100, max_fc = 2,
                             id_sel = NULL, sample_sel = NULL, order_col = TRUE,
                             text.size = 8) {
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
  
  smeta <- set$metadata %>% 
    filter(sample %in% colnames(tab) & !bad) %>% 
    arrange(day, treatment, batch)  
  smpls <- smeta %>% pull(sample)
  
  tab <- tab[, smpls]
  
  if (is.null(order_col)) {
    divs <- smeta %>%
      group_by(day) %>% 
      tally() %>% 
      mutate(div = cumsum(n) + 0.5)
  } else {
    divs <- NULL
  }
  
  
  ggheatmap(tab, order.col = order_col, with.x.text = TRUE, with.y.text = FALSE,
            dendro.line.size = 0.2, text.size = text.size,
            max.fc = max_fc, legend.name = "logFC", divs = divs)
}


plot_protein <- function(set, pids, what = "abu_norm", colour_var = "batch", shape_var = "sex",
                         ncol = NULL, point_size = 1.5, filt = "completion") {
  info <- set$info %>% 
    filter(id %in% pids) %>% 
    mutate(protein_descriptions = str_remove(protein_descriptions, ";.+$")) %>% 
    mutate(prot = glue::glue("{gene_names}: {protein_descriptions}"))
  d <- set$dat %>%
    mutate(val = get(what)) %>% 
    filter(id %in% pids) %>%
    left_join(set$metadata, by = "sample") %>% 
    filter(!bad & !!rlang::parse_expr(filt)) %>% 
    rename(colvar = !!colour_var, shapevar = !!shape_var) %>% 
    arrange(treatment, day) %>% 
    unite(x, c(treatment, day)) %>% 
    mutate(x = as_factor(x), xi = as.integer(x)) %>% 
    left_join(info, by = "id") %>% 
    select(id, prot, x, xi, val, colvar, shapevar)
  dm <- d %>% 
    group_by(id, prot, xi) %>% 
    summarise(M = mean(val))
  rm(set)
  ggplot(d, aes(x = x, y = val)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    scale_colour_manual(values = okabe_ito_palette, name = colour_var) +
    scale_shape_discrete(name = shape_var) +
    ggbeeswarm::geom_quasirandom(aes(colour = colvar, shape = shapevar), width = 0.2, size = point_size, alpha = 0.8) +
    geom_segment(data = dm, aes(x = xi - 0.3, y = M, xend = xi + 0.3, yend = M), size = 1, colour = "brown") +
    facet_wrap(~ prot, labeller = label_wrap_gen(), ncol = ncol) +
    labs(x = NULL, y = what)
}


plot_lograt_protein <- function(df, pids, ncol = NULL, colour_var = "age_group", shape_var = "sex", filt = "completion") {
  info <- df$info %>% 
    filter(id %in% pids) %>% 
    mutate(protein_descriptions = str_remove(protein_descriptions, ";.+$")) %>% 
    mutate(prot = glue::glue("{gene_names}: {protein_descriptions}"))
  d <- df$dat %>% 
    filter(id %in% pids) %>% 
    left_join(df$metadata, by = "participant_id") %>% 
    filter(!bad & !!rlang::parse_expr(filt)) %>% 
    rename(colvar = !!colour_var, shapevar = !!shape_var) %>% 
    left_join(info, by = "id") %>% 
    mutate(x = as_factor(treatment), xi = as.integer(x))
  dm <- d %>% 
    group_by(id, prot, xi) %>% 
    summarise(M = mean(logFC))
  rm(df)
  ggplot(d, aes(x = x, y = logFC)) +
    theme_bw() +
    theme(
      panel.grid = element_blank()
    ) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_beeswarm(aes(colour = colvar, shape = shapevar)) +
    scale_colour_manual(values = okabe_ito_palette) +
    geom_segment(data = dm, aes(x = xi - 0.3, y = M, xend = xi + 0.3, yend = M), size = 1, colour = "brown") +
    facet_wrap(~ prot, labeller = label_wrap_gen(), ncol = ncol) +
    labs(x = NULL, y = expression(log[2]~I[29]/I[1]))
}


plot_protein_coverage <- function(qc) {
  cov_mean <- qc %>%
    filter(coverage > 0) %>% 
    group_by(id) %>%
    summarise(M = mean(coverage)) %>%
    arrange(M) %>%
    mutate(x = row_number())
  d <- qc %>% 
    filter(coverage > 0) %>% 
    left_join(cov_mean, by = "id") %>% 
    select(x, coverage)
  rm(qc, cov_mean)
  
  d %>% 
    ggplot(aes(x = x, y = coverage)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    geom_scattermore(pointsize = 0.2, pixels = c(1000, 1000)) +
    labs(x = "Proteins", y = "Coverage")
}

plot_participant_1_29 <- function(meta) {
  pp <- get_full_participants(meta)
  meta %>%
    left_join(pp, by = "participant_id") %>% 
    mutate(`both 1 and 29` = !is.na(first)) %>% 
  ggplot(aes(x = participant_id, y = day, group = participant_id, colour = `both 1 and 29`, shape = completion)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      #legend.position = "none"
    ) +
    geom_line() +
    geom_point() +
    scale_colour_manual(values = c("grey", "black")) +
    scale_shape_manual(values = c(4, 16)) +
    facet_wrap(~ treatment, ncol = 1, scales = "free_x") +
    labs(x = "Participant ID", y = "Day")
}

plot_meta_numbers <- function(meta, fill_var = "batch") {
  meta %>% 
    mutate(fillvar = get(fill_var)) %>% 
    group_by(treatment, day, fillvar) %>% 
    tally() %>% 
  ggplot(aes(x = day, y = n, fill = fillvar)) +
    theme_bw() +
    theme(
      panel.grid = element_blank()
    ) +
    geom_col(position = "stack", colour = "grey50") +
    scale_fill_manual(values = okabe_ito_palette) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    facet_wrap(~treatment, nrow = 1) +
    labs(x = "Day", y = "Count", fill = fill_var)
}

plot_paricipant_stats <- function(part) {
  th <- theme_bw() +
    theme(panel.grid = element_blank())
  p_d <- function(x, name) {
    tibble(x = x) %>% 
      group_by(x) %>% 
      tally() %>% 
      ggplot(aes(x = x, y = n)) +
      th +
      geom_col() +
      scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
      labs(x = name, y = NULL)
  }
  p_c <- function(x, name) {
    tibble(x = x) %>% 
      group_by(x) %>% 
      tally() %>% 
      ggplot(aes(x = x, xend = x, y = n, yend = 0)) +
      th +
      geom_segment(colour = "grey60") +
      geom_point() +
      scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
      labs(x = name, y = NULL)
  }
  
  plot_grid(
    p_c(part$age, "Age"),
    p_c(part$time_from_symptoms, "Time from symtoms (days)"),
    p_c(part$days_of_treatment, "Days of treatment"),
    p_d(part$sex, "Sex"),
    p_d(part$treatment, "Treatment"),
    p_d(part$completion, "Completion")
  )
}

plot_fc_fc <- function(df1, df2, info, labx = "x", laby = "y", logfc_limit = 0.01, dif_limit = 1) {
  d <- df1 %>%
    full_join(df2, by = "id") %>%
    left_join(info, by = "id") %>%
    filter(FDR.x < logfc_limit | FDR.y < logfc_limit) %>% 
    mutate(dif = abs(logFC.x - logFC.y))
  
  ggplot(d, aes(x = logFC.x, y = logFC.y)) +
    theme_bw() +
    geom_point(colour = "grey50") +
    geom_abline(slope = 1, intercept = 0, colour = "red") +
    geom_text_repel(data = d %>% filter(dif > dif_limit), aes(label = gene_names)) +
    labs(x = labx, y = laby)
}


make_upset_batch_run <- function(meta) {
  batch_list <- meta %>% 
    group_split(batch) %>% 
    map(\(x) pull(x, sample)) %>% 
    set_names(str_c("batch_", levels(meta$batch)))
  run_list <- meta %>% 
    group_split(run_index) %>% 
    map(\(x) pull(x, sample)) %>% 
    set_names(str_c("run_", levels(meta$run_index)))
  c(batch_list, run_list)
}

