download_uniprot_mapping <- function(uri) {
  read_tsv(UNIPROT_MAPPING_FILE, col_names = c("uniprot", "what", "gid"), show_col_types = FALSE) %>% 
    filter(what == "Gene_Name") %>% 
    select(uniprot, gene_name = gid)
}


read_spectronaut_long_data <- function(data_file, info_file, uni_gene, bad_samples = BAD_SAMPLES) {
  meta <- readxl::read_excel(info_file, na = "n/a") %>% 
    set_names(c(
      "batch",
      "run_index",
      "raw_file_name",
      "participant_id",
      "day",
      "study",
      "age",
      "sex",
      "time_from_symptoms",
      "treatment",
      "completion",
      "day_of_treatment",
      "on_drug"
    )) %>% 
    mutate(raw_sample = str_remove(raw_file_name, ".raw")) %>% 
    select(-raw_file_name) %>% 
    mutate(
      treatment = recode(treatment, "Brensocatib 25mg once daily" = "drug", "Placebo 25mg once daily" = "placebo"),
      completion = as.logical(completion),
      sex = toupper(sex),
      on_drug = str_remove(on_drug, "\\s\\(discontinued\\)"),
      delay = if_else(time_from_symptoms > 10, "large", "small"),
      age_group = case_when(
        age < 50 ~ "<50",
        age >= 50 & age < 65 ~ "50-65",
        age >= 65 & age < 80 ~ "65-80",
        age >= 80 ~ ">80"
      )
    ) %>%
    arrange(age) %>% 
    mutate(age_group = as_factor(age_group)) %>% 
    arrange(treatment, day, batch) %>% 
    unite(treat_day, c(treatment, day), sep = "_", remove = FALSE) %>% 
    group_by(treatment, day) %>% 
    mutate(trep = row_number()) %>% 
    ungroup()
  
  raw <- read_tsv(data_file, show_col_types = FALSE) %>% 
    rename_with(.fn = ~str_remove(., "^R.|^PG.")) %>% 
    janitor::clean_names() %>% 
    rename(
      raw_sample = file_name
    ) %>% 
    filter(raw_sample %in% meta$raw_sample)
  
  samrep <- raw %>% 
    select(raw_sample, replicate) %>% 
    distinct()
  
  metadata <- meta %>% 
    left_join(samrep, by = "raw_sample") %>%
    mutate(
      s1 = recode(treatment, "drug" = "D", "placebo" = "P"),
      s2 = as.character(batch),
      s3 = sprintf("%02d", day),
      s4 = sprintf("%02d", trep)
    ) %>% 
    mutate(sample = glue::glue("{s1}{s2}_{s3}-{s4}") %>% as.character()) %>% 
    select(-c(s1, s2, s3, s4, trep)) %>% 
    mutate(across(c(batch, participant_id, day, sex, treatment), as.factor)) %>% 
    mutate(
      treatment = fct_relevel(treatment, "placebo"),
      sex = fct_relevel(sex, "F"),
      bad = sample %in% bad_samples
    )
  
  info <- raw %>% 
    select(protein_accessions, protein_descriptions, protein_names) %>% 
    distinct() %>% 
    mutate(id = row_number(), .before = 1)
  
  id_prot_gene <- info %>%
    select(id, uniprot = protein_accessions) %>%
    separate_rows(uniprot, sep = ";") %>%
    distinct() %>%
    mutate(uniprot = str_remove(uniprot, "\\-\\d+$")) %>% 
    left_join(uni_gene, by = c("uniprot"))
  gene_names <- id_prot_gene %>%
    group_by(id) %>%
    summarise(gene_names = str_c(unique(gene_name), collapse = ";"))
  info <- info %>% 
    left_join(gene_names, by = "id")
  
  d <- raw %>% 
    left_join(select(info, id, protein_accessions), by = "protein_accessions") %>% 
    left_join(select(metadata, raw_sample, sample), by = "raw_sample")
  
  dat <- d %>% 
    select(id, sample, quantity) %>% 
    drop_na() %>% 
    normalise_proteins()
  
  qc <- d %>% 
    select(id, sample, coverage, is_single_hit, qvalue)
  
  list(
    metadata = metadata,
    info = info,
    id_prot_gene = id_prot_gene,
    qc = qc,
    dat = dat
  )
}

read_spectronaut_wide_data <- function(file, meta, min_pep) {
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


normalise_proteins <- function(dat) {
  md <- dat %>% 
    mutate(abu = log10(quantity)) %>% 
    group_by(sample) %>% 
    summarise(med = median(abu)) %>% 
    mutate(med = med / mean(med))
  dat %>% 
    left_join(md, by = "sample") %>% 
    mutate(
      abu = log10(quantity),
      abu_norm = abu / med
    ) %>% 
    select(-c(quantity, med))
}


select_detected_proteins <- function(set) {
  set$dat %>% 
    select(id) %>% 
    distinct() %>% 
    left_join(set$id_prot_gene, by = "id")
}


dat2mat <- function(dat, what = "abu_norm", names = "sample") {
  dat %>% 
    pivot_wider(id_cols = id, names_from = !!names, values_from = !!what) %>% 
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

get_full_participants <- function(meta, days = c(1, 29)) {
  meta %>% 
    filter(day %in% days & !bad) %>% 
    pivot_wider(id_cols = c(participant_id), names_from = day, values_from = c(sample, batch)) %>% 
    drop_na() %>% 
    set_names(c("participant_id", "first", "last", "first_batch", "last_batch"))
}

logfc_days <- function(set, days = c(1, 29)) {
  pd <- get_full_participants(set$metadata, days)
  X <- dat2mat(set$dat, what = "abu_norm")
  X <- X / log10(2)   # we are going to find log2_FC
  
  dat <- map_dfr(1:nrow(pd), function(i) {
    r <- pd[i, ]
    tibble(
      id = rownames(X) %>% as.integer(),
      participant_id = r$participant_id,
      logFC = X[, r$last] - X[, r$first]
    )
  }) %>% 
    drop_na() %>% 
    droplevels()
  
  meta <- set$metadata %>% 
    select(participant_id, age, sex, treatment, completion, age_group) %>% 
    distinct() %>% 
    filter(participant_id %in% unique(pd$participant_id)) %>% 
    left_join(select(pd, participant_id, first_batch, last_batch), by = "participant_id") %>% 
    unite(batches, c(first_batch, last_batch), sep = "-") %>% 
    add_column(bad = FALSE)
  
  list(
    dat = dat,
    metadata = meta,
    info = set$info,
    id_prot_gene = set$id_prot_gene
  )
}