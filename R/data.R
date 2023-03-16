download_uniprot_mapping <- function(uri) {
  read_tsv(UNIPROT_MAPPING_FILE, col_names = c("uniprot", "what", "gid"), show_col_types = FALSE) |> 
    filter(what == "Gene_Name") |> 
    select(uniprot, gene_symbol = gid)
}

#' Read COVID-19 metadata from an Excel file
#'
#' This function reads a metadata Excel file containing information on
#' participants, samples, and treatments for a COVID-19 study, and performs data
#' cleaning and variable recoding.
#'
#' @param metadata_file A character string specifying the path to the Excel file
#'   containing the metadata
#' @param bad_samples A character vector specifying the samples to be marked as
#'   'bad'
#'
#' @return A cleaned data frame containing metadata for the samples, with
#'   recoded variables and a new 'sample' column created from relevant metadata
#'   columns
read_covid_metadata <- function(metadata_file, bad_samples) {
  readxl::read_excel(metadata_file, na = "n/a") |> 
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
      "days_of_treatment",
      "on_drug"
    )) |> 
    mutate(data_col = str_remove(raw_file_name, ".raw")) |> 
    select(-raw_file_name) |> 
    mutate(
      treatment = recode(treatment, "Brensocatib 25mg once daily" = "Brensocatib", "Placebo 25mg once daily" = "Placebo"),
      completion = as.logical(completion),
      sex = toupper(sex),
      on_drug = str_remove(on_drug, "\\s\\(discontinued\\)"),
      delay = if_else(time_from_symptoms > 20, "large", "small"),
      age_group = case_when(
        age < 50 ~ "<50",
        age >= 50 & age < 65 ~ "50-65",
        age >= 65 & age < 80 ~ "65-80",
        age >= 80 ~ ">80"
      )
    ) |>
    arrange(age) |> 
    mutate(age_group = as_factor(age_group)) |> 
    arrange(treatment, day, batch) |> 
    unite(treat_day, c(treatment, day), sep = "_", remove = FALSE) |> 
    group_by(treatment, day) |> 
    mutate(trep = row_number()) |> 
    ungroup() |> 
    mutate(
      s1 = recode(treatment, "Brensocatib" = "B", "Placebo" = "P"),
      s2 = as.character(batch),
      s3 = sprintf("%02d", day),
      s4 = sprintf("%02d", trep)
    ) |> 
    mutate(sample = glue::glue("{s1}{s2}_{s3}-{s4}") |> as.character()) |> 
    select(-c(s1, s2, s3, s4, trep)) |> 
    mutate(across(c(batch, run_index, participant_id, day, sex, treatment), as_factor)) |> 
    mutate(
      treatment = fct_relevel(treatment, "Placebo"),
      sex = fct_relevel(sex, "M"),
      run_index = fct_relevel(run_index, "278"),
      bad = sample %in% bad_samples
    )
}

read_raw_file <- function(data_file) {
  read_tsv(data_file, show_col_types = FALSE) |> 
    rename_with(.fn = ~str_remove(., "^R.|^PG.")) |> 
    janitor::clean_names() |> 
    rename(
      data_col = file_name
    ) 
}

#' Read Spectronaut long format data
#'
#' This function reads Spectronaut long format data, filters it based on
#' metadata, removes contaminants, and normalizes protein quantities.
#'
#' @param data_file The raw data file in Spectronaut long format
#' @param meta A data frame containing metadata for the samples
#' @param uni_gene A data frame containing UniProt ID to gene symbol mappings
#' @param contaminants A character vector containing contaminant protein
#'   accessions
#'
#' @return A list containing the following elements:
#' * metadata: The input metadata
#' * info: A data frame with protein information and gene names
#' * id_prot_gene: A data frame with protein ID, UniProt accession, and gene symbol information
#' * qc: A data frame containing quality control information
#' * dat: A data frame with normalized protein quantities
read_spectronaut_long_data <- function(data_file, meta, uni_gene, contaminants) {
  raw <- read_raw_file(data_file) |>
    filter(data_col %in% meta$data_col)
  
  info <- raw |>
    select(protein_accessions, protein_descriptions, protein_names) |>
    distinct() |>
    mutate(id = row_number(), .before = 1)
  
  # just in case, use all isoforms
  contm <- contaminants |> 
    str_remove("\\-\\d+$")
  
  id_prot_gene <- info |>
    select(id, uniprot = protein_accessions) |>
    separate_rows(uniprot, sep = ";") |>
    mutate(uniprot = str_remove(uniprot, "\\-\\d+$")) |>
    filter(!(uniprot %in% contm)) |>
    left_join(uni_gene, by = c("uniprot"), multiple = "all") |>
    mutate(gene_symbol = if_else(is.na(gene_symbol), uniprot, gene_symbol))
  good_ids <- unique(id_prot_gene$id)
  
  gene_names <- id_prot_gene |>
    group_by(id) |>
    summarise(gene_names = str_c(unique(gene_symbol), collapse = ";"))
  
  info <- info |>
    filter(id %in% good_ids) |>
    left_join(gene_names, by = "id")
  
  d <- raw |>
    # rigth_join will remove contaminants
    right_join(select(info, id, protein_accessions), by = "protein_accessions") |>
    left_join(select(meta, data_col, sample), by = "data_col")
  
  dat <- d |>
    select(id, sample, quantity) |>
    drop_na() |>
    normalise_proteins()
  
  qc <- d |>
    select(id, sample, coverage, is_single_hit, qvalue)
  
  list(
    metadata = meta,
    info = info,
    id_prot_gene = id_prot_gene,
    qc = qc,
    dat = dat
  )
}

read_spectronaut_wide_data <- function(file, meta, min_pep) {
  d <- read_csv(file, show_col_types = FALSE) |> 
    rename(id = 1) |> 
    mutate(id = row_number())  # prefer starting from 1
  info <- d |> 
    select(id, starts_with("PG")) |> 
    rename_with(~str_remove(.x, "PG.")) |> 
    janitor::clean_names()
  # columns are complicated, need to simplify naming
  cols <- tibble(name = colnames(d)) |> 
    filter(name != "id" & !(str_detect(name, "^PG"))) |>
    mutate(new = name) |> 
    mutate(new = str_replace(new, "Razor", "PG.Razor")) |> 
    mutate(new = str_replace(new, ".PG.", ":"))
  
  d_ren <- d |> 
    select(c(id, all_of(cols$name))) |> 
    set_names(c("id", cols$new))
  
  # Coverage is a string, cannot pivot it with numerical values
  cover <- d_ren |> 
    select(c(id, contains("Coverage"))) |> 
    pivot_longer(-id, names_to = c("data_col", "quantity"), names_pattern = "(.*):(.*)") |> 
    pivot_wider(id_cols = c(id, data_col), names_from = quantity, values_from = value) |> 
    rename(coverage = Coverage) |> 
    mutate(coverage = str_remove_all(coverage, "%")) |> 
    mutate(coverage = str_remove(coverage, ";.+$")) |> # with multiple hits select the first
    mutate(coverage = as.numeric(coverage) / 100)

  full <- d_ren |> 
    select(-contains("Coverage")) |> 
    pivot_longer(-id, names_to = c("data_col", "quantity"), names_pattern = "(.*):(.*)") |> 
    pivot_wider(id_cols = c(id, data_col), names_from = quantity, values_from = value) |> 
    rename(abu = Quantity, abu_norm = Quantity.Normalised, unique_pep = Razor.Unique.Peptides, q_value = Qvalue, single_hit = IsSingleHit) |> 
    left_join(cover, by = c("id", "data_col")) |> 
    mutate(single_hit = as.logical(single_hit)) |> 
    mutate(abu = log10(abu), abu_norm = log10(abu_norm)) |> 
    inner_join(select(meta, sample, data_col), by = "data_col") |> 
    select(-data_col) |> 
    mutate(sample = as_factor(sample)) |> 
    relocate(sample, .after = id)
  
  
  qc <- full |> 
    select(id, sample, unique_pep, q_value, coverage, single_hit)
  dat <- full |> 
    filter(unique_pep >= min_pep) |> 
    select(-c(unique_pep, q_value, coverage, single_hit)) |> 
    drop_na()
  id_prot_gene <- info |> 
    select(id, gene_symbol = genes, uniprot = protein_accessions) |> 
    separate_rows(c(gene_symbol, uniprot), sep = ";") |> 
    mutate(uniprot = str_remove(uniprot, "-\\d+$")) |>   # redundant isoforms
    distinct()
  
  list(
    info = info,
    id_prot_gene = id_prot_gene,
    qc = qc,
    dat = dat,
    metadata = meta
  )
}


#' Normalize protein quantities
#'
#' This function normalizes protein quantities in a data frame by calculating
#' the log10 abundance, dividing by the median abundance of each sample, and
#' scaling by the overall mean median abundance.
#'
#' @param dat A data frame containing protein quantities, with columns 'id',
#'   'sample', and 'quantity'
#'
#' @return A data frame with normalized protein abundances, containing columns
#'   'id', 'sample', 'abu', and 'abu_norm'
#' 
normalise_proteins <- function(dat) {
  md <- dat |> 
    mutate(abu = log10(quantity)) |> 
    group_by(sample) |> 
    summarise(med = median(abu)) |> 
    mutate(med = med / mean(med))
  dat |> 
    left_join(md, by = "sample") |> 
    mutate(
      abu = log10(quantity),
      abu_norm = abu / med
    ) |> 
    select(-c(quantity, med))
}


select_detected_proteins <- function(set) {
  set$dat |> 
    select(id) |> 
    distinct() |> 
    left_join(set$id_prot_gene, by = "id")
}


#' Convert a data frame to a matrix
#'
#' This function converts a data frame containing protein abundances to a
#' matrix, with rows representing protein IDs and columns representing samples.
#'
#' @param dat A data frame containing protein abundances
#' @param what A character string specifying which column should be used as
#'   values in the matrix (default: "abu_norm")
#' @param names A character string specifying which column should be used as
#'   column names in the matrix (default: "sample")
#'
#' @return A matrix with protein IDs as row names, sample names as column names,
#'   and the specified values in the cells
dat2mat <- function(dat, what = "abu_norm", names = "sample") {
  dat |> 
    pivot_wider(id_cols = id, names_from = !!names, values_from = !!what) |> 
    column_to_rownames("id") |> 
    as.matrix()
}


detection_samples <- function(set) {
  set$dat |> 
    group_by(id) |> 
    summarise(n_det = n()) |> 
    group_by(n_det) |> 
    tally() |> 
    arrange(desc(n_det))
}



make_prot_umap <- function(set, what = "abu_norm", seed = 4654,
                           n_neighbours = 15, min_dist = 0.01, n_dim = 10) {
  d <- set$dat |> 
    mutate(val = get(what)) |> 
    group_by(id) |> 
    mutate(fc = val - mean(val)) |> 
    ungroup()
  
  tab <- dat2mat(d, "fc")
  tab <- tab[apply(tab, 1, function(v) sum(is.na(v)) == 0), ]
  
  set.seed(seed)
  u <- uwot::umap(tab, n_neighbors = n_neighbours, min_dist = min_dist, n_components = n_dim)
  colnames(u) <- paste0("X", 1:n_dim)
  u |> 
    as_tibble(rownames = "id") |> 
    mutate(id = as.integer(id))
}

cluster_dim <- function(d, n_clust) { 
  X <- d |> 
    select(starts_with("X")) |>
    as.matrix()
  km <- kmeans(X, centers = n_clust)
  d |> 
    mutate(clust = km$cluster |> as_factor())
}

#' Get participants with complete data for specified days
#'
#' This function filters a metadata data frame to include only participants with
#' complete data for the specified days and without any bad samples.
#'
#' @param meta A data frame containing metadata for the samples, with columns
#'   'participant_id', 'day', 'sample', 'batch', and 'bad'
#' @param days A numeric vector specifying the days for which participants
#'   should have complete data (default: c(1, 29))
get_full_participants <- function(meta, days = c(1, 29)) {
  meta |> 
    filter(day %in% days & !bad) |> 
    pivot_wider(id_cols = c(participant_id), names_from = day, values_from = c(sample, batch)) |> 
    drop_na() |> 
    set_names(c("participant_id", "first", "last", "first_batch", "last_batch"))
}

#' Calculate log2 fold change between specified days
#'
#' This function calculates the log2 fold change (logFC) of protein abundances
#' between the specified days for a given set of data and metadata.
#'
#' @param set A list containing the following elements:
#' * dat: A data frame with normalized protein quantities
#' * metadata: A data frame containing metadata for the samples
#' * info: A data frame with protein information and gene names
#' * id_prot_gene: A data frame with protein ID, UniProt accession, and gene symbol information
#' @param days A numeric vector specifying the days between which to calculate
#'   logFC (default: c(1, 29))
#' @param min_det A numeric value specifying the minimum number of detections
#'   per treatment to include a protein in the analysis (default: 3)
#'
#' @return A list containing the following elements:
#' * dat: A data frame with logFC values, median-normalized logFC values, and quantile-normalized logFC values
#' * metadata: A data frame containing metadata for the samples
#' * info: A data frame with protein information and gene names
#' * id_prot_gene: A data frame with protein ID, UniProt accession, and gene symbol information
logfc_days <- function(set, days = c(1, 29), min_det = 3) {
  pd <- get_full_participants(set$metadata, days)
  X <- dat2mat(set$dat, what = "abu_norm")
  X <- X / log10(2)   # we are going to find log2_FC
  
  dat <- map_dfr(1:nrow(pd), function(i) {
    r <- pd[i, ]
    tibble(
      id = rownames(X) |> as.integer(),
      participant_id = r$participant_id,
      logFC = X[, r$last] - X[, r$first]
    )
  }) |> 
    drop_na() |> 
    droplevels()
  
  meta <- set$metadata |> 
    select(participant_id, age, sex, treatment, completion, age_group) |> 
    distinct() |> 
    filter(participant_id %in% unique(pd$participant_id)) |> 
    left_join(select(pd, participant_id, first_batch, last_batch), by = "participant_id") |> 
    unite(batches, c(first_batch, last_batch), sep = "-") |> 
    add_column(bad = FALSE)
  
  # remove proteins with less than min_det detection per treatment
  good_prots <- dat |> 
    left_join(meta, by = "participant_id") |> 
    group_by(id, treatment) |> 
    tally() |> 
    filter(n >= min_det) |>
    pull(id) |> 
    unique()
  
  # Median normalisation
  dat <- dat |> 
    filter(id %in% good_prots) |> 
    group_by(participant_id) |> 
    mutate(norm = median(logFC)) |> 
    mutate(logFC_med = logFC - norm) |> 
    select(-norm) |> 
    ungroup()
  
  # Quantile normalisation
  dat_norm <- dat |>
    pivot_wider(id_cols = id, names_from = participant_id, values_from = logFC) |> 
    column_to_rownames("id") |> 
    as.matrix() |> 
    preprocessCore::normalize.quantiles(keep.names = TRUE) |> 
    as_tibble(rownames = "id") |> 
    pivot_longer(-id, names_to = "participant_id", values_to = "logFC_quant") |> 
    mutate(id = as.integer(id)) |> 
    drop_na()
  
  dat <- dat |> 
    left_join(dat_norm, by = c("id", "participant_id"))
  
  list(
    dat = dat,
    metadata = meta,
    info = set$info,
    id_prot_gene = set$id_prot_gene
  )
}



get_file_columns <- function(file) {
  f <- read_tsv(file, n_max = 1)
  tibble(
    column = colnames(f)
  )
}

get_file_nrow <- function(file) {
  read_tsv(file, show_col_types = FALSE) |> 
    nrow()
}

get_file_head <- function(file, n = 10) {
  read_tsv(file, n_max = n, show_col_types = FALSE)
}


add_genes <- function(res, info) {
  g <- info |> 
    mutate(gene_symbol = str_remove(gene_names, ";.+$")) |> 
    select(id, gene_symbol)
  res |> 
    left_join(g, by = "id")
}


save_data_csv <- function(set, file) {
  set$dat |>
    mutate(
      sample = factor(sample, levels = set$metadata$sample),
      abu_norm = signif(abu_norm, 5)
    ) |> 
    arrange(sample, id) |> 
    pivot_wider(id_cols = id, names_from = sample, values_from = abu_norm) |>
    left_join(set$info |> select(id, gene_names, protein_accessions, protein_descriptions), by = "id") |> 
    relocate(c(gene_names, protein_accessions, protein_descriptions), .after = "id") |> 
    write_csv(file)
  return(1)
}

collect_participant_stats <- function(meta) {
  # used to collapse day and on_drug into string
  days <- levels(meta$day) |> as.integer()
  dots <- rep(".", length(days)) |> 
    set_names(days)
  
  dd_ <- function(day, on_drug) {
    d <- dots
    d[day] <- on_drug
    str_c(d, collapse = "")
  }
  
  dd <- meta |> 
    select(participant_id, day, on_drug) |> 
    mutate(on_drug = recode(on_drug, "yes" = "Y", "no" = "N", "baseline" = "B")) |>
    group_by(participant_id) |> 
    summarise(on_drug = dd_(day, on_drug))
  
  meta |> 
    select(participant_id, age, sex, time_from_symptoms, treatment, completion, days_of_treatment) |> 
    distinct() |> 
    left_join(dd, by = "participant_id")
}


read_contaminants <- function(fasta_file) {
  read_lines(fasta_file) |>
    str_subset("^>") |>
    str_remove("\\s.+$") |>
    str_remove(">")
}

