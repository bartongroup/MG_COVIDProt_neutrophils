TAXONOMY_ID <- 9606
KEGG_SPECIES <- "hsa"
ENSEMBL_DATASET <- "hsapiens_gene_ensembl"
ENSEMBL_VERSION <- "106"

UNIPROT_MAPPING_FILE <- "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz"

SPECTRONAUT_FILE <-  "data/20211103_115953_Protein Quant (Normal).tsv"
METADATA_FILE <- "info/20211103_115953_Identifications_STOP-COVID_proteomics_21.4.22.xlsx"

MIN_PEPTIDES <- 3
FDR_LIMIT <- 0.01

BAD_SAMPLES <- c("D-01-267")

make_metadata_ <- function(s_file) {
  read_tsv(s_file, show_col_types = FALSE) %>% 
    mutate(
      sample = raw_sample %>% str_remove("-N"),
      treatment = if_else(str_detect(raw_sample, "ST"), "drug", "placebo"),
      batch = 1
    ) %>% 
    mutate(batch = as_factor(batch)) %>% 
    mutate(day = sample %>%
             str_remove("-\\D+$") %>%
             str_extract("\\d+$") %>%
             na_if("006") %>% na_if("007") %>% 
             as.numeric() %>%
             as_factor()
    )
}


ctr_compare_days <- function(des) {
  limma::makeContrasts(
    treat_daydrug_1-treat_dayplacebo_1,
    treat_daydrug_8-treat_dayplacebo_8,
    treat_daydrug_15-treat_dayplacebo_15,
    treat_daydrug_29-treat_dayplacebo_29,
    levels = des
  )
}

