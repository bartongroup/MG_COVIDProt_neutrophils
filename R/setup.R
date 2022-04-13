TAXONOMY_ID <- 9606
KEGG_SPECIES <- "hsa"
ENSEMBL_DATASET <- "hsapiens_gene_ensembl"
ENSEMBL_VERSION <- "106"

SP_FILE <-  "data/march_report.csv"

MIN_PEPTIDES <- 3
FDR_LIMIT <- 0.01


make_metadata <- function(s_file) {
  read_tsv(s_file, show_col_types = FALSE) %>% 
    mutate(
      sample = raw_sample %>% str_remove("-N"),
      group = if_else(str_detect(raw_sample, "ST"), "drug", "placebo"),
      batch = 1
    ) %>% 
    mutate(batch = as_factor(batch)) %>% 
    mutate(day = sample %>%
             str_remove("-\\D+$") %>%
             str_extract("\\d+$") %>%
             na_if("006") %>% na_if("007") %>%
             as.numeric() %>%
             recode(`7` = 8, `14` = 15) %>%
             as_factor()
    )
}