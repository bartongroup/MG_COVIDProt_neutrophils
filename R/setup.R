TAXONOMY_ID <- 9606
SPECIES <- "human"
ENSEMBL_DATASET <- "hsapiens_gene_ensembl"
ENSEMBL_VERSION <- "106"

CHROMOSOMES <- c(1:22, "X", "Y", "MT")

# Contaminant file downloaded from
# http://lotus1.gwdg.de/mpg/mmbc/maxquant_input.nsf/7994124a4298328fc125748d0048fee2/$FILE/contaminants.fasta
# on 10 Nove 2022
CONTAMINANT_FILE <- "info/contaminants.fasta"

UNIPROT_MAPPING_FILE <- "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz"

SPECTRONAUT_FILE <-  "data/20211103_115953_Protein Quant (Normal).tsv"
COVID_METADATA_FILE <- "info/20211103_115953_Identifications_STOP-COVID_proteomics_21.4.22.xlsx"

MIN_PEPTIDES <- 3
FDR_LIMIT <- 0.01

COVID_BAD_SAMPLES <- c("B5_01-45")

make_metadata_ <- function(s_file) {
  read_tsv(s_file, show_col_types = FALSE) |> 
    mutate(
      sample = data_col |> str_remove("-N"),
      treatment = if_else(str_detect(data_col, "ST"), "Brensocatib", "lacebo"),
      batch = 1
    ) |> 
    mutate(batch = as_factor(batch)) |> 
    mutate(day = sample |>
             str_remove("-\\D+$") |>
             str_extract("\\d+$") |>
             na_if("006") |> na_if("007") |> 
             as.numeric() |>
             as_factor()
    )
}

BATCH_EXAMPLES <- c(419, 4291)
BASE_FILTER <- "completion & batch %in% c(3, 4, 5)"

ENRICHMENT_EXAMPLES <- tibble::tribble(
  ~contrast, ~ontology, ~term_id,
  "treatmentBrensocatib", "go", "GO:0008236",
  "treatmentBrensocatib", "go", "GO:0035578",
  "treatmentBrensocatib", "kg", "hsa04613",
  "treatmentBrensocatib", "re", "R-HSA-6798695",
  "age_group>80", "go", "GO:0072562",
  "age_group>80", "go", "GO:0031720",
  "age_group50-65", "re", "R-HSA-1247673"
)


