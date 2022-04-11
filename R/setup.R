SP_FILE = "data/march_report.csv"


make_metadata <- function(s_file) {
  read_tsv(s_file) %>% 
    mutate(
      sample = raw_sample %>% str_remove("-N"),
      group = if_else(str_detect(raw_sample, "ST"), "treat", "ctrl"),
      batch = 1
    ) %>% 
    mutate(batch = as_factor(batch)) %>% 
    mutate(day = sample %>%
             str_remove("-\\D+$") %>%
             str_extract("\\d+$") %>%
             as.numeric() %>%
             as_factor()
    )
}