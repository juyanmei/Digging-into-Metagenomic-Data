library(dplyr)
library(readr)
library(here)
library(stringr)

sample_df <- 
  read_csv(here("data", "rat_245S_baseinfo.csv")) %>% 
  select(-sample_id)

sample_df_ <- 
  read_delim(here("data", "SE_245_ebi.txt"), delim = "\t", col_names = F) %>%
  rename(fastq_id = X3) %>%
  group_by(fastq_id) %>%
  mutate(sample_id = str_split(X4, "\\.")[[1]][1]) %>%
  select(sample_id, fastq_id) %>%
  inner_join(., sample_df) %>%
  write_csv(., here("data", "rat_245S_baseinfo_v2.csv"))
