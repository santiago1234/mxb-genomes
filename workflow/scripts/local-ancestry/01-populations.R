library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# parameters

oneTGP_pops <-  snakemake@params$oneTGP_pops
path_oneT_meta <- snakemake@input$oneT_meta
path_mxb_meta <- snakemake@input$mxb_meta

# outfiles
sample_map_file <- snakemake@output$sample_map_file
query_pops_file <- snakemake@output$query_pops



# filter populations ------------------------------------------------------


oneT_meta <- read_tsv(path_oneT_meta)

oneTGP_pops <- 
  oneT_meta %>%
  select(`Sample name`, `Population code`, `Superpopulation code`) %>% 
  filter(str_detect(`Population code`, oneTGP_pops)) %>% 
  rename(Sample = `Sample name`, Population = `Superpopulation code`) 

mxb_meta <- read_tsv(path_mxb_meta) %>% 
  mutate(
    Sample = str_replace(Sample_ID, "MXB", "MXB_"),
    Population = "NAT",
    `Population code` = "MXB"
  ) %>% 
  select(Sample, Population, `Population code`)


populations <- 
  bind_rows(oneTGP_pops, mxb_meta) %>% 
  arrange(Sample, Population, `Population code`)




sample_map <- filter(populations, Population != "AMR")
write_tsv(sample_map, sample_map_file)


query_pop <- filter(populations, Population == "AMR")
write_tsv(query_pop, query_pops_file)