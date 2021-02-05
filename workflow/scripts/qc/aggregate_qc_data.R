library(readr)
library(dplyr)
library(purrr)

# input parameters --------------------------------------------------------

seq_depth_files <- snakemake@input$seq_deps
nvars_file <- snakemake@input$n_vars
output_nvras <- snakemake@output$vars_per_genome
output_depth <- snakemake@output$seqs_deps
# agg nvars file ----------------------------------------------------------


load_nvars_file <- function(x) {
  header <- c("Individual", "n_variants", "chr")
  read_delim(nvars_file[1], col_names = header, delim = " ")

}

nvars <- map_df(nvars_file, load_nvars_file)

nvars <-
  nvars %>%
  group_by(Individual) %>%
  summarise(n_variants = sum(n_variants))

write_csv(nvars, output_nvras)


# merge sequence depth ----------------------------------------------------

map_df(seq_depth_files, read_csv) %>% 
  write_csv(output_depth)

