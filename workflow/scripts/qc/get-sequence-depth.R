library(vcfR)
library(readr)
library(dplyr)
library(tidyr)

# this script extracts the sequence depth
# it only selects a random fraction 
# missing data can be determined from here
# since missing depth means missing genotype
# see: https://grunwaldlab.github.io/Population_Genetics_in_R/qc.html
# input parameters --------------------------------------------------------

# write_rds(snakemake, "snake.rds")
vcf_file <- snakemake@input[[1]]
depth_outfile <- snakemake@output[[1]]
chrn <- snakemake@wildcards$chrn

# fraction of variants to sample
fraction <- snakemake@params$fraction

set.seed(42)

# load vcf ----------------------------------------------------------------


vcf <- read.vcfR(vcf_file)

dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE) %>% 
  as_tibble(rownames = "varid")

dp <- sample_frac(dp, size = fraction)

# make tidy frame

dp <- 
  dp %>% 
  pivot_longer(cols = -c(varid), values_to = "DP", names_to = "Individual")

dp %>%
  mutate(chrn = chrn) %>%
  write_csv(file = depth_outfile)