library(tidyverse)

# uncoment to test the code
# input_tracts_file <- "results/3-pops/tract-distribution-by-population.csv"
input_tracts_file <- snakemake@input[[1]] 
output_plot <- snakemake@output[[1]]

tracts <- read_csv(input_tracts_file)

get_position <- function(bin) {
  bin %>% 
    str_replace("-\\d*$", "") %>% 
    as.numeric()
}
tracts$tract_len_cm <- get_position(tracts$tract_length)


tracts %>% 
  filter(freq_mean > 0,) %>% 
  ggplot(aes(x = tract_len_cm, y = freq_mean, fill=Ancestry)) +
  geom_line(linetype = 2, size = .1) +
  geom_point(shape = 21, size = 2, alpha = .9) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000)) +
  scale_x_continuous(breaks = c(0, 100, 200)) +
  scale_fill_viridis_d(option = "C") +
  facet_grid(~Subpopulation) +
  theme_bw() +
  labs(
    x = "Tract Length (cM)",
    y = "Relative frequency"
  )
  
ggsave(output_plot, height = 2, width = 8)
ggsave(str_replace(output_plot, ".png", ".pdf"), height = 2, width = 8)
