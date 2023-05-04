library(tidyverse)


load_params <- function(params_file) {
  
  quart <- str_extract(params_file, "q\\d") %>% 
    str_replace('q', '') %>% as.integer()
  
  vcat <- str_extract(params_file, "v_(intergenic|intronic)") %>% 
    str_replace('v_', '')
  
  
  read_tsv(params_file) %>% 
    mutate(
      quartile = quart,
      vcat = vcat
    )
}


params <- list.files('results/ConfidenceIntervals/',
           pattern = 'q*NAT-EXPANSION-v_(intergenic|intronic).tsv',
           full.names = T) %>% 
  map_df(load_params)

params <- 
  params %>% 
  rename(parameter = `#param`) %>% 
  mutate(
    ci_h = opt_value + 2 * std_err,
    ci_l = opt_value - 2 * std_err,
  )

params %>% 
  filter(parameter == 'TN') %>% 
  ggplot(aes(x = opt_value, y = quartile)) +
  geom_point(
    shape = 21
  ) +
  geom_errorbarh(aes(xmin=ci_l, xmax=ci_h), height = 0.3) +
  facet_grid(vcat~.) +
  scale_y_continuous(breaks = 1:4) +
  scale_x_continuous(labels = scales::label_number(scale =  1/1e3,suffix = 'K')) +
  theme_bw() +
  theme(
    panel.grid.minor.y  = element_blank(),
    legend.position = 'none'
  )
