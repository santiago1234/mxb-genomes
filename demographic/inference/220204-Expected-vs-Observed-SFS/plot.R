library(tidyverse)
library(scales)

d <- read_csv('results/model.csv')
d %>% 
  ggplot(aes(Frq1, Frq2, fill= log10(Density))) +
  geom_raster() +
  scale_fill_viridis_c(limits=c(2, 4), oob=squish) 


load_data <- function(df) {
  
  datype <- df %>% str_extract('(model|data|residuals)')
  
  read_csv(df) %>% 
    select(-Pops) %>% 
    pivot_longer(cols=-c(Pop1), values_to = 'Density', names_to='Pop2') %>% 
    separate(Pop1, into=c('Pop_1', 'Frq1'), sep='_') %>% 
    separate(Pop2, into=c('Pop_2', 'Frq2'), sep='_') %>% 
    mutate(Comparison = paste0(Pop_1, '-',Pop_2)) %>% 
    mutate(
      Frq1 = as.integer(Frq1),
      Frq2 = as.integer(Frq2),
    ) %>% 
    filter(
      !(Frq1 == 0 & Frq2 == 0),
      !(Frq1 == 20 & Frq2 == 20)
    ) %>% 
    mutate(
      D_type = datype
    )
  
}


all_d <- list.files('results/', pattern = '*csv', full.names = T) %>% 
  map_df(load_data)


all_d %>% 
  filter(
    D_type != 'residuals'
  ) %>% 
  ggplot(aes(Frq1, Frq2, fill= log10(Density))) +
  geom_raster() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_c(limits=c(2, 4), oob=squish, option = 'E') +
  facet_grid(D_type~Comparison) +
  theme(
    panel.background = element_blank()
  )
  

all_d %>% 
  filter(
    D_type == 'residuals'
  ) %>% 
  ggplot(aes(Frq1, Frq2, fill= Density)) +
  geom_raster() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradient2(low = 'blue', high = 'red',mid = 'white', limits=c(-110, 110), oob=squish) +
  facet_grid(D_type~Comparison) +
  theme(
    panel.background = element_blank()
  )

all_d %>% 
  filter(
    D_type == 'residuals'
  ) %>% 
  ggplot(
    aes(x = Density)
  ) +
  geom_histogram(binwidth = 5) +
  facet_wrap(.~Comparison, scales = 'free_y', nrow = 1) +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank()
  )
