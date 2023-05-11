library(tidyverse)
library(ggthemes)
# load data and combine data sets -----------------------------------------

d <- read_csv('results/pca-real.csv')

d_genomes <- 
  d %>% 
  select(PC_1:PC_3, Samplename, Subpopulation) %>% 
  mutate(
    source = 'Data\n(genomes)'
  )

d_browning <- read_csv('results/pca-BrowningEtAl2011.csv')
d_browning <- 
  d_browning %>% 
  select(PC_1:PC_3, Samplename) %>% 
  mutate(
    Subpopulation = str_sub(Samplename, 1, 3),
    source = 'Browning et al 2018\n(simulation)'
  )


d_mdl <- read_csv('results/pca-Medina2022.csv')

d_mdl <- 
  d_mdl %>% 
  select(PC_1:PC_3, Samplename) %>% 
  mutate(
    Subpopulation = str_sub(Samplename, 1, 3),
    source = 'Inferred model\n(simulation)',
    PC_2 = -PC_2,
    PC_1 = -PC_1
  )

d <- bind_rows(d_genomes, d_mdl, d_browning)
d$source <- factor(d$source, levels = 
                     c(
                       'Data\n(genomes)',
                       'Inferred model\n(simulation)',
                       'Browning et al 2018\n(simulation)'
                     )
                     )

d_amix <- filter(d, Subpopulation == 'MXL')

d %>% 
  filter(Subpopulation != 'MXL') %>% 
  ggplot(aes(x = PC_1, y = PC_2)) +
  geom_point(aes(fill = Subpopulation), shape = 21, size = 2) +
  geom_point(
    data = d_amix, shape = 4, color = '#440154'
  ) +
  facet_wrap(~source, scales='free') +
  scale_fill_manual(values = c('#21918c', '#5ec962', '#fde725')) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),legend.position = 'bottom'
  )
ggsave('plots/compare-pc1-pc2.pdf', height = 2.5, width = 5)

d %>% 
  filter(Subpopulation != 'MXL') %>% 
  ggplot(aes(x = PC_2, y = PC_3)) +
  geom_point(aes(fill = Subpopulation), shape = 21, size = 2) +
  geom_point(
    data = d_amix, shape = 4, color = '#440154'
  ) +
  facet_wrap(~source, scales='free') +
  scale_fill_manual(values = c('#21918c', '#5ec962', '#fde725')) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = 'bottom'
  )
ggsave('plots/compare-pc2-pc3.pdf', height = 2.5, width = 5)