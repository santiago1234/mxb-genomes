library(tidyverse)

fst_data <- read_csv('../220910-Fst-Comparison/results/data-fst.csv') %>% 
  mutate(
    sim_source = 'data'
  )

fst_fwdpy11 <- list.files('../220910-Fst-Comparison/results/sim-fst/', full.names = T) %>% 
  map_df(read_csv)  %>% 
  mutate(
    sim_source = 'fwdpy11'
  )

fst_msprime <- read_csv('results/fst-msprime-sim.csv') %>% 
  mutate(
    pop_pair = str_c(Pop1, 'x', Pop2)
  ) %>% 
  select(-Pop1, -Pop2) %>% 
  mutate(
    sim_source = 'msprime'
  )


data <- bind_rows(fst_data, fst_fwdpy11, fst_msprime)


data <- data %>% 
  filter(
    sim_id <= 50,
    category %in% c('synonymous', 'noncoding')
  )



order_pops <- data %>% 
  filter(
    sim_source == 'data',
    category == 'noncoding'
  ) %>% 
  group_by(pop_pair) %>% 
  summarise(median_fst = median(fst)) %>% 
  arrange(-median_fst) %>% 
  pull(pop_pair)

data$pop_pair <- factor(data$pop_pair, levels = order_pops)
data$sim_source <- factor(data$sim_source, levels = c('fwdpy11', 'msprime', 'data'))


data %>% 
  filter(sim_source != 'msprime') %>% 
  ggplot(aes(x = pop_pair, y = fst, fill = sim_source)) +
  geom_boxplot(size = 1/3, outlier.shape = '.', alpha = 0.9) +
  facet_grid(~category) +
  coord_flip(ylim = c(0, 0.35)) +
  scale_fill_manual(values = c('#a6cee3', '#b2df8a')) +
  theme_bw() +
  labs(
    y = 'Fst',
    x = 'Pop1 x Pop2'
  )

ggsave('plots/fst-comparison1.pdf', height = 4, width = 6)


data %>% 
  ggplot(aes(x = pop_pair, y = fst, fill = sim_source)) +
  geom_boxplot(size = 1/3, outlier.shape = '.', alpha = 0.9) +
  facet_grid(~category) +
  coord_flip(ylim = c(0, 0.35)) +
  scale_fill_manual(values = c('#a6cee3', '#1f78b4', '#b2df8a')) +
  theme_bw() +
  labs(
    y = 'Fst',
    x = 'Pop1 x Pop2'
  )

ggsave('plots/fst-comparison2.pdf', height = 4, width = 6)

