library(tidyverse)
library(ggforce)

data_fst <- read_csv('results/data-fst.csv') %>% 
  rename(
    fst_data = fst
  )

sim_fst <- list.files('results/sim-fst/', full.names = T) %>% 
  map_df(read_csv) %>% 
  rename(
    fst_sim = fst
  )

fst <- inner_join(data_fst, sim_fst, by = c('category', 'sim_id', 'pop_pair'))

fst <- 
  fst %>% 
  mutate(
    pop_pair = str_replace(pop_pair, 'x', '-')
  ) %>% 
  pivot_longer(cols = c(fst_data, fst_sim), names_to = 'sim_or_data', values_to = 'fst')

order_pair <- fst %>% 
  filter(category == 'noncoding', sim_or_data == 'fst_data') %>% 
  group_by(pop_pair) %>% 
  summarise(
    m_val = mean(fst)
  ) %>% 
  arrange(m_val) %>% 
  pull(pop_pair) %>% 
  rev()

fst$pop_pair <- factor(fst$pop_pair, levels = order_pair)

fst$category <- factor(fst$category, levels = c('noncoding', 'synonymous', 'nonsynonymous'))

fst %>% 
  ggplot(aes(x = fst, y = pop_pair, fill = sim_or_data)) +
  geom_boxplot(outlier.shape = NA, size = 0.2, alpha = 0.8) +
  scale_fill_manual(values = c('darkorange', 'dodgerblue4')) +
  coord_cartesian(xlim = c(0.01/2, 0.4)) +
  scale_x_sqrt(breaks = c(0.01, 0.1, 0.2, 0.3)) +
  facet_grid(~category) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave('plots/fst-comparison.pdf', height = 2.5, width = 8)
