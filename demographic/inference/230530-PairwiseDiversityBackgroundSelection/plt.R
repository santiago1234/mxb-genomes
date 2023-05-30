library(tidyverse)

pis <- read_csv('results/pi.csv')

order_pops <- c('YRI', 'IBS', 'CHB', 'MXB')
  
  
pis$pi <- pis$pi / pis$ml 
pis$pop <- factor(pis$pop, levels = order_pops)

pis %>% 
  ggplot(aes(x = pi, y = quartile, fill = mut_type)) +
  geom_point(shape = 21) +
  facet_grid(pop~.) +
  scale_fill_manual(values = c('#1b9e77', '#e7298a')) +
  labs(
    x = 'pi / mL',
    y = 'B-statistic (quartile)'
  ) +
  theme_bw() +
  theme(
    panel.grid.minor.y = element_blank()
  )

ggsave('plots/pi.pdf', height = 4, width = 5)

pis %>% 
  ggplot(aes(x = ml, y = quartile)) +
  geom_point() +
  facet_grid(pop~mut_type, scales = 'free_x')
ggsave('plots/ml.pdf', height = 4, width = 5)
