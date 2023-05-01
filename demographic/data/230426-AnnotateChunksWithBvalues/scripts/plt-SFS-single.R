library(tidyverse)

sfs <- read_csv('data/whole-genome/sfs-single-pop.csv')
fst <- read_csv('data/whole-genome/fst.csv')

sfs <- sfs %>% 
  mutate(
    Freq_scaled = Frequency / mL
  )

sfs %>% 
  ggplot(aes(x = Minor_allel_freq, y = Freq_scaled, color = Quantile, group = Quantile)) +
  geom_line(size = 1) +
  geom_line(size = 0.3, color = 'black') +
  geom_point(shape='.') +
  facet_grid(~Population) +
  scale_color_viridis_c(option = 'A') +
  scale_y_log10() +
  theme_bw()
ggsave('plots/sfs-by-B-quartile.pdf', height = 3.5, width = 9)

fst %>% 
  mutate(
    pair = paste(pop1, pop2, sep = '-')
  ) %>% 
  ggplot(aes(x = Fst, y = reorder(pair, Fst, median), fill = Quantile)) +
  scale_fill_viridis_c(option = 'A') +
  geom_point(shape = 21)
