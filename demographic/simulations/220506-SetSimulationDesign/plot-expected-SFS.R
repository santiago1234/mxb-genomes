library(tidyverse)
library(gridExtra)

sfs_expected <- read_csv('results/expected-sfs/SFS.csv')


p1 <- sfs_expected %>% 
  filter(
    !DerivedFreq %in% c(min(DerivedFreq), max(DerivedFreq))
  ) %>% 
  group_by(MutType) %>% 
  mutate(
    Frequency = Frequency / sum(Frequency)
  ) %>% 
  ggplot(aes(x = DerivedFreq, y = Frequency, color = MutType)) +
  geom_point(aes(shape = MutType)) +
  scale_y_log10() +
  scale_color_manual(values = c('#66a61e', '#7570b3', '#e7298a', '#d95f02')) +
  labs(
    title = 'Proportion in each bin'
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )


p2 <- sfs_expected %>% 
  filter(
    !DerivedFreq %in% c(min(DerivedFreq), max(DerivedFreq))
  ) %>% 
  ggplot(aes(x = DerivedFreq, y = Frequency, color = MutType)) +
  geom_point(aes(shape = MutType)) +
  scale_y_log10() +
  scale_color_manual(values = c('#66a61e', '#7570b3', '#e7298a', '#d95f02')) +
  labs(
    title = 'Theta = 1'
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )

pdf('plots/Expected-SFS-under-selection.pdf', height = 5.5, width = 6)
grid.arrange(p1, p2)
dev.off()
