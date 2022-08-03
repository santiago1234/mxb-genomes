library(tidyverse)


sfs_sim <- read_csv('results/SFSs_simulation.csv') %>% 
  mutate(
    SFS_FROM = 'simualtion'
  ) %>% 
  rename(
    Frequency_sim = Frequency
  )

sfs_expected <- read_csv('results/expected-sfs/SFS_scaled_to_data.csv') %>% 
  mutate(
    SFS_FROM = 'expected'
  ) %>% 
  rename(
    Frequency_expct = Frequency
  )



# combine results ---------------------------------------------------------

# sfs_expected <- sfs_expected %>% 
#   mutate(
#     MutType = if_else(MutType == 'noncoding', 'neutral', MutType)
#   )

sfs <- inner_join(sfs_expected, sfs_sim, by = c('DerivedFreq', 'MutType')) %>% 
  filter(
    !DerivedFreq %in% c(min(DerivedFreq), max(DerivedFreq))
  )


sfs$MutType <- factor(sfs$MutType, levels = c('noncoding', 'synonymous', 'missense', 'LOF'))


sfs %>% 
  ggplot(aes(x = DerivedFreq)) +
  geom_point(aes(y = Frequency_sim), color='dodgerblue4', size = 1.2) +
  geom_line(aes(y = Frequency_expct),
            color =  'darkorange',
            linetype = 2,
            size = 1 / 2
            ) +
  geom_point(
    aes(y = Frequency_expct),
    color =  'darkorange',
    size = 0.1
  ) +
  scale_y_log10() +
  facet_grid(~MutType) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  ) +
  labs(
    y = 'Frequency',
    title = 'Nsim = 100, L = 1Mb'
  )
ggsave('plots/Expected-vs-observed-SFS.pdf',height = 3, width = 6.5)

