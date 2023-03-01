library(tidyverse)


sfs_sim <- read_csv('../220728-Simulation-DFE-Demography/results/all-SFSs.csv')

sfs_data <- read_csv('../220901-Data-Genomes-SFS/results/sfs-data.csv')


ORDER_MUT_CATS <- c('noncoding', 'synonymous', 'missense', 'LOF')
sfs_sim$Mut_Type <- factor(sfs_sim$Mut_Type, levels = ORDER_MUT_CATS)
sfs_data$Mut_Type <- factor(sfs_data$Mut_Type, levels = ORDER_MUT_CATS)

pop_order <- c('YRI', 'CHB', 'IBS', 'MXB', 'MXL')
sfs_data$Pop <- factor(sfs_data$Pop, levels = pop_order)
sfs_sim$Pop <- factor(sfs_sim$Pop, levels = pop_order)

# turn of on with comments ------------------------------------------------
# This normalizes the spectrums to proportions

cols_to_mutate <- colnames(sfs_sim)[4:10]
sfs_sim <- sfs_sim %>%
  group_by(Pop, Mut_Type) %>%
  mutate_at(
    .vars = cols_to_mutate,
    .funs = function(x) x / sum(x)

  )

sfs_data <- sfs_data %>%
  group_by(Pop, Mut_Type) %>%
  mutate(
    Freq = Freq / sum(Freq)
  )
  
# plot --------------------------------------------------------------------


sfs_sim %>% 
  ggplot(aes(x = minor_allel_freq)) +
  geom_ribbon(
    aes(ymin = min, ymax = max),
    alpha = 0.2,
    fill =  'darkorange'
  ) +
  geom_point( 
    data = sfs_data,
    aes(y = Freq),
    color='dodgerblue4',
    shape = 19,
    size = 1.7
  ) +
  geom_point(
    data = sfs_data,
    aes(y = Freq),
    color='white',
    shape = 19,
    size = 0.3
  ) +
  geom_line(
    aes(y = Freq_q_5),
    size = 0.2,
    color = 'darkorange'
  ) +
  geom_point(
    aes(y = Freq_q_5),
    size = 0.7,
    alpha = 0.85,
    color = 'darkorange'
  ) +
  scale_y_log10() +
  facet_grid(Mut_Type ~ Pop, scales = 'free_y') +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  ) +
  labs(
    x = 'Minor Allele Frequency',
    y = 'Frequency'
  )

ggsave('plots/sfs-data-model.pdf', height = 4, width = 6)



# plot only MXL for main figure -------------------------------------------

sfs_sim <- sfs_sim %>% 
  filter(Pop == 'MXL')
  
sfs_data <- sfs_data %>% 
  filter(Pop == 'MXL')
  
sfs_sim %>% 
  ggplot(aes(x = minor_allel_freq)) +
  geom_ribbon(
    aes(ymin = min, ymax = max),
    alpha = 0.2,
    fill =  'darkorange'
  ) +
  geom_point( 
    data = sfs_data,
    aes(y = Freq),
    color='dodgerblue4',
    shape = 19,
    size = 1.7
  ) +
  geom_point(
    data = sfs_data,
    aes(y = Freq),
    color='white',
    shape = 19,
    size = 0.3
  ) +
  geom_line(
    aes(y = Freq_q_5),
    size = 0.2,
    color = 'darkorange'
  ) +
  geom_point(
    aes(y = Freq_q_5),
    size = 0.7,
    alpha = 0.85,
    color = 'darkorange'
  ) +
  scale_y_log10() +
  facet_grid(.~Mut_Type, scales = 'free_y') +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  ) +
  labs(
    x = 'Minor Allele Frequency',
    y = 'Frequency'
  )

ggsave('plots/sfs-data-model-MXL.pdf', height = 2, width = 6)
