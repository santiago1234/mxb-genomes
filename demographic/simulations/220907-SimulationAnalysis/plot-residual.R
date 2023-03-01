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


residual <- sfs_sim %>% 
  select(Pop, Mut_Type, minor_allel_freq, Freq_q_5) %>% 
  inner_join(sfs_data) %>% 
  mutate(
    poisson_residual = (Freq_q_5 - Freq) / sqrt(Freq_q_5)
  )


residual %>% 
  ggplot(aes(x = minor_allel_freq, y = poisson_residual)) +
  geom_point(fill = 'grey', shape=21) +
  geom_hline(yintercept = 0, linetype=2) +
  facet_grid(Mut_Type ~ Pop) +
  coord_cartesian(ylim = c(-0.15, 0.15)) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  ) +
  labs(
    x = 'Minor Allele Frequency',
    y = 'Residual'
  )

ggsave('plots/residuals.pdf', width = 6, height = 4)
