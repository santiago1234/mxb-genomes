library(tidyverse)

source('uitls.R')

sfs_sim <- read_csv('../220728-Simulation-DFE-Demography/results/all-SFSs.csv')

sfs_data <- read_csv('../220901-Data-Genomes-SFS/results/sfs-data.csv')

sfs_sim$Mut_Type <- factor(sfs_sim$Mut_Type, levels = ORDER_MUT_CATS)
sfs_data$Mut_Type <- factor(sfs_data$Mut_Type, levels = ORDER_MUT_CATS)


# turn of on with comments ------------------------------------------------

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

ggsave('plots/sfs-data-model.pdf', height = 6, width = 7.5)



# yet another way ---------------------------------------------------------

sfs_sim %>% 
  filter(Mut_Type != 'LOF') %>% 
  select(Pop:minor_allel_freq, Freq_q_5) %>% 
  group_by(Pop, Mut_Type) %>% 
  mutate(
    Freq = Freq_q_5 / sum(Freq_q_5)
  ) %>% 
  filter(Pop == 'MXL') %>% 
  ggplot(aes(x = minor_allel_freq, y = Freq, color = Mut_Type)) +
  geom_line() +
  scale_color_manual(values = MUT_COLORS) +
  geom_point() +
  scale_y_log10()

# more --------------------------------------------------------------------



sfs_sim %>% 
  ggplot(aes(x = minor_allel_freq)) +
  geom_ribbon(
    aes(ymin = Freq_q_25, ymax = Freq_q_75),
    alpha = 0.2,
    fill =  'darkorange'
  ) +
  geom_ribbon(
    aes(ymin = Freq_q_05, ymax = Freq_q_95),
    alpha = 0.2,
    fill =  'darkorange'
  ) +
  geom_ribbon(
    aes(ymin = min, ymax = max),
    alpha = 0.2,
    fill =  'darkorange'
  ) +
  # geom_point(
  #   aes(y = Freq_q_5),
  #   color =  'darkorange'
  # ) +
  geom_point(
    data = sfs_data,
    aes(y = Freq),
    color='dodgerblue4',
    shape = 1
  ) +
  # geom_point(
  #   data = sfs_data,
  #   aes(y = Freq),
  #   color='white',
  #   size = 0.5
  # ) +
  geom_line(
    aes(y = Freq_q_5),
    color = 'darkorange'
  ) +
  geom_line(
    aes(y = Freq_q_5),
    color = 'black',
    size = 0.1
  ) +
  scale_y_log10() +
  facet_grid(Mut_Type ~ Pop, scales = 'free_y') +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )



sfs_sim %>% 
  ggplot(aes(x = minor_allel_freq)) +
  geom_ribbon(
    aes(ymin = Freq_q_25, ymax = Freq_q_75),
    alpha = 0.2,
    fill =  'darkorange'
  ) +
  geom_ribbon(
    aes(ymin = Freq_q_05, ymax = Freq_q_95),
    alpha = 0.2,
    fill =  'darkorange'
  ) +
  geom_ribbon(
    aes(ymin = min, ymax = max),
    alpha = 0.2,
    fill =  'darkorange'
  ) +
  # geom_point(
  #   aes(y = Freq_q_5),
  #   color =  'darkorange'
  # ) +
  geom_point(
    data = sfs_data,
    aes(y = Freq),
    color='dodgerblue4',
    shape = 1
  ) +
  # geom_point(
  #   data = sfs_data,
  #   aes(y = Freq),
  #   color='white',
  #   size = 0.5
  # ) +
  geom_line(
    aes(y = Freq_q_5),
    color = 'darkorange'
  ) +
  geom_line(
    aes(y = Freq_q_5),
    color = 'black',
    size = 0.1
  ) +
  scale_y_log10() +
  facet_grid(Mut_Type ~ Pop, scales = 'free_y') +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )


# syn VS miss -------------------------------------------------------------

sfs_sim %>% 
  group_by(Mut_Type, Pop) %>% 
  mutate(
    Freq = Freq_q_5 / sum(Freq_q_5)
  ) %>% 
  ggplot(aes(x = minor_allel_freq, y = Freq, color = Mut_Type)) +
  scale_y_log10() +
  geom_line(size = 0.4) +
  geom_point() +
  facet_grid(Pop ~ ., scales = 'free_y') +
  scale_color_manual(values = MUT_COLORS) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )
  
sfs_data %>% 
  filter(Pop == 'MXL') %>% 
  group_by(Mut_Type, Pop) %>% 
  mutate(
    Freq = Freq / sum(Freq)
  ) %>% 
  ggplot(aes(x = minor_allel_freq, y = Freq, color = Mut_Type)) +
  scale_y_log10() +
  geom_line(size = 0.4) +
  geom_point() +
  facet_grid(. ~ Pop, scales = 'free_y') +
  scale_color_manual(values = MUT_COLORS) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )


# neutral between pops ----------------------------------------------------

sfs_sim %>% 
  filter(Mut_Type == 'noncoding') %>% 
  ggplot(aes(x = minor_allel_freq, y = Freq_q_5, fill = Pop)) +
  scale_y_log10() +
  geom_line(aes(color = Pop),size = 0.4) +
  geom_point(shape=21) +
  scale_fill_manual(values = POP_COLORS) +
  scale_color_manual(values = POP_COLORS) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )



