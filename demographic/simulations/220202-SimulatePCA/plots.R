library(tidyverse)
library(ggthemes)
library(gridExtra)

theme_set(theme_tufte(base_family = 'Helvetica'))

pcd <- read_csv('results/PCA-real-data.csv')

continental <- c('YRI', 'IBS', 'MXB', 'CHB')
continental_pc <- pcd %>% filter(Subpopulation %in% continental)
admix_pc <- pcd %>% filter(!Subpopulation %in% continental)

admix_pc_tmp <- select(admix_pc, -Subpopulation)

admix_pc <- 
  admix_pc %>% 
  mutate(
    ADMIX = Subpopulation,
    Subpopulation = 'ADMIX'
  )


plot_theme <-
  theme(
    legend.position = 'None',
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_rect(fill = NA)
  )

pr_12 <- admix_pc %>% 
  ggplot(aes(x = PC_1, y = PC_2)) +
  geom_point(data = admix_pc_tmp, color = 'grey', shape = 4, size = 1.0) +
  geom_point(data = continental_pc, aes(fill = Subpopulation), shape = 21) +
  geom_point(shape = 4, color = '#440154', size = 1.0) +
  scale_fill_manual(values = c('#21918c', '#5ec962', '#3b528b', '#fde725')) +
  facet_wrap(~ADMIX) +
  labs(
    title = 'Real data'
  ) +
  plot_theme
  

pr_23 <- admix_pc %>% 
  ggplot(aes(x = PC_2, y = PC_3)) +
  geom_point(data = admix_pc_tmp, color = 'grey', shape = 4, size = 1.0) +
  geom_point(data = continental_pc, aes(fill = Subpopulation), shape = 21) +
  geom_point(shape = 4, color = '#440154', size = 1.0) +
  scale_fill_manual(values = c('#21918c', '#5ec962', '#3b528b', '#fde725')) +
  facet_wrap(~ADMIX) +
  labs(
    title = 'Real data'
  ) +
  plot_theme



pcd <- read_csv('results/PCA-simulated-data.csv')

pcd <- 
  pcd %>% 
  separate(Samplename, into = c('Subpopulation', 'Samplename'), sep = '_')


continental <- c('YRI', 'IBS', 'MXB', 'CHB')
continental_pc <- pcd %>% filter(Subpopulation %in% continental)
admix_pc <- pcd %>% filter(!Subpopulation %in% continental)

admix_pc_tmp <- select(admix_pc, -Subpopulation)

admix_pc <- 
  admix_pc %>% 
  mutate(
    ADMIX = Subpopulation,
    Subpopulation = 'ADMIX'
  )


ps_12 <- admix_pc %>% 
  ggplot(aes(x = -PC_1, y = -PC_2)) +
  geom_point(data = admix_pc_tmp, color = 'grey', shape = 4, size = 1.0) +
  geom_point(data = continental_pc, aes(fill = Subpopulation), shape = 21) +
  geom_point(shape = 4, color = '#440154', size = 1.0) +
  scale_fill_manual(values = c('#21918c', '#5ec962', '#3b528b', '#fde725')) +
  facet_wrap(~ADMIX) +
  labs(
    title = 'Simulated data'
  ) +
  plot_theme


ps_23 <- admix_pc %>% 
  ggplot(aes(x = -PC_2, y = PC_3)) +
  geom_point(data = admix_pc_tmp, color = 'grey', shape = 4, size = 1.0) +
  geom_point(data = continental_pc, aes(fill = Subpopulation), shape = 21) +
  geom_point(shape = 4, color = '#440154', size = 1.0) +
  scale_fill_manual(values = c('#21918c', '#5ec962', '#3b528b', '#fde725')) +
  facet_wrap(~ADMIX) +
  labs(
    title = 'Simulated data'
  ) +
  plot_theme

pdf('plots/PCAs.pdf', height = 5, width = 5)
grid.arrange(pr_12, pr_23, ps_12, ps_23, ncol = 2)
dev.off()
