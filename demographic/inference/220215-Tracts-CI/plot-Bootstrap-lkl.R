library(tidyverse)
library(ggthemes)
library(ggridges)

load_info_file <- function(fp) {
  read_tsv(fp) %>% 
    select(loglik:BOOT) %>% 
    distinct()
  
}

infos <- 
  list.files("results/inference/", pattern = '*info.tsv', full.names = T) %>% 
  map_df(load_info_file)

order_mdls <- c("ppx_xxp", "ppx_xxp_pxx", "ccx_xxp", "ppx_ccx_xxp") %>% rev()
infos$mdl <- factor(infos$mdl, levels = order_mdls)

## get means to show in plot

# draw this line for comparison
best_lkl <- 
  infos %>% 
  filter(loglik > -1000) %>% 
  group_by(Population, mdl) %>% 
  summarise(mlkl = mean(loglik)) %>% 
  ungroup() %>% 
  group_by(Population) %>% 
  arrange(-mlkl) %>% 
  slice(1:1) %>% 
  mutate(fw = mdl) %>% 
  select(-mdl)



infos %>% 
  mutate(fw = mdl) %>% 
  #  filter(Population == 'MXL') %>% 
  filter(loglik > -1000) %>% 
  ggplot(aes(x = loglik, y = mdl)) + 
  geom_density_ridges(size = 0.4, fill = 'steelblue', alpha = 0.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_wrap(. ~ Population, scales = 'free_x', nrow = 1) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave('plots/fs-lkl-mdls.pdf', height = 1.5, width = 6)
ggsave('plots/fs-lkl-mdls.svg', height = 1.5, width = 6)



# 4pops -------------------------------------------------------------------

infos <- 
  list.files("results/inference-MXL-4pops/", pattern = '*info.tsv', full.names = T) %>% 
  map_df(load_info_file)

infos %>% 
  mutate(fw = mdl) %>% 
  #  filter(Population == 'MXL') %>% 
  filter(loglik > -1000) %>% 
  ggplot(aes(x = loglik, y = mdl)) + 
  geom_density_ridges(size = 0.4, fill = 'steelblue', alpha = 0.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_wrap(. ~ Population, scales = 'free_x', nrow = 1) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave('plots/4pops-fs-lkl-mdls.pdf', height = 1.5, width = 3)
ggsave('plots/4pops-fs-lkl-mdls.svg', height = 1.5, width = 3)
