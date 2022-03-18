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

order_mdls <- c("ppx_xxp", "ppx_xxp_pxx", "ccx_xxp", "ppx_ccx_xxp", 'ppp_pxp') %>% rev()
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



# compute Bayesian Information Criteria -----------------------------------

n_points <- 
  read_tsv('../../../resources/1TGP-samples-meta-data/integrated_call_samples_v3.20130502.ALL.panel') %>% 
  select(sample, pop) %>% 
  filter(pop %in% unique(infos$Population)) %>% 
  count(pop) %>% 
  rename(Population = pop)

k_params <- tribble(
  ~mdl, ~k,
  'ppx_xxp', 2,
  'ppx_xxp_pxx', 4,
  'ccx_xxp', 5,
  'ppx_ccx_xxp', 6,
  'ppp_pxp', 4
)


bic_d <- crossing(n_points, k_params)


## Compute a CI for the lilelihood

bic_d <- 
  infos %>% 
  group_by(Population, mdl) %>% 
  summarise(
    mean_lkl = mean(loglik),
    sd_lkl = sd(loglik),
  ) %>% 
  ungroup() %>% 
  mutate(
    ci_h_lkl = mean_lkl + 2 * sd_lkl,
    ci_l_lkl = mean_lkl - 2 * sd_lkl
  ) %>% 
  inner_join(bic_d)


bic <- function(k, n, log_likelihood) {
  
  k * log(n) - 2 * log_likelihood
}


bic_d <- bic_d %>% 
  mutate(
    bic = bic(k, n, mean_lkl),
    bic_ci_h = bic(k, n, ci_h_lkl), 
    bic_ci_l = bic(k, n, ci_l_lkl), 
  )

order_mdls <- c("ppx_xxp", "ppx_xxp_pxx", "ccx_xxp", "ppx_ccx_xxp", 'ppp_pxp') %>% rev()
bic_d$mdl <- factor(bic_d$mdl, levels = order_mdls)

bic_d <- bic_d %>% 
  group_by(Population) %>% 
  mutate(
    Best = bic == min(bic)
  )


bic_d %>% 
  ggplot(aes(x = bic, y = mdl, color = Best)) +
  geom_point() +
  geom_errorbarh(aes(xmin = bic_ci_l, xmax = bic_ci_h), height = 0) +
  facet_grid(~Population) +
  coord_cartesian(xlim = c(500, 1500)) +
  scale_color_manual(values = c('grey', 'steelblue')) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = 'none'
  ) +
  labs(
    x = 'BIC'
  )

ggsave('plots/fs-BIC-mdls.pdf', height = 1.5, width = 6)


# BIC done ----------------------------------------------------------------



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
