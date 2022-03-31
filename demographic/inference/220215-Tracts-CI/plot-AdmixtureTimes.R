library(tidyverse)


get_time_pars <- function(mdl, parameters) {
  switch(
    mdl,
    ppx_xxp = {parameters},
    ppx_xxp_pxx = {parameters[c(2, 3)]},
    ccx_xxp = {parameters[c(3, 4)]},
    ppx_ccx_xxp = {parameters[c(4, 6)]},
    ppc = {parameters[c(3)]}
  )
  
}


load_time_pars <- function(fp) {
  patern <- '(ppx_ccx_xxp|ppx_xxp_pxx|ccx_xxp|ppx_xxp|ppc)'
  
  pop <- str_extract(fp, '(CLM|PEL|PUR|MXL)')
  
  boostrap <- str_extract(fp, 'boot\\d\\d?') %>% 
    str_replace('boot', '') %>% 
    as.numeric()
  mdl_name <- str_extract(fp, patern)
  
  parameters <- 
    read_lines(fp) %>% 
    str_split('\t') %>% 
    .[[1]] %>% 
    as.numeric() %>% 
    (function(x) x*100)()
  
  
  pars <- get_time_pars(mdl_name, parameters)
  
  if (length(pars) == 1) {
    Time = c('T1')
  } else {
    Time = c('T1', 'T2')
  }
  
  tibble(Time_pars = pars, Time = Time) %>% 
    mutate(
      mdl = mdl_name,
      boot = boostrap,
      Population = pop
    )
  
  
}


params_files <- list.files(path = 'results/inference/', pattern = '*pars', full.names = T)


d <- map_df(params_files, load_time_pars) %>% 
  group_by(mdl, Time, Population)

d <- d %>% 
  mutate(Time_pars = -Time_pars) %>% 
  summarise(
    mean_par = mean(Time_pars),
    std = sd(Time_pars)
  ) %>% 
  mutate(
    ci_l = mean_par - 2*std,
    ci_u = mean_par + 2*std,
  ) %>% 
  ungroup()

order_mdls <- c("ppx_xxp", "ppx_xxp_pxx", "ccx_xxp", "ppx_ccx_xxp", 'ppc') %>% rev()
d$mdl <- factor(d$mdl, levels = order_mdls)

d %>% 
  ggplot(aes(y = mdl, color = Time)) +
  geom_point(aes(x = mean_par)) +
  geom_errorbar(aes(xmin = ci_l, xmax=ci_u), width = 0.1) +
  scale_color_manual(values = c('#440154', '#21918c')) +
  facet_grid(~Population) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = 'bottom'
  ) +
  labs(
    x = 'Generations ago',
    y = 'Model'
  )
ggsave('plots/fs-TimePars-mdls.pdf', height = 2, width = 6)
ggsave('plots/fs-TimePars-mdls.svg', height = 2, width = 6)
