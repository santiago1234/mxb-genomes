library(tidyverse)
library(scales)
library(ggthemes)

theme_set(theme_tufte(base_family = 'Helvetica'))

# plot fits for best mdls -------------------------------------------------

fits3pops <-
  read_csv('results/fits-data.csv') %>% 
  filter(
    mdl == 'ppx_ccx_xxp',
    Population != 'MXL',
    bootstrap == 0
  )

fits4pops <- 
  read_csv('results/fits-data-4pops-MXL.csv') %>% 
  filter(
    mdl == 'ppxx_ccxx_xxpp',
    bootstrap == 0
  ) %>% 
  mutate(
    Population = 'MXL'
  )

fits <- 
  bind_rows(fits3pops, fits4pops)

fits$cM <- fits$bins * 100

fits %>% 
  ggplot(aes(x = cM, y = pred, color = Ancestry, fill = Ancestry)) +
  geom_ribbon(aes(ymin = ci_l, ymax = ci_u), alpha = 0.3, color = "black", show.legend=FALSE, size= 0.2) +
  geom_line(show.legend=FALSE) +
  geom_line(show.legend=FALSE, color = "black", size = 0.2) +
  geom_point(aes(y=dat),shape = 21, size = 2, alpha = .7, color = "black") +
  scale_y_log10(
    oob = scales::squish_infinite,
    expand = c(0, 0),
    breaks = 10^(-1:5),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_x_continuous(expand = c(.03, 0), breaks = c(0, 100, 200)) +
  coord_cartesian(ylim = c(1, 14000), xlim = c(-1, 250)) +
  facet_grid(.~Population) +
  scale_color_manual(
    values = c('#fde725', '#21918c', '#5ec962', '#3b528b')
  ) +
  scale_fill_manual(
    values = c('#fde725', '#21918c', '#5ec962', '#3b528b')
  ) +
  theme_bw() +
  labs(
    y = 'Relative Frequency',
    x = 'Tract Length (cM)'
  )
ggsave('plots/mf-best-fit.pdf', height = 2, width = 6)


# draw ancestry fracts ----------------------------------------------------

ancfracs <- read_csv('results/anc-props-best-mdls.csv')

ancfracs <- ancfracs %>% 
  replace_na(list(EAS = 0)) %>% 
  pivot_longer(cols = c(EUR, NAT, AFR, EAS), names_to = 'Ancestry', values_to = 'Frac')


ancfracs %>% 
  ggplot(aes(x = -ga, y = Frac, fill = Ancestry)) +
  geom_bar(stat = "identity", width = 1.1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(~Population) +
  scale_fill_manual(
    values = c('#fde725', '#21918c', '#5ec962', '#3b528b')
  ) +
  labs(
    x = 'Generations ago',
    y = 'Ancestry proportion'
  )
ggsave('plots/mf-anctry-props-over-time.pdf', height = 1.3, width = 6)


# draw fist pulse ---------------------------------------------------------



ancfracs %>% 
  group_by(Population) %>% 
  filter(ga == max(ga)) %>% 
  filter(!Frac == 0) %>% 
  ggplot(aes(x = "", y = Frac, fill = Ancestry)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start=0) +
  facet_grid(~Population) +
  scale_fill_manual(
    values = c('#fde725', '#5ec962', '#3b528b')
  ) +
  theme_void()
ggsave('plots/mf-1st-pulses-pies.pdf', height = 1, width = 3)

### continuous pulses


tibble(x
  
)


# get cont mig ------------------------------------------------------------

ancfracs %>% 
  group_by(Population) %>% 
  filter(!ga == max(ga)) %>% 
  filter(ga == max(ga)) %>% 
  filter(!Frac == 0) %>% 
  mutate(
  )

cont_pulses <- 
  tribble(
  ~Population, ~Ancestry, ~Frac,
  'MXL', 'NAT', 0.0573,
  'MXL', 'EUR', 0.0502,
  
  'CLM', 'NAT',  0.012,
  'CLM', 'EUR', 0.0617,
  
  'PEL', 'NAT',  0.118,
  'PEL', 'EUR', 0.014,
  
  'PUR', 'AFR',  0.005,
  'PUR', 'EUR', 0.03
  
) %>% 
  group_by(Population) %>% 
  mutate(
    Frac = Frac / sum(Frac)
  )





cont_pulses %>% 
  ggplot(aes(x = "", y = Frac, fill = Ancestry)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start=0) +
  facet_grid(~Population) +
  scale_fill_manual(
    values = c('#fde725', '#5ec962', '#3b528b')
  ) +
  theme_void()
ggsave('plots/mf-continous-pulses-pies.pdf', height = 1, width = 3)
