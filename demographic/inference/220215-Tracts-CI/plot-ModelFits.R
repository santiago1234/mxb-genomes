library(tidyverse)
library(scales)

fits <- read_csv('results/fits-data.csv')

model_order <- c("ppx_xxp", "ppx_xxp_pxx", "ccx_xxp", "ppx_ccx_xxp")

fits$cM <- fits$bins * 100
fits$mdl <- factor(fits$mdl, levels = model_order)


fits %>% 
  sample_frac(size = 1) %>% 
  ggplot(aes(x = cM, y = pred, color = Ancestry, fill = Ancestry)) +
  geom_ribbon(aes(ymin = ci_l, ymax = ci_u), alpha = 0.3, color = "black", show.legend=FALSE, size= 0.2) +
  geom_line(show.legend=FALSE) +
  geom_line(show.legend=FALSE, color = "black", size = 0.2) +
  geom_point(aes(y=dat),shape = 21, size = 2, alpha = .9, color = "black") +
  scale_y_log10(
    oob = scales::squish_infinite,
    expand = c(0, 0),
    breaks = 10^(-1:5),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_x_continuous(expand = c(.03, 0), breaks = c(0, 100, 200)) +
  coord_cartesian(ylim = c(1, 14000), xlim = c(-1, 250)) +
  facet_grid(mdl~Population) +
  scale_color_manual(
    values = c('#fde725', '#5ec962', '#3b528b')
  ) +
  scale_fill_manual(
    values = c('#fde725', '#5ec962', '#3b528b')
  ) +
  theme_bw() +
  labs(
    y = 'Relative Frequency',
    x = 'Tract Length (cM)'
  )
ggsave('plots/fs-fits-mdls.pdf', height = 6, width = 6.5)
ggsave('plots/fs-fits-mdls.svg', height = 6, width = 6.5)



# 4 pops model ------------------------------------------------------------


fits <- read_csv('results/fits-data-4pops-MXL.csv')

fits$cM <- fits$bins * 100

fits %>% 
  sample_frac(size = 1) %>% 
  ggplot(aes(x = cM, y = pred, color = Ancestry, fill = Ancestry)) +
  geom_ribbon(aes(ymin = ci_l, ymax = ci_u), alpha = 0.3, color = "black", show.legend=FALSE, size= 0.2) +
  geom_line(show.legend=FALSE) +
  geom_line(show.legend=FALSE, color = "black", size = 0.2) +
  geom_point(aes(y=dat),shape = 21, size = 2, alpha = .9, color = "black") +
  scale_y_log10(
    oob = scales::squish_infinite,
    expand = c(0, 0),
    breaks = 10^(-1:5),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_x_continuous(expand = c(.03, 0), breaks = c(0, 100, 200)) +
  coord_cartesian(ylim = c(1, 14000), xlim = c(-1, 250)) +
  facet_grid(~mdl) +
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
ggsave('plots/4pops-fs-fits-mdls.pdf', height = 2, width = 5)
ggsave('plots/4pops-fs-fits-mdls.svg', height = 2, width = 5)
 