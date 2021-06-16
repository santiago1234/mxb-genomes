library(tidyverse)
library(scales)

tractsfiles <- list.files("results", pattern = "res.csv", full.names = T)

tracts <- map_df(tractsfiles, read_csv)

tracts <- 
  tracts %>% 
  mutate(
    id = paste0(population, "\n", parameters)
  )

tracts$cM <- tracts$bins * 100

tracts %>% 
  ggplot(aes(x = cM, y = pred, color = Ancestry, fill = Ancestry)) +
  geom_ribbon(aes(ymin = ci_l, ymax = ci_u), alpha = 0.3, color = NA, show.legend=FALSE) +
  geom_line(show.legend=FALSE) +
  geom_point(aes(y=dat),shape = 21, size = 2, alpha = .9, color = "black") +
  scale_y_log10(
    oob = scales::squish_infinite,
    expand = c(0, 0),
    breaks = 10^(-1:5),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_x_continuous() +
  coord_cartesian(ylim = c(1, 14000), xlim = c(-1, 250)) +
  scale_fill_viridis_d(option = "C") +
  scale_color_viridis_d(option = "C") +
  labs(
    x = "Tract Length (cM)",
    y = "Relative frequency"
  ) +
  facet_grid(~id) +
  theme(legend.position = "bottom")

ggsave("plots/taino.pdf", height = 3, width = 9)
ggsave("plots/taino.png", height = 3, width = 9)
