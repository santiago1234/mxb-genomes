library(tidyverse)
library(ggthemes)
theme_set(theme_tufte(base_family = "Helvetica"))

track_length <- read_csv("data/tractlen-HG01893.csv")

track_length %>% 
  ggplot(aes(x = bins_cm, y = x, color = assignment)) +
  geom_smooth(
    size = 0.5, se = FALSE
  ) +
  geom_point() +
  geom_point(shape = 1, color = "black") +
  scale_color_viridis_d(option = "C") +
  geom_rangeframe(color = "black") +
  scale_y_log10() +
  labs(
    x = "Tract Length (cM)",
    y = "Number of tracts",
    subtitle = "Tract length distribution"
  )


ggsave("plots/tractlen-HG01893.pdf", height = 4, width = 6)
ggsave("plots/tractlen-HG01893.png", height = 4, width = 6)
