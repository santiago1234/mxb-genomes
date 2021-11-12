library(tidyverse)
library(scales)


source("../ploting.R")


newborn_id <- "id3126_NA19771_X_NA19648"
mom_id <- "NA19771"
dad_id <- "NA19648"
parents_path_tract <- "../../210514-ProcessXGMixOutForTracts/data/3-pops/tracts/MXL/"

newborn <- load_individual("./", newborn_id) %>% 
  mutate(relation = "newborn")

mom <- load_individual(parents_path_tract, mom_id) %>% 
  mutate(relation = "parent 1")

dad <- load_individual(parents_path_tract, dad_id) %>% 
  mutate(relation = "parent 2")


tracts <- bind_rows(mom, dad, newborn)


tracts %>% 
  ggplot(aes(x = spos, fill = Ancestry)) +
  geom_rect(
    aes(xmin = spos, xmax=epos, ymin = ymin_haplo, ymax = ymax_haplo)
  ) +
  geom_hline(yintercept = 0.5, color ="white", size = 0.1) +
  facet_grid(chrn ~ relation,switch = "y") +
  scale_x_continuous(
    expand = c(0,0),
    breaks = c(10, 50, 100, 200) * 1e6,
    labels = label_number(scale = 1e-6, suffix = "Mb")
  ) +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    strip.text.y = element_text(angle = 0),
    panel.background = element_blank(),
    legend.position = c(.9, .2)
  ) +
  labs(
    y = "Chromosome"
  ) +
  scale_fill_brewer(type="qual", palette = "Accent") +
  #scale_fill_viridis_d(option = "C") +
  labs(
    x = "Position relative to chromosome start"
  )
ggsave("plot.pdf", height = 5, width = 8)
