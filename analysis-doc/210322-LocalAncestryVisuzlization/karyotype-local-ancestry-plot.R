library(tidyverse)
library(scales)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))
bed <- read_tsv("data/bed/allchrn-HG01893.tsv")
individual <- "HG01893"

# aling chromosomes so the 1st position is zero ---------------------------

bed <- 
  bed %>% 
  group_by(chm) %>% 
  mutate(
    spos = spos - min(spos),
    epos = epos - min(spos)
  ) %>% 
  ungroup() %>% 
  select(-c(sgpos, egpos))

# draw a indicator variable for Hapl  otypes limits -------------------------

limit_y_max <- function(haplo) {
  if (haplo == 1) {
    return(Inf)
  }
  if (haplo == 0) {
    return(0.5) #0.49
  }
}

limit_y_min <- function(haplo) {
  if (haplo == 1) {
    return(0.5)
  }
  if (haplo == 0) {
    return(-Inf)
  }
}

bed <- bed %>% 
  mutate(
    ymin_haplo = map_dbl(Haplotype, limit_y_min),
    ymax_haplo = map_dbl(Haplotype, limit_y_max),
  )

bed %>% 
  ggplot(aes(x = spos, fill = assignment)) +
  geom_rect(
    aes(xmin = spos, xmax=epos, ymin = ymin_haplo, ymax = ymax_haplo)
  ) +
  geom_hline(yintercept = 0.5, color ="white", size = 0.1) +
  facet_grid(chm ~.,switch = "y") +
  scale_x_continuous(
    expand = c(0,0),
    labels = label_number(scale = 1e-6, suffix = "Mb")
  ) +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    strip.text.y = element_text(angle = 0),
    panel.background = element_blank(),
    legend.position = c(.8, .2)
  ) +
  labs(
    y = "Chromosome"
  ) +
  scale_fill_viridis_d(option = "C") +
  labs(
    x = "Position relative to chromosome start",
    title = individual
  )
ggsave("test.pdf", height = 5, width = 5)
ggsave("test.png", height = 5, width = 5)


# this function is to visualy check that points are placed correct --------


bed %>% 
  ggplot(aes(x = spos, fill = assignment)) +
  geom_point(aes(y = Haplotype),shape = 21, position = "jitter") +
  geom_hline(yintercept = 0.5, color ="white", size = 1) +
  facet_grid(chm ~., switch = "y") +
  scale_x_continuous(labels = label_number(scale = 1e-6, suffix = "Mb")) +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.background = element_blank(),
    legend.position = c(1, 0)
  ) +
  labs(
    y = NULL
  ) +
  scale_fill_viridis_d(option = "C") +
  labs(
    x = "Position relative to chromosome start",
    subtitle = individual
  )
ggsave("test-points.pdf", height = 3, width = 8)
