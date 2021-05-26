library(tidyverse)
library(scales)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))

# define supplied parameters ----------------------------------------------

# input_bed <- "data/3-pops/bed/allchrn-HG01893.tsv"
# output_plot <- "plots/karyo-HG01893.png"
# bed <- read_tsv(input_bed)

input_bed <- snakemake@input[[1]]
output_plot <- snakemake@output[[1]]
  
individual <- basename(input_bed) %>%
  str_replace(".bed", "")

bed <- read_tsv(input_bed)

# aling chromosomes so the 1st position is zero ---------------------------

# This creates a bug when drawing the tracts
# bed <- 
#   bed %>% 
#   group_by(chrn) %>% 
#   mutate(
#     spos = spos - min(spos),
#     epos = epos - min(spos)
#   ) %>% 
#   ungroup() %>% 
#   select(-c(sgpos, egpos))

# draw a indicator variable for Hapl  otypes limits -------------------------

limit_y_max <- function(haplo) {
  if (haplo == "A") {
    return(Inf)
  }
  if (haplo == "B") {
    return(0.5) #0.49
  }
}

limit_y_min <- function(haplo) {
  if (haplo == "A") {
    return(0.5)
  }
  if (haplo == "B") {
    return(-Inf)
  }
}

bed <- bed %>% 
  mutate(
    ymin_haplo = map_dbl(Haplotype, limit_y_min),
    ymax_haplo = map_dbl(Haplotype, limit_y_max),
  )

bed %>% 
  ggplot(aes(x = spos, fill = Ancestry)) +
  geom_rect(
    aes(xmin = spos, xmax=epos, ymin = ymin_haplo, ymax = ymax_haplo)
  ) +
  geom_hline(yintercept = 0.5, color ="white", size = 0.1) +
  facet_grid(chrn ~.,switch = "y") +
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

ggsave(output_plot, height = 5, width = 4)
ggsave(
  str_replace(output_plot, ".png", ".pdf"),
  height = 5, width = 4)



# code to plot start position points --------------------------------------

# bed$y <- if_else(bed$ymin_haplo == 0.5, true = 0, false = 1)
# 
# bed %>% 
#   ggplot(aes(x = spos, color = Ancestry)) +
#   geom_point(aes(y = y), position = "jitter", size = 0.3) +
#   geom_hline(yintercept = 0.5, color ="white", size = 0.1) +
#   facet_grid(chrn ~.,switch = "y") +
#   scale_x_continuous(
#     expand = c(0,0),
#     labels = label_number(scale = 1e-6, suffix = "Mb")
#   ) +
#   theme(
#     axis.ticks.y = element_blank(),
#     axis.text.y = element_blank(),
#     strip.text.y = element_text(angle = 0),
#     panel.background = element_blank(),
#     legend.position = c(.8, .2)
#   ) +
#   labs(
#     y = "Chromosome"
#   ) +
#   scale_color_viridis_d(option = "C") +
#   labs(
#     x = "Position relative to chromosome start",
#     title = individual
#   )
