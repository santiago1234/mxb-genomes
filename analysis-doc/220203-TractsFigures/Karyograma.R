library(tidyverse)
library(scales)
library(ggthemes)

theme_set(theme_tufte(base_family = 'Helvetica'))

cols_n <- c('chrn', 'spos', 'epos', 'Ancestry', 's_gp', 'e_gp')


bed_a <- read_tsv(
  '../210514-ProcessXGMixOutForTracts/data/3-pops/tracts/MXL/NA19651_anc_A_cM.bed',
  col_names = cols_n
) %>% 
  mutate(
    Haplotype = "A"
  )


bed_b <- read_tsv(
  '../210514-ProcessXGMixOutForTracts/data/3-pops/tracts/MXL/NA19651_anc_B_cM.bed',,
  col_names = cols_n
) %>% 
  mutate(
    Haplotype = "B"
  )

bed <- bind_rows(bed_a, bed_b)
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


# chromosome size file ----------------------------------------------------

chr_sizes <- read_tsv("~/Research/proyectos/22-01-ManilaGaleon-Juan/maas-MDS-analysis/chrom-sizes", col_names = c('chrn', 'epos')) %>% 
  mutate(spos = 1)

haplo <- tibble(Haplotype = c('A', 'B'))

chr_sizes <- 
  crossing(chr_sizes, haplo) %>% 
  mutate(
    ymin_haplo = map_dbl(Haplotype, limit_y_min),
    ymax_haplo = map_dbl(Haplotype, limit_y_max),
  )



# the plot ----------------------------------------------------------------

anc <- tribble(
  ~ANC_N, ~Ancestry,
  'European', 'European',
  'Native American', 'Native American',
  'African', 'African',
  'East Asian-Melanesian', 'East Asian-Melanesian'
)


bed2 <- inner_join(bed, anc)

bed %>% 
  ggplot(aes(x = spos)) +
  geom_rect(
    data = chr_sizes,
    aes(xmin = spos, xmax=epos, ymin = ymin_haplo, ymax = ymax_haplo),
    fill = NA,
    color = 'black',
    size = 1/7
  ) +
  geom_rect(
    aes(xmin = spos, xmax=epos, ymin = ymin_haplo, ymax = ymax_haplo, fill = Ancestry)
  ) +
  geom_hline(yintercept = 0.5, color ="white", size = 0.1) +
  facet_grid(chrn ~ ., switch = "y") +
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
    legend.position = c(.8, .2),
    panel.spacing = unit(0.15, "lines")
  ) +
  labs(
    y = "Chromosome"
  ) +
  labs(
    x = "Position relative to chromosome start"
  ) +
  scale_fill_manual(
    values = c('#fde725', '#5ec962', '#3b528b')
  )
ggsave('plots/karyogram-by-tracts.pdf', height = 3, width = 2.5)
