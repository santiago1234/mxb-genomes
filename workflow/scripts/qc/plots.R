library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggridges)
library(ggrepel)
library(ggthemes)
library(scales)
library(stringr)

theme_set(theme_tufte(base_family = "Helvetica"))

# snakemake params --------------------------------------------------------
# input parameters
var_per_genome <- snakemake@input$vars_per_genome
seqs_dep <- snakemake@input$seqs_deps


# output plots
vars_per_genome_plt <- snakemake@output$vars_per_genome_plt
depth_per_sample_plt <- snakemake@output$depth_per_sample_plt
depth_in_chr22_plt <- snakemake@output$depth_in_chr22_plt
miss_ind_plt <- snakemake@output$miss_ind_plt
miss_var_plt <- snakemake@output$miss_var_plt

# helper function to change the name to pdf
pdf_name <- function(x) str_replace(x, "png", "pdf")
# plot n vars per genome histogram ----------------------------------------

nvars <- read_csv(var_per_genome)

## identify outliers

lowerq <- quantile(nvars$n_variants)[2]
upperq <- quantile(nvars$n_variants)[4]
iqr <- upperq - lowerq

threshold_upper = (iqr * 1.5) + upperq
threshold_lower = lowerq - (iqr * 1.5)

outlier_points <- nvars %>% 
  filter(n_variants > threshold_upper | n_variants < threshold_lower)


nvars %>% 
  ggplot(aes(x = n_variants)) +
  geom_density(fill = "grey90", color = NA) +
  geom_point(
    aes(y = 0, fill = Individual %in% outlier_points$Individual)
    ,size = 3, shape = 21, color = "white"
  ) +
  geom_text_repel(data = outlier_points, aes(y = 0, label = Individual), color = "grey60") +
  scale_color_manual(values = c("orange", "grey60")) +
  scale_fill_manual(values = c("orange", "grey60")) +
  scale_x_continuous(
    labels = unit_format(unit = "M", scale = 1e-6, accuracy = 0.001)
  )  + 
  geom_rangeframe() +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(
    x = "# variants per genome",
    y = "Density",
    subtitle = "Outliers are highlighted in grey color",
    title = "Distribution of variants per genome"
  )

ggsave(filename = vars_per_genome_plt, height = 3, width = 5)
ggsave(filename =pdf_name(vars_per_genome_plt), 
       height = 3, width = 5)

# depth plot --------------------------------------------------------------

sdepth <- read_csv(seqs_dep)

mean_depth <- mean(sdepth$DP, na.rm = T)

sdepth %>% 
  filter(DP < 100) %>% 
  ggplot(aes(x = DP, y = reorder(Individual, DP, mean))) +
  geom_density_ridges(
    color = "white",
    fill = "grey45"
  ) +
  geom_vline(xintercept = mean_depth, size = .1) +
  theme(
    axis.text.y = element_text(size = 5, color = "grey"),
    plot.caption = element_text(colour = "grey")
  ) +
  labs(
    y = NULL,
    x = "Sequencing depth",
    caption = "points with depth > 100 were removed for\nvisualization"
  )

ggsave(depth_per_sample_plt, height = 7, width = 3)
ggsave(pdf_name(depth_per_sample_plt), height = 7, width = 3)


# depth in chromosome region ----------------------------------------------

sdepth_22 <- 
  sdepth %>% 
  filter(chrn == 22)


sdepth_22 <- 
  sdepth_22 %>% 
  select(-chrn) %>% 
  separate(varid, sep = "_", into = c("chrn", "position")) %>% 
  mutate(position = as.integer(position))

set.seed(42)
# show 10 individuals
my_individuals <- sample(sdepth_22$Individual, 10) %>% 
  c(., "MXB_15210")

sdepth_22 %>% 
  filter(Individual %in% my_individuals, DP < 100) %>% 
  ggplot(aes(x = position, y = DP)) +
  geom_point(shape = ".", alpha = .1, color = "steelblue") +
  geom_rangeframe(sides = "l") +
  geom_hline(yintercept = mean_depth, size = .25, color = "grey20", linetype = 2) +
  facet_grid(Individual~.) +
  theme(
    strip.text.y = element_text(size = 5),
    plot.caption = element_text(colour = "grey")
  ) +
  labs(
    title = "Sequencing depth in chromosome 22",
    x = "Position (bp)",
    y = "Sequencing depth",
    caption = "dashed lines are the mean depth across all samples"
  )
  
ggsave(depth_in_chr22_plt, height = 7, width = 4)
ggsave(pdf_name(depth_in_chr22_plt), height = 7, width = 4)


# missing data ------------------------------------------------------------
# If the depth is missing that means the the variant is also missing

# missing by individual

miss_ind <- 
  sdepth %>% 
  group_by(Individual) %>% 
  summarise(p_miss = mean(is.na(DP)))


miss_ind %>% 
  ggplot(aes(x = p_miss, y = reorder(Individual, p_miss))) +
  geom_point() +
  geom_rangeframe() +
  labs(
    x = "Missingness %",
    y = NULL,
    title = "Missing data by sample"
  )

ggsave(miss_ind_plt, height = 7, width = 5)
ggsave(pdf_name(miss_ind_plt), height = 7, width = 5)

# missing by variant

miss_var <- 
  sdepth %>% 
  group_by(varid) %>% 
  summarise(p_miss = mean(is.na(DP)))

miss_var %>% 
  ggplot(aes(x = p_miss)) +
  geom_histogram(bins = 30) +
  scale_y_log10() +
  labs(
    x = "Missingness %",
    y = "Number of variants",
    title = "Distribution of missing data by variant",
    caption = "only a fraction of all variants are shown"
  ) +
  theme(
    plot.caption = element_text(colour = "grey")
  )

ggsave(miss_var_plt, height = 3, width = 4)
ggsave(pdf_name(miss_var_plt), height = 3, width = 4)
