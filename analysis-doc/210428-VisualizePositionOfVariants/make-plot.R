library(tidyverse)
library(scales)
library(gridExtra)


var_d <- read_csv("data/densities-all.csv")
var_d$chromosome <- as.integer(var_d$chromosome)


is_outlier_upper <- function(a_vector) {
  
  # Check for outliers, I just check the outliers that
  # are in the upper range
  Q <- quantile(a_vector, probs=c(.25, .75), na.rm = FALSE)
  
  iqr <- IQR(a_vector)
  
  up <-  Q[2] + 1.5 * iqr # Upper Range  
  
  a_vector > up

}


# remove upper outliers

var_d <- 
  var_d %>% 
  group_by(
    chromosome, DataSource, mask_filter
  ) %>% 
  mutate(
    upper_out = is_outlier_upper(variant_density)
  ) %>% 
  filter(!upper_out) %>% 
  ungroup()


# scale the data

var_d <- 
  var_d %>% 
  group_by(
    chromosome, DataSource, mask_filter
  ) %>% 
  mutate(
    vd_scaled = rescale(variant_density, to=c(0, 1))
)


chromosomes_to_plt <- c(1, 5, 15, 22)
x <- filter(var_d, chromosome %in% chromosomes_to_plt)

x <- 
  x %>% 
  mutate(
    id = paste0(mask_filter, " masking |", DataSource)
  )



# plots -------------------------------------------------------------------


low_density <- 0.1 / 3
x$chromosome <- paste0("chr", x$chromosome)
x$chromosome <- factor(x$chromosome, levels = paste0("chr", chromosomes_to_plt))

p_before <- x %>% 
  filter(mask_filter == "before") %>% 
  ggplot(aes(x = position_bp, y = vd_scaled)) +
  geom_point(aes(color = vd_scaled < low_density), size = 1/2) +
  facet_grid(DataSource~chromosome, scales = "free", space = "free_x") +
  scale_x_continuous(
    expand = c(0,0),
    labels = label_number(scale = 1e-6, suffix = "Mb")
  ) +
  scale_color_manual(name = "Low density",values = c("grey", "orange")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 40, hjust=1),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    text = element_text(size=20)
  ) +
  labs(
    title = "Variant density, BEFORE mask filter",
    y = "Density",
    x = "Genomic position"
  )


p_after <- x %>% 
  filter(mask_filter == "after") %>% 
  ggplot(aes(x = position_bp, y = vd_scaled)) +
  geom_point(aes(color = vd_scaled < low_density), size = 1/2) +
  facet_grid(DataSource~chromosome, scales = "free", space = "free_x") +
  scale_x_continuous(
    expand = c(0,0),
    labels = label_number(scale = 1e-6, suffix = "Mb")
  ) +
  scale_color_manual(name = "Low density",values = c("grey", "orange")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 40, hjust=1),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    text = element_text(size=20)
  ) +
  labs(
    title = "Variant density, AFTER mask filter",
    y = "Density",
    x = "Genomic position"
  )


png("plots/variant_density.png", width=1250, height=700)
grid.arrange(p_before, p_after)
dev.off()
