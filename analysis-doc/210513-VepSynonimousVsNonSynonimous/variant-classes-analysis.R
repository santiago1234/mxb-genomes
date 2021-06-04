library(tidyverse)
library(scales)

# genaral preprocessing ---------------------------------------------------
# collapse the data

# load the frame will all the data
read_csv("results/sfs-impact_HIGH.csv")

load_cat <- function(filepath) {
  # loads the given category
  # get the category name
  cat <- basename(filepath) %>% 
    str_replace_all("(sfs-|.csv)", "")
  
  read_csv(filepath) %>% 
    mutate(
      Variant_Category = cat
    )
}

all_cats <- 
  list.files("results", pattern = "sfs-", full.names = T) %>% 
  map_df(load_cat)

# Drop the corners
all_cats <- 
  all_cats %>% 
  filter(n != min(n), n != max(n))

# population to focus on
populations <- c("YRI", "IBS", "MXL", "MXB")

all_cats <- 
  all_cats %>% 
  filter(Population %in% populations) %>% 
  mutate(
    Population = factor(Population, levels = populations)
  )

# normalize spectrum

all_cats <- 
  all_cats %>% 
  group_by(Population, Variant_Category) %>% 
  mutate(
    spectrum_norm = Freq / sum(Freq)
  ) %>% 
  ungroup()

# general categories ---------------------------------------------------------

general_pattern_cats <- c("synonymous", "missense", "loss_of_function", "intergenic", "intron", "utr_5_or_3")

general <- all_cats %>% 
  filter(Variant_Category %in% general_pattern_cats)


# The SFS for populations and Variant Categories
general %>% 
  ggplot(aes(x = n, y = Freq)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Variant_Category ~ Population, scale = "free_y") +
  scale_y_sqrt() +
  labs(
    title = "SFS",
    y = "Freq",
    caption = "y scale was transformed with sqrt"
  )

ggsave("plots/SFS-all.pdf", height = 7, width = 8)
ggsave("plots/SFS-all.png", height = 7, width = 8)

# COOL ANALYSIS -----------------------------------------------------------


cross_compare <- function(df_long, x_var, y_var) {
  # Args:
  # df_long: 
  cruza <- 
    df_long %>% 
    filter(Variant_Category %in% c(x_var, y_var)) %>% 
    pivot_wider(
      values_from = spectrum_norm,
      names_from = Variant_Category,
      id_cols = c(n, Population)
    )
  
  order_vars <- c("n", "Population", x_var, y_var)
  cruza <- cruza[, order_vars]
  
  colnames(cruza) <- c("n", "Population", "x", "y")
  cruza %>% 
    mutate(
      xv = paste0("xv: ", x_var),
      yv = paste0("yv: ", y_var)
    )
  
}


my_vars <- unique(general$Variant_Category)

d_cruza <- 
  expand_grid(v1 = my_vars, v2 = my_vars) %>% 
  mutate(
    fr = map2(v1, v2, function(x, y) cross_compare(general, x, y))
  ) %>% 
  unnest(fr) %>% 
  select(-c(v1, v2))


d_cruza %>% 
  ggplot(aes(x = x, y = y, color = n)) +
  geom_abline(size = 0.2) +
  geom_point(size = 0.5) +
  scale_x_log10(breaks = breaks_log(n = 4, base = 10)) +
  scale_y_log10() +
  scale_color_viridis_c(breaks = c(1, 20, 40, 60, 79), option = "B") + #B,C
  facet_grid(yv~xv) +
  theme(panel.grid = element_blank()) +
  labs(
    x = "log10  Proportion",
    color = "Derived\nAllele\nFreq\n",
    y = "log10  Proportion",
    title = "SFS Variant Categories Comparison",
    subtitle = "xv: x-var, yv: y-var"
  )

ggsave("plots/pairwise-comparison-Categories.pdf", height = 8, width = 10)
ggsave("plots/pairwise-comparison-Categories.png", height = 8, width = 10)

# focus on synonymous variants

general %>% 
  filter(Variant_Category %in% c("synonymous", "intergenic", "missense")) %>% 
  ggplot(aes(x = n, y = spectrum_norm)) +
  geom_line(aes(linetype = Variant_Category)) +
  scale_y_log10() +
  facet_wrap(~Population) +
  labs(
    y = "log10  Proportion",
    x = "Derived allele frequency",
    linetype = "Variant Class: ",
    title = "SFS for different variant classes"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("plots/SynVsIntergenic.pdf", height = 4, width = 5)
ggsave("plots/SynVsIntergenic.png", height = 4, width = 5)

  
  
  
