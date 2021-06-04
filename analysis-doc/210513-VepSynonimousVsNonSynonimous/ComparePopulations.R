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
populations <- c("YRI", "IBS", "MXL", "MXB", "CHB")

all_cats <- 
  all_cats %>% 
  filter(Population %in% populations) %>% 
  mutate(
    Population = factor(Population, levels = populations)
  )


all_cats <- 
  all_cats %>% 
  group_by(Population, Variant_Category) %>% 
  mutate(
    spectrum_norm = Freq / sum(Freq)
  ) %>% 
  ungroup()

# general pattern ---------------------------------------------------------

general_pattern_cats <- c("synonymous", "missense", "loss_of_function", "intergenic", "intron", "utr_5_or_3")

general <- all_cats %>% 
  filter(Variant_Category %in% general_pattern_cats)


cross_compare2 <- function(df_long, x_var, y_var) {
  # Args:
  # df_long: 
  cruza <- 
    df_long %>% 
    filter(Population %in% c(x_var, y_var)) %>% 
    #select(-Freq) %>% 
    pivot_wider(
      values_from = spectrum_norm,
      names_from = Population,
      id_cols = c(n, Variant_Category)
    )
  
  order_vars <- c("n", "Variant_Category", x_var, y_var)
  cruza <- cruza[, order_vars]
  
  colnames(cruza) <- c("n", "Variant_Category", "x", "y")
  cruza %>% 
    mutate(
      xv = paste0("xv: ", x_var),
      yv = paste0("yv: ", y_var)
    )
  
}

general$Population <- as.character(general$Population)
my_vars <- unique(general$Population)

d_cruza <- 
  expand_grid(v1 = my_vars, v2 = my_vars) %>% 
  mutate(
    fr = map2(v1, v2, function(x, y) cross_compare2(general, x, y))
  ) %>% 
  unnest(fr) %>% 
  select(-c(v1, v2))




d_cruza %>% 
  filter(Variant_Category == "missense") %>% 
  ggplot(aes(x = x, y = y, color = n)) +
  geom_abline(size = 0.2) +
  geom_point(size = 0.65, color = "black") +
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
    title = "SFS Missense variants",
    subtitle = "xv: x-var, yv: y-var"
  )
ggsave("plots/sfs-missense.pdf", height = 7, width = 8)
ggsave("plots/sfs-missense.png", height = 7, width = 8)


d_cruza %>% 
  filter(Variant_Category == "synonymous") %>% 
  ggplot(aes(x = x, y = y, color = n)) +
  geom_abline(size = 0.2) +
  geom_point(size = 0.65, color = "black") +
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
    title = "SFS Synonymous variants",
    subtitle = "xv: x-var, yv: y-var"
  )
ggsave("plots/sfs-synonymous.pdf", height = 7, width = 8)
ggsave("plots/sfs-synonymous.png", height = 7, width = 8)



d_cruza %>% 
  filter(Variant_Category == "loss_of_function") %>% 
  ggplot(aes(x = x, y = y, color = n)) +
  geom_abline(size = 0.2) +
  geom_point(size = 0.65, color = "black") +
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
    title = "SFS Loss of Functions variants",
    subtitle = "xv: x-var, yv: y-var"
  )
ggsave("plots/sfs-lof.pdf", height = 7, width = 8)
ggsave("plots/sfs-lof.png", height = 7, width = 8)


# SFS plot ----------------------------------------------------------------


pops <- c("YRI", "IBS", "MXL", "MXB")
variant_cats <- c("intergenic", "synonymous", "missense", "loss_of_function")

general %>% 
  filter(
    Population %in% pops,
    Variant_Category %in% variant_cats
  ) %>% 
  mutate(
    Variant_Category = factor(Variant_Category, levels = variant_cats)
  ) %>% 
  ggplot(aes(x = n, y = spectrum_norm, color = Population)) +
  geom_line(size = 0.2) +
  geom_point(size = 0.4) +
  scale_color_viridis_d(option = "D", direction = -1)  +
  scale_y_log10() +
  facet_grid(~Variant_Category) +
  labs(
    y = "log10  Proportion",
    title = "SFS shape in different populations",
    x = "Derived allel frequency"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )
ggsave("plots/sfs-Pops.pdf", height = 4, width = 7.5)
ggsave("plots/sfs-Pops.png", height = 4, width = 7.5)
 