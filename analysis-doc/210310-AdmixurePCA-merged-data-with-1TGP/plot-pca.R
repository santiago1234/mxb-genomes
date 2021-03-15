library(tidyverse)
library(gghighlight)
library(patchwork)

mxb_meta <- read_tsv("../../resources/genomes-metadata/50Genomes_info.txt") %>% 
  select(Sample_ID, Region) %>% 
  mutate(
    Superpopulation = "MXB",
    Sample = str_replace(Sample_ID, "MXB", "MXB_")
  ) %>% 
  mutate(
    Population = "MXB"
  ) %>% 
  select(Sample, Population, Superpopulation)

oneT_meta <- read_tsv("../../resources/1TGP-samples-meta-data/igsr-1000genomes.tsv") %>% 
  select(`Sample name`, `Population code`, `Population name`, `Superpopulation code`) %>% 
  rename(Sample = `Sample name`, Population = `Population code`, Superpopulation = `Superpopulation code`) %>% 
  select(Sample, Population, Superpopulation)

pop_info <- bind_rows(mxb_meta, oneT_meta)

# plink removes the MXB_ from the sample name
pop_info <- 
  pop_info %>% 
  mutate(
    Sample = str_replace(Sample, "MXB_", "")
  )

write_csv(pop_info, "data/pop_info.csv")

n_com <- 20
header <- paste0("PC_", 1:n_com)
header <- c("M", "ID", header)
pca_res <- read_delim(
  "results/pca/pca.eigenvec",
  delim = " ", col_names = header
) %>% 
  rename(Sample = ID)

eigen_vals <- read_csv("results/pca/pca.eigenval", col_names = "eigen")

eigen_vals <- 
  eigen_vals %>% 
  mutate(
    component = header[c(-1, -2)],
    var_exp = (eigen / sum(eigen)) * 100
  )

# add pop-labels
pca_res <- inner_join(pca_res, pop_info)

# get explained variance

get_exp_var <- function(component_val) {
  filter(eigen_vals, component == component_val) %>% 
    pull(var_exp) %>% 
    round(2)
}

pc1_var <- get_exp_var("PC_1")
pc2_var <- get_exp_var("PC_2")
pc3_var <- get_exp_var("PC_3")

p1 <- pca_res %>% 
  mutate(id = paste0(Superpopulation, " | ", Population)) %>% 
  ggplot(aes(PC_1, PC_2, fill = id)) +
  geom_point(shape = 21) +
  gghighlight() +
  scale_fill_viridis_d() +
  facet_wrap(~id) +
  theme_bw() +
  labs(
    title = "PCA",
    subtitle = "Comp 1 & 2",
    x = paste0("PC_1 (", pc1_var, "%)"),
    y = paste0("PC_2 (", pc2_var, "%)")
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) 



p2 <- pca_res %>% 
  mutate(id = paste0(Superpopulation, " | ", Population)) %>% 
  ggplot(aes(PC_1, PC_3, fill = id)) +
  geom_point(shape = 21) +
  gghighlight() +
  scale_fill_viridis_d() +
  facet_wrap(~id) +
  theme_bw() +
  labs(
    title = "PCA",
    subtitle = "Comp 1 & 3",
    x = paste0("PC_1 (", pc1_var, "%)"),
    y = paste0("PC_3 (", pc3_var, "%)")
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) 

p1 + p2
ggsave("plots/pca.pdf", height = 6, width = 12)
ggsave("plots/pca.png", height = 6, width = 12)

