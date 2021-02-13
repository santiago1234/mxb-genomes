library(dplyr)
library(readr)
library(tidyr)
library(ggthemes)
library(gghighlight)

theme_set(theme_tufte(base_family = "Helvetica"))

meta_dat <- read_tsv("resources/genomes-metadata/50Genomes_info.txt")
eigen_vals <- read_csv("results/QC/pca_results.eigenval", col_names = "eigen")

n_com <- 20
header <- paste0("PC_", 1:n_com)
header <- c("M", "ID", header)
pca_res <- read_delim(
  "results/QC/pca_results.eigenvec",
  delim = " ", col_names = header
)



eigen_vals <- 
  eigen_vals %>% 
  mutate(
    component = header[c(-1, -2)],
    var_exp = (eigen / sum(eigen)) * 100
  )

pca_res <- 
  pca_res %>% 
  mutate(
    Sample_ID = paste0(M, ID)
    ) %>% 
  select(Sample_ID, everything()) %>% 
  select(-M, -ID) %>% 
  inner_join(meta_dat)



# basic PCA showing comp 1 and 2 ------------------------------------------

low_cov_samples <- c("MXB20306", "MXB20305", "MXB15210")

pca_res_d1 <- 
  pca_res %>% 
  select(Sample_ID:PC_3) %>% 
  pivot_longer(cols = c(PC_2, PC_3), names_to = "PC", values_to = "comp")

  
pca_res_d1 %>% 
  ggplot(aes(x = PC_1, y = comp)) +
  geom_point(aes(color = Sample_ID %in% low_cov_samples)) +
  facet_wrap(~PC, scales = "free") +
  scale_color_manual(name = "Low coverage", values = c("grey", "steelblue")) +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    legend.position = "bottom"
  ) +
  labs(
    caption = "PC1 6%, PC2 5.7%, and PC3 5.4%",
    y = "PC_X"
  )

ggsave("results/plots/qc/pca_comp1_2_3.png", height = 3, width = 5)
ggsave("results/plots/qc/pca_comp1_2_3.pdf", height = 3, width = 5)

set.seed(12)
etnics <- 
  pca_res$Ethnicity_Geography %>% 
  unique() %>% 
  sample(length(.))

pca_res$Ethnicity_Geography_color <- factor(
  pca_res$Ethnicity_Geography, levels = etnics
)

pca_res %>% 
  ggplot(aes(x = PC_1, y = PC_2, fill = Ethnicity_Geography_color)) +
  geom_point(shape=21, size=2) +
  scale_fill_viridis_d(option = "A") +
  gghighlight() +
  labs(
    title = "Proyection of samples in the two principal components"
  ) +
  facet_wrap(~Ethnicity_Geography) +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line()
  )

ggsave("results/plots/qc/pca_withEthnicity_Geography.png", height = 7, width = 9)
ggsave("results/plots/qc/pca_withEthnicity_Geography.pdf", height = 7, width = 9)
