library(tidyverse)
library(gridExtra)
library(cowplot)
# get the native american ancestry ----------------------------------------

admxt <- read_csv("../210310-AdmixurePCA-merged-data-with-1TGP/results/admixture.csv")

nat_p <- 
  admxt %>% 
  filter(
    K == "K = 3"
  ) %>% 
  filter(cluster_grp == "MXB", Superpopulation == "AMR") %>% 
  select(Sample, p, Population) %>% 
  rename(Samplename = Sample, NAT_p = p)

pops <- c("PEL", "CLM", "PUR", "MXL")
nat_p$Population <- factor(nat_p$Population, levels = pops)

dev_counts <- read_csv("results/derived_counts.csv")
dev_counts <- inner_join(dev_counts, nat_p)


# plots -------------------------------------------------------------------


make_plot <- function(CAT) {
  
  mxl <- dev_counts %>% 
    filter(Population == "MXL")
  mxl %>% 
    filter(variant == CAT) %>% 
    ggplot(aes(x = NAT_p, y = derived_count)) +
    ggpubr::stat_cor(color = "black", size = 2) +
    geom_smooth(method = "lm") +
    geom_point(color="grey30") +
    facet_wrap(Region~VarFreq, scales = "free_y") +
    theme_bw() +
    labs(
      title = CAT,
      x = "NAT ancestry proportion",
      y = "# Derived variants per individual"
    ) +
    theme(legend.position = "none") 
  
}


plots <- map(c("INTERGENIC", "SYNONYMOUS", "DELETERIOUS"), make_plot)

pdf("plots/mutation-burden-mxl.pdf", height = 5, width = 14)
plot_grid(plots[[1]], plots[[2]], plots[[3]], nrow = 1)
dev.off()