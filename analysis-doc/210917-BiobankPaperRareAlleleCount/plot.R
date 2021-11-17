library(tidyverse)
library(gridExtra)
library(cowplot)
## make plot for MXBiobank paper

## there are two parameters the variant category and the ancestry

admxt <- read_csv("../210310-AdmixurePCA-merged-data-with-1TGP/results/admixture.csv") %>% 
  rename(Samplename = Sample)


ancestry_data <- function(ancestry) {
  # we are using the admixture results
  # with K = 3
  # we only want the MXL samples
  admxt %>% 
    filter(
      K == 'K = 3',
      cluster_grp == ancestry,
      Population == 'MXL'
    )
  
}


## Variant counts data

dev_counts <- read_csv("results/derived_counts.csv")
#dev_counts <- inner_join(dev_counts, nat_p)

dev_counts$VarFreq <- factor(dev_counts$VarFreq, levels = c("Rare (DAF <= 5%)", "Common (DAF <= 100%)"))


# plotting function -------------------------------------------------------

make_plot <- function(CAT, ancestry) {
  
  mxl <- ancestry_data(ancestry) %>% 
    inner_join(dev_counts)
  
  
  x_pos <- if_else(ancestry == "MXB", 0.01, 0.25)
  mxl %>% 
    filter(variant == CAT) %>% 
    ggplot(aes(x = p, y = derived_count)) +
    geom_smooth(method = "lm", color = "grey50") +
    geom_point(color="grey70") +
    ggpubr::stat_cor(
      color = "firebrick",
      size = 3.5,
      label.x.npc  = x_pos,
      label.y.npc = 0.1,
      p.accuracy = 0.001,
      p.digits = 3,
      r.accuracy = 0.01
    ) +
    facet_wrap(Region~VarFreq, scales = "free_y") +
    theme_bw() +
    labs(
      title = CAT,
      x = paste(str_replace(ancestry, "MXB", "NAT"), "ancestry proportion"),
      y = "# Derived variants per individual"
    ) +
    theme(legend.position = "none") 
  
}

plot_ancestry <- function(ancestry) {
  map(c("INTERGENIC", "SYNONYMOUS", "MISSENSE"), function(x) make_plot(x, ancestry))
  
}


pdf("plots/mutation-burden-mxl.pdf", height = 12, width = 14)
plt_mxb <- plot_ancestry('MXB')
plt_eur <- plot_ancestry('EUR')
plt_afr <- plot_ancestry('AFR')
plot_grid(
  plt_mxb[[1]], plt_mxb[[2]], plt_mxb[[3]],
  plt_eur[[1]], plt_eur[[2]], plt_eur[[3]],
  plt_afr[[1]], plt_afr[[2]], plt_afr[[3]],
  nrow = 3
)
dev.off()

