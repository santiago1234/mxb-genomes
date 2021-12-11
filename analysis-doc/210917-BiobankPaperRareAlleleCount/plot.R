library(tidyverse)
library(gridExtra)
library(cowplot)
library(ggthemes)
## make plot for MXBiobank paper

## there are two parameters the variant category and the ancestry
theme_set(theme_base(base_family = "Helvetica"))

admxt <- read_csv("../210310-AdmixurePCA-merged-data-with-1TGP/results/admixture.csv") %>% 
  rename(Samplename = Sample) %>% 
  filter(
    K == 'K = 3', # Using 3 clusters with ADMIXTURE
    Population == 'MXL'
  )

# dev counts contains data for AMR samples (i.e MXL, CLM, PEL, PUR)
dev_counts <- read_csv("results/derived_counts.csv")
dev_counts$VarFreq <- factor(dev_counts$VarFreq, levels = c("Rare (DAF <= 5%)", "Common (DAF <= 100%)"))


dat <- dev_counts %>% 
  inner_join(admxt)

# order variables
dat$variant <- factor(dat$variant, levels = c("INTERGENIC", "SYNONYMOUS", "MISSENSE", "DELETERIOUS"))
dat$cluster_grp <- factor(dat$cluster_grp, levels = c("MXB", "EUR", "AFR"))


dat %>% 
  filter(Region == 'GENOME') %>% 
  ggplot(aes(x = p, y = derived_count, color = cluster_grp)) +
  geom_smooth(aes(fill = cluster_grp), method = "lm", size = 1/3, alpha = 0.2) +
  geom_point(alpha = 0.8) +
  ggpubr::stat_cor(
    size = 3,
    label.x.npc = 0.45,
    label.y.npc = 1,
    method = "spearman",
    p.accuracy = 0.001,
    p.digits = 3,
    r.accuracy = 0.01
  ) +
  facet_wrap(VarFreq ~ variant, scales = "free", ncol = 4) +
  scale_color_manual(values = c("#E04B4B", "#63BC6A", "#6094C3")) +
  scale_fill_manual(values = c("#E04B4B", "#63BC6A", "#6094C3")) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 1/5),
    strip.background = element_blank(),
    legend.position = "none"
  ) +
  labs(
    title = "MEGA array SNPs",
    y = "Number of derived alleles carried by an individual",
    x = "Ancestry from Americas, Europe or Africa quantified using ADMIXTURE"
  )
ggsave("plots/array-variants.pdf", height = 5, width = 11)


dat %>% 
  filter(Region == 'GENOME') %>% 
  ggplot(aes(x = p, y = derived_count, color = cluster_grp)) +
  geom_smooth(aes(fill = cluster_grp), method = "lm", size = 1/3, alpha = 0.2) +
  geom_point(alpha = 0.8) +
  ggpubr::stat_cor(
    size = 3,
    label.x.npc = 0.45,
    label.y.npc = 1,
    method = "spearman",
    p.accuracy = 0.001,
    p.digits = 3,
    r.accuracy = 0.01
  ) +
  facet_wrap(VarFreq ~ variant, scales = "free", ncol = 4) +
  scale_color_manual(values = c("#E04B4B", "#63BC6A", "#6094C3")) +
  scale_fill_manual(values = c("#E04B4B", "#63BC6A", "#6094C3")) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 1/5),
    strip.background = element_blank(),
    legend.position = "none"
  ) +
  labs(
    title = "All SNPs (Genome)",
    y = "Number of derived alleles carried by an individual",
    x = "Ancestry from Americas, Europe or Africa quantified using ADMIXTURE"
  )
ggsave("plots/genome-variants.pdf", height = 5, width = 11)

