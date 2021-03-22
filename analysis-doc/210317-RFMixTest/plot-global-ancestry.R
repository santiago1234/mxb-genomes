library(tidyverse)


glob_anc <- read_tsv('results/rfmix/rf.chr22.rfmix.Q', skip = 1) %>% 
  rename(
    sample = `#sample`
  )

panel <- read_tsv("../../resources/1TGP-samples-meta-data/igsr-1000genomes.tsv") %>% 
  select(`Sample name`, `Population code`) %>% 
  rename(sample = `Sample name`, Population = `Population code`)


# 3 add panel information

glob_anc <- inner_join(glob_anc, panel) %>% 
  pivot_longer(cols = c(AFR, EAS, EUR, NAT), names_to = "ancestry", values_to = 'p')


glob_anc %>% 
  ggplot(aes(x = ancestry, y = p, color = Population)) +
  geom_jitter()


# reorder samples so that the ancestry with the highest average  --------
# proportion within a population goes in either ascending -----------------
# order among the individuals ---------------------------------------------

order_x <- 
  glob_anc %>% 
  filter(ancestry == "EUR") %>% 
  arrange(-p) %>% 
  pull(sample)

glob_anc$sample <- factor(glob_anc$sample, levels = order_x)

order_pops <- c("EUR", "NAT","AFR", "EAS") %>% rev()
glob_anc$ancestry <- factor(glob_anc$ancestry, levels = order_pops)

# bar-type plots ----------------------------------------------------------

glob_anc %>% 
  ggplot(aes(x = sample, y = p, fill = ancestry)) +
  geom_bar(stat = 'identity', width = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_d(name = "assigment",option = "A", direction = -1) +
  facet_grid(~Population, space = 'free_x', scales = 'free') +
  labs(
    x = "individual",
    title = 'Global diploid ancestry estimates for chromosome 22 with RFMix'
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_line(size = 1/10)
  )

ggsave("plot.pdf", height = 2, width = 7)
ggsave("plot.png", height = 2, width = 7)
