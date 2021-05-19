library(tidyverse)



sfs_syn <- read_csv("results/sfs-synonymous.csv") %>% 
  mutate(
    cat = "Synonymous"
  )

sfs_nonsyn <- read_csv("results/sfs-nonsynonymous.csv") %>% 
  mutate(
    cat = "Non-synonymous"
  )


sfs <- bind_rows(sfs_nonsyn, sfs_syn) %>% 
  filter(!n %in% c(min(.$n), max(.$n))) # Drop fixated allels

# focus on important pops -------------------------------------------------

pops <- c("YRI", "IBS", "MXL", "MXB")


sfs %>% 
  filter(
    Population %in% pops
  ) %>%
  mutate(
    Population = factor(Population, levels = pops)
  ) %>% 
  ggplot(aes(x = n, y = Freq, fill=Population)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_sqrt() +
  scale_fill_viridis_d(option = "D") +
  facet_grid(cat ~ .) +
  theme_bw()

ggsave("plots/sfs.pdf", height = 4, width = 9)
ggsave("plots/sfs.png", height = 4, width = 9)

# construc the plot in Henna et al 2015 -----------------------------------

maf_low <- sfs %>% 
  filter(n / max(n) < 0.18) %>% 
  mutate(
    MAF = "MAF < 0.18"
  )

maf_high <- sfs %>% 
  filter(n / max(n) > 0.82) %>% 
  mutate(
    MAF = "MAF > 0.82"
  )

maf <- bind_rows(maf_low, maf_high) %>% 
  filter(
    Population %in% pops
  ) %>%
  mutate(
    Population = factor(Population, levels = pops),
    id = paste0(cat, "\n", MAF)
  )

maf %>% 
  ggplot(aes(x = n, y = Freq, fill=Population)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
  scale_y_sqrt() +
  facet_wrap(~id, scales = "free") +
  scale_fill_viridis_d(option = "D") +
  labs(
    y = "Number of SNPs",
    x = "Derived allele frequency"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )

ggsave("plots/sfs-maf.pdf", height = 4, width = 9)
ggsave("plots/sfs-maf.png", height = 4, width = 9)


# Syn VS non-Syn ----------------------------------------------------------

sfs_dNdS <- 
  sfs %>% 
  pivot_wider(names_from = cat, values_from = Freq) %>% 
  mutate(
    logdNdS = log10(`Non-synonymous` / Synonymous)
  ) %>% 
  filter(
    Population %in% pops
  ) %>%
  mutate(
    Population = factor(Population, levels = pops)
  )


sfs_dNdS %>% 
  ggplot(aes(x = n, y = logdNdS)) +
  geom_bar(stat = "identity") +
  facet_grid(Population ~.) +
  theme_minimal() +
  labs(
    y = "log10 Non-Syn/Syn",
    x = "Derived allele frequency"
  )


ggsave("plots/NtoNS.pdf", height = 4, width = 6)
ggsave("plots/NtoNS.png", height = 4, width = 6)


