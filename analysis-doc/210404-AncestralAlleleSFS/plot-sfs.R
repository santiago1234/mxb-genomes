library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(gghighlight)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))

sfs <- read_csv("data/sfs-22.csv")

corners <- c(0, 100) # this are the invariant site frquencies

sfs <- 
  sfs %>% 
  filter(!n %in% corners)


sfs %>% 
  ggplot(aes(x = n, y = Freq, color = Population)) +
  geom_line(color = "black", size = 0.5) +
  facet_wrap(~Population, scales = "free") +
  scale_y_log10() +
  labs(
    x = "Derived allele frequency",
    y = "log10 Frequency",
    subtitle = "SFS unfolded proyected (chr 22)"
  ) +
  gghighlight(use_direct_label = FALSE) +
  geom_rangeframe(color = "grey30")

ggsave("plots/sfs-chr22.png", height = 5, width = 6)
