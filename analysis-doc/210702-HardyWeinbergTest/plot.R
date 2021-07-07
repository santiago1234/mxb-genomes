library(tidyverse)
library(gridExtra)
library(UpSetR)

spectrums <- read_csv("results/all-spectrums-tidy.csv")

spectrums %>% 
  ggplot(aes(x = n1, y = n2, fill = log10(freq))) +
  geom_tile() +
  scale_fill_viridis_c(na.value = NA, option = "B") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(hw_pass~Sample) +
  theme(
    panel.background = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    x = "MXB, derived allele freq",
    y = "MXL, derived allele freq"
  )


# upset R plot ------------------------------------------------------------

droped_snps <- list(
  MX = read_lines("data/hw-nopasan-MX.txt"),
  MXB = read_lines("data/hw-nopasan-MXB.txt"),
  NAT = read_lines("data/hw-nopasan-NAT.txt")
)

upset(fromList(droped_snps))

