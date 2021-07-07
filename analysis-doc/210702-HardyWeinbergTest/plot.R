library(tidyverse)
library(gridExtra)
library(UpSetR)



spectrums <- read_csv("results/all-spectrums-tidy.csv") %>% 
  mutate(
    hw_pass = if_else(hw_pass == "pasan", "snps kept", "snps removed")
  )

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

ggsave("plots/hw-sfs.pdf", height = 4, width = 7)
ggsave("plots/hw-sfs.png", height = 4, width = 7)

# upset R plot ------------------------------------------------------------

droped_snps <- list(
  MX = read_lines("data/hw-nopasan-MX.txt"),
  MXB = read_lines("data/hw-nopasan-MXB.txt"),
  NAT = read_lines("data/hw-nopasan-NAT.txt")
)

upset(fromList(droped_snps))


# the final spectrum ------------------------------------------------------


spect_clean <- read_csv("results/final-sfs-tidy.csv")

spect_clean %>% 
  ggplot(aes(x = n1, y = n2, fill = log10(freq))) +
  geom_tile() +
  scale_fill_viridis_c(na.value = NA, option = "B") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    panel.background = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    x = "MXB, derived allele freq",
    y = "MXL, derived allele freq"
  )
ggsave("plots/filterd-HW-MXL&PEL.pdf", height = 4, width = 4)
ggsave("plots/filterd-HW-MXL&PEL.png", height = 4, width = 4)
