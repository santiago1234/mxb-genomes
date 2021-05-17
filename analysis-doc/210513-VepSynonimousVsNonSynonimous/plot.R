library(tidyverse)
library(gghighlight)
library(ggthemes)
library(gridExtra)
theme_set(theme_tufte(base_family = "Helvetica"))

sfs_syn <- read_csv("results/sfs-synonymous.csv") %>% 
  mutate(
    cat = "Synonymous"
  )

sfs_nonsyn <- read_csv("results/sfs-nonsynonymous.csv") %>% 
  mutate(
    cat = "Non-synonymous"
  )


sfs <- bind_rows(sfs_nonsyn, sfs_syn) %>% 
  filter(!n %in% c(0, 100)) # Drop fixated allels


sfs %>% 
  filter(Population == "MXL") %>% 
  ggplot(aes(x = n, y = Freq)) +
  geom_line() +
  scale_y_log10() +
  facet_grid(.~cat)



sfs %>% 
  ggplot(aes(x = n, y = Freq, group = Population)) +
  geom_area() +
  scale_y_log10() +
  facet_grid(cat ~ Population)



sfs %>% 
  filter(Population == "MXB") %>% 
  ggplot(aes(x = n, y = Freq)) +
  geom_bar(aes(fill = cat), stat = "identity", position = position_dodge(width = 0.7)) +
  scale_y_log10() +
  scale_fill_manual(values = c("#d95f02", "#7570b3"))



sfs %>% 
  ggplot(aes(x = n, y = Freq, color = cat, fill = cat)) +
  geom_point(size = 0.5) +
  geom_line() +
  geom_area(position = "dodge", alpha = 0.1) +
  scale_y_log10() +
  facet_grid(~Population) +
  scale_fill_manual(values = c("#d95f02", "#7570b3")) +
  scale_color_manual(values = c("#d95f02", "#7570b3"))


# non-syn case ------------------------------------------------------------

p_nsyn <- sfs %>% 
  filter(cat == "Non-synonymous") %>% 
  ggplot(aes(x = n, y = Freq, color = Population)) +
  geom_line(color= "black") +
  geom_rangeframe(color = "black", size = 0.5) +
  gghighlight(use_direct_label = FALSE)  +
  scale_y_log10() +
  facet_wrap(~Population) +
  labs(
    title = "Non-synonymous Variants"
  )

p_syn <- sfs %>% 
  filter(cat == "Synonymous") %>% 
  ggplot(aes(x = n, y = Freq, color = Population)) +
  geom_line(color= "black") +
  geom_rangeframe(color = "black", size = 0.5) +
  gghighlight(use_direct_label = FALSE)  +
  scale_y_log10() +
  facet_wrap(~Population) +
  labs(
    title = "Synonymous Variants"
  )

grid.arrange(p_syn, p_nsyn, nrow = 1)


# make another plot -------------------------------------------------------

sfs_l <- 
  sfs %>% 
  pivot_wider(names_from = cat, values_from = Freq) %>% 
  mutate(
    foldChange = log10(`Non-synonymous` / Synonymous)
  )

sfs_l %>% 
  #filter(Population == "CHB") %>% 
  ggplot(aes(x = n, y = foldChange)) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~Population)
  

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
  #scale_fill_brewer(type = "qual", palette = 2) +
  scale_fill_viridis_d(option = "C") +
  facet_grid(cat ~ .) +
  theme_bw()

