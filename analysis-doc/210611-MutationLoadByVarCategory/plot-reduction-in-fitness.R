library(tidyverse)

load <- read_csv("data/exploreload.py")

load %>% 
  ggplot(aes(x = q, y = l, color = s)) +
  geom_line(aes(group=s)) +
  scale_color_viridis_c() +
  geom_abline(linetype = 2, size = 0.2) +
  facet_grid(~h) +
  labs(
    x = "q: Alternative allele frequency",
    y = "l = s(2hq + (1-2h)*q^2)",
    title = "reduction in fitness contributed by a locus"
  )

ggsave("plots/fitness-reduction.png", height = 2.5, width = 7)


# laod computed for missense ----------------------------------------------

l <- read_csv("results/load-missense.csv")

populations <- c("YRI", "IBS", "MXL", "MXB")

l <- 
  l %>% 
  filter(
    !n %in% c(min(n), max(n)),
    Population %in% populations
  )


l %>% 
  ggplot(aes(x = n, y = Freq, color = Population)) +
  geom_point(size = 0.4) +
  geom_line() +
  scale_color_viridis_d(option = "D", direction = -1) +
  scale_y_log10() +
  labs(
    title = "SFS"
  )
ggsave("plots/sfs.png", height = 2.5, width = 4)
  
l %>% 
  ggplot(aes(x = n, y = 1 - fitness_reduction)) +
  geom_point() +
  geom_line() +
  labs(
    title = "s = 0.08, h = 0.25"
  )
ggsave("plots/NloadContribution.png", height = 2.5, width = 4)

l %>% 
  ggplot(aes(x = n, y = s_log_load, color = Population)) +
  geom_point() +
  geom_line() +
  scale_color_viridis_d(option = "D", direction = -1)
ggsave("plots/LoadSpectrum.png", height = 2.5, width = 4)

l %>% 
  group_by(Population) %>% 
  summarise(
    Log_L = sum(s_log_load)
  )
 