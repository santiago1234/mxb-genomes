library(tidyverse)
library(scales)

res <- read_csv("results/results.csv")


res$cM <- res$bins * 100

res %>% 
  ggplot(aes(x = cM, y = pred, color = Ancestry, fill = Ancestry)) +
  geom_ribbon(aes(ymin = ci_l, ymax = ci_u), alpha = 0.3, color = "black", show.legend=FALSE, size= 0.2) +
  geom_line(show.legend=FALSE) +
  geom_line(show.legend=FALSE, color = "black", size = 0.2) +
  geom_point(aes(y=dat),shape = 21, size = 2, alpha = .9, color = "black") +
  scale_y_log10(
    oob = scales::squish_infinite,
    expand = c(0, 0),
    breaks = 10^(-1:5),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_x_continuous(expand = c(.01, 0)) +
  coord_cartesian(ylim = c(1, 14000), xlim = c(-1, 250)) +
  scale_fill_viridis_d(option = "C") +
  scale_color_viridis_d(option = "C") +
  labs(
    title = "MXL tract length distribution",
    x = "Tract Length (cM)",
    y = "Relative frequency"
  ) +
  theme(
    legend.position = c(.9, .8)
  )
ggsave("plots/tracts.png", height = 3, width = 4.4)
ggsave("plots/tracts.pdf", height = 3, width = 4.4)



ancp <- read_csv("results/ancprop.csv")
pulses <- read_csv("results/migpulses.csv")

ancp <- 
  ancp %>% 
  mutate(
    gen_past = max(ga) - ga
  ) %>% 
  pivot_longer(cols = EUR:AFR, names_to = "anc", values_to = "anc_prop")


max_gen <- max(ancp$ga)

pulses <- 
  pulses %>% 
  mutate(
    gen_past = max(ga) - ga
  )
  

ancp %>%
  ggplot(aes(x = gen_past, y = anc_prop, fill = anc)) +
  geom_bar(stat = "identity", width = 1) +
  geom_vline(xintercept = pulses$gen_past - 0.5, linetype = 2, size = 0.3) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 1)) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = c(0, max_gen),
    labels = c(paste(max_gen, " GA"), "today")
  ) +
  scale_fill_viridis_d(option = "C") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    x = "Time",
    y = "Ancestry\npropotion",
    subtitle = "Magnitud and origin of migrants"
  )

ggsave("plots/migmat.pdf", height = 2, width = 4)



# draw pie charts ---------------------------------------------------------

pulses %>% 
  filter(gen_past == 0) %>% 
  select(EUR:AFR) %>% 
  pivot_longer(cols = EUR:AFR) %>% 
  ggplot(aes(x = "", y = value, fill = name)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_viridis_d(option = "C") +
  theme_void() +
  theme(legend.position = "none")
ggsave("plots/pulse1.pdf", height = 1, width = 1)
