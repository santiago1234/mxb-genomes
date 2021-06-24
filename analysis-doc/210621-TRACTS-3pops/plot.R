library(tidyverse)
library(scales)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))

fits <- read_csv("data/fits.csv")

model_order <- c("ppx_xxp", "ppx_xxp_pxx", "ccx_xxp", "ppx_ccx_xxp")

fits$cM <- fits$bins * 100
fits$model <- factor(fits$model, levels = model_order)

fits %>% 
  ggplot(aes(x = cM, y = pred, color = Ancestry, fill = Ancestry)) +
  geom_ribbon(aes(ymin = ci_l, ymax = ci_u), alpha = 0.3, color = "black", show.legend=FALSE, size= 0.2) +
  geom_line(show.legend=FALSE) +
  geom_line(show.legend=FALSE, color = "black", size = 0.2) +
  geom_point(aes(y=dat),shape = 21, size = 2, alpha = .9, color = "black") +
  geom_rangeframe(color = "black") +
  scale_y_log10(
    oob = scales::squish_infinite,
    expand = c(0, 0),
    breaks = 10^(-1:5),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_x_continuous(expand = c(.03, 0)) +
  coord_cartesian(ylim = c(1, 14000), xlim = c(-1, 250)) +
  scale_fill_viridis_d(option = "C") +
  scale_color_viridis_d(option = "C") +
  facet_grid(~model) +
  labs(
    x = "Tract Length (cM)",
    y = "Relative frequency"
  ) +
  theme(
    legend.position = c(.92, .7),
    legend.box.background = element_rect(colour = "black", size = 0.1)
  )

ggsave("plots/tract-length.pdf", height = 2.5, width = 7)

# mig mat -----------------------------------------------------------------

ancp <- read_csv("data/ancp.csv")
pulses <- read_csv("data/pulses.csv")

ancp <- 
  ancp %>%
  group_by(model) %>% 
  mutate(
    gen_past = max(ga) - ga
  ) %>% 
  ungroup() %>% 
  pivot_longer(cols = c(NAT, EUR, AFR), names_to = "anc", values_to = "anc_prop")



ancp$model <- factor(ancp$model, levels = model_order)
ancp %>%
  ggplot(aes(x = gen_past, y = anc_prop, fill = anc)) +
  geom_bar(stat = "identity", width = 1) +
  #geom_vline(xintercept = pulses$gen_past - 0.5, linetype = 2, size = 0.3) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 1)) +
  # scale_x_continuous(
  #   expand = c(0, 0),
  #   breaks = c(0, 3, 5, 6, 7, 11, 15),
  #   labels = c(15, 11, 7, 6, 5, 3, 0)
  # ) +
  scale_fill_viridis_d(option = "C") +
  facet_grid(~model, scales = "free_x") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.spacing.x = unit(1.5, "lines")
  ) +
  labs(
    x = "Time",
    y = "Ancestry\npropotion",
    subtitle = "Magnitud and origin of migrants"
  )

ggsave("plots/models-anc.pdf", height = 1.5, width = 7)


# draw pie charts ---------------------------------------------------------

# pie ppx_xxp

pulses %>% 
  filter(model == "ppx_xxp") %>% 
  filter(ga == 9) %>% 
  select(c(EUR, NAT, AFR)) %>% 
  pivot_longer(cols = EUR:AFR) %>% 
  ggplot(aes(x = "", y = value, fill = name)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_viridis_d(option = "C") +
  theme_void() +
  theme(legend.position = "none")
 
ggsave("plots/ppx_xppPulse.pdf", height = 1, width = 1)



pulses %>% 
  filter(model == "ppx_ccx_xxp") %>% 
  filter(ga == 17) %>% 
  select(c(EUR, NAT, AFR)) %>% 
  pivot_longer(cols = EUR:AFR) %>% 
  ggplot(aes(x = "", y = value, fill = name)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_viridis_d(option = "C") +
  theme_void() +
  theme(legend.position = "none")

ggsave("plots/ppx_ccx1stPulse.pdf", height = 1, width = 1)


pulses %>% 
  filter(model == "ppx_ccx_xxp") %>% 
  filter(ga == 15) %>% 
  select(c(EUR, NAT, AFR)) %>% 
  pivot_longer(cols = EUR:AFR) %>% 
  ggplot(aes(x = "", y = value, fill = name)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_viridis_d(option = "C") +
  theme_void() +
  theme(legend.position = "none")

ggsave("plots/ppx_ccxContinuousPulse.pdf", height = 1, width = 1)
