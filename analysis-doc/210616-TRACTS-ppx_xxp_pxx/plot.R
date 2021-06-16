library(tidyverse)
library(scales)

res <- read_csv("res.csv")


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
    title = "MXL",
    x = "Tract Length (cM)",
    y = "Relative frequency"
  )
  # theme(
  #   legend.position = c(.8, .8),
  #   axis.line.x = element_line(color="black", size = 0.2),
  #   axis.line.y = element_line(color="black", size = 0.2)
  # ) 
ggsave("result.pdf", height = 3, width = 4.4)
