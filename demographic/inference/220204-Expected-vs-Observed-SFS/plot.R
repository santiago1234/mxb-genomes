library(tidyverse)
library(scales)

d <- read_csv('results/expected-observed-1dSFS.csv')
d_background <- select(d, -Population)

order_panel <- c('YRI', 'IBS', 'CHB', 'MXB')
d$Population <- factor(d$Population, levels = order_panel)

d_data <- filter(d, SF_from == 'Data')
d_model <- filter(d, SF_from == 'Expected')

d_background <- rename(d_model, p = Population)

N <- 1.5
d_data %>% 
  ggplot(aes(x = Minor_allel_freq, y =  Frequency)) +
  geom_line(
    data = d_background,
    aes(group = p),
    color = 'grey',
    size = 1 / 4
    ) +
  geom_point(size = 3.5 / N, color='dodgerblue4') +
  geom_point(size = 1.2 / N, color='white') +
  geom_line(
    data = d_model,
    color =  'darkorange',
    linetype = 2,
    size = 1 / 2
  ) +
  geom_point(
    data = d_model,
    color =  'darkorange',
    size = 1.3 / N
  ) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = 'l',size = 1 / 5) +
  facet_wrap(~Population) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )
ggsave('plots/obs-vs-exp.pdf', height = 3, width = 5)
 