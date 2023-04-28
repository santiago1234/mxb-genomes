library(tidyverse)

cols_n <- c('chunk_ID', 'Bscore')
bscores <- read_delim('results/bscores.txt', delim = ' ', col_names = cols_n)

b_quants <- quantile(bscores$Bscore)

bscores %>% 
  ggplot(aes(x = Bscore)) +
  geom_histogram(bins = 50, fill = 'grey60') +
  geom_vline(
    xintercept = b_quants,
    color = 'steelblue'
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_bw() +
  coord_cartesian(
    xlim = c(min(bscores$Bscore), max(bscores$Bscore))
  ) +
  theme(
    panel.grid = element_blank()
  ) +
  labs(
    subtitle = 'Histogram of B-scores'
  )
ggsave('plots/bscores-hist.pdf', height = 2, width = 4)
