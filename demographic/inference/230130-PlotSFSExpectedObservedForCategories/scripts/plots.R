library(tidyverse)
library(scales)

dfiles <- list.files('results/', pattern = 'sfs*', full.names = TRUE)
sfs <- map_df(dfiles, read_csv)


order_panel <- c('YRI', 'IBS', 'CHB', 'MXB')
sfs$Population <- factor(sfs$Population, levels = order_panel)


d_data <- filter(sfs, SF_from == 'Data')
d_model <- filter(sfs, SF_from == 'Expected')

d_data %>% 
  ggplot(aes(x = Minor_allel_freq, y =  Frequency)) +
  geom_point(size = 2, color='dodgerblue4') +
  geom_point(size = 1, color='white') +
  geom_line(
    data = d_model,
    color =  'darkorange',
    linetype = 2,
    size = 1 / 2
  ) +
  geom_point(
    data = d_model,
    color =  'darkorange',
    size = 1
  ) +
  facet_grid(category~Population, scales = 'free_y') +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = 'l',size = 1 / 5) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank()
    
  )
ggsave('plots/sfs-expected-observed.pdf', height = 4, width = 6)
