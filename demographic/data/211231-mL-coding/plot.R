library(tidyverse)

mls <- list.files('results/', pattern = 'mLs-chr*', full.names = T)

mLs <- mls %>% 
  map_df(read_csv)

mLs <- mLs %>% 
  group_by(Q) %>% 
  summarise(
    mL = sum(mL)
  )

mLs$Q <- factor(mLs$Q, levels = c('missense', 'synonymous', 'LOF'))

mLs %>% 
  ggplot(aes(x = Q, y = mL)) +
  geom_point() +
  geom_linerange(aes(ymin = 0, ymax=mL)) +
  ggrepel::geom_text_repel(aes(label=round(mL, 4))) +
  labs(
    y = 'mL',
    x = 'Consequence'
  )
ggsave('plots/mLs.pdf', height = 2.5, width = 4)



# counts ------------------------------------------------------------------

counts <- list.files('results/', pattern = 'counts-chr*', full.names = T)

counts <- counts %>% 
  map_df(read_csv)

counts <- counts %>% 
  group_by(Q) %>% 
  summarise(
    n = sum(n)
  )

counts$Q <- factor(mLs$Q, levels = c('missense', 'synonymous', 'LOF'))

counts %>% 
  ggplot(aes(x = Q, y = n)) +
  geom_point() +
  geom_linerange(aes(ymin = 0, ymax=n)) +
  labs(
    y = 'mL',
    x = 'Consequence'
  )
ggsave('plots/counts.pdf', height = 2.5, width = 4)


# ratios ------------------------------------------------------------------

mLs <- mLs %>% 
  pivot_wider(names_from = Q, values_from = mL) %>% 
  mutate(
    miss_to_syn = log2(missense / synonymous),
    lof_to_syn = log2(LOF / synonymous),
    lof_to_miss = log2(LOF / missense),
  ) %>% 
  select(contains('_to_')) %>%
  pivot_longer(cols = everything(), names_to = 'comparison', values_to = 'log2_ratio') %>% 
  mutate(
    data = 'scaled by u'
  )

counts <- counts %>% 
  pivot_wider(names_from = Q, values_from = n) %>% 
  mutate(
    miss_to_syn = log2(missense / synonymous),
    lof_to_syn = log2(LOF / synonymous),
    lof_to_miss = log2(LOF / missense),
  ) %>% 
  select(contains('_to_')) %>%
  pivot_longer(cols = everything(), names_to = 'comparison', values_to = 'log2_ratio') %>% 
  mutate(
    data = 'counts'
  )

bind_rows(counts, mLs) %>% 
  ggplot(aes(x = data, y = log2_ratio)) +
  geom_linerange(aes(ymin = 0, ymax = log2_ratio)) +
  geom_point() +
  facet_grid(~comparison) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave('plots/comparison-log2-ratio.pdf', height = 3, width = 4)
