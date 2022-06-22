library(tidyverse)

sfs <- read_csv('results/SFSs_convoluted.csv')

sfs <- 
  sfs %>% 
  filter(
    !Derived_Freq %in% c(min(Derived_Freq), max(Derived_Freq))
  ) %>% 
  pivot_longer(cols = -Derived_Freq, names_to = 'SFS', values_to = 'Count')


sfs %>% 
  ggplot(aes(x = Derived_Freq, y = Count, color=SFS)) +
  geom_point() +
  scale_y_log10()
ggsave('SFS-counts.pdf')

sfs <- sfs %>% 
  group_by(SFS) %>% 
  mutate(
    Scaled = Count / sum(Count)
  )


sfs %>% 
  ggplot(aes(x = Derived_Freq, y = Scaled, color=SFS)) +
  geom_point() +
  scale_y_log10() +
  geom_line()
ggsave('SFS-porportions.pdf')