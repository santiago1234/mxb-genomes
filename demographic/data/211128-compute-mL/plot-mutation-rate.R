library(tidyverse)
library(ggrepel)


# not coding mLs ----------------------------------------------------------

intergenic_ml <- read_lines('results/mL-intergenic.txt') %>% 
  str_split(" ") %>% 
  unlist() %>% 
  .[2] %>% 
  as.numeric()


intronic_ml <- read_lines('results/mL-introns.txt') %>% 
  str_split(" ") %>% 
  unlist() %>% 
  .[2] %>% 
  as.numeric()

# coding mLs --------------------------------------------------------------

load_ml_coding <- function(fp) {
  fp %>% 
    read_lines() %>% 
    as.numeric()
}

lof_ml <- load_ml_coding('../211231-mL-coding/results/mL-lof.txt')
missense_ml <- load_ml_coding('../211231-mL-coding/results/mL-missense.txt')
synonymous_ml <- load_ml_coding('../211231-mL-coding/results/mL-synonymous.txt')


# make a tibble with values -----------------------------------------------

mLs <- 
  tribble(
  ~Region,~Type,~mL,
  'non-coding', 'intergenic', intergenic_ml,
  'non-coding', 'intronic', intronic_ml,
  'coding', 'lof', lof_ml,
  'coding', 'missense', missense_ml,
  'coding', 'synonymous', synonymous_ml,
)


mLs %>% 
  ggplot(aes(y=reorder(Type, mL), x=mL)) +
  geom_point() +
  geom_segment(aes(x=0, xend=mL, y=reorder(Type, mL), yend=reorder(Type, mL)), size = 0.3) +
  geom_text_repel(aes(label=round(mL, 4)), color='darkblue') +
  theme_bw() +
  theme(panel.grid = element_blank())  +
  labs(
    y = NULL
  )
ggsave("mutation-rate.pdf", height = 2.5, width = 4)
