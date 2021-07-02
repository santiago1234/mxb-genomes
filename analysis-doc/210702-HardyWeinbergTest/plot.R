library(tidyverse)


load_spectrum <- function(sf) {
   hw_cut <- str_extract(sf, pattern = "e-\\d*")
   hw_cut <- paste0("HW-p > ", 1, hw_cut)
   spectrum <- read_csv(sf)
   
   n_vars <- spectrum %>% 
     select(-index) %>% 
     as.matrix() %>% 
     sum()
   
  id_name <- paste0(hw_cut, "\nSNPs=", n_vars) 
  
  spectrum %>% 
    rename(pop1 = index) %>% 
    pivot_longer(cols = -pop1, values_to = "freq") %>% 
    separate(col = pop1, into = c("pop1", "n1"), sep = "-") %>% 
    separate(col = name, into = c("pop2", "n2"), sep = "-") %>% 
    mutate(
      n1 = as.numeric(n1),
      n2 = as.numeric(n2),
      id = id_name
    )
    
  
}


spectrums <- 
  list.files("results", pattern = "sfs", full.names = T) %>% 
  map_df(load_spectrum)


ids <- unique(spectrums$id)
# order by snps left wich correlates with hw
the_order <- str_extract(ids, "\\d*$") %>% 
  as.numeric() %>% 
  order(decreasing = T)

spectrums$id <- factor(spectrums$id, levels = ids[the_order])


spectrums %>%
  ggplot(aes(x = n2, y = n1, fill = log10(freq))) +
  geom_tile() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_c(option = "C", na.value = "white") +
  labs(
    x = "MXL allele freq",
    y = "MXB allele freq",
    title = "Joint SFS"
  ) +
  facet_grid(~id) +
  theme(
    legend.position = "bottom"
  )

ggsave("plots/jsf-by-hwe-cuttoff.pdf", height = 3, width = 8)
ggsave("plots/jsf-by-hwe-cuttoff.png", height = 3, width = 8)
