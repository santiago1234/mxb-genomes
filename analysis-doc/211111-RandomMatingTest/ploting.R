library(tidyverse)
library(scales)



# load tracts data --------------------------------------------------------

load_tracts <- function(tracts_file) {
  
  cnames <- c('chrn', 'spos', 'epos', 'Ancestry', 'sgpos', 'egpos')
  
  haplo <- basename(tracts_file) %>% 
    str_extract(pattern = "anc_[AB]_") %>% 
    str_extract("[AB]")
  
  read_tsv(tracts_file, col_names = cnames) %>% 
    mutate(
      haplo = haplo
    )
  
}


load_individual <- function(path_to_tracts, id, end = "_anc_A_cM.bed") {
  # tf (tracts file)
  
  tf_A <- paste0(id, end) %>% 
    file.path(path_to_tracts, .)
  
  tf_B <- paste0(id, end) %>% 
    str_replace("_A_", "_B_") %>% 
    file.path(path_to_tracts, .)
  
  tracts <- bind_rows(
    load_tracts(tf_A),
    load_tracts(tf_B),
  )
  

  # add this for plotting ---------------------------------------------------
  limit_y_max <- function(haplo) {
    if (haplo == "A") {
      return(Inf)
    }
    if (haplo == "B") {
      return(0.5) #0.49
    }
  }
  
  limit_y_min <- function(haplo) {
    if (haplo == "A") {
      return(0.5)
    }
    if (haplo == "B") {
      return(-Inf)
    }
  }
  
  tracts %>% 
    mutate(
      ymin_haplo = map_dbl(haplo, limit_y_min),
      ymax_haplo = map_dbl(haplo, limit_y_max),
    )
  
}



# plot --------------------------------------------------------------------

plot_tracts <- function(tracts) {
  
  tracts %>% 
    ggplot(aes(x = spos, fill = Ancestry)) +
    geom_rect(
      aes(xmin = spos, xmax=epos, ymin = ymin_haplo, ymax = ymax_haplo)
    ) +
    geom_hline(yintercept = 0.5, color ="white", size = 0.1) +
    facet_grid(chrn ~.,switch = "y") +
    scale_x_continuous(
      expand = c(0,0),
      breaks = c(10, 50, 100, 200) * 1e6,
      labels = label_number(scale = 1e-6, suffix = "Mb")
    ) +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      strip.text.y = element_text(angle = 0),
      panel.background = element_blank(),
      legend.position = c(.8, .2)
    ) +
    labs(
      y = "Chromosome"
    ) +
    scale_fill_brewer(type="qual", palette = "Accent") +
    #scale_fill_viridis_d(option = "C") +
    labs(
      x = "Position relative to chromosome start"
    )
  
}
