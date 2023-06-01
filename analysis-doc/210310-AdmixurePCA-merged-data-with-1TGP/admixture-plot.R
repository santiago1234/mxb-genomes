library(tidyverse)


# population information --------------------------------------------------

ids <- read.table("data/1TGP_and_50MXB-ALL-ldpruned_data.fam") %>% 
  select(V2) %>% 
  rename(Sample = V2) %>% 
  as_tibble()

pop_info <- read_csv("data/pop_info.csv")


# process admixture results -----------------------------------------------

load_admixuture_res <- function(path_to_res) {
  
  admix <- read.table(path_to_res) %>% 
    as_tibble() %>% 
    rename_all(.funs = function(x) str_replace(x, "V", ""))
  
  # the number of ancestral populations
  k_val <- ncol(admix)
  
  # add population information
  # the order of rows is the same as in the plink input file
  admix$Sample <- ids$Sample
  
  admix <- 
    admix %>% 
    pivot_longer(
      cols = -Sample,
      names_to = "k",
      values_to = "p"
    ) %>% 
    mutate(
      K = paste0("K = ", k_val)
    )
  
  # add population metadata
  admix <- inner_join(admix, pop_info)
  
  
  # now, i will add a cluster name
  # the idea of this is to align the cluster for colors
  
  cl_aligment <- admix %>% 
    group_by(k, Superpopulation) %>% 
    summarise(
      ancestry = mean(p)
    ) %>% 
    ungroup() %>% 
    group_by(k) %>% 
    filter(ancestry == max(ancestry)) %>% 
    ungroup() %>% 
    mutate(
      cluster_grp = Superpopulation
    ) %>% 
    select(k, cluster_grp)
  
  res_k <- length(unique(cl_aligment$cluster_grp))
  stopifnot(res_k == as.numeric(k_val))
  
  # add the cluster grp
  admix %>% 
    inner_join(cl_aligment)
  
  
}


# load data ---------------------------------------------------------------


admixture <- 
  list.files("results/", pattern = "*.Q", full.names = T) %>% 
  map_df(load_admixuture_res)

order_x <- admixture %>% 
  filter(
    K == "K = 4",
    cluster_grp == "EUR"
  ) %>% 
  arrange(p) %>% 
  pull(Sample)


admixture$Sample <- factor(admixture$Sample, levels = order_x)

# order for populations ---------------------------------------------------
pop_order <- c("YRI", "IBS", "GBR", "PUR", "CLM", "MXL", "PEL", "MXB", "CHB")
admixture$Population <- factor(admixture$Population, levels = pop_order)
write_csv(admixture, "results/admixture.csv")

color_order <- c("EUR","MXB", "AFR", "EAS", "AMR") %>% rev()
admixture$cluster_grp <- factor(admixture$cluster_grp, levels = color_order)
admixture %>% 
  #mutate(id = paste0(Superpopulation, "/", Population)) %>% 
  ggplot(aes(x = Sample, y = p, fill = cluster_grp)) +
  geom_col(width = 1) +
  facet_grid(K~Population, space = "free_x", scales = "free_x") +
  scale_fill_viridis_d(direction = -1) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, .25,.5, .75)) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing.x = unit(0.05, "lines"),
    legend.position = "none"
  ) +
  labs(
    x = NULL
  )

ggsave("plots/admixture.pdf", height = 2.5, width = 6)
ggsave("plots/admixture.png", height = 2.5, width = 6)



# add plot for paper fig1B ------------------------------------------------

my_pops <- c('YRI', 'IBS', 'CHB', 'MXB', 'CLM', 'MXL', 'PEL', 'PUR')

admixture2 <- 
  admixture %>% 
  filter(K == 'K = 4') %>% 
  filter(Population %in% my_pops)

admixture2$Population <- factor(admixture2$Population, levels = my_pops)

admixture2 %>% 
  filter(p > 0.0001) %>% 
  ggplot(aes(x = Sample, y = p, color = cluster_grp)) +
  geom_col(width = 1) +
  facet_grid(.~Population, space = "free_x", scales = "free_x") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, .25,.5, .75)) +
  scale_color_manual(
    values = c('#800080', '#fde725', '#3b528b', '#5ec962')
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing.x = unit(0.05, "lines"),
    legend.position = 'none'
  ) +
  labs(
    x = NULL
  )

ggsave("plots/admixtureK4.pdf", height = 1.2, width = 6)
