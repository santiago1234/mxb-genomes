library(tidyverse)
library(scales)

theme_set(theme_bw())
# functions to load the data ----------------------------------------------


params_type <- tribble(
  ~Param, ~Type,
  # Effective Popsize,
  "Ne", "Effective population size",
  "N_A", "Effective population size",
  "NB", "Effective population size",
  "NEu0", "Effective population size",
  "NEuF", "Effective population size",
  "NAs0", "Effective population size",
  "NAsF", "Effective population size",
  "NCmxI", "Effective population size",
  "NCmxF", "Effective population size",
  
  # Time in years
  "TA", "Time years ago",
  "TB", "Time years ago",
  "TF", "Time years ago",
  "TN", "Time years ago",
  
  # Migration rates
  "mAfB", "Migration rate",
  "mAfEu", "Migration rate",
  "mAfAs", "Migration rate",
  "mEuAs", "Migration rate",
  
)

make_fp <- function(varcat, q) {
  # Make relative path to the confidence table
  # for the given varcat
  
  c(
    "ooa" = paste0('results/ConfidenceIntervals/q', q, '-mdl_OOA-v_', varcat, ".tsv"),
    "nat" = paste0('results/ConfidenceIntervals/q', q, '-mdl_NAT-EXPANSION-v_', varcat, ".tsv")
  )
  
}

load_data <- function(varcat, q) {
  df <- make_fp(varcat, q)
  params <- map_df(df, read_tsv) %>% 
    replace_na(list(`#param` = 'N_A')) #One parameter is named NA which R consider a missing value
  
  params %>% 
    mutate(
      Category = varcat,
      quartile = q
    )
  
}

dpars <- expand_grid(quartile = 1:4, vcat = c('intergenic', 'intronic')) %>% 
  mutate(
    d = map2(vcat, quartile, load_data)
  ) %>% 
  select(d) %>% 
  unnest(d) %>% 
  rename(Param = `#param`) %>% 
  inner_join(params_type)

dpars$Param <- factor(dpars$Param, levels = params_type$Param)

# make nice visualization -------------------------------------------------

## => => => => => => => => => => <= <= <= <= <= <= <= <= <= <= <= <= 
## Time
## => => => => => => => => => => <= <= <= <= <= <= <= <= <= <= <= <= 

my_labs_for_times <- label_number(accuracy = 1, scale = 1/1e3, 
                                  prefix = "", suffix = " Kya",
                                  big.mark = ",", decimal.mark = ".")


y_labs <- paste0('q', 1:4)

position_quartile <- function(q, category) {
  # a function to position better the quartiles
  offset <- 0.2
  if (category == 'intronic') return(q - offset)
  else return(q + offset)
}


dpars <- dpars %>% 
  mutate(
    quartile_position = map2_dbl(quartile, Category, position_quartile)
  )


dpars %>% 
  filter(Type == 'Time years ago') %>% 
  ggplot(aes(x = opt_value, y = quartile_position, color = Category)) +
  geom_point(shape = 19) +
  geom_errorbarh(
    aes(xmax = opt_value + 2*std_err, xmin = opt_value - 2*std_err),
    height = 0
  ) +
  facet_grid(Param~.) +
  scale_x_log10(labels = my_labs_for_times) +
  scale_y_continuous(breaks = 1:4, labels = y_labs) +
  theme(
    panel.spacing = unit(0.2, "lines"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  coord_cartesian(
    xlim = c(25*1e3, 900*1e3)
  ) +
  scale_color_manual(values = c('#1b9e77', '#e7298a'))

ggsave('plots/PARS-time.pdf', height = 4.25, width = 5)

## => => => => => => => => => => <= <= <= <= <= <= <= <= <= <= <= <= 
## Effective size
## => => => => => => => => => => <= <= <= <= <= <= <= <= <= <= <= <= 

my_labs_for_sizes <- label_number(accuracy = 1, scale = 1/1e3, 
                                  prefix = "", suffix = " K",
                                  big.mark = ",", decimal.mark = ".")
dpars %>% 
  filter(Type == 'Effective population size') %>% 
  mutate(
    ci_l = opt_value - 2*std_err,
    ci_h = opt_value + 2*std_err,
  ) %>% 
  mutate(
    ci_l = if_else(ci_l < 0, 700, ci_l)
  ) %>% 
  ggplot(aes(x = opt_value, y = quartile_position, color = Category)) +
  geom_point(shape = 19) +
  geom_errorbarh(
    aes(xmax = ci_l, xmin = ci_h),
    height = 0
  ) +
  facet_grid(Param~.) +
  scale_x_log10(labels = my_labs_for_times) +
  scale_y_continuous(breaks = 1:4, labels = y_labs) +
  theme(
    panel.spacing = unit(0.2, "lines"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  coord_cartesian(
    xlim = c(0.5*1e3, 70*1e3)
  ) +
  scale_color_manual(values = c('#1b9e77', '#e7298a'))

ggsave('plots/PARS-size.pdf', height = 9, width = 5)


## => => => => => => => => => => <= <= <= <= <= <= <= <= <= <= <= <= 
## Migration Rate
## => => => => => => => => => => <= <= <= <= <= <= <= <= <= <= <= <= 

my_custom_rates <- label_scientific(accuracy = 1, scale = 1, 
                                    prefix = "", suffix = "",
                                    big.mark = ",", decimal.mark = ".")



dpars %>% 
  filter(Type == 'Migration rate') %>% 
  mutate(
    ci_l = opt_value - 2*std_err,
    ci_h = opt_value + 2*std_err,
  ) %>% 
  mutate(
    ci_l = if_else(ci_l < 0, 10**(-6), ci_l)
  ) %>% 
  ggplot(aes(x = opt_value, y = quartile_position, color = Category)) +
  geom_point(shape = 19) +
  geom_errorbarh(
    aes(xmax = ci_l, xmin = ci_h),
    height = 0
  ) +
  facet_grid(Param~.) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(breaks = 1:4, labels = y_labs) +
  theme(
    panel.spacing = unit(0.2, "lines"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  scale_color_manual(values = c('#1b9e77', '#e7298a')) +
  coord_cartesian(
    xlim = c(10**(-5.9), 10**(-3.7))
  )

ggsave('plots/PARS-mig.pdf', height = 4.25, width = 5)

