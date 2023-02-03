library(gt)
library(tidyverse)
library(scales)
#######################################################################
#                            Load the data                            #
#######################################################################

pars_tab <- list.files('data/', pattern = 'para-table-*', full.names = TRUE) %>% 
  map_df(read_csv)

descriptions <- read_csv('data/descriptions.csv')
#######################################################################
#           Format numeric values for a nice representation           #
#######################################################################

format_values <- function(unit, val) {
  # Make pretty names for parameter values  
  my_custom_years <-  label_number(accuracy = 0.001, scale = 1/1e3, 
                                   prefix = "", suffix = "",
                                   big.mark = ",", decimal.mark = ".")
  
  my_custom_effective <- label_number(accuracy = 1, scale = 1,
                                      prefix = "", suffix = "",
                                      big.mark = ", ", decimal.mark = ".")
  
  my_custom_rates <- label_scientific(accuracy = 1, scale = 1, 
                                      prefix = "", suffix = "",
                                      big.mark = ",", decimal.mark = ".")
  

  formater <- 
    list(
      `Time event  (Thousands of years)` = my_custom_years,
      `Effective population size (# of individuals)` = my_custom_effective,
      `Migration rate (Fraction of individuals per generation moving between populations)` = my_custom_rates
    )
  
  
  if (is.na(val)) return('-')
  formater[[unit]](val)
  
}


pars_tab$opt_value_fmt <- map2_chr(pars_tab$parma_type, pars_tab$opt_value, format_values) %>% 
  paste0('$', ., '$')

pars_tab$std_err_fmt <- map2_chr(pars_tab$parma_type, pars_tab$std_err, format_values) %>% 
  paste0('$', ., '$')


# -------------------------------------------------------------------------
tab_std_errs <- 
  pars_tab %>% 
  select(param_name, std_err_fmt, snp_category, parma_type) %>% 
  pivot_wider(
    id_cols = c(param_name, parma_type),
    values_from = std_err_fmt, names_from = snp_category,
    names_prefix = 'Standard error '
  )

tab_opt_vals <- 
  pars_tab %>% 
  select(param_name, opt_value_fmt, snp_category, parma_type) %>% 
  pivot_wider(
    id_cols = c(param_name, parma_type),
    values_from = opt_value_fmt, names_from = snp_category,
    names_prefix = 'Optimal Value '
  )

tab_gt <- inner_join(tab_opt_vals, tab_std_errs) %>% 
  group_by(parma_type) %>% 
  inner_join(descriptions) %>% 
  select(param_code, description, everything())


# surround paras colum with $ to use math mode for numbers ----------------


# code to put table in latex format ---------------------------------------





# code to put table in latex format ---------------------------------------

tab_latex <- tab_gt %>% 
  select(-param_name) %>% 
  gt(rowname_col = 'param_code') %>% 
  tab_spanner(
    label = md('**Synonymous**'),
    columns = c('Optimal Value synonymous', 'Standard error synonymous')
  ) %>% 
  tab_spanner(
    label = md('**Intronic**'),
    columns = c('Optimal Value intronic', 'Standard error intronic')
  ) %>% 
  tab_spanner(
    label = md('**Intergenic**'),
    columns = c('Optimal Value intergenic', 'Standard error intergenic')
  ) %>% 
  tab_header(
    title = 'Inferred parameters',
    subtitle = md("_moments_")
  ) %>% 
  tab_caption(
    md('**Table 1**. This table displays the inferred parameter values  and the corresponding standard errors (SEs) calculated from allele frequencies across different SNP categories.')
  )


gt_latex_dependencies()

tab_latex %>% 
  as_latex() %>% 
  cat()



