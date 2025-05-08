
##TODO - sort out converting modelling to faster variant than MCMC sampling

library(tidyverse)
library(readxl)
library(janitor)
library(brms)
library(cmdstanr)
library(tidybayes)

#0 - Take user inputs
sample_metadata_file <- '../example_data/rbd_edna-extraction_plate-map_original.xlsx'
sheet_metadata <- "Plate-map-tidy"
sample_quant_file <- list.files('../example_data/processed_quant_plates', 
                                pattern = '*csv$',
                                full.names = TRUE)
y_var <- 'ng_per_ul'
minimum_pipettable_volume <- 0.75 #uL
maximum_low_volume <- 4
target_dna <- 2 #multiplier_on mean

# sample_quant_file <- list.files('../example_data/processed_requant_plates/', 
#                                 pattern = '*csv$',
#                                 full.names = TRUE)


read_plate_map <- function(filepath, sheet_name){
  read_xlsx(filepath,
            sheet = sheet_name, 
            na = c('na', 'NA', 'N/A')) %>%
    clean_names() %>% 
    filter(!is.na(plate_row), !is.na(plate_row), !is.na(plate_id)) %>%
    mutate(across(where(is.character), str_to_lower)) %>%
    rename(dna_plate_id = plate_id,
           dna_plate_col = plate_col,
           dna_plate_row = plate_row)
}

#1 - read in sample info
dna_plate_map <- read_plate_map(sample_metadata_file, sheet_metadata)

#2 - read in suite of quant outcomes
quant_plates <- read_csv(sample_quant_file, 
                         show_col_types = FALSE) 

#3 - merge together
merged_dna_quants <- inner_join(dna_plate_map,
                                quant_plates,
                                by = c('dna_plate_id',
                                       'dna_plate_col',
                                       'dna_plate_row')) %>%
  select(dna_extract_tube_id,
         sample_type,
         dna_plate_id,
         dna_plate_col,
         dna_plate_row,
         
         quant_row,
         quant_column,
         sample_volume,
         rfu,
         all_of(y_var)) %>%
  mutate(across(ends_with('row'), str_to_upper),
         sample_id = str_c(dna_plate_id, '_', dna_plate_row, dna_plate_col),
         is_control = str_detect(sample_type, 'control'))

#### Model DNA Quants ####
model_formula <- str_c(y_var,  '~ is_control + (1 | sample_id)')

# get_prior(bf(as.formula(model_formula),
#              shape ~ is_control + (1 | sample_id)),
#           family = Gamma(link = 'log'),
#           data = merged_dna_quants)

#Initialize empty BRMS model
dna_model_bayes <- brm(bf(as.formula(model_formula),
                          shape ~ is_control + (1 | sample_id)),

                       prior = prior(student_t(3, 1.5, 2.5),
                                     class = 'Intercept') +
                         prior(student_t(3, 0, 2.5),
                               class = 'Intercept',
                               dpar = 'shape') +

                         prior(student_t(3, 0, 2.5),
                               class = 'b') +
                         prior(student_t(3, 0, 2.5),
                               class = 'b',
                               dpar = 'shape') +

                         prior(gamma(2, 1),
                               class = 'sd') +
                         prior(gamma(2, 1),
                               class = 'sd',
                               dpar = 'shape'),

                       family = Gamma(link = 'log'),
                       data = filter(merged_dna_quants, 
                                     !!sym(y_var) > 0),

                       #save_model = 'dna_concentration.stan'

                       # chains = 0,
                       # iter = 5000,
                       # warmup = 2000,
                       # thin = 10,

                       # backend = 'cmdstanr',
                       # cores = parallelly::availableCores() - 1,
                       empty = TRUE)

dna_model <- cmdstan_model('dna_concentration.stan')

data_stan <- c(make_standata(bf(as.formula(model_formula),
                                shape ~ is_control + (1 | sample_id)),
                             family = Gamma(link = 'log'),
                             data = filter(merged_dna_quants, 
                                           !!sym(y_var) > 0)),
               list(intercept_prior = c(1.5, 2.5),
                    beta_prior = c(0, 2.5),
                    var_prior = c(2, 1),
                    inteceptShape_prior = c(0, 2.5),
                    betaShape_prior = c(0, 2.5),
                    varShape_prior = c(2, 1)))

dna_fit_prefit <- dna_model$laplace(
  data = data_stan
)




dna_fit <- dna_model$sample(
  data = data_stan,
  chains = 4,
  parallel_chains = parallelly::availableCores(),
  iter_sampling = 2000,
  iter_warmup = 1000,
  thin = 10,
  refresh = 500#, # print update every 500 iters
  # init = dna_fit_prefit
)

# dna_fit <- dna_model$laplace(
#   data = data_stan
# )

dna_model_bayes$fit <- read_csv_as_stanfit(dna_fit$output_files(), 
                                           model = dna_model)
dna_model_bayes <- rename_pars(dna_model_bayes)

#### Calculate DNA amounts ####
quant_files_summarized <- select(merged_dna_quants, 
                                 -starts_with('quant'),
                                 -rfu, -all_of(y_var)) %>%
  distinct() %>%
  add_epred_draws(dna_model_bayes, allow_new_levels = TRUE) %>%
  point_interval(.width = 0.95) %>%
  rename(!!sym(str_c(y_var, '_mean')) := .epred,
         !!sym(str_c(y_var, '_lwr95')) := .lower,
         !!sym(str_c(y_var, '_upr95')) := .upper) %>%
  select(-.width, -.point, -.interval, -sample_id, -.row) %>%
  
  mutate(!!sym(str_c(y_var, '_normspread')) := (!!sym(str_c(y_var, '_upr95')) - !!sym(str_c(y_var, '_lwr95'))) / !!sym(str_c(y_var, '_mean')),
         dna_plate_well_id = as.factor(sprintf("%s%02d", dna_plate_row, dna_plate_col))) %>%
  relocate(dna_plate_well_id, 
           .after = dna_plate_id) %>%
  arrange(dna_plate_id, dna_plate_col, dna_plate_row)

#### Calculate useful constants ####
mean_dna_concentration <- add_epred_draws(newdata = data.frame(is_control = FALSE),
                            dna_model_bayes,
                            re.form = NA) %>%
  point_interval() %>%
  pull(.epred)

dna_quantity_interval <- lm(as.formula(str_c('log(', y_var, '_mean) ~ is_control')), 
                            quant_files_summarized) %>%
  predict(newdata = data.frame(is_control = c(TRUE, FALSE)), 
          interval = 'prediction', level = 0.95) %>%
  as_tibble() %>%
  mutate(is_control = c(TRUE, FALSE),
         across(where(is.numeric), exp)) %>%
  select(is_control, mean_dna = fit, lwr_limit = lwr, upr_limit = upr)


dna_variability_interval <- lm(as.formula(str_c('log(', y_var, '_normspread) ~ is_control')), 
                               data = quant_files_summarized) %>%
  predict(newdata = data.frame(is_control = c(TRUE, FALSE)), 
          interval = 'prediction', level = 0.95) %>%
  as_tibble() %>%
  mutate(is_control = c(TRUE, FALSE),
         across(where(is.numeric), exp)) %>%
  select(is_control, lwr_limit = lwr, upr_limit = upr)

#### Flag Samples for Dilutions etc. ####
quant_files_flagged <- quant_files_summarized %>%
  
  mutate(ul_per_rxn = (target_dna * mean_dna_concentration) / !!sym(str_c(y_var,'_mean')),
         ul_per_rxn = case_when(ul_per_rxn > maximum_low_volume & is_control ~ maximum_low_volume,
                                ul_per_rxn > maximum_low_volume & !is_control ~ maximum_low_volume,
                                ul_per_rxn < minimum_pipettable_volume ~ NA_real_,
                                TRUE ~ ul_per_rxn),
         ul_per_rxn = round(ul_per_rxn),
         
         dilution_factor = case_when(is.na(ul_per_rxn) ~ !!sym(str_c(y_var,'_mean')) / mean_dna_concentration,
                                     TRUE ~ 0),
         dilution_factor = round(dilution_factor),
         
         postDilution_concentration = if_else(dilution_factor > 0,
                                          !!sym(str_c(y_var,'_mean')) / dilution_factor,
                                          NA_real_),
         postDilution_ul_per_rxn = (target_dna * mean_dna_concentration) / postDilution_concentration,
         postDilution_ul_per_rxn = round(postDilution_ul_per_rxn),
         
         rxn_ng = if_else(is.na(ul_per_rxn),
                          postDilution_concentration * postDilution_ul_per_rxn,
                          !!sym(str_c(y_var,'_mean')) * ul_per_rxn)) %>%
  
  left_join(select(dna_variability_interval, is_control, upr_limit),
            by = 'is_control') %>%
  
  mutate(flags = case_when(is_control & !!sym(str_c(y_var,'_upr95')) > mean_dna_concentration ~ "Contaminated Control",
                           !is.na(postDilution_concentration) & !!sym(str_c(y_var,'_normspread')) > upr_limit ~ 'Excess & Variable DNA',
                           !is.na(postDilution_concentration) ~ 'Excess DNA',
                           !!sym(str_c(y_var,'_normspread')) > upr_limit ~ 'Variable DNA',
                           TRUE ~ 'Good Sample'),
         .after = sample_type)  %>%
  
  select(dna_plate_id, dna_plate_col, dna_plate_row,
         dna_extract_tube_id, dna_plate_well_id,
         sample_type,
         all_of(str_c(y_var, c('_mean', '_lwr95', '_upr95'))),
         flags,
         ul_per_rxn, dilution_factor, 
         starts_with('postDilution'), 
         rxn_ng)


#### Make Plot ####
colour_fill_choices <- distinct(quant_files_flagged, 
       sample_type, flags) %>%
  mutate(int_name = str_c(sample_type, flags, sep = '.'),
         colour = case_when(flags == 'Good Sample' ~ 'black',
                            flags == 'Excess DNA' ~ '#F8766D',
                            flags == 'Variable DNA' ~ '#00BFC4',
                            flags == 'Contaminated Control' ~ 'purple',
                            flags == 'Excess & Variable DNA' ~ 'orange',
                            TRUE ~ 'pink'),
         fill = if_else(sample_type == 'sample', 'white', colour))

colour_values <- distinct(colour_fill_choices, flags, colour) %>%
  with(., set_names(colour, flags))

fill_values <- distinct(colour_fill_choices, int_name, fill) %>%
  with(., set_names(fill, int_name))

quant_plot <- quant_files_flagged %>%
  mutate(dna_plate_well_id = fct_reorder(dna_plate_well_id, 
                                         desc(dna_plate_well_id))) %>%
   ggplot(aes(y = dna_plate_well_id,
          x = !!sym(str_c(y_var, '_mean')),
          shape = sample_type,
          col = flags, 
          fill = interaction(sample_type, flags)))  +
   geom_vline(data = dna_quantity_interval,
              aes(xintercept = mean_dna,
                  linetype = is_control)) +
   
  
   geom_errorbar(aes(xmin = !!sym(str_c(y_var, '_lwr95')),
                     xmax = !!sym(str_c(y_var, '_upr95'))),
                 show.legend = FALSE) +
  geom_point(size = 1.5) +
  
  scale_colour_manual(values = colour_values) +
  
  scale_fill_manual(values = fill_values) +
  
   scale_shape_manual(values = c('sample' = 'circle filled', 
                                 'extraction control' = 'triangle filled', 
                                 'field control' = 'square filled', 
                                 'filter control' = 'triangle down filled'),
                      breaks = c('sample',
                                 'field control', 
                                 'filter control',
                                 'extraction control'),
                      labels = str_to_title) +
   scale_linetype_manual(values = c('TRUE' = 'dashed',
                                    'FALSE' = 'solid'),
                         labels = c('TRUE' = 'Control',
                                    'FALSE' = 'Sample')) +
  guides(shape = guide_legend(override.aes = list(fill = c('white',
                                                           rep('black', 3)), 
                                                  size = 3)),
         fill = 'none',
         colour = guide_legend(override.aes = list(size = 3))) +
  
   facet_wrap(~dna_plate_id,
               nrow = 1,
               scales =  "free_y") +
   theme_bw() +
   labs(x = "DNA Concentration (ng/uL) \u00b1 95% CI",
     y = "Sample Well ID",
     color = "Sample Flag",
     shape = "Sample Type",
     linetype = 'Sample Type') +
   theme(axis.text.x = element_text(size = 7))

quant_plot +
  scale_x_continuous(transform = scales::log10_trans(),
                     labels = scales::comma_format())

#### Output Files ####
write_csv(quant_files_flagged, 
          '../example_data/concentration_and_dilution/dna_concentration_and_dilution_r1.csv')




