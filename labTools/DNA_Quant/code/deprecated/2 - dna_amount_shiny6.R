#Run on: http://10.5.146.65/DNA_Quantification/
#Copy as "app.R" to /srv/shiny-server/DNA_Quantification/2-DNA_concentration
#Logs go here: /var/log/shiny-server/

##TODO - init with VB
##TODO - use gamlss? - or make bayes vs freq a choice? freq for empirical priors?
##TODO - sample_id_store - back to sample_id when writing CSV

library(shiny) |> suppressMessages() |> suppressWarnings()
library(shinyBS) |> suppressMessages() |> suppressWarnings()         # For collapsible panels
library(shinyFiles)  |> suppressMessages() |> suppressWarnings()     # For file saving dialogs on the server side
library(readxl) |> suppressMessages() |> suppressWarnings()
library(readr) |> suppressMessages() |> suppressWarnings()
library(dplyr) |> suppressMessages() |> suppressWarnings()
library(janitor) |> suppressMessages() |> suppressWarnings()
library(stringr) |> suppressMessages() |> suppressWarnings()
library(purrr) |> suppressMessages() |> suppressWarnings()
library(brms) |> suppressMessages() |> suppressWarnings()
library(cmdstanr) |> suppressMessages() |> suppressWarnings()
library(callr) |> suppressMessages() |> suppressWarnings()
library(ggplot2) |> suppressMessages() |> suppressWarnings()
library(tidybayes) |> suppressMessages() |> suppressWarnings()       # for add_epred_draws and point_interval
library(forcats) |> suppressMessages() |> suppressWarnings()         # for fct_reorder
library(scales) |> suppressMessages() |> suppressWarnings()          # for log10_trans and comma_format

if(Sys.info()["nodename"] == 'gawain'){
  set_cmdstan_path('/home/shiny/.cmdstan/cmdstan-2.36.0')
}

# rebuild_cmdstan(
#   dir = cmdstan_path(),
#   cores = 20,
#   quiet = FALSE,
#   timeout = 600
# )

# Helper function to compute common prefix (ignoring case)
common_prefix <- function(s1, s2) {
  s1 <- tolower(s1)
  s2 <- tolower(s2)
  cp <- ""
  for (i in seq_len(min(nchar(s1), nchar(s2)))) {
    if (substr(s1, i, i) == substr(s2, i, i)) {
      cp <- paste0(cp, substr(s1, i, i))
    } else {
      break
    }
  }
  cp <- sub("[ _-]+$", "", cp)
  return(cp)
}

# Function to read in the plate map from Excel
read_plate_map <- function(filepath, sheet_name) {
  read_xlsx(filepath,
            sheet = sheet_name, 
            na = c('na', 'NA', 'N/A'), 
            .name_repair = ~ vctrs::vec_as_names(..., 
                                                 repair = "unique",
                                                 quiet = TRUE)) %>%
    clean_names() %>% 
    # filter(!is.na(plate_row), !is.na(plate_row), !is.na(plate_id)) %>%
    mutate(across(where(is.character), str_to_lower)) %>%
    rename(plate_id = any_of('sample_plate')) %>% #colnames
    rename_with(~str_replace(., 'sample', 'plate'),
                .cols = c(starts_with('sample'), 
                          -any_of('sample_id'))) %>%
    rename_with(~str_replace(., 'column$', 'col')) %>%
    rename(dna_plate_id = plate_id,
           dna_plate_col = plate_col,
           dna_plate_row = plate_row)
}

jaccard_similarity <- function(vec1, vec2) {
  # Compute unique values in each vector
  uniq1 <- unique(vec1)
  uniq2 <- unique(vec2)
  
  # Intersection and union of the unique values
  inter <- intersect(uniq1, uniq2)
  union_vals <- union(uniq1, uniq2)
  
  # Jaccard index = size of intersection / size of union
  if (length(union_vals) == 0) return(0)
  length(inter) / length(union_vals)
}

get_subsequent_row_col <- function(plate_col, df1){
  if (length(plate_col) == 0) {
    stop("The best candidate column was not found in df1.")
  }
  
  best_candidate_index <- which(colnames(df1) == plate_col)
  # Look at the columns that come after the best candidate column
  subsequent_columns <- colnames(df1)[(best_candidate_index + 1):length(colnames(df1))]
  
  # Identify the first column with "row" in its name (case-insensitive)
  row_candidate <- subsequent_columns[grep("row", subsequent_columns, ignore.case = TRUE)][1]
  col_candidate <- subsequent_columns[grep("col", subsequent_columns, ignore.case = TRUE)][1]
  
  if (is.na(row_candidate)) {
    stop("No column after the best candidate column contains 'row' in its name.")
  } 
  c(row_candidate, col_candidate)
}

identify_join_cols <- function(df1, df2){
  candidate_cols <- df1 %>% select(where(is.character)) %>% colnames()
  
  # Compute the overlap score between df2$dna_plate_id and each candidate column in df1
  overlap_scores <- sapply(candidate_cols, function(col_name) {
    jaccard_similarity(df1[[col_name]], df2$dna_plate_id)
  })
  
  
  col_pair <- set_names('dna_plate_id',
                        names(which.max(overlap_scores)))
  
  row_col_match <- set_names(c('dna_plate_row', 'dna_plate_col'), 
                             get_subsequent_row_col(names(col_pair), df1))
  
  out <- c(col_pair, row_col_match)
  
  message('joining: dna_plate_map - ', str_c(names(out), collapse = '; '))
  message('joining: quant_plates - ', str_c(out, collapse = '; '))
  
  out
}

join_quants_map <- function(map_data = dna_plate_map, quant_data = quant_plates, quant_type = 'original'){
  # #quant_type: "original" or "requant"
  # join_vars <- c('dna_plate_id',
  #                'dna_plate_col',
  #                'dna_plate_row')
  # 
  # if(quant_type == 'requant'){
  #   names(join_vars) <- c('zymo_plate_id',
  #                         'zymo_plate_col',
  #                         'zymo_plate_row')
  # }
  join_vars <- identify_join_cols(map_data, quant_data)
  message('DEBUG: join_vars: ', str_c(join_vars, collapse = '; '))
  message('DEBUG: names(join_vars): ', str_c(names(join_vars), collapse = '; '))
  
  
  out <- inner_join(map_data,
                    quant_data,
                    by = join_vars)
  attr(out, 'join_cols') <- join_vars
  out
}

process_merged_data <- function(df){
  the_cols <- colnames(df)
  join_cols <- names(attr(df, 'join_cols'))
  
  if(!any(str_detect(the_cols, 'sample_type'))){
    df <- mutate(df, sample_type = 'sample')
  }
  # write_csv(df, 'tmp1.csv')
  if(any(str_detect(the_cols, 'sample_id'))){
    df <- rename(df, sample_id_store = 'sample_id')
  }
  # write_csv(df, 'tmp2.csv')
  df %>%
    mutate(across(ends_with('row'), str_to_upper),
           sample_id = str_c(!!sym(join_cols[1]), '_', !!sym(join_cols[2]), !!sym(join_cols[3])),
           is_control = str_detect(sample_type, 'control'))
}

ui <- fluidPage(
  titlePanel("Data Input Shiny App"),
  tabsetPanel(
    tabPanel("Data Input",
             sidebarLayout(
               sidebarPanel(
                 fileInput("excel_file", "Upload Plate Map File", 
                           accept = c(".xlsx", ".xls")),
                 selectInput("sheet_name", "Excel Sheet Name", choices = NULL),
                 fileInput("csv_files", "Upload Quant Output Files (select all)", 
                           multiple = TRUE, accept = ".csv"),
                 selectInput("quant_type", "Quant Stage", 
                             choices = c('Original' = 'original', 
                                         "Requantification" = 'requant'), 
                             selected = 'original'),
                 actionButton("load_data", "Load Data"),
               ),
               mainPanel(
                 h4("Merged Data"),
                 DT::DTOutput("merged_table")
               )
             )
    ),
    tabPanel("DNA Concentrations",
             sidebarLayout(
               sidebarPanel(
                 # 1. Select DNA Concentration input
                 selectInput("y_var", "Select DNA Concentration", choices = NULL),
                 
                 # 2. Model Settings (collapsible)
                 bsCollapse(
                   bsCollapsePanel("Model Settings", 
                                   numericInput("num_chains", "Number of Chains", value = 4, min = 1, step = 1),
                                   numericInput("iter_sampling", "Number of Sampling Iterations", value = 2000, min = 1, step = 1),
                                   numericInput("iter_warmup", "Number of Warmup Iterations", value = 1000, min = 1, step = 1),
                                   numericInput("thin", "Thinning Interval", value = 10, min = 1, step = 1),
                                   style = "primary"
                   ),
                   open = NULL
                 ),
                 
                 # 3. Priors Settings (collapsible, with paired inputs)
                 bsCollapse(
                   bsCollapsePanel("Priors Settings",
                                   h4("Main Model Priors"),
                                   fluidRow(
                                     column(6, numericInput("intercept_prior_mean", "Intercept Prior Mean", value = 1.5)),
                                     column(6, numericInput("intercept_prior_sd", "Intercept Prior SD", value = 2.5))
                                   ),
                                   fluidRow(
                                     column(6, numericInput("beta_prior_mean", "Beta Prior Mean", value = 0)),
                                     column(6, numericInput("beta_prior_sd", "Beta Prior SD", value = 2.5))
                                   ),
                                   fluidRow(
                                     column(6, numericInput("var_prior_param1", "Variance Prior Parameter 1", value = 2)),
                                     column(6, numericInput("var_prior_param2", "Variance Prior Parameter 2", value = 1))
                                   ),
                                   h4("Shape Priors"),
                                   fluidRow(
                                     column(6, numericInput("interceptShape_prior_mean", "Intercept Shape Prior Mean", value = 0)),
                                     column(6, numericInput("interceptShape_prior_sd", "Intercept Shape Prior SD", value = 2.5))
                                   ),
                                   fluidRow(
                                     column(6, numericInput("betaShape_prior_mean", "Beta Shape Prior Mean", value = 0)),
                                     column(6, numericInput("betaShape_prior_sd", "Beta Shape Prior SD", value = 2.5))
                                   ),
                                   fluidRow(
                                     column(6, numericInput("varShape_prior_param1", "Variance Shape Prior Parameter 1", value = 2)),
                                     column(6, numericInput("varShape_prior_param2", "Variance Shape Prior Parameter 2", value = 1))
                                   ),
                                   style = "primary"
                   ),
                   open = NULL
                 ),
                 
                 # 4. Flag Settings (grouped but not collapsible)
                 bsCollapse(
                   bsCollapsePanel("Flag Settings",
                                   numericInput("min_volume", "Minimum Pipettable Volume", value = 0.75, min = 0),
                                   numericInput("max_low_volume", "Maximum Low Volume", value = 4, min = 0),
                                   numericInput("target_dna", "Target Amount of DNA", value = 2, min = 0),
                                   style = "primary"
                   ),
                   open = "Flag Settings"
                 ),
                 
                 # Other controls
                 checkboxInput("log_transform", "Plot Log10 Transform", value = TRUE), 
                 actionButton("fit_model", "Fit Model"),
                 
                 # Download control using shinyFiles save button
                 # shinyFiles::shinySaveButton("download_csv", "Download CSV", "Select download location", filetype = list(csv = "csv"))
                 downloadButton("download_csv", "Download Sample Concentrations (CSV)")
               ),
               mainPanel(
                 plotOutput("dna_plot")
               )
             )
    )
  )
)

server <- function(input, output, session) {
  
  # Set up roots for shinyFiles (here, we allow navigation of the home directory)
  roots <- c(home = "~")
  shinyFileSave(input, "download_csv", roots = roots, session = session, filetypes = c("csv"))
  
  # Update the sheet selection when an Excel file is uploaded
  observeEvent(input$excel_file, {
    req(input$excel_file)
    sheets <- excel_sheets(input$excel_file$datapath)
    updateSelectInput(session, "sheet_name", choices = sheets)
  })
  
  # Update y_var choices based on numeric columns in the combined CSV files
  observe({
    req(input$csv_files)
    
    # Combine all CSV files into one tibble
    quant_plates <- map_dfr(input$csv_files$datapath, ~ read_csv(.x, show_col_types = FALSE))
    
    # Identify numeric columns in the CSV data, excluding names containing 'col', 'column', or 'volume'
    numeric_cols <- str_subset(
      names(quant_plates)[sapply(quant_plates, is.numeric)],
      pattern = 'col|column|volume',
      negate = TRUE
    )
    
    updateSelectInput(session, "y_var", choices = numeric_cols,
                      selected = str_subset(numeric_cols, '_per_ul'))
  })
  
  # Reactive expression for joined data, triggered by the Load Data button
  joined_data <- eventReactive(input$load_data, {
    req(input$excel_file, input$sheet_name, input$csv_files, input$y_var)
    
    # Read the plate map using the provided function
    dna_plate_map <- read_plate_map(input$excel_file$datapath, input$sheet_name)
    
    # Read and combine the CSV files
    quant_plates <- map_dfr(input$csv_files$datapath, ~ read_csv(.x, show_col_types = FALSE)) %>%
      rename_with(~str_replace(., 'column', 'col'))
    
    message('DEBUG: quant_type: ', input$quant_type)
    # Merge the plate map with the quant data using an inner join
    merged_dna_quants <- join_quants_map(dna_plate_map,
                                         quant_plates,
                                         input$quant_type) %>%
      process_merged_data() 
    
    merged_dna_quants
  })
  
  # Render the merged table when joined_data is updated
  output$merged_table <- DT::renderDT({
    joined_data()  # Ensure the merged data is re-rendered
  })
  
  # # Render the merged table in the Data Input tab
  # output$merged_table <- DT::renderDT({
  #   joined_data()
  # })
  
  # Modified model_fit using callr::r_bg for asynchronous Stan sampling
  model_fit <- eventReactive(input$fit_model, {
    req(joined_data(), input$y_var)
    
    # message("DEBUG: Local Files: ", str_c(list.files(), collapse = ', '))
    message("DEBUG: Chains: ", input$num_chains)
    message("DEBUG: Cores: ", min(c(parallelly::availableCores(),
                                    input$num_chains)))
    message("DEBUG: Threads: ", parallelly::availableCores() %/% min(c(parallelly::availableCores(),
                                                                       input$num_chains)))
    message("DEBUG: Available Cores: ", parallelly::availableCores())
    
    
    model_formula <- str_c(input$y_var, '~ is_control + (1 | sample_id)')
    dna_form <- bf(as.formula(model_formula),
                   shape ~ is_control + (1 | sample_id))
    dna_model_bayes <- brm(dna_form,
                           family = Gamma(link = 'log'),
                           data = filter(joined_data(),
                                         !!sym(input$y_var) > 0),
                           empty = TRUE)
    
    if(Sys.info()["nodename"] == 'gawain'){
      dna_model <- cmdstan_model('model/dna_concentration_threaded.stan',
                                 cpp_options = list(stan_threads = TRUE))
    } else {
      dna_model <- cmdstan_model('dna_concentration_threaded.stan',
                                 cpp_options = list(stan_threads = TRUE))
    }

    message("DEBUG: Read Model")
    # Note: Stan expects "inteceptShape_prior" (without the "r") as in the original code.
    data_stan <- c(
      make_standata(bf(as.formula(model_formula),
                       shape ~ is_control + (1 | sample_id)),
                    family = Gamma(link = 'log'),
                    data = filter(joined_data(), 
                                  !!sym(input$y_var) > 0),
                    threads = parallelly::availableCores() %/% min(c(parallelly::availableCores(),
                                                                     input$num_chains))),
      list(intercept_prior = c(input$intercept_prior_mean, input$intercept_prior_sd),
           beta_prior = c(input$beta_prior_mean, input$beta_prior_sd),
           var_prior = c(input$var_prior_param1, input$var_prior_param2),
           inteceptShape_prior = c(input$interceptShape_prior_mean, input$interceptShape_prior_sd),
           betaShape_prior = c(input$betaShape_prior_mean, input$betaShape_prior_sd),
           varShape_prior = c(input$varShape_prior_param1, input$varShape_prior_param2))
    )
    
    # Create a persistent output directory for the CmdStanR CSV files.
    
    
    if(Sys.info()["nodename"] == 'gawain'){
      output_dir <- "model/cmdstan_output"
    } else {
      output_dir <- "./cmdstan_output"
    }
    # output_dir <- tempdir()
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    message("DEBUG: Made TEMP-DIR")
    
    # Define the function to run sampling in the background using callr::r_bg.
    sampling_function <- function(data_stan, output_dir, num_chains, iter_sampling, iter_warmup, thin) {
      library(cmdstanr)
      if(Sys.info()["nodename"] == 'gawain'){
        set_cmdstan_path('/home/shiny/.cmdstan/cmdstan-2.36.0')
      }
      
      message("DEBUG: Inside Sub-process")
      if(Sys.info()["nodename"] == 'gawain'){
        dna_model <- cmdstan_model('model/dna_concentration_threaded.stan',
                                   cpp_options = list(stan_threads = TRUE))
      } else {
        dna_model <- cmdstan_model('dna_concentration_threaded.stan',
                                   cpp_options = list(stan_threads = TRUE))
      }
      
      fit <- dna_model$sample(
        data = data_stan,
        chains = num_chains,
        parallel_chains = min(c(parallelly::availableCores(), num_chains)),
        threads_per_chain = parallelly::availableCores() %/% min(c(parallelly::availableCores(), num_chains)),
        iter_sampling = iter_sampling + iter_warmup,
        iter_warmup = iter_warmup,
        thin = thin,
        refresh = 100,
        show_messages = TRUE,
        show_exceptions = FALSE,
        output_dir = output_dir
      )
      message("DEBUG: Chain 1:", stringr::str_c(fit$output(1), collapse = '\n'))
      
      fit$output_files()
    }
    
    # Launch sampling in the background with user-defined settings.
    message("DEBUG: Start background R")
    bg_process <- callr::r_bg(
      func = sampling_function,
      args = list(data_stan = data_stan,
                  output_dir = output_dir,
                  num_chains = input$num_chains,
                  iter_sampling = input$iter_sampling,
                  iter_warmup = input$iter_warmup,
                  thin = input$thin),
      stdout = "|",
      stderr = "2>&1"
    )
    
    # Use withProgress to update the UI while sampling runs.
    withProgress(message = "Fitting Stan model", value = 0, {
      progress <- 0
      while (bg_process$is_alive()) {
        Sys.sleep(1)
        new_lines <- bg_process$read_output_lines()
        # message("DEBUG: Background R:", str_c(new_lines, collapse = '\n'))
        if (length(new_lines) > 0) {
          percents <- sapply(new_lines, function(line) {
            m <- regmatches(line, regexpr("[0-9]+%", line))
            if (length(m) > 0 && nzchar(m)) as.numeric(gsub("%", "", m)) else NA
          })
          percents <- percents[!is.na(percents)]
          if (length(percents) > 0) {
            max_pct <- max(percents) / 100
            if (max_pct > progress) {
              progress <- max_pct
              setProgress(value = progress,
                          detail = paste0(round(progress * 100), "% complete"))
            }
          }
        }
      }
      setProgress(value = 1, detail = "100% complete")
    })
    
    # Wait for the background process to finish and retrieve the output file paths.
    output_files <- bg_process$get_result()
    message("DEBUG: End Sub-process")
    
    # Reconstruct the Stan fit using the output files.
    dna_model_bayes$fit <- read_csv_as_stanfit(output_files, model = dna_model)
    dna_model_bayes <- rename_pars(dna_model_bayes)
    
    # Delete the temporary output directory and its files.
    unlink(output_dir, recursive = TRUE)
    
    dna_model_bayes
  })
  
  # Display the model summary in the DNA Concentrations tab
  output$model_summary <- renderPrint({
    model_fit()
  })
  
  #### DNA Summary Calculations ####
  
  dna_summary <- eventReactive(input$fit_model, {
    req(joined_data(), input$y_var, model_fit())
    dna_model_bayes <- model_fit()
    
    # write_csv(joined_data(), 'tmp.csv')
    # write_rds(dna_model_bayes, 'tmp.rds')
    
    quant_files_summarized <- joined_data() %>%
      select(sample_id, sample_type, is_control) %>%
      distinct() %>%
      add_epred_draws(dna_model_bayes, 
                      allow_new_levels = TRUE) %>%
      point_interval(.width = 0.95) %>%
      rename(!!sym(str_c(input$y_var, "_mean")) := .epred,
             !!sym(str_c(input$y_var, "_lwr95")) := .lower,
             !!sym(str_c(input$y_var, "_upr95")) := .upper) %>%
      select(-.width, -.point, -.interval, -.row) %>%
      mutate(!!sym(str_c(input$y_var, "_normspread")) := 
               (!!sym(str_c(input$y_var, "_upr95")) - !!sym(str_c(input$y_var, "_lwr95"))) / !!sym(str_c(input$y_var, "_mean")))
    
    mean_dna_concentration <- add_epred_draws(newdata = data.frame(is_control = FALSE),
                                              dna_model_bayes,
                                              re.form = NA) %>%
      point_interval() %>%
      pull(.epred)
    
    dna_quantity_interval <- lm(as.formula(str_c("log(", input$y_var, "_mean) ~ is_control")), 
                                quant_files_summarized) %>%
      predict(newdata = data.frame(is_control = c(TRUE, FALSE)), 
              interval = "prediction", level = 0.95) %>%
      as_tibble() %>%
      mutate(is_control = c(TRUE, FALSE),
             across(where(is.numeric), exp)) %>%
      select(is_control, mean_dna = fit, lwr_limit = lwr, upr_limit = upr)
    
    dna_variability_interval <- lm(as.formula(str_c("log(", input$y_var, "_normspread) ~ is_control")), 
                                   data = quant_files_summarized) %>%
      predict(newdata = data.frame(is_control = c(TRUE, FALSE)), 
              interval = "prediction", level = 0.95) %>%
      as_tibble() %>%
      mutate(is_control = c(TRUE, FALSE),
             across(where(is.numeric), exp)) %>%
      select(is_control, lwr_limit = lwr, upr_limit = upr)
    
    list(
      quant_files_summarized = quant_files_summarized,
      mean_dna_concentration = mean_dna_concentration,
      dna_quantity_interval = dna_quantity_interval,
      dna_variability_interval = dna_variability_interval
    )
  })
  
  # Reactive expression for the flagged table that will be used for both plotting and download.
  flagged_table <- reactive({
    req(joined_data(), dna_summary())
    quant_files_summarized <- dna_summary()$quant_files_summarized
    mean_dna_concentration <- dna_summary()$mean_dna_concentration
    dna_quantity_interval <- dna_summary()$dna_quantity_interval
    dna_variability_interval <- dna_summary()$dna_variability_interval
    
    # write_csv(quant_files_summarized, 'model/before_flagging.csv')
    
    quant_files_flagged <- quant_files_summarized %>%
      mutate(ul_per_rxn = (input$target_dna * mean_dna_concentration) / !!sym(str_c(input$y_var, "_mean")),
             ul_per_rxn = case_when(ul_per_rxn > input$max_low_volume & is_control ~ input$max_low_volume,
                                    ul_per_rxn > input$max_low_volume & !is_control ~ input$max_low_volume,
                                    ul_per_rxn < input$min_volume ~ NA_real_,
                                    TRUE ~ ul_per_rxn),
             ul_per_rxn = round(ul_per_rxn),
             dilution_factor = case_when(is.na(ul_per_rxn) ~ !!sym(str_c(input$y_var, "_mean")) / mean_dna_concentration,
                                         TRUE ~ 0),
             dilution_factor = round(dilution_factor),
             postDilution_concentration = if_else(dilution_factor > 0,
                                                  !!sym(str_c(input$y_var, "_mean")) / dilution_factor,
                                                  NA_real_),
             postDilution_ul_per_rxn = (input$target_dna * mean_dna_concentration) / postDilution_concentration,
             postDilution_ul_per_rxn = round(postDilution_ul_per_rxn),
             rxn_ng = if_else(is.na(ul_per_rxn),
                              postDilution_concentration * postDilution_ul_per_rxn,
                              !!sym(str_c(input$y_var, "_mean")) * ul_per_rxn)) %>%
      left_join(select(dna_variability_interval, is_control, upr_limit),
                by = "is_control") %>%
      mutate(flags = case_when(
        is_control & !!sym(str_c(input$y_var, "_upr95")) > mean_dna_concentration ~ "Contaminated Control",
        !is.na(postDilution_concentration) & !!sym(str_c(input$y_var, "_normspread")) > upr_limit ~ "Excess & Variable DNA",
        !is.na(postDilution_concentration) ~ "Excess DNA",
        !!sym(str_c(input$y_var, "_normspread")) > upr_limit ~ "Variable DNA",
        TRUE ~ "Good Sample"),
        .after = sample_type) %>%

      select(sample_id,
             sample_type,
             all_of(str_c(input$y_var, c("_mean", "_lwr95", "_upr95"))),
             flags,
             ul_per_rxn, dilution_factor,
             starts_with("postDilution"),
             rxn_ng)
    
    # write_csv(quant_files_flagged, 'model/after_flagging.csv')
    
    quant_files_flagged
  })
  
  #### Plot Generation Reactive ####
  
  plot_obj <- reactive({
    req(flagged_table())
    req(dna_summary())
    quant_files_flagged <- flagged_table()
    mean_dna_concentration <- dna_summary()$mean_dna_concentration
    dna_quantity_interval <- dna_summary()$dna_quantity_interval
    dna_variability_interval <- dna_summary()$dna_variability_interval
    
    colour_fill_choices <- distinct(quant_files_flagged, sample_type, flags) %>%
      mutate(int_name = str_c(sample_type, flags, sep = "."),
             colour = case_when(
               flags == "Good Sample" ~ "black",
               flags == "Excess DNA" ~ "#F8766D",
               flags == "Variable DNA" ~ "#00BFC4",
               flags == "Contaminated Control" ~ "purple",
               flags == "Excess & Variable DNA" ~ "orange",
               TRUE ~ "pink"),
             fill = if_else(sample_type == "sample", "white", colour))
    
    colour_values <- distinct(colour_fill_choices, flags, colour) %>%
      with(., set_names(colour, flags))
    
    fill_values <- distinct(colour_fill_choices, int_name, fill) %>%
      with(., set_names(fill, int_name))
    
    quant_plot <- quant_files_flagged %>%
      # mutate(dna_plate_well_id = str_c(coalesce(dna_plate_row, zymo_plate_row), 
      #                                  sprintf('%02d', coalesce(dna_plate_col, zymo_plate_col))),
      #        facet_plate = coalesce(dna_plate_id, zymo_plate_id)) %>%
      mutate(dna_plate_well_id = str_extract(sample_id, '_[a-zA-Z][0-9]+$') %>% str_remove('_'),
             facet_plate = str_remove(sample_id, str_c('_', dna_plate_well_id))) %>%
      ggplot(aes(y = dna_plate_well_id,
                 x = !!sym(str_c(input$y_var, "_mean")),
                 shape = sample_type,
                 col = flags, 
                 fill = interaction(sample_type, flags))) +
      geom_vline(data = dna_quantity_interval,
                 aes(xintercept = mean_dna,
                     linetype = is_control)) +
      geom_errorbar(aes(xmin = !!sym(str_c(input$y_var, "_lwr95")),
                        xmax = !!sym(str_c(input$y_var, "_upr95"))),
                    show.legend = FALSE) +
      geom_point(size = 1.5) +
      scale_colour_manual(values = colour_values) +
      scale_fill_manual(values = fill_values) +
      scale_shape_manual(values = c("sample" = 'circle filled', 
                                    "extraction control" = 'triangle filled', 
                                    "field control" = 'triangle down filled', 
                                    "filter control" = 'square filled'),
                         breaks = c("sample", "field control", "filter control", "extraction control"),
                         labels = str_to_title) +
      scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"),
                            labels = c("TRUE" = "Control", "FALSE" = "Sample")) +
      guides(shape = guide_legend(override.aes = list(fill = if_else(str_detect(unique(quant_files_flagged$sample_type), 'control'), 
                                                                     'black', 'white'),
                                                      size = 3)),
             fill = "none",
             colour = guide_legend(override.aes = list(size = 3))) +
      facet_wrap(~facet_plate, nrow = 1, scales = "free_y") +
      theme_bw() +
      labs(x = "DNA Concentration (ng/uL) \u00b1 95% CI",
           y = "Sample Well ID",
           color = "Sample Flag",
           shape = "Sample Type",
           linetype = "Sample Type") +
      theme(axis.text.x = element_text(size = 7))
    
    if (input$log_transform) {
      quant_plot <- quant_plot + 
        scale_x_continuous(transform = scales::log10_trans(),
                           labels = scales::comma_format())
    }
    
    quant_plot
  })
  
  output$dna_plot <- renderPlot({
    plot_obj()
  })
  
  # Download handler using shinyFiles save button selection.
  output$download_csv <- downloadHandler(
    filename = function() {
      # Attempt to extract the chosen file path using parseSavePath().
      fileinfo <- shinyFiles::parseSavePath(roots, input$download_csv)
      if(nrow(fileinfo) > 0) {
        return(fileinfo$datapath)
      } else {
        req(input$excel_file, input$csv_files)
        prefix <- common_prefix(input$excel_file$name, input$csv_files$name[1])
        prefix <- tolower(prefix)
        paste0(prefix, "_", input$quant_type, "_DNA_concentrations_", Sys.Date(), ".csv")
      }
    },
    content = function(file) {
      out_table <- select(joined_data(), -rfu, -input$y_var, 
                          -is_control,
                          -starts_with('quant')) %>% 
        full_join(flagged_table(),
                  by = 'sample_id') %>%
        mutate(flags = na_if(flags, 'Good Sample'),
               quant_stage = input$quant_type) %>%
        relocate(quant_stage, 
                 .before = everything())
      
      if(any(str_detect(colnames(out_table), 'sample_id_store'))){
        out_table <- mutate(out_table, 
                            sample_id = sample_id_store,
                            .keep = 'unused')
      }
      
      write_csv(out_table,
                file = file,
                na = '')
    },
    contentType = "text/csv"
  )
  
}

shinyApp(ui = ui, server = server)
