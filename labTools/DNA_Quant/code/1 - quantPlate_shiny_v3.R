#Run on: http://10.5.146.65/DNA_Quantification/
#Copy as "app.R" to /srv/shiny-server/DNA_Quantification/1-quant_plate
#Logs go here: /var/log/shiny-server/

##TODO - take as inputs sample & standard volumes. default = 1,5

#### Libraries ####
library(shiny) |> suppressMessages() |> suppressWarnings()
library(shinyFiles) |> suppressMessages() |> suppressWarnings()
library(readr) |> suppressMessages() |> suppressWarnings()
library(readxl) |> suppressMessages() |> suppressWarnings()
library(dplyr) |> suppressMessages() |> suppressWarnings()
library(stringr) |> suppressMessages() |> suppressWarnings()
library(tidyr) |> suppressMessages() |> suppressWarnings()
library(DT) |> suppressMessages() |> suppressWarnings()
library(janitor) |> suppressMessages() |> suppressWarnings()
library(ggplot2) |> suppressMessages() |> suppressWarnings()
library(patchwork) |> suppressMessages() |> suppressWarnings()
library(jsonlite) |> suppressMessages() |> suppressWarnings()
library(yardstick) |> suppressMessages() |> suppressWarnings()
library(outliers) |> suppressMessages() |> suppressWarnings()
library(purrr) |> suppressMessages() |> suppressWarnings()
library(waiter) |> suppressMessages() |> suppressWarnings()

#### Functions ####
# Helper to read CSV or Excel files
read_file <- function(path) {
  tryCatch({
    if (is.null(path) || !file.exists(path)) {
      stop("File does not exist or no file selected")
    }
    
    ext <- tools::file_ext(path)
    if(tolower(ext) %in% c("csv")) {
      data <- read_csv(path, show_col_types = FALSE)
    } else if(tolower(ext) %in% c("xlsx", "xls")) {
      data <- read_excel(path)
    } else {
      stop("Unsupported file type. Please upload a CSV, XLS, or XLSX file.")
    }
    
    if (nrow(data) == 0) {
      stop("The file is empty. Please check your file.")
    }
    
    return(data)
  }, error = function(e) {
    stop(paste("Error reading file:", e$message))
  })
}

read_map_file <- function(path){
  tryCatch({
    raw_input <- read_file(path) %>%
      clean_names() %>%
      rename(plate_id = any_of('sample_plate')) %>%
      rename_with(~str_replace(., 'col$', 'column')) %>%
      rename_with(~str_replace(., 'sample', 'plate'), 
                  .cols = any_of(c('sample_row', 'sample_column')))
    
    the_cols <- colnames(raw_input)
    
    # Basic validation - check we have some key columns
    if (!any(c("plate_id", "sample_id") %in% the_cols)) {
      stop("Plate map file must contain 'plate_id' or 'sample_id' column")
    }
    
    # Continue with original logic
    if(any(str_detect(the_cols, '(sample)_id')) & 
       !(any(str_detect(the_cols, '(sample|plate)_row')) & 
         any(str_detect(the_cols, '(sample|plate)_column')))){
      
      out <- raw_input %>%
        mutate(plate_row = if_else(!str_detect(plate_id, '[sS]tandard'),
                                   str_extract(sample_id, '[A-Za-z]+'),
                                   NA_character_),
               plate_column = if_else(!str_detect(plate_id, '[sS]tandard'),
                                      str_extract(sample_id, '[0-9]+'),
                                      NA_character_) %>% 
                 as.integer(),
               dna_concentration = if_else(!str_detect(plate_id, '[sS]tandard'),
                                           NA_real_,
                                           parse_number(sample_id) %>% suppressWarnings()))
      
    } else if (!any(str_detect(the_cols, '(sample)_id')) & 
               (any(str_detect(the_cols, '(sample|plate)_row')) & 
                any(str_detect(the_cols, '(sample|plate)_column')))) {
      
      standards <- filter(raw_input, str_detect(plate_id, '[sS]tandard')) %>%
        select(plate_column, plate_row) %>%
        filter(plate_column > 0, plate_row > 0) %>%
        mutate(across(where(is.character), parse_number)) 
      
      if (nrow(standards) == 0) {
        stop("No standards found in plate map. Please ensure you have rows with 'standard' in the plate_id.")
      }
      
      if (all(standards$plate_column < standards$plate_row)) {
        dna_col <- "plate_column"
      } else if (all(standards$plate_row < standards$plate_column)) {
        dna_col <- "plate_row"
      } else {
        stop("Cannot determine DNA concentration column from standards layout.")
      }
      
      out <- raw_input %>%
        mutate(dna_concentration = if_else(!str_detect(plate_id, '[sS]tandard'),
                                           NA_real_,
                                           !!sym(dna_col)),
               plate_row = if_else(str_detect(plate_id, '[sS]tandard'),
                                   NA_character_, 
                                   plate_row),
               plate_column = if_else(str_detect(plate_id, '[sS]tandard'),
                                      NA_integer_, 
                                      plate_column))
      if(!any(class(out$dna_concentration) %in% 'numeric')){
        out$dna_concentration <- parse_number(out$dna_concentration)
      }
      
    } else if(any(str_detect(the_cols, '(sample)_id')) & 
              (any(str_detect(the_cols, '(sample|plate)_row')) & 
               any(str_detect(the_cols, '(sample|plate)_column')))) {
      
      out <- raw_input %>%
        mutate(dna_concentration = if_else(!str_detect(plate_id, '[sS]tandard'),
                                           NA_real_,
                                           parse_number(sample_id) %>% suppressWarnings()))
      
    } else {
      stop(paste("Unexpected plate map file format. Please check that your file contains the required columns.",
                 "Available columns:", paste(the_cols, collapse = ", ")))
    }
    
    if(!any(str_detect(the_cols, 'volume'))){
      out <- mutate(out, 
                    sample_volume = if_else(!str_detect(plate_id, '[sS]tandard'),
                                            1, 5))
    } else {
      volume_column <- str_subset(the_cols, 'volume')
      out <- rename(out, 
                    sample_volume = all_of(volume_column)) %>%
        mutate(sample_volume = case_when(!is.na(sample_volume) ~ sample_volume,
                                         is.na(sample_volume) & str_detect(plate_id, '[sS]tandard') ~ 5,
                                         TRUE ~ 1))
    }
    
    out %>%
      mutate(across(where(is.character), str_to_lower))
    
  }, error = function(e) {
    stop(paste("Error processing plate map file:", e$message))
  })
}

read_machine_file <- function(path){
  tryCatch({
    raw_data <- read_file(path)
    
    # Check for required columns (before cleaning names)
    raw_cols <- tolower(colnames(raw_data))
    if (!("wells" %in% raw_cols && "value" %in% raw_cols)) {
      stop(paste("Raw RFU data file must contain 'Wells' and 'Value' columns.",
                 "Found columns:", paste(colnames(raw_data), collapse = ", ")))
    }
    
    data <- raw_data %>%
      clean_names() %>%
      select(wells, value)
    
    # Check wells format
    well_pattern <- "^[A-Za-z]+[0-9]+$"
    if (!any(str_detect(data$wells, well_pattern))) {
      stop("Wells column does not contain valid well identifiers (e.g., A1, B2, J10)")
    }
    
    data %>%
      mutate(quant_row = str_extract(wells, '[a-zA-Z]+') %>% str_to_lower(),
             quant_column = str_extract(wells, '[0-9]+') %>%
               as.integer(),
             rfu = value,
             .keep = 'unused')
    
  }, error = function(e) {
    stop(paste("Error processing raw RFU data file:", e$message))
  })
}


# Read and merge quant data
read_quant_data <- function(path_raw_data, path_plate_map, quant_kit) {
  tryCatch({
    if (is.null(path_raw_data) || is.null(path_plate_map)) {
      stop("Both raw data file and plate map file must be selected")
    }
    
    quant_raw_data <- read_machine_file(path_raw_data)
    quant_plate_map <- read_map_file(path_plate_map)
    
    standard_unit <- case_when(quant_kit == "accublue-nextgen" ~ "pg",
                               quant_kit == "accublue" ~ "ng",
                               quant_kit == "accuclear" ~ 'ng')
    
    quant_data <- left_join(quant_plate_map,
                            quant_raw_data,
                            by = c("quant_row", 'quant_column'))
    
    # Check if join was successful
    if (nrow(quant_data) == 0) {
      stop("No data resulted from joining plate map and raw data files")
    }
    
    # Check if we have RFU values
    rfu_count <- sum(!is.na(quant_data$rfu))
    if (rfu_count == 0) {
      stop("No matching data found between plate map and raw data. Please verify that the well positions (quant_row/quant_column) match between files.")
    }
    
    quant_data <- quant_data %>%
      mutate(!!sym(str_c(standard_unit, '_per_ul')) := dna_concentration,
             !!sym(str_c(standard_unit, '_per_well')) := !!sym(str_c(standard_unit, '_per_ul')) * sample_volume) %>%
      select(any_of(c('plate_id', 'plate_row', 'plate_column',
                      'quant_row', 'quant_column',
                      'sample_volume', 'rfu',
                      str_c(standard_unit, '_per_ul'),
                      str_c(standard_unit, '_per_well'))))
    
    return(quant_data)
    
  }, error = function(e) {
    stop(e$message)
  })
}


# Helper to compute common prefix (ignoring capitalization) and return lowercase.
common_prefix <- function(s1, s2) {
  s1 <- tolower(s1)
  s2 <- tolower(s2)
  min_len <- min(nchar(s1), nchar(s2))
  cp <- ""
  for(i in seq_len(min_len)) {
    if(substr(s1, i, i) == substr(s2, i, i)) {
      cp <- paste0(cp, substr(s1, i, i))
    } else {
      break
    }
  }
  cp <- sub("[ _-]+$", "", cp)
  return(cp)
}

loocv <- function(fit){
  #https://datacadamia.com/lang/r/cross_validation
  h <- lm.influence(fit)$h
  mean((residuals(fit)/(1-h))^2)
}

fit_model <- function(df_model, model_name, input_values){
  x_var <- input_values$x_var
  y_var <- input_values$y_var
  log_x <- input_values[[paste0("log_x_", model_name)]]
  log_y <- input_values[[paste0("log_y_", model_name)]]
  x_term <- ifelse(log_x, paste0("log10(", x_var, ")"), x_var)
  y_term <- ifelse(log_y, paste0("log10(", y_var, ")"), y_var)
  
  if(is.character(df_model[['included_in_model']])){
    df_model[['included_in_model']] <- df_model[['included_in_model']] == 'Yes'
  }
  
  if(model_name == "linear") {
    formula <- as.formula(paste(y_term, "~", x_term))
  } else if(model_name == "quadratic") {
    formula <- as.formula(paste(y_term, "~ poly(", x_term, ", 2)"))
  } else if(model_name == "cubic") {
    formula <- as.formula(paste(y_term, "~ poly(", x_term, ", 3)"))
  } else if(model_name == "power") {
    formula <- as.formula(paste("log10(", y_var, ") ~ log10(", x_var, ")"))
  }
  model_obj <- tryCatch(
    lm(formula, data = filter(df_model, included_in_model)),
    error = function(e) {
      showNotification(paste("Error fitting", model_name, "model:", e$message), type = "error")
      NULL
    }
  )
  
  model_obj
}

post_process_model <- function(sub_model, data, model_name, input_values){
  y_var <- input_values$y_var
  log_y <- input_values[[paste0("log_y_", model_name)]]
  y_term <- ifelse(log_y, paste0("log10(", y_var, ")"), y_var)
  
  # predictions 
  obs_y <- data[[y_var]]
  if(str_detect(y_term, 'log')){
    pred_y <- 10^predict(sub_model, data)
  } else {
    pred_y <- predict(sub_model, data)
  }
  
  tibble(fit_model = list(sub_model),
         fitted = list(pred_y),
         rmse = possibly(rmse_vec, otherwise = NA_real_)(truth = obs_y, 
                                                         estimate = pred_y),
         r2 = possibly(rsq_vec, otherwise = NA_real_)(truth = obs_y, 
                                                      estimate = pred_y),
         loocv = loocv(sub_model))
}

identify_outliers <- function(x, test = dixon.test, alpha = 0.05, ...){
  # Returns a logical vector indicating which values are outliers found by the specified test (either dixon.test or grubbs.test) at the specified alpha.
  # Specifically testing for outliers above the average
  
  p <- 1/Inf
  outlier_standards <- logical(length(x))
  
  while(p < alpha & sum(!outlier_standards) >= 3){
    test_out <- test(x[!outlier_standards], two.sided = FALSE)
    
    p <- test_out$p.value
    if(p < alpha){
      outlier_standards[!outlier_standards][which.max(x[!outlier_standards])] <- TRUE
    }
    
  }
  
  outlier_standards
}

calc_standard_metrics <- function(data, model_choices, input_values){
  #number_over_mean_standards = 1 #make argument later if implementing additional metrics
  #number_over_mean_standards - if there are 9 standards will fit all combinations above half + 1 standards
  
  y_var <- input_values$y_var
  
  data <- filter(data, included_in_model == 'Yes') %>%
    select(-included_in_model)
  
  n_standard <- nrow(data)
  
  ## Make all combinations of standards
  # Use this range if implementing number_over_mean_standards
  #(ceiling(n_standard/2) + number_over_mean_standards) : n_standard
  standard_combos <- map((n_standard - 1): n_standard, 
                         ~combn(data$standard_index, .x, simplify = FALSE)) %>%
    do.call(c, .) %>%
    map_dfr(~mutate(data, 
                    included_in_model = (standard_index %in% .x)),
            .id = 'standard_combination') %>%
    nest(standard_data = -c(standard_combination)) %>%
    mutate(standard_combination = as.integer(standard_combination)) %>%
    distinct(standard_data, 
             .keep_all = TRUE) 
  
  ## Fit all models to each combination of standards
  fit_standards <- standard_combos %>%
    expand_grid(model_name = model_choices) %>%
    # sample_n(10) %>%
    rowwise() %>%
    mutate(n_standards = sum(standard_data$included_in_model),
           fit_model(df_model = standard_data,
                     input_values = input_values,
                     model_name = model_name) %>%
             post_process_model(standard_data, model_name, input_values = input_values)) %>%
    ungroup %>%
    summarise(across(where(is.numeric), mean, .names = '{col}_mean'),
              # across(c(where(is.numeric), -ends_with('mean')), list),
              model_stats = tibble(model_name, across(c(where(is.numeric), -ends_with('mean')))) %>%
                list(),
              individual_models = list(tibble(model_name, fit_model, fitted)),
              .by = c(standard_combination, standard_data,
                      n_standards)) %>%
    filter(!if_any(where(is.numeric), is.infinite))
  
  ## For each standard test to see if there is a significant difference in the rmse of models which do/don't include that standard
  # mean_rmse_with.without <- fit_standards %>%
  #   select(standard_combination, standard_data,
  #          rmse_mean) %>%
  #   unnest(standard_data) %>%
  #   summarise(mean_true = mean(rmse_mean[included_in_model]),
  #             mean_false = mean(rmse_mean[!included_in_model]),
  #             # mean_true = 1 / mean(1 / rmse_mean[included_in_model]),
  #             # mean_false = 1 / mean(1 / rmse_mean[!included_in_model]),
  #             
  #             #Test is FALSE - TRUE (i.e. if greater in false then the model is better without that term.)
  #             t.test(log(rmse_mean) ~ included_in_model, alternative = 'less') %>% 
  #               broom::tidy() %>% select(statistic, parameter, p.value, alternative),
  #             .by = c(standard_index))
  
  ## ID Outliers from LOO-CV
  full_model_prediction <- filter(fit_standards, 
                                  n_standards == n_standard) %>%
    select(standard_data, ends_with('individual_models')) %>%
    unnest(individual_models) %>%
    rowwise(model_name) %>%
    reframe(tibble(standard_index = pull(standard_data, standard_index),
                   obs = pull(standard_data, y_var),
                   full_fitted = fitted))
  
  jackknife_improvements <- filter(fit_standards, n_standards == n_standard - 1) %>%
    select(standard_combination, standard_data,
           ends_with('individual_models')) %>%
    unnest(individual_models) %>%
    rowwise(standard_combination, model_name) %>%
    reframe(tibble(standard_index = pull(standard_data, standard_index),
                   obs = pull(standard_data, y_var),
                   loo_fitted = fitted,
                   included_in_model = pull(standard_data, included_in_model))) %>%
    filter(!included_in_model) %>%
    select(-included_in_model) %>%
    left_join(full_model_prediction,
              .,
              by = c('standard_index',
                     'model_name',
                     'obs')) %>%
    
    mutate(rel_improvement = sqrt((sqrt((full_fitted - obs)^2) - sqrt((loo_fitted - obs)^2))^2) / obs,
           .keep = 'unused') %>%
    mutate(outlier = identify_outliers(rel_improvement),
           .by = model_name)
  
  ## Make Output
  #standard_index; rfu; ng_per_ul; included_in_model; linear_Cook; linear_Hat; quadratic_Cook; quadratic_Hat; cubic_Cook; cubic_Hat; power_Cook; power_Hat
  out <- jackknife_improvements %>%
    select(-standard_combination) %>%
    summarise(avg_relImprovement = mean(rel_improvement),
              prop_outlier = str_c(sum(outlier), '/', n()),
              model_metrics = tibble(model_name, rel_improvement, outlier) %>%
                                     pivot_wider(names_from = model_name,
                                                 values_from = c(rel_improvement, outlier),
                                                 names_glue = '{model_name}_{.value}'),
              .by = standard_index) 
  
  out
}

parse_outlier <- function(x){
  numerator <- str_extract(x, '^[0-9]+') %>% as.integer()
  denominator <- str_extract(x, '[0-9]+$') %>% as.integer()
  
  numerator / denominator
}

parse_number_range <- function(input_string) {
  if (is.null(input_string) || length(input_string) == 0 || nchar(trimws(input_string)) == 0) {
    return(integer(0))
  }
  
  # Start with trimmed input
  cleaned_input <- trimws(input_string)
  
  # Step-by-step cleaning to handle edge cases
  # Remove trailing separators (commas, dashes, spaces)
  while (grepl("[,\\s-]$", cleaned_input)) {
    cleaned_input <- gsub("[,\\s-]+$", "", cleaned_input)
  }
  
  # Remove leading separators  
  while (grepl("^[,\\s-]", cleaned_input)) {
    cleaned_input <- gsub("^[,\\s-]+", "", cleaned_input)
  }
  
  # If empty after cleaning, return empty vector
  cleaned_input <- trimws(cleaned_input)
  if (nchar(cleaned_input) == 0) {
    return(integer(0))
  }
  
  # Split by comma first, then clean each part
  parts <- unlist(strsplit(cleaned_input, ","))
  
  # Clean each part individually
  parts <- sapply(parts, function(p) {
    p <- trimws(p)
    # Remove trailing dashes from individual parts
    p <- gsub("-+$", "", p)
    return(trimws(p))
  })
  
  # Remove empty parts
  parts <- parts[parts != "" & !is.na(parts)]
  
  if (length(parts) == 0) {
    return(integer(0))
  }
  
  # Process each part to extract numbers
  result <- c()
  for (part in parts) {
    part <- trimws(part)
    if (nchar(part) == 0) next
    
    if (grepl("-", part) && !grepl("^-", part) && !grepl("-$", part)) {
      # Handle range like "48-52" (but not "-48" or "48-")
      range_parts <- unlist(strsplit(part, "-"))
      range_parts <- trimws(range_parts)
      range_parts <- range_parts[range_parts != ""]
      
      if (length(range_parts) == 2) {
        start_num <- suppressWarnings(as.numeric(range_parts[1]))
        end_num <- suppressWarnings(as.numeric(range_parts[2]))
        
        if (!is.na(start_num) && !is.na(end_num) && start_num <= end_num) {
          result <- c(result, seq(start_num, end_num))
          next
        }
      }
      
      # If range parsing failed, try to get the first valid number
      for (rp in range_parts) {
        num_val <- suppressWarnings(as.numeric(rp))
        if (!is.na(num_val)) {
          result <- c(result, num_val)
          break
        }
      }
    } else {
      # Handle single number
      num_val <- suppressWarnings(as.numeric(part))
      if (!is.na(num_val)) {
        result <- c(result, num_val)
      }
    }
  }
  
  # Remove duplicates and convert to integer
  result <- unique(result)
  result <- result[!is.na(result)]
  return(as.integer(result))
}

#### UI ####
#### UI ####
ui <- fluidPage(
  waiter::useWaiter(),
  
  tags$head(tags$style(HTML("
    .shiny-notification-error {
      background-color: #f8d7da !important;
      border-left: 5px solid #dc3545 !important;
      color: #721c24 !important;
      border-radius: 5px;
      font-weight: bold;
      padding: 15px !important;
      margin: 10px !important;
    }
    
    .shiny-notification-message {
      background-color: #d4edda !important;
      border-left: 5px solid #28a745 !important;
      color: #155724 !important;
      border-radius: 5px;
      font-weight: bold;
      padding: 15px !important;
      margin: 10px !important;
    }
    
    /* Instructions styling */
    .instructions-box {
      background-color: #f8f9fa;
      border: 1px solid #dee2e6;
      border-radius: 5px;
      padding: 15px;
      margin-bottom: 20px;
    }
    
    .instructions-box h4 {
      margin-top: 0;
      color: #495057;
    }
    
    /* Mobile/smaller screens - instructions above sidebar */
    @media (max-width: 1199px) {
      .instructions-mobile {
        display: block;
        order: -1; /* Ensures it appears first */
      }
      
      .instructions-desktop {
        display: none;
      }
      
      .tab-content-mobile {
        display: flex;
        flex-direction: column;
      }
    }
    
    /* Desktop/larger screens - instructions in main panel */
    @media (min-width: 1200px) {
      .instructions-mobile {
        display: none;
      }
      
      .instructions-desktop {
        display: block;
      }
    }
  "))),
  
  titlePanel("Step 1: Model DNA Concentration"),
  
  div(class = "tab-content-mobile",
      sidebarLayout(
        sidebarPanel(
          # Mobile instructions for Data Input tab
          conditionalPanel(
            condition = "input.main_tab == 'Data Input'",
            div(class = "instructions-mobile instructions-box",
                h4("Instructions:"),
                tags$ul(
                  tags$li("Select the RFU Data and Plate Map Files"),
                  tags$li("Click ", tags$code("Load Data")),
                  tags$li("Confirm that the Quant Kit and Standard Rows were correctly detected"),
                  tags$li("Go to ", tags$code("Model Results"), " Tab")
                )
            )
          ),
          
          # Mobile instructions for Model Results tab
          conditionalPanel(
            condition = "input.main_tab == 'Model Results'",
            div(class = "instructions-mobile instructions-box",
                h4("Instructions:"),
                tags$ul(
                  tags$li("Evaluate how many data points are within the range of the standards, between the dashed lines (don't include 0 ng_per_well). If many data points are outside the range of the standards, then these estimates are unreliable and the quant protocol should be reevaluated."),
                  tags$li("Scrutinize the standards. Sometimes a standard is off and it should be dropped if it doesn't follow the trend of the other standards."),
                  tags$li("Scrutinize the top ranked model. Is it missing many data points? Does it look like a good fit? If something is off, consult with Sharon and Chris."),
                  tags$li("Go to ", tags$code("Finalize"), " Tab")
                )
            )
          ),
          
          # Mobile instructions for Finalize tab
          conditionalPanel(
            condition = "input.main_tab == 'Finalize'",
            div(class = "instructions-mobile instructions-box",
                h4("Instructions:"),
                tags$ul(
                  tags$li("Confirm model and standards to use"),
                  tags$li("Select ", tags$code("Generate Zip File")),
                  tags$li("Select ", tags$code("Download Zip File"))
                )
            )
          ),
          
          # Your existing sidebar content
          # Shown only on Data Input tab
          conditionalPanel(
            condition = "input.main_tab == 'Data Input'",
            fileInput("file_raw", "Upload Raw RFU Data File", accept = c(".csv", ".xls", ".xlsx")),
            fileInput("file_plate", "Upload Plate Map File", accept = c(".csv", ".xls", ".xlsx")),
            actionButton("load_data", "Load Data"),
            uiOutput("data_status"),
            br(),  
            br(),  
            selectInput("quant_kit", "Quant Kit Used",
                        choices = c("accublue-nextgen", "accublue", "accuclear"),
                        selected = "accuclear"),
            
            # Advanced settings - collapsible section
            tags$details(
              tags$summary(tags$b("Advanced Settings (click to expand/collapse)")),
              br(),
              uiOutput("x_var_selector"),
              uiOutput("y_var_selector"),
              style = "margin-top: 15px;"
            ),
            br(),  
            br(),  
            textInput("selected_standards", "Enter Standard Rows (comma-separated or range, e.g. '193-196, 199'):", "")
          ),
          # "Standards to Drop" appears on multiple tabs
          conditionalPanel(
            condition = "input.main_tab == 'Standards Preview' || input.main_tab == 'Model Results'",
            textInput("drop_standards", "Standards to Drop (comma-separated or range of standard_index, e.g. '1-3, 5'):", "")
          ),
          # Shown only on Model Results: Modelling Options and average influence option
          conditionalPanel(
            condition = "input.main_tab == 'Model Results'",
            tags$details(
              tags$summary(tags$b("Select Models and Log Transform Options (click to expand/collapse)")),
              fluidRow(
                column(6,
                       wellPanel(
                         h5("Linear Model"),
                         checkboxInput("model_linear", "Enable", TRUE),
                         checkboxInput("log_x_linear", "Log X", FALSE),
                         checkboxInput("log_y_linear", "Log Y", FALSE),
                         checkboxInput("show_fit_line_linear", "Show Best-Fit Line", TRUE),
                         checkboxInput("show_ci_linear", "Show CI", FALSE),
                         checkboxInput("show_pi_linear", "Show PI", FALSE)
                       )
                ),
                column(6,
                       wellPanel(
                         h5("Quadratic Model"),
                         checkboxInput("model_quadratic", "Enable", TRUE),
                         checkboxInput("log_x_quadratic", "Log X", FALSE),
                         checkboxInput("log_y_quadratic", "Log Y", FALSE),
                         checkboxInput("show_fit_line_quadratic", "Show Best-Fit Line", TRUE),
                         checkboxInput("show_ci_quadratic", "Show CI", FALSE),
                         checkboxInput("show_pi_quadratic", "Show PI", FALSE)
                       )
                )
              ),
              fluidRow(
                column(6,
                       wellPanel(
                         h5("Cubic Model"),
                         checkboxInput("model_cubic", "Enable", TRUE),
                         checkboxInput("log_x_cubic", "Log X", FALSE),
                         checkboxInput("log_y_cubic", "Log Y", FALSE),
                         checkboxInput("show_fit_line_cubic", "Show Best-Fit Line", TRUE),
                         checkboxInput("show_ci_cubic", "Show CI", FALSE),
                         checkboxInput("show_pi_cubic", "Show PI", FALSE)
                       )
                ),
                column(6,
                       wellPanel(
                         h5("Power Model"),
                         checkboxInput("model_power", "Enable", TRUE),
                         checkboxInput("log_x_power", "Log X", TRUE),
                         checkboxInput("log_y_power", "Log Y", TRUE),
                         checkboxInput("show_fit_line_power", "Show Best-Fit Line", TRUE),
                         checkboxInput("show_ci_power", "Show CI", FALSE),
                         checkboxInput("show_pi_power", "Show PI", FALSE)
                       )
                )
              )
            ),
            h4("Plot Scaling Options"),
            checkboxGroupInput("log_axis", "Log10 Scale Plot Axes:",
                               choices = c("X Axis" = "log_x", "Y Axis" = "log_y"),
                               selected = c("log_x", "log_y")),
            checkboxInput("overlay_samples", "Overlay fitted sample points", TRUE),
            selectInput("influence_option", "Influence Display Option:",
                        choices = c("Standard Suggestions", "Standard Stats - Average", "Standard Stats - Model"),
                        selected = "Standard Suggestions")
            
          )
        ),
        mainPanel(
          # Custom navigation header that appears above tab content on all tabs
          div(
            style = "border-bottom: 1px solid #ddd; margin-bottom: 20px; padding-bottom: 10px;",
            div(
              style = "display: flex; justify-content: space-between; align-items: center;",
              # Left side - could add additional navigation elements here if needed
              div(),
              # Right side - Return to Menu button
              tags$a(href = "http://10.5.146.65/DNA_Quantification/", 
                     target = "_blank",  # Remove this line if you want same tab
                     "Return to Menu",
                     style = paste0("color: #337ab7; text-decoration: none; font-weight: bold; ",
                                    "padding: 8px 16px; border: 1px solid #337ab7; ",
                                    "border-radius: 4px; background-color: #f8f9fa; ",
                                    "transition: background-color 0.2s;"),
                     onmouseover = "this.style.backgroundColor='#e9ecef'",
                     onmouseout = "this.style.backgroundColor='#f8f9fa'")
            )
          ),
          
          # Your existing tabsetPanel
          tabsetPanel(id = "main_tab",
                      tabPanel("Data Input", 
                               # Desktop instructions for Data Input tab
                               div(class = "instructions-desktop instructions-box",
                                   h4("Instructions:"),
                                   tags$ul(
                                     tags$li("Select the RFU Data and Plate Map Files"),
                                     tags$li("Click ", tags$code("Load Data")),
                                     tags$li("Confirm that the Quant Kit and Standard Rows were correctly detected"),
                                     tags$li("Go to ", tags$code("Model Results"), " Tab")
                                   )
                               ),
                               DTOutput("data_table")
                      ),
                      tabPanel("Model Results", 
                               # Desktop instructions for Model Results tab
                               div(class = "instructions-desktop instructions-box",
                                   h4("Instructions:"),
                                   tags$ul(
                                     tags$li("Evaluate how many data points are within the range of the standards, between the dashed lines (don't include 0 ng_per_well). If many data points are outside the range of the standards, then these estimates are unreliable and the quant protocol should be reevaluated."),
                                     tags$li("Scrutinize the standards. Sometimes a standard is off and it should be dropped if it doesn't follow the trend of the other standards."),
                                     tags$li("Scrutinize the top ranked model. Is it missing many data points? Does it look like a good fit? If something is off, consult with Sharon and Chris."),
                                     tags$li("Go to ", tags$code("Finalize"), " Tab")
                                   )
                               ),
                               tagList(
                                 plotOutput("model_grid"),
                                 DTOutput("standards_influence_table"),
                                 DTOutput("model_comp_table"),
                                 htmlOutput("model_comp_desc")
                               )
                      ),
                      tabPanel("Finalize",
                               # Desktop instructions for Finalize tab
                               div(class = "instructions-desktop instructions-box",
                                   h4("Instructions:"),
                                   tags$ul(
                                     tags$li("Confirm model and standards to use"),
                                     tags$li("Select ", tags$code("Generate Zip File")),
                                     tags$li("Select ", tags$code("Download Zip File"))
                                   )
                               ),
                               fluidRow(
                                 column(6,
                                        selectInput("chosen_models", "Choose Model(s) for Prediction:",
                                                    choices = c("linear", "quadratic", "cubic", "power"),
                                                    selected = NULL, multiple = TRUE)
                                 ),
                                 column(6,
                                        selectInput("final_standards", "Choose Standards to Use:",
                                                    choices = NULL, multiple = TRUE)
                                 )
                               ),
                               br(),
                               actionButton("generate_zip", "Generate Zip File"),
                               uiOutput("download_zip_ui")
                      )
          )
        )
      )
  )
)

#### Functionality ####
server <- function(input, output, session) {
  
  data_all <- reactiveVal(NULL)
  first_model_fit_done <- reactiveVal(FALSE)
  auto_populated_once <- reactiveVal(FALSE)
  
  # A reactive that computes the default quant kit based on the raw file's name.
  quantKitDefault <- reactive({
    req(input$file_raw)  # Ensure a file is uploaded.
    filename <- tolower(input$file_raw$name)
    if (grepl("accublue", filename) & grepl("nextgen", filename)) {
      "accublue-nextgen"
    } else if (grepl("accuclear", filename)) {
      "accuclear"
    } else if (grepl("accublue", filename) & !grepl("nextgen", filename)){
      "accublue"
    } else {
      "accuclear"  # Fallback default.
    }
  })
  
  # Whenever a new raw file is uploaded, update the select input for quant_kit.
  observeEvent(input$file_raw, {
    updateSelectInput(session, "quant_kit", selected = quantKitDefault())
  })
  
  observeEvent(input$load_data, {
    # Clear any existing notifications
    removeNotification("file_error")
    
    tryCatch({
      # Check file inputs first
      if (is.null(input$file_raw)) {
        stop("Please select a raw RFU data file")
      }
      if (is.null(input$file_plate)) {
        stop("Please select a plate map file")
      }
      
      # Reset state variables as if the app were restarted:
      first_model_fit_done(FALSE)
      auto_populated_once(FALSE)
      
      data_combined <- read_quant_data(input$file_raw$datapath, 
                                       input$file_plate$datapath, 
                                       input$quant_kit)
      data_all(data_combined)
      
      # Auto-populate "Enter Standard Rows" with rows where plate_id == 'standard'
      std_rows <- which(tolower(data_combined$plate_id) == "standard")
      if (length(std_rows) > 0) {
        default_std <- if (min(std_rows) == max(std_rows)) {
          as.character(min(std_rows))
        } else if ((max(std_rows) - min(std_rows)) == (length(std_rows) - 1)) {
          paste0(min(std_rows), "-", max(std_rows))
        } else {
          stringr::str_c(std_rows, collapse = ', ')
        }
        updateTextInput(session, "selected_standards", value = default_std)
      }
      
      showNotification(
        paste("Data loaded successfully!", nrow(data_combined), "rows loaded with", 
              sum(!is.na(data_combined$rfu)), "RFU measurements"), 
        type = "message", 
        duration = 5
      )
      
    }, error = function(e) {
      showNotification(
        paste("File Loading Error:", e$message), 
        type = "error", 
        duration = 12,
        id = "file_error"
      )
      data_all(NULL)  # Clear any existing data
    })
  })
  
  
  output$data_table <- renderDT({
    req(data_all())
    data_all() %>%
      mutate(Row = row_number(),
             .before = 1) %>%
      datatable(options = list(pageLength = 10),
                rownames = FALSE)
  })
  
  output$data_status <- renderUI({
    if (!is.null(data_all()) && nrow(data_all()) > 0) {
      tags$div(
        style = "color: #28a745; font-weight: bold; padding: 10px; margin: 5px 0;",
        icon("check-circle"), 
        paste("Data loaded:", nrow(data_all()), "rows")
      )
    } else {
      tags$div(
        style = "color: #6c757d; font-style: italic; padding: 10px; margin: 5px 0;",
        icon("info-circle"), 
        "No data loaded. Please select files and click 'Load Data'."
      )
    }
  })
  
  output$x_var_selector <- renderUI({
    req(data_all())
    numeric_cols <- str_subset(names(data_all())[sapply(data_all(), is.numeric)], 
                               'replicate|quant|volume',
                               negate = TRUE)
    
    # Add a placeholder option with an empty string as value.
    selectInput("x_var", "Select X Variable (Independent Variable):", 
                choices = c("Please choose" = "", numeric_cols), 
                selected = str_subset(numeric_cols, '[rR][fF][uU]'))
  })
  
  output$y_var_selector <- renderUI({
    req(data_all())
    numeric_cols <- str_subset(names(data_all())[sapply(data_all(), is.numeric)], 
                               'replicate|quant|volume',
                               negate = TRUE)
    
    selectInput("y_var", "Select Y Variable (Dependent Variable):", 
                choices = c("Please choose" = "", numeric_cols), 
                selected = str_subset(numeric_cols, '_per_well'))
  })
  
  # Compute standards based on "selected_standards" and remove those specified in "drop_standards"
  standards_data <- reactive({
    req(data_all(), input$selected_standards)
    req(input$x_var != "", input$y_var != "")
    
    # Parse selected standards
    expanded_rows <- parse_number_range(input$selected_standards)
    
    # Check which parsed numbers are actually valid
    valid_mask <- expanded_rows %in% seq_len(nrow(data_all()))
    
    # Filter to valid row indices
    expanded_rows <- expanded_rows[valid_mask]
    
    # If no valid rows, return empty data frame with correct structure
    if (length(expanded_rows) == 0) {
      
      # Show what rows DO exist that contain standards
      if (!is.null(data_all()) && nrow(data_all()) > 0) {
        std_rows <- which(tolower(data_all()$plate_id) == "standard")
        
        if (length(std_rows) > 0) {
          showNotification(
            paste("Row", paste(expanded_rows[!valid_mask], collapse = ", "), 
                  "does not exist in your dataset. Available standard rows are:", 
                  paste(std_rows, collapse = ", ")), 
            type = "warning", duration = 10
          )
        }
      }
      
      empty_df <- data.frame(
        standard_index = integer(0),
        stringsAsFactors = FALSE
      )
      empty_df[[input$x_var]] <- numeric(0)
      empty_df[[input$y_var]] <- numeric(0)
      empty_df$included_in_model <- character(0)
      return(empty_df)
    }
    
    df <- data_all()[expanded_rows, , drop = FALSE]
    
    # Parse standards to drop
    drop_indices <- parse_number_range(input$drop_standards)
    drop_indices <- drop_indices[drop_indices %in% seq_len(nrow(df))]
    
    df <- df %>%
      mutate(standard_index = row_number(),
             included_in_model = if_else(standard_index %in% drop_indices, "No", "Yes")
      ) %>%
      select(standard_index, all_of(c(input$x_var, input$y_var)), included_in_model)
    
    return(df)
  })
  
  # Create a reactiveVal to track if the drop_standards field has been manually modified.
  drop_standards_modified <- reactiveVal(FALSE)
  
  # Observer to detect user modifications to drop_standards.
  observeEvent(input$drop_standards, {
    # If the user types anything (i.e. non-empty), mark as modified.
    if (nchar(input$drop_standards) > 0) {
      drop_standards_modified(TRUE)
    } else {
      drop_standards_modified(FALSE)
    }
  }, ignoreInit = TRUE)
  
  # Create a flag to ensure auto-population happens only once.
  auto_populated_once <- reactiveVal(FALSE)
  
  observe({
    req(standards_data(), input$y_var)
    # Only auto-populate if it hasn't been done yet
    if (!auto_populated_once()) {
      df_standards <- standards_data()
      # Identify rows (using the standards table indices) where the selected y variable equals 0
      zero_rows <- which(df_standards[[input$y_var]] == 0)
      if (length(zero_rows) > 0) {
        default_drop <- paste(zero_rows, collapse = ", ")
        updateTextInput(session, "drop_standards", value = default_drop)
      }
      auto_populated_once(TRUE)  # Mark that auto-population has occurred
    }
  })
  
  # Update "final_standards" choices in the Finalize tab.
  observe({
    req(standards_data())
    updateSelectInput(session, "final_standards",
                      choices = standards_data()$standard_index,
                      selected = standards_data() %>% filter(included_in_model == "Yes") %>% pull(standard_index))
  })
  
  # ---------------------------------------------------------------------------
  # Model fitting with prediction data moved into this reactive:
  # ---------------------------------------------------------------------------
  model_fit <- reactive({
    req(standards_data(), input$x_var != "", input$y_var != "")
    req(auto_populated_once() == TRUE)
    
    # On the first run, require that drop_standards is not empty.
    if (!first_model_fit_done()) {
      req(input$drop_standards)
    }
    
    df_model <- standards_data() %>% filter(included_in_model == "Yes")
    df_model <- as.data.frame(df_model)
    if(nrow(df_model) == 0) {
      showNotification("No standards available for model fitting.", type = "error")
      return(list())
    }
    rownames(df_model) <- df_model$standard_index
    results <- list()
    
    x_var <- input$x_var
    y_var <- input$y_var
    
    x_range_min <- min(df_model[[x_var]], na.rm = TRUE)
    x_range_max <- max(as.data.frame(standards_data())[[x_var]], na.rm = TRUE)
    grid_x <- seq(x_range_min, x_range_max, length.out = 1000)
    base_grid <- tibble(!!sym(x_var) := grid_x)
    
    # Compute sample data as rows NOT in the standards (based on selected_standards)
    selected_rows <- unlist(strsplit(input$selected_standards, ","))
    selected_rows <- trimws(selected_rows)
    expanded_rows <- unlist(lapply(selected_rows, function(x) {
      if (grepl("-", x)) {
        range_vals <- as.numeric(unlist(strsplit(x, "-")))
        if(length(range_vals)==2) return(seq(range_vals[1], range_vals[2]))
      } else return(as.numeric(x))
    }))
    expanded_rows <- expanded_rows[expanded_rows %in% seq_len(nrow(data_all()))]
    all_rows <- seq_len(nrow(data_all()))
    sample_rows <- setdiff(all_rows, expanded_rows)
    sample_data <- data_all()[sample_rows, ]
    # message("DEBUG: Models: ", str_c(str_subset(names(input), 'model_'), collapse = '; '))
    for(model_name in c("linear", "quadratic", "cubic", "power")) {
      if(input[[paste0("model_", model_name)]]) {
        log_x <- input[[paste0("log_x_", model_name)]]
        log_y <- input[[paste0("log_y_", model_name)]]
        x_term <- ifelse(log_x, paste0("log10(", x_var, ")"), x_var)
        y_term <- ifelse(log_y, paste0("log10(", y_var, ")"), y_var)
       
        # if(model_name == "linear") {
        #   formula <- as.formula(paste(y_term, "~", x_term))
        # } else if(model_name == "quadratic") {
        #   formula <- as.formula(paste(y_term, "~ poly(", x_term, ", 2)"))
        # } else if(model_name == "cubic") {
        #   formula <- as.formula(paste(y_term, "~ poly(", x_term, ", 3)"))
        # } else if(model_name == "power") {
        #   formula <- as.formula(paste("log10(", y_var, ") ~ log10(", x_var, ")"))
        # }
        # model_obj <- tryCatch(
        #   lm(formula, data = df_model),
        #   error = function(e) {
        #     showNotification(paste("Error fitting", model_name, "model:", e$message), type = "error")
        #     NULL
        #   }
        # )
        
        model_obj <- fit_model(df_model = df_model, model_name = model_name, input_values = input)
        
        if (!is.null(model_obj)) {
          # Compute grid predictions
          pred_int <- predict(model_obj, newdata = base_grid, interval = "prediction")
          conf_int <- predict(model_obj, newdata = base_grid, interval = "confidence")
          if(log_y) {
            pred_int <- 10^(pred_int)
            conf_int <- 10^(conf_int)
          }
          grid_df <- base_grid %>% mutate(
            predicted_y = pred_int[,1],
            ci_lower = pred_int[,2],
            ci_upper = pred_int[,3],
            pi_lower = conf_int[,2],
            pi_upper = conf_int[,3]
          )
          
          # Compute sample predictions (for all samples not in standards)
          sample_pred <- predict(model_obj, newdata = sample_data)
          if(log_y) {
            sample_pred <- 10^(sample_pred)
          }
          sample_df <- sample_data %>% mutate(predicted = sample_pred)
          
          results[[model_name]] <- list(
            fit = model_obj,
            grid = grid_df,
            sample = sample_df
          )
        }
      }
    }
    first_model_fit_done(TRUE)
    return(results)
  })
  
  model_metrics <- reactive({
    req(model_fit(), standards_data(), 
        input$x_var != '', input$y_var != '')
    
    # write_rds(model_fit(), 'model_fit.rds')
    # write_rds(standards_data(), 'standards_data.rds')
    # write_rds(input, 'input.rds')
    
    df_std <- standards_data()
    req(nrow(df_std) > 0)
    
    # Check if model_fit() returned any models
    models <- model_fit()
    if (length(models) == 0) {
      # Return empty data frame with correct structure
      return(data.frame(
        Model = character(0),
        Adj_R2 = numeric(0),
        AICc = numeric(0),
        RMSE = numeric(0),
        RSE = numeric(0),
        LOOCV = numeric(0),
        RANK = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    
    df <- df_std %>% filter(included_in_model == "Yes")
    n <- nrow(df)
    
    # Only process models that actually exist and have valid fits
    valid_models <- names(models)[sapply(models, function(x) !is.null(x$fit))]
    
    if (length(valid_models) == 0) {
      # Return empty data frame if no valid models
      return(data.frame(
        Model = character(0),
        Adj_R2 = numeric(0),
        AICc = numeric(0),
        RMSE = numeric(0),
        RSE = numeric(0),
        LOOCV = numeric(0),
        RANK = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    
    metrics <- lapply(valid_models, function(model_name) {
      model <- models[[model_name]]$fit
      if (is.null(model)) return(NULL)
      
      tryCatch({
        summ <- summary(model)
        r2 <- summ$r.squared
        adj_r2 <- summ$adj.r.squared
        
        aic_val <- AIC(model)
        bic_val <- BIC(model)
        
        log_y_flag <- input[[paste0("log_y_", model_name)]]
        observed <- df[[input$y_var]]
        
        # Add the Jacobian adjustment for models with log10 transformed y.
        # https://stats.stackexchange.com/questions/61332/comparing-aic-of-a-model-and-its-log-transformed-version
        if (log_y_flag) {
          adjustment <- sum(2 * log(observed)) / log(10)
          aic_val <- aic_val + adjustment
        }
        
        k <- length(coef(model))
        aicc_val <- aic_val + (2 * k * (k + 1)) / (n - k - 1)
        
        fitted_vals <- if (log_y_flag) { 10^(fitted(model)) } else { fitted(model) }
        rmse_val <- sqrt(mean((observed - fitted_vals)^2))
        rse_val <- summary(model)$sigma
        loocv_val <- loocv(model)
        
        data.frame(Model = toupper(model_name),
                   R2 = round(r2, 4),
                   Adj_R2 = round(adj_r2, 4),
                   AIC = round(aic_val, 4),
                   AICc = round(aicc_val, 4),
                   BIC = round(bic_val, 4),
                   RMSE = round(rmse_val, 4),
                   RSE = round(rse_val, 4),
                   LOOCV = round(loocv_val, 4),
                   stringsAsFactors = FALSE)
      }, error = function(e) {
        # If any error occurs in metric calculation, return NULL
        showNotification(paste("Error calculating metrics for", model_name, ":", e$message), 
                         type = "warning", duration = 3)
        return(NULL)
      })
    })
    
    # Remove NULL entries
    metrics <- metrics[!sapply(metrics, is.null)]
    
    if (length(metrics) == 0) {
      # Return empty data frame if no valid metrics
      return(data.frame(
        Model = character(0),
        Adj_R2 = numeric(0),
        AICc = numeric(0),
        RMSE = numeric(0),
        RSE = numeric(0),
        LOOCV = numeric(0),
        RANK = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    
    # Safely combine metrics
    combined_metrics <- tryCatch({
      do.call(rbind, metrics)
    }, error = function(e) {
      showNotification(paste("Error combining model metrics:", e$message), 
                       type = "error", duration = 5)
      return(data.frame(
        Model = character(0),
        Adj_R2 = numeric(0),
        AICc = numeric(0),
        RMSE = numeric(0),
        RSE = numeric(0),
        LOOCV = numeric(0),
        RANK = numeric(0),
        stringsAsFactors = FALSE
      ))
    })
    
    if (is.null(combined_metrics) || nrow(combined_metrics) == 0) {
      return(data.frame(
        Model = character(0),
        Adj_R2 = numeric(0),
        AICc = numeric(0),
        RMSE = numeric(0),
        RSE = numeric(0),
        LOOCV = numeric(0),
        RANK = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    
    # Safely process the combined metrics
    result <- tryCatch({
      combined_metrics %>%
        select(-R2, -BIC, -AIC) %>%
        mutate(RANK = (dense_rank(AICc) +
                         dense_rank(LOOCV) +
                         dense_rank(RMSE) +
                         dense_rank(RSE))) %>%
        mutate(RANK = dense_rank(RANK))
    }, error = function(e) {
      showNotification(paste("Error processing model metrics:", e$message), 
                       type = "error", duration = 5)
      return(data.frame(
        Model = character(0),
        Adj_R2 = numeric(0),
        AICc = numeric(0),
        RMSE = numeric(0),
        RSE = numeric(0),
        LOOCV = numeric(0),
        RANK = numeric(0),
        stringsAsFactors = FALSE
      ))
    })
    
    return(result)
  })
  
  output$model_comp_table <- renderDT({
    req(model_metrics())
    dt <- model_metrics()
    # best_R2    <- max(dt$R2, na.rm = TRUE)
    best_Adj_R2<- max(dt$Adj_R2, na.rm = TRUE)
    # best_AIC   <- min(dt$AIC, na.rm = TRUE)
    best_AICc  <- min(dt$AICc, na.rm = TRUE)
    # best_BIC   <- min(dt$BIC, na.rm = TRUE)
    best_RMSE  <- min(dt$RMSE, na.rm = TRUE)
    best_RSE  <- min(dt$RSE, na.rm = TRUE)
    best_RANK  <- min(dt$RANK, na.rm = TRUE)
    best_LOOCV <- min(dt$LOOCV, na.rm = TRUE)
    
    datatable(arrange(dt, RANK), options = list(pageLength = 5, orderClasses = TRUE)) %>%
      # formatStyle('R2', backgroundColor = styleEqual(best_R2, 'lightgreen')) %>%
      formatStyle('Adj_R2', backgroundColor = styleEqual(best_Adj_R2, 'lightgreen')) %>%
      # formatStyle('AIC', backgroundColor = styleEqual(best_AIC, 'lightgreen')) %>%
      formatStyle('AICc', backgroundColor = styleEqual(best_AICc, 'lightgreen')) %>%
      # formatStyle('BIC', backgroundColor = styleEqual(best_BIC, 'lightgreen')) %>%
      formatStyle('RMSE', backgroundColor = styleEqual(best_RMSE, 'lightgreen')) %>%
      formatStyle('RSE', backgroundColor = styleEqual(best_RSE, 'lightgreen')) %>%
      formatStyle('RANK', backgroundColor = styleEqual(best_RANK, 'lightgreen')) %>%
      formatStyle('LOOCV', backgroundColor = styleEqual(best_LOOCV, 'lightgreen'))
  })
  
  output$model_comp_desc <- renderUI({
    HTML(paste0(
      "<b>Model Comparison Metrics:</b><br/>",
      "<ul>",
      # "<li><b>R2</b>: Proportion of variance explained (higher is better).</li>",
      "<li><b>Adj_R2</b>: Adjusted R2, penalizes for number of predictors (higher is better).</li>",
      "<li><b>AIC</b>: Akaike Information Criterion (lower is better).</li>",
      "<li><b>AICc</b>: AIC corrected for small sample sizes (lower is better).</li>",
      # "<li><b>BIC</b>: Bayesian Information Criterion (lower is better).</li>",
      "<li><b>RMSE</b>: Root Mean Squared Error (lower is better).</li>",
      "<li><b>RSE</b>: Residual Standard Error (lower is better).</li>",
      "<li><b>LOOCV</b>: Leave-one-out Error (lower is better, higher means more overfit).</li>",
      "</ul>",
      "Default model selection is based on the best AICc, RMSE, RSE, and LOOCV. If multiple models are tied then predictions are averaged"
    ))
  })
  
  # Set default chosen_models based on lowest overall Ranking.
  observe({
    req(model_metrics())
    dt <- model_metrics()
    min_rank <- min(dt$RANK, na.rm = TRUE)
    default_models <- tolower(dt$Model[dt$RANK == min_rank])
    updateSelectInput(session, "chosen_models", selected = default_models)
  })
  
  # Compute influence values: Cook's D and hat values.
  influence_table <- reactive({
    req(standards_data(), model_fit())
    
    # Get the base standards data.
    base_df <- standards_data()
    
    # Calculate model-specific influence metrics.
    for(model_name in names(model_fit())) {
      model <- model_fit()[[model_name]]$fit
      cd <- cooks.distance(model)
      hv <- hatvalues(model)
      base_df[[paste0(model_name, "_Cook")]] <- ifelse(base_df$included_in_model == "Yes",
                                                       round(cd[as.character(base_df$standard_index)], 4),
                                                       NA)
      base_df[[paste0(model_name, "_Hat")]] <- ifelse(base_df$included_in_model == "Yes",
                                                      round(hv[as.character(base_df$standard_index)], 4),
                                                      NA)
    }
    
    #Calculate LOO Improvement from removing standards
    # message("DEBUG: standards columns: ", str_c(colnames(base_df), collapse = '; '))
    
    loo_metrics <- calc_standard_metrics(base_df, names(model_fit()), input)
    
    # Choose the display option based on the dropdown.
    option <- input$influence_option
    if(option == "Standard Suggestions") {
      # Just show the standards data.
      base_df <- standards_data()
      base_df <- base_df %>%
        left_join(select(loo_metrics, standard_index, prop_outlier),
                  by = 'standard_index')
      return(base_df)
    } else if(option == "Standard Stats - Average") {
      # Compute average influence metrics across models.
      cook_matrix <- sapply(names(model_fit()), function(model_name) {
        as.numeric(base_df[[paste0(model_name, "_Cook")]])
      })
      hat_matrix <- sapply(names(model_fit()), function(model_name) {
        as.numeric(base_df[[paste0(model_name, "_Hat")]])
      })
      base_df$Avg_Cook <- round(rowMeans(cook_matrix, na.rm = TRUE), 4)
      base_df$Avg_Hat <- round(rowMeans(hat_matrix, na.rm = TRUE), 4)
      
      # Remove the individual model columns.
      cols_to_remove <- c(sapply(names(model_fit()), function(model_name) paste0(model_name, "_Cook")),
                          sapply(names(model_fit()), function(model_name) paste0(model_name, "_Hat")))
      base_df <- base_df[, !(names(base_df) %in% cols_to_remove)]
      
      base_df <- left_join(base_df, 
                           select(loo_metrics, -model_metrics),
                           by = 'standard_index') %>%
        mutate(across(where(is.numeric), ~round(., 4))) %>%
        relocate(.after = last_col(), ends_with('Cook'), ends_with('Hat')) %>%
        select(-ends_with('Cook'), -ends_with('Hat'))
      
      return(base_df)
    } else if(option == "Standard Stats - Model") {
      # Show model-specific influence values.
      base_df <- left_join(base_df, 
                select(loo_metrics, standard_index, model_metrics) %>%
                  unnest(model_metrics),
                by = 'standard_index') %>%
        mutate(across(where(is.numeric), ~round(., 4))) %>%
        select(-ends_with('Cook'), -ends_with('Hat'))
      return(base_df)
    } else {
      # Fallback: just return the standards data.
      return(standards_data())
    }
  })
  
  output$standards_influence_table <- renderDT({
    req(influence_table())
    option <- input$influence_option
    
    if(option == "Standard Suggestions"){
      out_tab <- influence_table() %>%
        mutate(prop_outlier = parse_outlier(prop_outlier))
      
      message('DEBUG: ', str_c(colnames(out_tab), collapse = '; '))
      datatable(out_tab, 
                options = list(pageLength = 10,
                               columnDefs = list(list(visible = FALSE, targets = which(colnames(out_tab) == 'prop_outlier') - 1))),
                caption = "Standards Influence on Models (Cook's D and Hat Values)") %>%
        formatStyle(
          "prop_outlier",
          target = "row",
          # Use styleInterval to assign colors:
          # - "lightgreen" for scores up to 0.4 (good),
          # - "yellow" for scores between 0.4 and 0.7 (moderate),
          # - "tomato" for scores above 0.7 (bad)
          backgroundColor = styleInterval(c(0.4, 0.7), c("lightgreen", "yellow", "tomato"))
        )
      #%>%
        # formatStyle('LOOCV', backgroundColor = styleEqual(best_LOOCV, 'lightgreen'))
    } else {
      datatable(influence_table(), options = list(pageLength = 10),
                caption = "Standards Influence on Models (Cook's D and Hat Values)")
    }

  })
  
  model_grid_plot <- reactive({
    req(model_fit(), standards_data(), data_all())
    results <- model_fit()
    df_std <- standards_data()
    x_var <- input$x_var
    y_var <- input$y_var
    
    log_x_plot <- "log_x" %in% input$log_axis  
    log_y_plot <- "log_y" %in% input$log_axis
    
    range_rfu <- filter(df_std, included_in_model == 'Yes') %>%
      pull(!!sym(x_var)) %>%
      range()
    
    plots <- list()
    
    for (model_name in names(results)) {
      log_x_model <- input[[paste0("log_x_", model_name)]]
      log_y_model <- input[[paste0("log_y_", model_name)]]
      
      x_label <- ifelse(log_x_plot, paste0("log10(", x_var, ")"), x_var)
      y_label <- ifelse(log_y_plot, paste0("log10(", y_var, ")"), y_var)
      
      grid_df <- results[[model_name]]$grid
      
      p <- ggplot() +
        (if (input[[paste0("show_ci_", model_name)]]) {
          grid_df_ci <- if (log_y_plot) {
            filter(grid_df, ci_lower > 0)
          } else {
            grid_df
          }
          geom_ribbon(data = grid_df_ci,
                      aes(x = !!sym(x_var), ymin = ci_lower, ymax = ci_upper),
                      fill = "gray40", alpha = 0.5, inherit.aes = FALSE)
        } else { 
          NULL 
        }) +
        (if (input[[paste0("show_pi_", model_name)]]) {
          grid_df_pi <- if (log_y_plot) {
            filter(grid_df, pi_lower > 0)
          } else {
            grid_df
          }
          geom_ribbon(data = grid_df_pi,
                      aes(x = !!sym(x_var), ymin = pi_lower, ymax = pi_upper),
                      fill = "gray40", alpha = 0.4, inherit.aes = FALSE)
        } else { 
          NULL 
        }) +
        (if (input[[paste0("show_fit_line_", model_name)]]) {
          geom_line(data = grid_df,
                    aes(x = !!sym(x_var), y = predicted_y),
                    color = "black", linewidth = 1)
        } else { 
          NULL 
        }) +
        geom_point(data = df_std,
                   aes(x = !!sym(x_var), y = !!sym(y_var), color = included_in_model),
                   size = 2, alpha = 0.25) +
        geom_text(data = df_std,
                  aes(x = !!sym(x_var), y = !!sym(y_var), color = included_in_model,
                      label = standard_index),
                  size = 3, hjust = 'inward', vjust = 'inward',
                  show.legend = FALSE) +
        geom_vline(xintercept = range_rfu,
                   linetype = 2) +
        annotate("text",
                 x = range_rfu,
                 y = c(Inf, 0),
                 label = c("below lod", 'above lod'),
                 angle = c(90, 270),
                 size = 3,
                 hjust = 'inward',
                 vjust = 'inward') +
        scale_color_manual(values = c("Yes" = "blue", "No" = "red")) +
        theme_minimal() +
        labs(title = paste(toupper(model_name), "Model Fit"),
             x = x_label,
             y = y_label,
             color = "Included in Model")
      
      if (input$overlay_samples) {
        sample_df <- results[[model_name]]$sample
        less_lod <- mean(pull(sample_df, !!sym(x_var)) < range_rfu[1])
        more_lod <- mean(pull(sample_df, !!sym(x_var)) > range_rfu[2])
        p <- p + geom_point(data = sample_df,
                            aes(x = !!sym(x_var), y = predicted),
                            shape = 16, color = "black", size = 2, alpha = 0.5) +
          labs(caption = stringr::str_c(scales::percent(less_lod), ' samples below LoD\n',
                                        scales::percent(more_lod), ' samples above LoD'))
      }
      
      if (log_x_plot) {
        p <- p + scale_x_continuous(trans = scales::log10_trans(), 
                                    labels = scales::comma_format())
      } else {
        p <- p + scale_x_continuous(labels = scales::comma_format())
      }
      
      if (log_y_plot) {
        p <- p + scale_y_continuous(trans = scales::log10_trans(), 
                                    labels = scales::comma_format())
      } else {
        p <- p + scale_y_continuous(labels = scales::comma_format())
      }
      
      plots[[model_name]] <- p 
    }
    
    updated_y_range <- purrr::map(plots, ~layer_scales(.x)$y$range$range) %>%
      unlist() %>%
      range() 
    
    if(log_y_plot){
      updated_y_range <- 10^updated_y_range
    }
    
    updated_x_range <- purrr::map(plots, ~layer_scales(.x)$x$range$range) %>%
      unlist() %>%
      range()
    if(log_x_plot){
      updated_x_range <- 10^updated_x_range
    }
    
    
    
    combined_plot <- wrap_plots(plots, ncol = 2, 
                                guides = 'collect', 
                                axes = 'collect', 
                                axis_titles = 'collect') &
      coord_cartesian(ylim = updated_y_range,
                      xlim = updated_x_range)
    combined_plot
  })
  
  output$model_grid <- renderPlot({
    req(model_grid_plot())
    model_grid_plot()
  })
  
  # Download handler for CSV output using sample predictions from model_fit():
  output$download_csv <- downloadHandler(
    filename = function() {
      req(input$file_raw, input$file_plate)
      prefix <- common_prefix(input$file_raw$name, input$file_plate$name)
      prefix <- tolower(prefix)
      paste0(prefix, "_sample_predictions_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(data_all(), input$chosen_models, input$selected_standards)
      chosen_models <- input$chosen_models
      models_results <- model_fit()
      if(length(chosen_models) == 0) {
        stop("No model chosen for prediction.")
      }
      
      pred_list <- lapply(chosen_models, function(model_name) {
        if(is.null(models_results[[model_name]])) {
          stop(paste("Model", model_name, "is not available."))
        }
        sample_df <- models_results[[model_name]]$sample
        return(sample_df$predicted)
      })
      pred_matrix <- do.call(cbind, pred_list)
      avg_pred <- rowMeans(pred_matrix)
      
      # Use sample_data from the first chosen model (the sample rows are identical)
      sample_data <- models_results[[chosen_models[1]]]$sample
      sample_data$predicted_y <- avg_pred
      
      # Remove columns that are entirely NA.
      sample_data <- sample_data[, colSums(is.na(sample_data)) < nrow(sample_data)]
      # Rename "predicted_y" to the name of the y variable.
      names(sample_data)[names(sample_data) == "predicted_y"] <- input$y_var
      
      # write_rds(sample_data, 'tmp.rds')
      sample_data <- select(sample_data,
                            starts_with('plate'),
                            starts_with('quant'),
                            sample_volume,
                            all_of(c(input$x_var, 
                                     input$y_var))) %>%
        rename_with(~str_c('dna_', .), 
                    .cols = starts_with('plate')) %>%
        mutate(!!sym(str_c(str_extract(input$y_var, '^[np]g'), '_per_ul')) := !!sym(input$y_var) / sample_volume,
               .keep = 'unused')
      
      write_csv(sample_data, file)
    },
    contentType = "text/csv"
  )
  
  # Download handler for JSON output.
  output$download_json <- downloadHandler(
    filename = function() {
      req(input$file_raw, input$file_plate)
      prefix <- common_prefix(input$file_raw$name, input$file_plate$name)
      prefix <- tolower(prefix)
      paste0(prefix, "_model_details_", Sys.Date(), ".json")
    },
    content = function(file) {
      req(input$chosen_models)
      model_details <- list(
        chosen_models = input$chosen_models,
        standards_selected = input$final_standards
      )
      json_text <- toJSON(model_details, pretty = TRUE, auto_unbox = TRUE)
      writeLines(json_text, con = file)
    },
    contentType = "application/json"
  )
  
  
  # Reactive value to store the generated ZIP file's path.
  zip_info <- reactiveVal(NULL)
  
  # Create a waiter object (you can customize the spinner as needed)
  w <- waiter::Waiter$new(html = waiter::spin_fading_circles(), color = "rgba(255,255,255,0.8)")
  
  # When the user clicks "Generate Zip File":
  observeEvent(input$generate_zip, {
    # Show the spinner
    w$show()
    
    # Ensure required reactives are available
    req(data_all(), model_fit(), model_grid_plot(), input$chosen_models, input$selected_standards)
    
    # Generate a prefix from file names (your existing helper)
    prefix <- tolower(common_prefix(input$file_raw$name, input$file_plate$name))
    
    # Create a temporary directory and file paths
    tmpdir <- tempdir()
    tmpdir <- file.path("outdir", basename(tmpdir))
    dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
    
    zip_path <- file.path(paste0(prefix, "_output_", Sys.Date(), ".zip"))
    prediction_path <- file.path(tmpdir, paste0(prefix, "_sample_predictions_", Sys.Date(), ".csv"))
    inputs_path <- file.path(tmpdir, paste0(prefix, "_inputs_", Sys.Date(), ".txt"))
    json_path <- file.path(tmpdir, paste0(prefix, "_settings_", Sys.Date(), ".json"))
    model_path <- file.path(tmpdir, paste0(prefix, "_models_", Sys.Date(), ".rds"))
    code_path <- file.path(tmpdir, paste0(prefix, "_sourceCode_", Sys.Date(), ".R"))
    plots_path <- file.path(tmpdir, paste0(prefix, "_plots_", Sys.Date(), ".png"))
    
    #### Output predictions ####
    chosen_models <- input$chosen_models
    models_results <- model_fit()
    if(length(chosen_models) == 0) {
      stop("No model chosen for prediction.")
    }
    
    pred_list <- lapply(chosen_models, function(model_name) {
      if(is.null(models_results[[model_name]])) {
        stop(paste("Model", model_name, "is not available."))
      }
      sample_df <- models_results[[model_name]]$sample
      return(sample_df$predicted)
    })
    pred_matrix <- do.call(cbind, pred_list)
    avg_pred <- rowMeans(pred_matrix)
    
    # Use sample_data from the first chosen model (the sample rows are identical)
    sample_data <- models_results[[chosen_models[1]]]$sample
    sample_data$predicted_y <- avg_pred
    
    # Remove columns that are entirely NA.
    sample_data <- sample_data[, colSums(is.na(sample_data)) < nrow(sample_data)]
    # Rename "predicted_y" to the name of the y variable.
    names(sample_data)[names(sample_data) == "predicted_y"] <- input$y_var
    
    # write_rds(sample_data, 'tmp.rds')
    sample_data <- select(sample_data,
                          starts_with('plate'),
                          starts_with('quant'),
                          sample_volume,
                          all_of(c(input$x_var, 
                                   input$y_var))) %>%
      rename_with(~str_c('dna_', .), 
                  .cols = starts_with('plate')) %>%
      mutate(!!sym(str_c(str_extract(input$y_var, '^[np]g'), '_per_ul')) := !!sym(input$y_var) / sample_volume,
             .keep = 'unused')
    
    write_csv(sample_data, prediction_path)
    
    #### Output inputs ####
    input_list <- lapply(names(input), function(name) {
      paste(name, ":", toString(input[[name]]))
    })
    write_lines(input_list, inputs_path)
    
    #### Output JSON ####
    model_details <- list(
      chosen_models = input$chosen_models,
      standards_selected = input$final_standards
    )
    json_text <- toJSON(model_details, pretty = TRUE, auto_unbox = TRUE)
    write_lines(json_text, json_path)
    
    #### Output Models ####
    write_rds(model_fit()[[input$chosen_models]]$fit, 
              model_path)
    
    #### Output Code ####
    source_code_file <- "app.R"
    if (file.exists(source_code_file)) {
      file.copy(source_code_file, code_path, overwrite = TRUE)
    } else {
      # If the source file is not found, write a placeholder message.
      write_lines("# Source code file not found. The app may be running in an interactive session.",
                  code_path)
    }
    
    
    #### Output Plots ####
    # Save the combined plot as a PNG file
    ggsave(filename = plots_path, plot = model_grid_plot(), 
           device = "png", width = 10, height = 8)
    
    
    #### Zip contents ####
    files_to_zip <- c(code_path, 
                      inputs_path,
                      json_path,
                      model_path,
                      prediction_path,
                      plots_path)
    
    
    # Change working directory to the temporary directory
    message('DEBUG: Original WorkingDIR: ', getwd())
    owd <- setwd(tmpdir)
    on.exit(setwd(owd), add = TRUE)
    message('DEBUG: New WorkingDIR: ', getwd())
    
    # Create the zip file (using base names for a cleaner archive)
    message('DEBUG: files_to_zip: ', 
            str_c(files_to_zip, collapse = '; '))
    message('DEBUG: Files to zip existance: ',
            str_c(file.exists(basename(files_to_zip)), 
                  collapse = '; '))
    
    zip(zipfile = zip_path, files = basename(files_to_zip))
    
    # zip_file(file.path(tmpdir, zip_path))
    message('DEBUG: zip: ', zip_path)
    message('DEBUG: zip existance: ', file.exists(zip_path))
    # Hide the spinner once done
    
    zip_info(list(zip = file.path(tmpdir, zip_path), tmpdir = tmpdir))
    w$hide()
  })
  
  # Render a download button once the ZIP file is generated
  output$download_zip_ui <- renderUI({
    req(zip_info())
    message('DEBUG: WD: ', getwd())
    message('DEBUG: zip: ', zip_info()$zip)
    downloadButton("download_zip", "Download Zip File")
  })
  
  # Download handler that serves the already generated ZIP file
  output$download_zip <- downloadHandler(
    filename = function() {
      req(zip_info())
      basename(zip_info()$zip)
    },
    content = function(file) {
      file.copy(zip_info()$zip, file)
      # Delete the temporary directory and its contents after download
      unlink(zip_info()$tmpdir, recursive = TRUE)
      # Reset the reactive value
      zip_info(NULL)
    },
    contentType = "application/zip"
  )
  
}

shinyApp(ui, server)

