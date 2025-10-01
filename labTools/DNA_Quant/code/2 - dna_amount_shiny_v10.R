#Run on: http://10.5.146.65/DNA_Quantification/
#Copy as "app.R" to /srv/shiny-server/DNA_Quantification/2-DNA_concentration
#Logs go here: /var/log/shiny-server/

##TODO - init with VB
##TODO - use gamlss? - or make bayes vs freq a choice? freq for empirical priors?

library(shiny) |> suppressMessages() |> suppressWarnings()
library(shinyBS) |> suppressMessages() |> suppressWarnings()         # For collapsible panels
library(shinyFiles)  |> suppressMessages() |> suppressWarnings()     # For file saving dialogs on the server side
library(readxl) |> suppressMessages() |> suppressWarnings()
library(writexl) |> suppressMessages() |> suppressWarnings()
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
library(waiter) |> suppressMessages() |> suppressWarnings() 
library(jsonlite) |> suppressMessages() |> suppressWarnings() 
library(DT) |> suppressMessages() |> suppressWarnings()

if(Sys.info()["nodename"] == 'gawain'){
  set_cmdstan_path('/home/shiny/.cmdstan/cmdstan-2.36.0')
}

#### Helper Functions ####
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
  tryCatch({
    read_xlsx(filepath,
              sheet = sheet_name, 
              na = c('na', 'NA', 'N/A'), 
              .name_repair = ~ vctrs::vec_as_names(..., 
                                                   repair = "unique",
                                                   quiet = TRUE)) %>%
      clean_names() %>% 
      mutate(across(where(is.character), str_to_lower)) %>%
      rename_with(~str_replace(., 'sample', 'plate'),
                  .cols = c(starts_with('sample'), 
                            -any_of(c('sample_id', 'sample_type')))) %>%
      rename_with(~str_replace(., 'column$', 'col')) %>%
      rename_with(~str_remove(., '^pcr[12]_'),
                  .cols = c(contains('plate_id'), 
                            contains('plate_col'), 
                            contains('plate_row'))) %>%
      rename(dna_plate_id = plate_id,
             dna_plate_col = plate_col,
             dna_plate_row = plate_row)
  }, error = function(e) {
    stop(paste("Error reading plate map:", e$message))
  })
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
  tryCatch({
    candidate_cols <- df1 %>% select(where(is.character)) %>% colnames()
    
    # Compute the overlap score between df2$dna_plate_id and each candidate column in df1
    overlap_scores <- sapply(candidate_cols, function(col_name) {
      jaccard_similarity(df1[[col_name]], df2$dna_plate_id)
    })
    
    if (all(overlap_scores == 0)) {
      stop("No matching plate IDs found between plate map and quantification data. Please check that the plate IDs match between files.")
    }
    
    col_pair <- set_names('dna_plate_id',
                          names(which.max(overlap_scores)))
    
    row_col_match <- set_names(c('dna_plate_row', 'dna_plate_col'), 
                               get_subsequent_row_col(names(col_pair), df1))
    
    out <- c(col_pair, row_col_match)
    
    message('joining: dna_plate_map - ', str_c(names(out), collapse = '; '))
    message('joining: quant_plates - ', str_c(out, collapse = '; '))
    
    out
  }, error = function(e) {
    stop(paste("Error identifying join columns:", e$message))
  })
}

join_quants_map <- function(map_data = dna_plate_map, quant_data = quant_plates, quant_type = 'original', session = NULL){
  tryCatch({
    # write_rds(map_data, 'map_data.rds')
    # write_rds(quant_data, 'quant_data.rds')
    
    join_vars <- identify_join_cols(map_data, quant_data)
    
    out <- inner_join(map_data,
                      quant_data,
                      by = join_vars)
    
    # Display informational message to user
    info_message <- paste0('Detected ', nrow(map_data), ' samples with ', 
                           nrow(quant_data) / nrow(map_data), ' replicates each.\n',
                           'Please check that the number of samples and replicates matches what is expected.')
    
    if (!is.null(session)) {
      showNotification(info_message, type = "message",
                       duration = NULL,
                       closeButton = TRUE, session = session)
    } else {
      message(info_message) # fallback for console
    }
    
    # Check for unique sample names
    if(nrow(map_data) * (nrow(quant_data) / nrow(map_data)) != nrow(out)){
      warning_message <- 'Check sample names are unique with plate ID, plate row, and plate column'
      if (!is.null(session)) {
        showNotification(warning_message, type = "warning",
                         duration = NULL,
                         closeButton = TRUE, session = session)
      } else {
        message(warning_message) # fallback for console
      }
    }
    
    if (nrow(out) == 0) {
      stop("No matching records found between plate map and quantification data. Please check that the plate IDs, rows, and columns match between files.")
    }
    
    attr(out, 'join_cols') <- join_vars
    out
  }, error = function(e) {
    stop(paste("Error joining data:", e$message))
  })
}

process_merged_data <- function(df){
  tryCatch({
    the_cols <- colnames(df)
    join_cols <- names(attr(df, 'join_cols'))
    
    if(!any(str_detect(the_cols, 'sample_type'))){
      df <- mutate(df, sample_type = 'sample')
    }
    if(any(str_detect(the_cols, '^sample_id$'))){
      df <- rename(df, sample_id_store = 'sample_id')
    }
    df %>%
      mutate(across(ends_with('row'), str_to_upper),
             sample_id = str_c(!!sym(join_cols[1]), '_', !!sym(join_cols[2]), !!sym(join_cols[3])),
             is_control = str_detect(sample_type, 'control'))
  }, error = function(e) {
    stop(paste("Error processing merged data:", e$message))
  })
}

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
      max-width: 500px !important;
      white-space: pre-wrap !important;
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
    
    .filter-controls {
      background-color: #f8f9fa;
      border: 1px solid #dee2e6;
      border-radius: 5px;
      padding: 15px;
      margin-bottom: 20px;
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
    
    /* Style for the output prefix section */
    .output-prefix-section {
      background-color: #f0f8ff;
      border: 1px solid #b0d4f1;
      border-radius: 5px;
      padding: 15px;
      margin: 15px 0;
    }
    
    .output-prefix-section h5 {
      margin-top: 0;
      color: #0066cc;
    }
  "))),
  
  titlePanel("Step 2: Calculate Mean DNA Concentration"),
  
  # Custom navigation header
  div(
    style = "border-bottom: 1px solid #ddd; margin-bottom: 20px; padding-bottom: 10px;",
    div(
      style = "display: flex; justify-content: space-between; align-items: center;",
      div(),
      tags$a(href = "http://10.5.146.65/DNA_Quantification/", 
             target = "_blank",
             "Return to Menu",
             style = paste0("color: #337ab7; text-decoration: none; font-weight: bold; ",
                            "padding: 8px 16px; border: 1px solid #337ab7; ",
                            "border-radius: 4px; background-color: #f8f9fa; ",
                            "transition: background-color 0.2s;"),
             onmouseover = "this.style.backgroundColor='#e9ecef'",
             onmouseout = "this.style.backgroundColor='#f8f9fa'")
    )
  ),
  
  tabsetPanel(
    tabPanel("Data Input",
             div(class = "tab-content-mobile",
                 # Mobile instructions (appears above controls on smaller screens)
                 div(class = "instructions-mobile instructions-box",
                     h4("Instructions:"),
                     tags$ul(
                       tags$li("Upload the DNA Plate Map File (not the quant plate map)"),
                       tags$li("Be sure the correct sheet within the excel file is selected"),
                       tags$li("Upload the DNA concentrations output by Step 1: Model DNA Concentrations"),
                       tags$li("Select ", tags$code("Load Data")),
                       tags$li("Go to ", tags$code("Filter Data"), " Tab to review and filter your data"),
                       tags$li("Then go to ", tags$code("Calculate Mean DNA Concentrations"), " Tab")
                     )
                 ),
                 
                 sidebarLayout(
                   sidebarPanel(
                     fileInput("excel_file", "Upload Source DNA Plate Map File(s)", 
                               accept = c(".xlsx", ".xls"), multiple = TRUE),
                     selectInput("sheet_name", "Excel Sheet Name", choices = NULL),
                     fileInput("csv_files", "Upload Modeled DNA Concentration Files (select all)", 
                               multiple = TRUE, accept = ".csv"),
                     selectInput("quant_type", "Quant Stage", 
                                 choices = c('Original' = 'original', 
                                             "Requantification" = 'requant'), 
                                 selected = 'original'),
                     actionButton("load_data", "Load Data"),
                     uiOutput("data_status")
                   ),
                   mainPanel(
                     # Desktop instructions (appears in main panel on larger screens)
                     div(class = "instructions-desktop instructions-box",
                         h4("Instructions:"),
                         tags$ul(
                           tags$li("Upload the DNA Plate Map File (not the quant plate map)"),
                           tags$li("Be sure the correct sheet within the excel file is selected"),
                           tags$li("Upload the DNA concentrations output by Step 1: Model DNA Concentrations"),
                           tags$li("Select ", tags$code("Load Data")),
                           tags$li("Go to ", tags$code("Filter Data"), " Tab to review and filter your data"),
                           tags$li("Then go to ", tags$code("Calculate Mean DNA Concentrations"), " Tab")
                         )
                     ),
                     h4("Merged Data"),
                     DT::DTOutput("merged_table")
                   )
                 )
             )
    ),
    
    # Filter Data Tab
    tabPanel("Filter Data",
             div(class = "tab-content-mobile",
                 # Mobile instructions
                 div(class = "instructions-mobile instructions-box",
                     h4("Instructions:"),
                     tags$ul(
                       tags$li("Review your merged data in the table below"),
                       tags$li("Use the filter controls on the left to filter your data"),
                       tags$li("The table will update automatically as you apply filters"),
                       tags$li("Only the filtered data will be used for model fitting"),
                       tags$li("Proceed to ", tags$code("Calculate Mean DNA Concentrations"), " when ready")
                     )
                 ),
                 
                 sidebarLayout(
                   sidebarPanel(
                     h4("Data Filtering Controls"),
                     
                     # Dynamic filter controls will be generated here
                     uiOutput("filter_controls"),
                     
                     br(),
                     actionButton("reset_filters", "Reset All Filters", 
                                  class = "btn-warning"),
                     br(), br(),
                     
                     # Filter summary
                     div(
                       style = "background-color: #e9ecef; border-radius: 5px; padding: 10px;",
                       h5("Filter Summary:"),
                       verbatimTextOutput("filter_summary")
                     )
                   ),
                   mainPanel(
                     # Desktop instructions
                     div(class = "instructions-desktop instructions-box",
                         h4("Instructions:"),
                         tags$ul(
                           tags$li("Review your merged data in the table below"),
                           tags$li("Use the filter controls on the left to filter your data"),
                           tags$li("The table will update automatically as you apply filters"),
                           tags$li("Only the filtered data will be used for model fitting"),
                           tags$li("Proceed to ", tags$code("Calculate Mean DNA Concentrations"), " when ready")
                         )
                     ),
                     
                     h4("Filtered Data"),
                     p("This is the data that will be used for model fitting:"),
                     DT::DTOutput("filtered_table")
                   )
                 )
             )
    ),
    
    tabPanel("Calculate Mean DNA Concentrations",
             div(class = "tab-content-mobile",
                 # Mobile instructions
                 div(class = "instructions-mobile instructions-box",
                     h4("Instructions:"),
                     tags$ul(
                       tags$li("Select ", tags$code("Fit Model"), " (can take a minute)"),
                       tags$li("Evaluate the resulting plot"),
                       tags$li("Review and optionally edit the output filename prefix"),
                       tags$li("Adjust settings as necessary:"),
                       tags$ul(
                         tags$li("Excess DNA: changes which samples are flagged for dilution"),
                         tags$li("Other flag settings can be adjusted based on your experimental needs")
                       ),
                       tags$li("Select ", tags$code("Generate Zip File"), " then ", tags$code("Download Zip File"))
                     )
                 ),
                 
                 sidebarLayout(
                   sidebarPanel(
                     # DNA Concentration input
                     selectInput("y_var", "Select DNA Concentration", choices = NULL),
                     
                     # Other controls
                     actionButton("fit_model", "Fit Model"),
                     
                     # Output Filename Prefix section - NEW
                     div(class = "output-prefix-section",
                         h5("Output File Settings"),
                         textInput("outnamePrefix", 
                                   "Output Filename Prefix:", 
                                   value = "",
                                   placeholder = "e.g., plate123_analysis"),
                         helpText("This prefix will be used for all generated output files. Leave blank to use auto-generated prefix from input filenames.")
                     ),
                     
                     # Download control
                     actionButton("generate_zip", "Generate Zip File"),
                     uiOutput("download_zip_ui"),
                     br(), br(),
                     checkboxInput("log_transform", "Plot Log10 Transform", value = TRUE), 
                     
                     # All your existing collapsible panels remain the same...
                     bsCollapse(
                       bsCollapsePanel("Flag Settings",
                                       numericInput("min_volume", "Minimum Pipettable Volume (μL)", value = 0.75, min = 0),
                                       numericInput("max_low_volume", "Max Vol (μL) to Pipette", value = 4, min = 0),
                                       numericInput("target_dna", "Target Amount of DNA to Transfer (ng)", value = 2, min = 0),
                                       numericInput("mean_multiple", "Excess DNA is 'X' times more than the upper 95% CI", value = 2, min = 0),
                                       style = "primary"
                       ),
                       open = "Flag Settings"
                     ),
                     
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
                     )
                   ),
                   mainPanel(
                     # Desktop instructions
                     div(class = "instructions-desktop instructions-box",
                         h4("Instructions:"),
                         tags$ul(
                           tags$li("Select ", tags$code("Fit Model"), " (can take a minute)"),
                           tags$li("Evaluate the resulting plot"),
                           tags$li("Review and optionally edit the output filename prefix"),
                           tags$li("Adjust settings as necessary:"),
                           tags$ul(
                             tags$li("Excess DNA: changes which samples are flagged for dilution"),
                             tags$li("Other flag settings can be adjusted based on your experimental needs")
                           ),
                           tags$li("Select ", tags$code("Generate Zip File"), " then ", tags$code("Download Zip File"))
                         )
                     ),
                     plotOutput("dna_plot", height = "calc(100vh - 250px)")
                   )
                 )
             )
    )
  )
)

#### Server ####
server <- function(input, output, session) {
  
  # Set up roots for shinyFiles (here, we allow navigation of the home directory)
  roots <- c(home = "~")
  shinyFileSave(input, "download_csv", roots = roots, session = session, filetypes = c("csv"))
  
  # Update the sheet selection when Excel file(s) are uploaded
  observeEvent(input$excel_file, {
    tryCatch({
      req(input$excel_file)
      
      # Get sheets from all Excel files and combine unique sheet names
      all_sheets <- character(0)
      for (i in 1:nrow(input$excel_file)) {
        if (file.exists(input$excel_file$datapath[i])) {
          sheets <- excel_sheets(input$excel_file$datapath[i])
          all_sheets <- unique(c(all_sheets, sheets))
        }
      }
      
      if (length(all_sheets) > 0) {
        updateSelectInput(session, "sheet_name", choices = all_sheets)
      } else {
        showNotification("No sheets found in Excel files.", type = "error",
                         duration = NULL,
                         closeButton = TRUE)
        updateSelectInput(session, "sheet_name", choices = NULL)
      }
    }, error = function(e) {
      showNotification(paste("Error reading Excel file(s):", e$message), type = "error",
                       duration = NULL,
                       closeButton = TRUE)
      updateSelectInput(session, "sheet_name", choices = NULL)
    })
  })
  
  # Update y_var choices based on numeric columns in the combined CSV files
  observe({
    if (is.null(input$csv_files)) return()
    
    tryCatch({
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
    }, error = function(e) {
      # Only show error if CSV files were actually selected
      if (!is.null(input$csv_files)) {
        showNotification(paste("Error processing CSV files:", e$message), type = "error",
                         duration = NULL,
                         closeButton = TRUE)
      }
    })
  })
  
  observeEvent(input$y_var, {
    # pull the unit from the selected column name: "ng" or "pg"
    unit <- stringr::str_extract(input$y_var, "[pn]g")   # returns "ng" or "pg"
    
    if (!is.na(unit)) {
      updateNumericInput(
        session,
        inputId = "target_dna",
        label   = paste0("Target Amount of DNA (", unit, ")")
      )
    }
  })
  
  # Reactive expression for joined data, triggered by the Load Data button
  joined_data <- eventReactive(input$load_data, {
    tryCatch({
      req(input$excel_file, input$sheet_name, input$csv_files, input$y_var)
      
      # Read and combine multiple plate map files
      dna_plate_map <- map_dfr(1:nrow(input$excel_file), function(i) {
        filepath <- input$excel_file$datapath[i]
        filename <- input$excel_file$name[i]
        
        tryCatch({
          read_plate_map(filepath, input$sheet_name)
        }, error = function(e) {
          # If the sheet doesn't exist in this file, return empty tibble
          # This allows some files to have the sheet and others not
          message("Sheet '", input$sheet_name, "' not found in ", filename, ": ", e$message)
          return(tibble())
        })
      })
      
      # Check if any data was loaded
      if (nrow(dna_plate_map) == 0) {
        stop(paste("No data found in sheet '", input$sheet_name, "' across all Excel files. Please check the sheet name."))
      }
      
      # Read and combine the CSV files
      quant_plates <- map_dfr(input$csv_files$datapath, ~ read_csv(.x, show_col_types = FALSE)) %>%
        rename_with(~str_replace(., 'column', 'col'))
      
      # Auto-populate the output prefix based on common prefix - NEW
      if (nchar(input$outnamePrefix) == 0) {  # Only update if empty
        prefix <- common_prefix(input$excel_file$name[1], input$csv_files$name[1])
        if (nchar(prefix) > 0) {
          updateTextInput(session, "outnamePrefix", value = prefix)
        }
      }
      
      # write_csv(dna_plate_map, 'dna_plate_map.csv')
      # write_csv(quant_plates, 'quant_plates.csv')
      
      # Merge the plate map with the quant data using an inner join
      # PASS SESSION TO THE FUNCTION
      merged_dna_quants <- join_quants_map(dna_plate_map,
                                           quant_plates,
                                           input$quant_type,
                                           session = session) %>%  # Add session parameter
        process_merged_data() 
      
      showNotification(paste("Data loaded successfully!", nrow(merged_dna_quants), "samples loaded from", nrow(input$excel_file), "Excel file(s)."), 
                       type = "message",
                       duration = NULL,
                       closeButton = TRUE)
      
      merged_dna_quants
    }, error = function(e) {
      showNotification(paste("Data Loading Error:", e$message),
                       duration = NULL,
                       closeButton = TRUE)
      return(NULL)
    })
  })
  
  
  # NEW: Generate dynamic filter controls
  output$filter_controls <- renderUI({
    req(joined_data())
    
    data <- joined_data()
    controls <- list()
    
    # Get column names and types
    char_cols <- names(data)[sapply(data, function(x) is.character(x) || is.factor(x))]
    num_cols <- names(data)[sapply(data, is.numeric)]
    
    # Create filter controls for character/factor columns
    for (col in char_cols) {
      unique_vals <- sort(unique(data[[col]], na.rm = TRUE))
      if (length(unique_vals) <= 20) {  # Only create filter if manageable number of options
        controls[[paste0(col, "_filter")]] <- checkboxGroupInput(
          inputId = paste0(col, "_filter"),
          label = paste("Filter", str_to_title(str_replace_all(col, "_", " ")), ":"),
          choices = unique_vals,
          selected = unique_vals,
          inline = TRUE
        )
      }
    }
    
    # Create filter controls for numeric columns
    for (col in num_cols) {
      if (sum(!is.na(data[[col]])) > 0) {  # Only if there are non-NA values
        col_range <- range(data[[col]], na.rm = TRUE)
        controls[[paste0(col, "_range")]] <- sliderInput(
          inputId = paste0(col, "_range"),
          label = paste("Filter", str_to_title(str_replace_all(col, "_", " ")), ":"),
          min = col_range[1],
          max = col_range[2],
          value = col_range,
          step = ifelse(diff(col_range) > 100, 1, 0.01)
        )
      }
    }
    
    controls
  })
  
  # NEW: Filtered data reactive
  filtered_data <- reactive({
    req(joined_data())
    
    data <- joined_data()
    
    # Start with all data  
    filtered <- data
    
    # Track if any filters have been applied
    any_filters_applied <- FALSE
    
    # Check if we're in initial load state (all filters are NULL)
    # by looking at all checkbox filters
    char_cols <- names(data)[sapply(data, function(x) is.character(x) || is.factor(x))]
    all_filters_null <- TRUE
    for (col in char_cols) {
      unique_vals <- sort(unique(data[[col]], na.rm = TRUE))
      if (length(unique_vals) <= 20) {
        if (!is.null(input[[paste0(col, "_filter")]])) {
          all_filters_null <- FALSE
          break
        }
      }
    }
    
    # Apply character/factor filters using dplyr
    for (col in char_cols) {
      unique_vals <- sort(unique(data[[col]], na.rm = TRUE))
      if (length(unique_vals) <= 20) {  # Only columns with checkbox controls
        filter_id <- paste0(col, "_filter")
        filter_input <- input[[filter_id]]
        
        # Skip if filter UI hasn't been created yet (NULL on first render)
        # BUT only skip if ALL filters are NULL (initial state)
        if (is.null(filter_input) && all_filters_null) {
          next
        }
        
        # Filter has been created or we're past initial load
        any_filters_applied <- TRUE
        
        # Treat NULL as empty selection if we're past initial load
        if (is.null(filter_input) || length(filter_input) == 0) {
          # Empty vector OR NULL after init - user unchecked all boxes, filter out everything
          filtered <- filtered %>% slice(0)
          break
        } else {
          # Normal filtering when some values are selected
          filtered <- filtered %>% 
            filter(!!sym(col) %in% filter_input | is.na(!!sym(col)))
        }
      }
    }
    
    # Apply numeric filters using dplyr (only if we still have rows)
    if (nrow(filtered) > 0) {
      num_cols <- names(data)[sapply(data, is.numeric)]
      for (col in num_cols) {
        if (sum(!is.na(data[[col]])) > 0) {
          range_input <- input[[paste0(col, "_range")]]
          if (!is.null(range_input) && length(range_input) == 2) {
            any_filters_applied <- TRUE
            filtered <- filtered %>% 
              filter(between(!!sym(col), range_input[1], range_input[2]) | is.na(!!sym(col)))
          }
        }
      }
    }
    
    return(filtered)
  })
  
  
  # NEW: Reset filters
  observeEvent(input$reset_filters, {
    req(joined_data())
    
    data <- joined_data()
    
    # Reset character filters
    char_cols <- names(data)[sapply(data, function(x) is.character(x) || is.factor(x))]
    for (col in char_cols) {
      unique_vals <- sort(unique(data[[col]], na.rm = TRUE))
      if (length(unique_vals) <= 20) {
        updateCheckboxGroupInput(
          session,
          paste0(col, "_filter"),
          selected = unique_vals
        )
      }
    }
    
    # Reset numeric filters
    num_cols <- names(data)[sapply(data, is.numeric)]
    for (col in num_cols) {
      if (sum(!is.na(data[[col]])) > 0) {
        col_range <- range(data[[col]], na.rm = TRUE)
        updateSliderInput(
          session,
          paste0(col, "_range"),
          value = col_range
        )
      }
    }
  })
  
  observe({
    req(filtered_data())
    
    filter_messages <- attr(filtered_data(), "filter_messages")
    
    if (length(filter_messages) > 0) {
      # Combine all messages
      combined_message <- paste(filter_messages, collapse = "\n")
      
      # Show notification with all filter effects
      showNotification(
        combined_message,
        type = "message",
        duration = 5,
        closeButton = TRUE
      )
    }
  })
  
  # NEW: Filter summary
  output$filter_summary <- renderText({
    req(joined_data())
    
    # Get filtered data
    filtered <- filtered_data()
    
    original_rows <- nrow(joined_data())
    filtered_rows <- nrow(filtered)
    
    # Get filter messages for detailed summary
    filter_messages <- attr(filtered, "filter_messages")
    
    summary_text <- paste0(
      "Original data: ", original_rows, " rows\n",
      "Filtered data: ", filtered_rows, " rows\n",
      "Rows removed: ", original_rows - filtered_rows, " (", 
      round((original_rows - filtered_rows) / original_rows * 100, 1), "%)"
    )
    
    # Add warning if all data filtered out
    if (filtered_rows == 0) {
      summary_text <- paste0(
        summary_text, "\n\n",
        "⚠️ WARNING: All data filtered out!\n",
        "Adjust filters to see data."
      )
    }
    
    summary_text
  })
  
  # Data status indicator
  output$data_status <- renderUI({
    data <- joined_data()
    if (!is.null(data) && nrow(data) > 0) {
      tags$div(
        style = "color: #28a745; font-weight: bold; padding: 10px; margin: 5px 0; border: 1px solid #28a745; border-radius: 5px; background-color: #f8fff9;",
        icon("check-circle"), 
        paste("✔ Data loaded:", nrow(data), "samples")
      )
    } else {
      tags$div(
        style = "color: #6c757d; font-style: italic; padding: 10px; margin: 5px 0;",
        icon("info-circle"), 
        "No data loaded. Please select files and click 'Load Data'."
      )
    }
  })
  
  # Render the merged table when joined_data is updated
  output$merged_table <- DT::renderDT({
    joined_data()  # Ensure the merged data is re-rendered
  })
  
  # NEW: Render filtered table
  output$filtered_table <- DT::renderDT({
    req(joined_data())
    
    data <- filtered_data()
    
    # Return the data as-is, even if empty
    # DT will show "No data available in table" automatically for empty data frames
    data
  }, options = list(scrollX = TRUE, pageLength = 25))
  
  # Modified model_fit using callr::r_bg for asynchronous Stan sampling
  model_fit <- eventReactive(input$fit_model, {
    tryCatch({
      req(filtered_data(), input$y_var)
      
      # NEW: Check if filtered data is empty
      if (nrow(filtered_data()) == 0) {
        showNotification(
          "Cannot fit model: No data remains after filtering. Please adjust your filters.", 
          type = "error",
          duration = NULL,
          closeButton = TRUE
        )
        return(NULL)
      }
      
      # NEW: Check if there's sufficient data for modeling after removing zeros
      model_data <- filter(filtered_data(), !!sym(input$y_var) > 0)
      if (nrow(model_data) < 2) {
        showNotification(
          paste("Cannot fit model: Insufficient data after filtering and removing zero values.",
                "Need at least 2 observations, but only have", nrow(model_data)), 
          type = "error",
          duration = NULL,
          closeButton = TRUE
        )
        return(NULL)
      }
      
      
      model_formula <- str_c(input$y_var, '~ is_control + (1 | sample_id)')
      dna_form <- bf(as.formula(model_formula),
                     shape ~ is_control + (1 | sample_id))
      dna_model_bayes <- brm(dna_form,
                             family = Gamma(link = 'log'),
                             data = filter(model_data,
                                           !!sym(input$y_var) > 0),
                             empty = TRUE)
      
      if(Sys.info()["nodename"] == 'gawain'){
        dna_model <- cmdstan_model('model/dna_concentration_threaded.stan',
                                   cpp_options = list(stan_threads = TRUE))
      } else {
        dna_model <- cmdstan_model('dna_concentration_threaded.stan',
                                   cpp_options = list(stan_threads = TRUE))
      }
      
      # Note: Stan expects "inteceptShape_prior" (without the "r") as in the original code.
      data_stan <- c(
        make_standata(bf(as.formula(model_formula),
                         shape ~ is_control + (1 | sample_id)),
                      family = Gamma(link = 'log'),
                      data = model_data,
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
      
      if (!dir.exists(output_dir)) {
        dir.create(output_dir)
      }
      
      # Define the function to run sampling in the background using callr::r_bg.
      sampling_function <- function(data_stan, output_dir, num_chains, iter_sampling, iter_warmup, thin) {
        library(cmdstanr)
        if(Sys.info()["nodename"] == 'gawain'){
          set_cmdstan_path('/home/shiny/.cmdstan/cmdstan-2.36.0')
        }
        
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
        
        fit$output_files()
      }
      
      # Launch sampling in the background with user-defined settings.
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
      
      # Reconstruct the Stan fit using the output files.
      dna_model_bayes$fit <- read_csv_as_stanfit(output_files, model = dna_model)
      dna_model_bayes <- rename_pars(dna_model_bayes)
      
      # Delete the temporary output directory and its files.
      unlink(output_dir, recursive = TRUE)
      
      dna_model_bayes
    }, error = function(e) {
      showNotification(paste("Model Fitting Error:", e$message), type = "error",
                       duration = NULL,
                       closeButton = TRUE)
      return(NULL)
    })
  })
  
  #### DNA Summary Calculations ####
  dna_summary <- eventReactive(input$fit_model, {
    tryCatch({
      req(filtered_data(), input$y_var, model_fit())
      dna_model_bayes <- model_fit()
      
      quant_files_summarized <- filtered_data() %>%
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
    }, error = function(e) {
      showNotification(paste("Error calculating DNA summary:", e$message), type = "error",
                       duration = NULL,
                       closeButton = TRUE)
      return(NULL)
    })
  })
  
  # Reactive expression for the flagged table that will be used for both plotting and download.
  flagged_table <- reactive({
    tryCatch({
      req(filtered_data(), dna_summary())
      
      # NEW: Check for valid summary
      if (is.null(dna_summary())) {
        return(NULL)
      }
      
      quant_files_summarized <- dna_summary()$quant_files_summarized
      mean_dna_concentration <- dna_summary()$mean_dna_concentration
      dna_quantity_interval <- dna_summary()$dna_quantity_interval
      dna_variability_interval <- dna_summary()$dna_variability_interval
      
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
               !!sym(str_c('rxn_', str_extract(input$y_var, '[pn]g'))) := if_else(is.na(ul_per_rxn),
                                                                                  postDilution_concentration * postDilution_ul_per_rxn,
                                                                                  !!sym(str_c(input$y_var, "_mean")) * ul_per_rxn)) %>%
        left_join(select(dna_variability_interval, is_control, var_upr_limit = upr_limit),
                  by = "is_control") %>%
        left_join(select(dna_quantity_interval, is_control, mean_upr_limit = upr_limit),
                  by = 'is_control') %>%
        mutate(flags = case_when(
          is_control & !!sym(str_c(input$y_var, "_upr95")) > mean_dna_concentration ~ "Contaminated Control",
          !!sym(str_c(input$y_var, "_mean")) > (input$mean_multiple * mean_upr_limit) & !!sym(str_c(input$y_var, "_normspread")) > var_upr_limit ~ "Excess & Variable DNA",
          !!sym(str_c(input$y_var, "_mean")) > (input$mean_multiple * mean_upr_limit) ~ "Excess DNA",
          !!sym(str_c(input$y_var, "_normspread")) > var_upr_limit ~ "Variable DNA",
          TRUE ~ "Good Sample"),
          .after = sample_type) %>%
        select(sample_id,
               sample_type,
               all_of(str_c(input$y_var, c("_mean", "_lwr95", "_upr95"))),
               flags,
               ul_per_rxn, dilution_factor,
               starts_with("postDilution"),
               any_of(c('rxn_ng', 'rxn_pg')))
      
      quant_files_flagged
    }, error = function(e) {
      showNotification(paste("Error generating flagged table:", e$message), type = "error",
                       duration = NULL,
                       closeButton = TRUE)
      return(NULL)
    })
  })
  
  #### Plot Generation Reactive ####
  
  plot_obj <- reactive({
    tryCatch({
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
        mutate(dna_plate_well_id = str_extract(sample_id, '_[a-zA-Z][0-9]+$') %>% str_remove('_'),
               facet_plate = str_remove(sample_id, str_c('_', dna_plate_well_id)),
               sample_type = fct_relevel(sample_type, 'sample', after = 0L)) %>%
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
                                      'control' = 'circle filled',
                                      "extraction control" = 'triangle filled',
                                      "field control" = 'triangle down filled',
                                      "filter control" = 'square filled',
                                      "pcr control" = "diamond filled"),
                           labels = str_to_title) +
        
        scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"),
                              labels = c("TRUE" = "Control", "FALSE" = "Sample")) +
        guides(shape = guide_legend(override.aes = list(fill = if_else(str_detect(unique(quant_files_flagged$sample_type), 
                                                                                  'control'),
                                                                       'black', 'white'),
                                                        size = 3)),
               fill = "none",
               colour = guide_legend(override.aes = list(size = 3))) +
        facet_wrap(~facet_plate, nrow = 1, scales = "free_y") +
        theme_bw() +
        labs(x = str_c("DNA Concentration (", str_extract(input$y_var, '[pn]g'), "/uL) ± 95% CI"),
             y = "Sample Well ID",
             color = "Sample Flag",
             shape = "Sample Type",
             linetype = "Sample Type (Mean)") +
        theme(axis.text.x = element_text(size = 7))
      
      if (input$log_transform) {
        quant_plot <- quant_plot + 
          scale_x_continuous(transform = scales::log10_trans(),
                             labels = scales::comma_format())
      }
      quant_plot
    }, error = function(e) {
      showNotification(paste("Error generating plot:", e$message), type = "error",
                       duration = NULL,
                       closeButton = TRUE)
      ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                 label = paste("Plot error:", e$message), 
                 size = 4, color = "red", hjust = 0.5) +
        theme_void()
    })
  })
  
  output$dna_plot <- renderPlot({
    if (is.null(filtered_data()) || (length(filtered_data()) > 0 && nrow(filtered_data()) == 0)) {
      # No data loaded yet
      ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                 label = "No data loaded.\nPlease go to 'Data Input' tab and load your files.", 
                 size = 6, color = "gray50", hjust = 0.5) +
        theme_void()
    } else if (is.null(input$y_var) || input$y_var == "") {
      # Y variable not selected
      ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                 label = "No concentration variable selected.\nPlease select a DNA concentration variable above.", 
                 size = 6, color = "gray50", hjust = 0.5) +
        theme_void()
    } else {
      # Try to check if model is fitted
      model_result <- tryCatch({
        model_fit()
      }, error = function(e) {
        NULL
      })
      
      if (is.null(model_result)) {
        # Model not fitted yet
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, 
                   label = "Model not fit yet.\nClick 'Fit Model' button to generate plot.", 
                   size = 6, color = "gray50", hjust = 0.5) +
          theme_void()
      } else {
        # Model is fitted, generate the plot
        plot_obj()
      }
    }
  })
  
  # Download handler using shinyFiles save button selection.
  output$download_csv <- downloadHandler(
    filename = function() {
      req(input$excel_file, input$csv_files)
      # Use user-provided prefix or fall back to common prefix
      prefix <- if (nchar(input$outnamePrefix) > 0) {
        input$outnamePrefix
      } else {
        common_prefix(input$excel_file$name[1], input$csv_files$name[1])
      }
      prefix <- tolower(prefix)
      paste0(prefix, "_", input$quant_type, "_DNA_concentrations_", Sys.Date(), ".csv")
    },
    content = function(file) {
      out_table <- select(filtered_data(), -rfu, -input$y_var, 
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
  
  # Reactive value to store the generated ZIP file's path.
  zip_info <- reactiveVal(NULL)
  
  # Create a waiter object (you can customize the spinner as needed)
  w <- waiter::Waiter$new(html = waiter::spin_fading_circles(), color = "rgba(255,255,255,0.8)")
  
  # When the user clicks "Generate Zip File":
  observeEvent(input$generate_zip, {
    # Show the spinner
    w$show()
    
    # Ensure required reactives are available
    req(filtered_data(), flagged_table())
    
    # Use the user-provided prefix or fall back to common prefix
    prefix <- if (nchar(input$outnamePrefix) > 0) {
      tolower(input$outnamePrefix)
    } else {
      tolower(common_prefix(input$excel_file$name[1], input$csv_files$name[1]))
    }
    
    # Create a temporary directory and file paths
    tmpdir <- tempdir()
    tmpdir <- file.path("outdir", basename(tmpdir))
    dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
    
    zip_path <- file.path(paste0(prefix, "_output_", Sys.Date(), ".zip"))
    prediction_path <- file.path(tmpdir, paste0(prefix, "_sample_concentrations_", Sys.Date(), ".csv"))
    excel_path <- file.path(tmpdir, paste0(prefix, "_sample_predictions_", Sys.Date(), ".xlsx"))
    inputs_path <- file.path(tmpdir, paste0(prefix, "_inputs_", Sys.Date(), ".txt"))
    json_path <- file.path(tmpdir, paste0(prefix, "_settings_", Sys.Date(), ".json"))
    model_code <- file.path(tmpdir, paste0(prefix, "_models_", Sys.Date(), ".stan"))
    model_path <- file.path(tmpdir, paste0(prefix, "_models_", Sys.Date(), ".rds"))
    code_path <- file.path(tmpdir, paste0(prefix, "_sourceCode_", Sys.Date(), ".R"))
    plots_path <- file.path(tmpdir, paste0(prefix, "_plots_", Sys.Date(), ".png"))
    
    #### Output predictions ####
    out_table <- select(filtered_data(), -rfu, -input$y_var, 
                        -is_control,
                        -starts_with('quant')) %>% 
      distinct() %>%
      full_join(select(flagged_table(), 
                       matches('sample_id'),
                       starts_with(input$y_var),
                       flags),
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
              file = prediction_path,
              na = '')
    
    #### Output as Excel ####
    write_xlsx(out_table, excel_path)
    
    #### Output inputs ####
    input_list <- lapply(names(input), function(name) {
      paste(name, ":", toString(input[[name]]))
    })
    write_lines(input_list, inputs_path)
    
    #### Output JSON ####
    model_details <- list(
      bayesian_fitting_settings = list('Number of Chains'= input$num_chains, 
                                       'Number of Sampling Iterations' = input$iter_sampling, 
                                       'Number of Warmup Iterations' = input$iter_warmup, 
                                       'Thinning' = input$thin),
      mean_priors = list('intercept mean' = input$intercept_prior_mean, 'intercept sd' = input$intercept_prior_sd, 
                         'beta (sample/control) mean' = input$beta_prior_mean, 'beta (sample/control) sd' = input$beta_prior_sd, 
                         'Random Effect Prior (shape)' = input$var_prior_param1, 'Random Effect Prior (rate)' = input$var_prior_param2),
      shape_priors = list('intercept mean' = input$interceptShape_prior_mean, 'intercept sd' = input$interceptShape_prior_sd, 
                          'beta (sample/control) mean' = input$betaShape_prior_mean, 'beta (sample/control) sd' = input$betaShape_prior_sd, 
                          'Random Effect Prior (shape)' = input$varShape_prior_param1, 'Random Effect Prior (rate)' = input$varShape_prior_param2),
      flag_settings = list("Minimum Pipettable Volume" = input$min_volume, 
                           "Maximum Low Volume" = input$max_low_volume,
                           "Target DNA Amount" = input$target_dna)
    )
    json_text <- toJSON(model_details, pretty = TRUE, auto_unbox = TRUE)
    write_lines(json_text, json_path)
    
    #### Output Models ####
    write_rds(model_fit(), 
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
    
    ### Output STAN Model Code ####
    if(Sys.info()["nodename"] == 'gawain'){
      file.copy('model/dna_concentration_threaded.stan', model_code, overwrite = TRUE)
    } else {
      file.copy('dna_concentration_threaded.stan', model_code, overwrite = TRUE)
    }
    
    #### Output Plots ####
    # Save the combined plot as a PNG file
    ggsave(filename = plots_path, plot = plot_obj(), 
           device = "png", width = 10, height = 8)
    
    #### Zip contents ####
    files_to_zip <- c(code_path, 
                      inputs_path,
                      json_path,
                      model_path,
                      model_code,
                      prediction_path,
                      excel_path,
                      plots_path)
    
    # Change working directory to the temporary directory
    owd <- setwd(tmpdir)
    on.exit(setwd(owd), add = TRUE)
    
    # Create the zip file (using base names for a cleaner archive)
    zip(zipfile = zip_path, files = basename(files_to_zip))
    
    zip_info(list(zip = file.path(tmpdir, zip_path), tmpdir = tmpdir))
    w$hide()
  })
  
  # Render a download button once the ZIP file is generated
  output$download_zip_ui <- renderUI({
    req(zip_info())
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

shinyApp(ui = ui, server = server)