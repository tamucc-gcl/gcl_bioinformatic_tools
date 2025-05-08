##TODO - priors & model fit specs to collapsible options
##TODO - precomplied model for shinyverse??
##TODO - move the "prediction table to the model fit not the plot object" - should speed up flag toggling
##TODO - sort out progress bar issue... - maybe https://search.r-project.org/CRAN/refmans/tfaddons/html/callback_tqdm_progress_bar.html


library(shiny)
library(readxl)
library(readr)
library(dplyr)
library(janitor)
library(stringr)
library(purrr)
library(brms)  # ensure brms is installed

# Function to read in the plate map from Excel
read_plate_map <- function(filepath, sheet_name) {
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
                           multiple = TRUE, accept = ".csv")
               ),
               mainPanel(
                 h4("Merged Data"),
                 dataTableOutput("merged_table")
               )
             )
    ),
    tabPanel("DNA Amounts",
             sidebarLayout(
               sidebarPanel(
                 selectInput("y_var", "Select DNA Concentration", choices = NULL),
                 numericInput("min_volume", "Minimum Pipettable Volume", 
                              value = 0.75, min = 0),
                 numericInput("max_low_volume", "Maximum Low Volume", 
                              value = 4, min = 0),
                 numericInput("target_dna", "Target Amount of DNA", 
                              value = 2, min = 0),
                 actionButton("fit_model", "Fit Model")
               ),
               mainPanel(
                 h4("Model Summary"),
                 verbatimTextOutput("model_summary")
               )
             )
    )
  )
)

server <- function(input, output, session) {
  
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
    
    updateSelectInput(session, "y_var", choices = numeric_cols)
  })
  
  # Reactive expression to join the plate map with the quant CSV data
  joined_data <- reactive({
    req(input$excel_file, input$sheet_name, input$csv_files, input$y_var)
    
    # Read the plate map using the provided function
    dna_plate_map <- read_plate_map(input$excel_file$datapath, input$sheet_name)
    
    # Read and combine the CSV files
    quant_plates <- map_dfr(input$csv_files$datapath, ~ read_csv(.x, show_col_types = FALSE))
    
    # Merge the plate map with the quant data using an inner join
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
             all_of(input$y_var)) %>%
      mutate(across(ends_with('row'), str_to_upper),
             sample_id = str_c(dna_plate_id, '_', dna_plate_row, dna_plate_col),
             is_control = str_detect(sample_type, 'control'))
    
    merged_dna_quants
  })
  
  # Render the merged table in the Data Input tab
  output$merged_table <- renderDataTable({
    joined_data()
  })
  
  # Reactive expression to fit the model using BRMS once data is loaded and y_var is selected
  model_fit <- eventReactive(input$fit_model, {
    req(joined_data(), input$y_var)
    
    # Insert your BRMS model specifics here.
    model_formula <- str_c(y_var,  '~ is_control + (1 | sample_id)')
    
    #Initialize empty BRMS model
    dna_model_bayes <- brm(bf(as.formula(model_formula),
                              shape ~ is_control + (1 | sample_id)),
                           family = Gamma(link = 'log'),
                           data = merged_dna_quants,
                           empty = TRUE)
    
    dna_model <- cmdstan_model('dna_concentration.stan')
    
    data_stan <- c(make_standata(bf(as.formula(model_formula),
                                    shape ~ is_control + (1 | sample_id)),
                                 family = Gamma(link = 'log'),
                                 data = merged_dna_quants),
                   list(intercept_prior = c(1.5, 2.5),
                        beta_prior = c(0, 2.5),
                        var_prior = c(2, 1),
                        
                        inteceptShape_prior = c(0, 2.5),
                        betaShape_prior = c(0, 2.5),
                        varShape_prior = c(2, 1)))
    
    dna_fit <- dna_model$sample(
      data = data_stan,
      chains = 4,
      parallel_chains = parallelly::availableCores(),
      iter_sampling = 2000,
      iter_warmup = 1000,
      thin = 10,
      refresh = 500 # print update every 500 iters
    )
    
    dna_model_bayes$fit <- read_csv_as_stanfit(dna_fit$output_files(), 
                                               model = dna_model)
    dna_model_bayes <- rename_pars(dna_model_bayes)
    
  })
  
  # Display the model summary (or placeholder text) in the DNA Amounts tab
  output$model_summary <- renderPrint({
    model_fit()
  })
  
}

shinyApp(ui = ui, server = server)
