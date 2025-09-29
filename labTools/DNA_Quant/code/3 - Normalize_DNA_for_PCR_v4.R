library(shiny) |> suppressMessages() |> suppressWarnings()
library(DT) |> suppressMessages() |> suppressWarnings()
library(readr) |> suppressMessages() |> suppressWarnings()
library(dplyr) |> suppressMessages() |> suppressWarnings()
library(stringr) |> suppressMessages() |> suppressWarnings()
library(ggplot2) |> suppressMessages() |> suppressWarnings()
library(writexl) |> suppressMessages() |> suppressWarnings()

ui <- fluidPage(
  tags$head(tags$style(HTML("
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
      background-color: #e8f4fd;
      border: 1px solid #bee5eb;
      border-radius: 5px;
      padding: 12px;
      margin: 15px 0;
    }
    
    .output-prefix-section label {
      font-weight: bold;
      color: #004085;
    }
    
    .output-prefix-section .form-control {
      margin-top: 5px;
    }
    
    .output-prefix-section .help-block {
      margin-top: 5px;
      font-size: 0.9em;
      color: #666;
    }
  "))),
  
  titlePanel("Step 3: DNA Normalization Calculator"),
  
  # Custom navigation header that appears above content
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
  
  div(class = "tab-content-mobile",
      # Mobile instructions (appears above controls on smaller screens)
      div(class = "instructions-mobile instructions-box",
          h4("Instructions:"),
          tags$ul(
            tags$li("Load the sample concentration file(s) from step 2."),
            tags$ul(
              tags$li("You may have more than one file if some samples were diluted and requanted")
            ),
            tags$li("If you are using this to calculate volumes for pooling, then check the box"),
            tags$li("Add the max volume to transfer per source well in your source plate"),
            tags$ul(
              tags$li("Underestimate the volume")
            ),
            tags$li("Fill in the remaining values which describe the characteristics of the sample after it is normalized"),
            tags$ul(
              tags$li("The amount of DNA in the destination well should be estimated based on the mean DNA concentration of the \"samples\" and the desired goals of the normalization")
            ),
            tags$br(),
            tags$li("Evaluate the columns: ", tags$code("transfer_volume"), ", ", tags$code("ul_to_add"), ", ", tags$code("actual_ng_dna"), " with respect to ", tags$code("sample_type")),
            tags$ul(
              tags$li("It is a good sign if the ", tags$code("transfer_volume"), " is less than the max for \"samples\" and the max for \"neg ctrls\""),
              tags$li("It is a good sign if the ", tags$code("actual_ng_dna"), " is equal to the target for \"samples\" and less than the target for \"neg ctrls\""),
              tags$li("It is desirable to maximize the difference between the ", tags$code("actual_ng_dna"), " for \"samples\" and \"neg ctrls\" while ensuring that there is not too much difference among \"samples\""),
              tags$li("The ", tags$code("ul_to_add"), " column indicates how much water to add or remove via centrivapping, but likely holds little meaning for pooling")
            ),
            tags$li("Adjust settings as necessary and repeat until desired results are achieved"),
            tags$br(),
            tags$li("Optionally edit the output filename prefix"),
            tags$li("Download the processed data")
          )
      ),
      
      sidebarLayout(
        sidebarPanel(
          fileInput("sample_files", "Upload Sample-concentration CSV From Step 2", accept = ".csv"),
          
          checkboxInput("pooling", "Pooling", value = FALSE),
          
          # Conditional panels based on pooling checkbox
          conditionalPanel(
            condition = "!input.pooling",
            # Section break
            hr(),
            h4("Characteristics of Samples Before Normalization"),
            numericInput("max_vol", "Max Vol to Transfer per Source Well (µL):", value = 90, min = 1),
            
            # Section break
            hr(),
            h4("Characteristics of Samples After Normalization"),
            numericInput("number_PCR_rxns", "Number of Reactions Needed per Destination Well:", value = 36, min = 1),
            numericInput("ul_per_PCR", "Targeted Final Volume per Rxn in Destination Well (µL):", value = 2, min = 0, step = 0.5),
            numericInput("DNA_per_PCR", "DNA per Rxn in Destination Well (ng):", value = 10, min = 0, step = 0.1)
          ),
          
          conditionalPanel(
            condition = "input.pooling",
            hr(),
            h4("Characteristics of Samples Before Normalization"),
            numericInput("max_vol", "Max Vol to Transfer per Source Well (µL):", value = 90, min = 1),
            hr(),
            h4("Characteristics of Samples After Normalization"),
            # Hidden input for number of transfers (set to 1)
            div(style = "display: none;",
                numericInput("number_PCR_rxns_pooling", "Number of Transfers per Source Well to Destination Well:", value = 1, min = 1)
            ),
            numericInput("ul_per_PCR_pooling", "Targeted Final Volume per Source Well in Destination Well (µL):", value = 0, min = 0, step = 0.5),
            numericInput("DNA_per_PCR_pooling", "Targeted Amount of DNA per Source Well in Destination Well (ng):", value = 10, min = 0, step = 0.1)
          ),
          
          # Output filename prefix section - NEW
          hr(),
          div(class = "output-prefix-section",
              tags$label("Output Filename Prefix:"),
              textInput("outnamePrefix", 
                        label = NULL,
                        value = "",
                        placeholder = "e.g., plate123_normalized"),
              tags$div(class = "help-block",
                       "Customize the output filename. Date and extension will be added automatically.")
          ),
          
          br(),
          # Error messages display
          conditionalPanel(
            condition = "output.show_error",
            div(
              style = "color: #d9534f; background-color: #f2dede; border: 1px solid #ebccd1; border-radius: 4px; padding: 10px; margin-bottom: 10px;",
              strong("Error: "),
              textOutput("error_message", inline = TRUE)
            )
          ),
          downloadButton("download", "Download processed data")
        ),
        mainPanel(
          # Desktop instructions (appears in main panel on larger screens)
          div(class = "instructions-desktop instructions-box",
              h4("Instructions:"),
              tags$ul(
                tags$li("Load the sample concentration file(s) from step 2."),
                tags$ul(
                  tags$li("You may have more than one file if some samples were diluted and requanted")
                ),
                tags$li("If you are using this to calculate volumes for pooling, then check the box"),
                tags$li("Add the max volume to transfer per source well in your source plate"),
                tags$ul(
                  tags$li("Underestimate the volume")
                ),
                tags$li("Fill in the remaining values which describe the characteristics of the sample after it is normalized"),
                tags$ul(
                  tags$li("The amount of DNA in the destination well should be estimated based on the mean DNA concentration of the \"samples\" and the desired goals of the normalization")
                ),
                tags$br(),
                tags$li("Evaluate the columns: ", tags$code("transfer_volume"), ", ", tags$code("ul_to_add"), ", ", tags$code("actual_ng_dna"), " with respect to ", tags$code("sample_type")),
                tags$ul(
                  tags$li("It is a good sign if the ", tags$code("transfer_volume"), " is less than the max for \"samples\" and the max for \"neg ctrls\""),
                  tags$li("It is a good sign if the ", tags$code("actual_ng_dna"), " is equal to the target for \"samples\" and less than the target for \"neg ctrls\""),
                  tags$li("It is desirable to maximize the difference between the ", tags$code("actual_ng_dna"), " for \"samples\" and \"neg ctrls\" while ensuring that there is not too much difference among \"samples\""),
                  tags$li("The ", tags$code("ul_to_add"), " column indicates how much water to add or remove via centrivapping, but likely holds little meaning for pooling")
                ),
                tags$li("Adjust settings as necessary and repeat until desired results are achieved"),
                tags$br(),
                tags$li("Optionally edit the output filename prefix"),
                tags$li("Download the processed data")
              )
          ),
          
          # Plot controls
          conditionalPanel(
            condition = "output.show_plot_controls",
            fluidRow(
              column(4,
                     selectInput("plot_y_axis", "Select Y-axis for boxplot:",
                                 choices = NULL,
                                 selected = NULL)
              ),
              column(8, br()) # spacing
            )
          ),
          
          # Boxplot
          conditionalPanel(
            condition = "output.show_plot",
            plotOutput("boxplot", height = "400px"),
            br()
          ),
          
          # Results table
          DTOutput("results_table")
        )
      )
  )
)

server <- function(input, output, session){
  
  # Reactive value to store error state and messages
  error_state <- reactiveValues(
    has_error = FALSE,
    message = ""
  )
  
  # Auto-populate the output prefix when a file is uploaded - NEW
  observeEvent(input$sample_files, {
    if (!is.null(input$sample_files) && nchar(input$outnamePrefix) == 0) {
      # Extract base filename without date and extension
      filename <- input$sample_files$name
      
      # Remove extension
      base_name <- tools::file_path_sans_ext(filename)
      
      # Try to remove common date patterns at the end of the filename
      # Patterns like _2024-01-15, _20240115, _2024_01_15, etc.
      base_name <- str_remove(base_name, "_\\d{4}[-_]?\\d{2}[-_]?\\d{2}$")
      base_name <- str_remove(base_name, "_\\d{8}$")
      
      # Also remove common suffixes that might be before the date
      base_name <- str_remove(base_name, "_(processed|final|output|export|data)$")
      
      updateTextInput(session, "outnamePrefix", value = base_name)
    }
  })
  
  # Reactive values that change based on pooling state
  reactive_inputs <- reactive({
    if (input$pooling) {
      list(
        ul_per_PCR = input$ul_per_PCR_pooling,
        number_PCR_rxns = 1,  # Fixed at 1 for pooling
        DNA_per_PCR = input$DNA_per_PCR_pooling  # Use the user input, not smart defaults
      )
    } else {
      list(
        ul_per_PCR = input$ul_per_PCR,
        number_PCR_rxns = input$number_PCR_rxns,
        DNA_per_PCR = input$DNA_per_PCR
      )
    }
  })
  
  # Update plot y-axis choices when data changes
  observe({
    df <- process_data()
    
    if (!is.null(df) && !error_state$has_error) {
      # Find relevant columns for plotting - more flexible pattern
      plot_cols <- colnames(df) %>% 
        str_subset("(ng_per_ul|pg_per_ul|transfer_volume|actual_ng_dna|actual_pg_dna)")
      
      if (length(plot_cols) > 0) {
        # Create nice labels for the columns
        choice_labels <- setNames(plot_cols, plot_cols %>%
                                    str_replace_all("_", " ") %>%
                                    str_to_title() %>%
                                    str_replace("Ul", "µL") %>%
                                    str_replace("Ng", "ng") %>%
                                    str_replace("Pg", "pg") %>%
                                    str_replace('Pcr', 'RXN'))
        
        # Preserve current selection if it's still valid, otherwise use first option
        current_selection <- input$plot_y_axis
        selected_value <- if (!is.null(current_selection) && current_selection %in% plot_cols) {
          current_selection
        } else {
          str_subset(plot_cols, 'actual')
        }
        
        updateSelectInput(session, "plot_y_axis", 
                          choices = choice_labels,
                          selected = selected_value)
      }
    }
  })
  
  # Show/hide plot controls
  output$show_plot_controls <- reactive({
    df <- process_data()
    !error_state$has_error && !is.null(df) && "sample_type" %in% colnames(df)
  })
  outputOptions(output, "show_plot_controls", suspendWhenHidden = FALSE)
  
  # Show/hide plot
  output$show_plot <- reactive({
    df <- process_data()
    !error_state$has_error && !is.null(df) && 
      "sample_type" %in% colnames(df) && 
      !is.null(input$plot_y_axis)
  })
  outputOptions(output, "show_plot", suspendWhenHidden = FALSE)
  
  # Generate boxplot
  output$boxplot <- renderPlot({
    df <- process_data()
    
    message('\nProcessed Data - preBoxPlot')
    message('DEBUG: Number data rows: ', nrow(df))
    
    if (is.null(df) || error_state$has_error || is.null(input$plot_y_axis)) {
      return(NULL)
    }
    
    if (!"sample_type" %in% colnames(df)) {
      return(NULL)
    }
    
    # Get the selected column name and create a nice label
    selected_col <- input$plot_y_axis
    y_label <- selected_col %>%
      str_replace_all("_", " ") %>%
      str_to_title() %>%
      str_replace("Ul", "µL") %>%
      str_replace("Ng", "ng") %>%
      str_replace("Pg", "pg") %>%
      str_replace('Pcr', 'RXN')
    
    # Create the plot
    message('\nProcessed Data - Just before boxplot')
    message('DEBUG: Number data rows: ', nrow(df))
    p <- ggplot(df, 
                aes(x = sample_type, 
                    y = !!sym(selected_col),
                    fill = sample_type)) +
      geom_boxplot(alpha = 0.7, outliers = FALSE) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      labs(x = "Sample Type", 
           y = y_label,
           title = paste("Distribution of", y_label, "by Sample Type"))
    
    print(p)
    # write_rds(df, 'df.rds')
    # write_rds(p, 'p.rds')
    # write_rds(y_label, 'y_label.rds')
    # write_rds(selected_col, 'selected_col.rds')
  })
  
  # Function to calculate smart defaults for pooling mode
  calculate_pooling_defaults <- reactive({
    if (is.null(input$sample_files)) {
      return(list(DNA_per_PCR = 10))  # fallback default
    }
    
    tryCatch({
      df <- read_csv(input$sample_files$datapath, show_col_types = FALSE) %>%
        group_by(across(any_of(c('dna_extract_id', 'dna_extract_tube_id')))) %>%
        filter(quant_stage == "requant" | !any(quant_stage == "requant")) %>%
        slice_head(n = 1) %>%
        ungroup()
      
      # Find concentration column
      conc_cols <- colnames(df) %>% str_subset("^[pn]g_per_ul_mean$")
      if (length(conc_cols) > 0) {
        conc_col <- conc_cols[1]
        mean_conc <- mean(df[[conc_col]], na.rm = TRUE)
        max_vol <- if (input$pooling) 90 else input$max_vol
        
        # Convert to ng if needed
        if (str_detect(conc_col, "^pg")) {
          mean_conc <- mean_conc / 1000  # convert pg to ng
        }
        
        smart_dna_per_pcr <- mean_conc * max_vol
        return(list(DNA_per_PCR = round(smart_dna_per_pcr, 1)))
      }
      
      return(list(DNA_per_PCR = 10))  # fallback
    }, error = function(e) {
      return(list(DNA_per_PCR = 10))  # fallback
    })
  })
  
  # Update pooling DNA input when file changes or pooling is toggled
  observe({
    if (input$pooling && !is.null(input$sample_files)) {
      defaults <- calculate_pooling_defaults()
      updateNumericInput(session, "DNA_per_PCR_pooling", value = defaults$DNA_per_PCR)
    }
  }) %>% 
    bindEvent(list(input$pooling, input$sample_files, input$max_vol), ignoreNULL = FALSE)
  
  # Function to validate file format and content
  validate_file <- function(file_path) {
    tryCatch({
      # Try to read the CSV file
      df <- read_csv(file_path, show_col_types = FALSE)
      
      # Check if file is empty
      if (nrow(df) == 0) {
        return(list(valid = FALSE, message = "Uploaded file is empty"))
      }
      
      # Check for required columns
      contained_columns <- str_detect(colnames(df), "(^quant_stage$)|(dna_extract(_tube)?_id)")
      if (sum(contained_columns) != 2) {
        required_cols <- c("dna_extract_id", "dna_extract_tube_id", "quant_stage")
        missing_cols <- setdiff(required_cols, colnames(df))
        missing_cols <- str_subset(missing_cols, 'tube', negate = TRUE)
        return(list(valid = FALSE, message = paste("Missing required columns:", paste(missing_cols, collapse = ", "))))
      }
      
      # required_cols <- c("dna_extract_id", "dna_extract_tube_id", "quant_stage")
      # missing_cols <- setdiff(required_cols, colnames(df))
      # if (length(missing_cols) > 0) {
      #   return(list(valid = FALSE, message = paste("Missing required columns:", paste(missing_cols, collapse = ", "))))
      # }
      
      # Check for concentration columns (pg_per_ul_mean or ng_per_ul_mean)
      conc_cols <- colnames(df) %>% str_subset("^[pn]g_per_ul_mean$")
      if (length(conc_cols) == 0) {
        return(list(valid = FALSE, message = "No concentration column found (expected 'pg_per_ul_mean' or 'ng_per_ul_mean')"))
      }
      
      # Check if concentration column has numeric values
      conc_col <- conc_cols[1]
      if (!is.numeric(df[[conc_col]])) {
        return(list(valid = FALSE, message = paste("Concentration column", conc_col, "must contain numeric values")))
      }
      
      # Check for negative or zero concentrations
      if (any(df[[conc_col]] <= 0, na.rm = TRUE)) {
        return(list(valid = FALSE, message = "Concentration values must be positive numbers"))
      }
      
      # Check for missing values in critical columns
      if(any(str_detect(colnames(df), 'dna_extract_id'))){
        if (any(is.na(df$dna_extract_id))) {
          return(list(valid = FALSE, message = "Missing values found in 'dna_extract_id' column"))
        }
      } else if(any(str_detect(colnames(df), 'dna_extract_tube_id'))){
        if (any(is.na(df$dna_extract_tube_id))) {
          return(list(valid = FALSE, message = "Missing values found in 'dna_extract_tube_id' column"))
        }
      }
      
      return(list(valid = TRUE, message = ""))
      
    }, error = function(e) {
      return(list(valid = FALSE, message = paste("File reading error:", e$message)))
    })
  }
  
  process_data <- reactive({
    # Reset error state
    error_state$has_error <- FALSE
    error_state$message <- ""
    
    # Check if file is uploaded - but don't show error on initial load
    if (is.null(input$sample_files)) {
      return(NULL)
    }
    
    # Get reactive input values
    inputs <- reactive_inputs()
    message('\n\n\n\nDEBUG: pooling toggle (TRUE = on, FALSE = off): ', input$pooling)
    
    message('\nJoint Settings')
    message('DEBUG: Max Vol: ', input$max_vol)
    
    message('\nNon-Pooling Settings')
    message('DEBUG: ul per PCR: ', input$ul_per_PCR)
    message('DEBUG: number PCR rxns: ', input$number_PCR_rxns)
    message('DEBUG: DNA per PCR: ', input$DNA_per_PCR)
    
    message('\nPooling Settings')
    message('DEBUG: ul per PCR: ', input$ul_per_PCR_pooling)
    message('DEBUG: number PCR rxns: ', input$number_PCR_rxns_pooling)
    message('DEBUG: DNA per PCR: ', input$DNA_per_PCR_pooling)
    
    # Validate file format and content
    validation <- validate_file(input$sample_files$datapath)
    if (!validation$valid) {
      error_state$has_error <- TRUE
      error_state$message <- validation$message
      return(NULL)
    }
    
    tryCatch({
      dna_concentration <- read_csv(input$sample_files$datapath, show_col_types = FALSE) %>%
        # keep only the original quant OR (if present) its requant
        group_by(across(any_of(c('dna_extract_id', 'dna_extract_tube_id')))) %>%
        filter(quant_stage == "requant" | !any(quant_stage == "requant")) %>%
        slice_head(n = 1) %>%
        ungroup()
      
      message('\nRaw Data')
      message('DEBUG: Number data rows: ', nrow(dna_concentration))
      
      ## 1. pick the concentration column
      conc_cols     <- colnames(dna_concentration) %>% str_subset("^[pn]g_per_ul_mean$")
      dna_conc_var  <- conc_cols[1]   # take the first match
      
      ## 2. convert requested DNA-per-PCR to pg if the file is in pg/µL
      dna_per_pcr <- if (str_detect(dna_conc_var, "^pg")) inputs$DNA_per_PCR * 1000 else inputs$DNA_per_PCR
      
      ## 3. do the volume maths
      # Get max_vol - use default if pooling is selected
      # max_vol_value <- if (input$pooling) 90 else input$max_vol
      # max_vol_value <- input$max_vol
      
      dna_concentration <- dna_concentration %>%
        mutate(goal_volume_ul = (inputs$number_PCR_rxns * dna_per_pcr) / .data[[dna_conc_var]],
               transfer_volume = pmin(goal_volume_ul, input$max_vol),
               ul_to_add       = inputs$number_PCR_rxns * inputs$ul_per_PCR - transfer_volume,
               !!str_glue("actual_{str_extract(dna_conc_var, '^[pn]g')}_dna_per_pcr") :=
                 .data[[dna_conc_var]] * transfer_volume / inputs$number_PCR_rxns) %>%
        select(-goal_volume_ul)
      
      # Additional validation on processed data
      if (nrow(dna_concentration) == 0) {
        error_state$has_error <- TRUE
        error_state$message <- "No valid samples found after processing"
        return(NULL)
      }
      
      return(dna_concentration)
      
    }, error = function(e) {
      error_state$has_error <- TRUE
      error_state$message <- paste("Data processing error:", e$message)
      return(NULL)
    })
  })
  
  output$results_table <- renderDT({
    df <- process_data()
    
    message('\nProcessed Data')
    message('DEBUG: Number data rows: ', nrow(df))
    
    # Don't show table if there's an error
    if (error_state$has_error || is.null(df)) {
      return(NULL)
    }
    
    # Round numeric columns for display
    num <- vapply(df, is.numeric, logical(1))
    df[num] <- lapply(df[num], round, 1)
    datatable(df, options = list(pageLength = 25))
  })
  
  # Show/hide error messages
  output$show_error <- reactive({
    error_state$has_error
  })
  outputOptions(output, "show_error", suspendWhenHidden = FALSE)
  
  output$error_message <- renderText({
    error_state$message
  })
  
  output$download <- downloadHandler(
    filename = function() {
      # Use user-provided prefix or fall back to default
      prefix <- if (nchar(input$outnamePrefix) > 0) {
        input$outnamePrefix
      } else {
        "processed_samples"
      }
      glue::glue("{prefix}_{Sys.Date()}.csv")
    },
    content = function(file) {
      # Check if file is uploaded first
      if (is.null(input$sample_files)) {
        showNotification(
          "Please upload a CSV file before downloading",
          type = "error",
          duration = NULL,
          closeButton = TRUE
        )
        return(NULL)
      }
      
      # Perform validation before download
      df <- process_data()
      
      if (error_state$has_error || is.null(df)) {
        # Show error notification
        showNotification(
          paste("Cannot download:", if(error_state$message != "") error_state$message else "No valid data to download"),
          type = "error",
          duration = NULL,
          closeButton = TRUE
        )
        return(NULL)
      }
      
      tryCatch({
        write_csv(df, file)
        showNotification("File downloaded successfully!", type = "message",
                         duration = NULL,
                         closeButton = TRUE)
      }, error = function(e) {
        showNotification(
          paste("Download failed:", e$message),
          type = "error",
          duration = NULL,
          closeButton = TRUE
        )
      })
    }
  )
}

shinyApp(ui, server)