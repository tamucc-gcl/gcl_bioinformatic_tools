library(shiny)
library(DT)
library(readr)
library(dplyr)

ui <- fluidPage(
  titlePanel("PCR Sample Concentration Processor"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("sample_files", "Upload Sample Concentration CSV(s)", 
                multiple = TRUE, accept = ".csv"),
      numericInput("ul_per_PCR", "µL per PCR:", value = 2, min = 0, step = 0.5),
      numericInput("number_PCR_rxns", "Number of PCR reactions:", value = 36, min = 1),
      numericInput("DNA_per_PCR", "DNA per PCR (ng):", value = 10, min = 0, step = 0.1),
      numericInput("max_vol", "Max volume per sample (µL):", value = 90, min = 1),
      downloadButton("download", "Download Processed Data")
    ),
    
    mainPanel(
      DTOutput("results_table")
    )
  )
)

server <- function(input, output, session) {
  
  process_data <- reactive({
    req(input$sample_files)
    
    # message('DEBUG: ', str_c(input$sample_files$datapath, collapse = '; '))
    
    dna_concentration <- read_csv(input$sample_files$datapath, 
                                  show_col_types = FALSE) %>%
      #If there is a requant remove the original otherwise keep original
      group_by(dna_extract_tube_id) %>%
      filter(quant_stage == "requant" | !("requant" %in% quant_stage)) %>%
      slice_head(n = 1) %>%
      ungroup()
    
    # message('DEBUG: ', str_c(colnames(dna_concentration), collapse = '; '))
    # 
    # message('DEBUG: ', input$DNA_per_PCR)
    # message('DEBUG: ', input$max_vol)
    # message('DEBUG: ', input$ul_per_PCR)
    # message('DEBUG: ', input$number_PCR_rxns)
    
    # Perform calculations
    dna_concentration <- dna_concentration %>%
      mutate(goal_volume_ul = (input$number_PCR_rxns * input$DNA_per_PCR) / ng_per_ul_mean) %>%
      mutate(transfer_volume = if_else(goal_volume_ul > input$max_vol,
                                       input$max_vol,
                                       goal_volume_ul),
             ul_to_add = (input$number_PCR_rxns * input$ul_per_PCR) - transfer_volume,
             actual_ng_dna_per_pcr = (ng_per_ul_mean * transfer_volume) / input$number_PCR_rxns) %>%
      select(-goal_volume_ul)
    
    dna_concentration
  })
  
  output$results_table <- renderDT({
    df <- process_data()
    numeric_cols <- sapply(df, is.numeric)
    df[numeric_cols] <- round(df[numeric_cols], 1)
    df
  })
  
  output$download <- downloadHandler(
    filename = function() { paste0("processed_samples_", Sys.Date(), ".csv") },
    content = function(file) {
      write_csv(process_data(), file)
    })
}

shinyApp(ui, server)
