library(shiny)
library(DT)
library(readr)
library(dplyr)
library(stringr)   # you only need stringr, not the full tidyverse

ui <- fluidPage(
  titlePanel("PCR Sample Concentration Processor"),
  
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
  
  sidebarLayout(
    sidebarPanel(
      fileInput("sample_files", "Upload Sample-concentration CSV", accept = ".csv"),
      numericInput("ul_per_PCR",     "µL per PCR:",          value = 2,  min = 0,  step = 0.5),
      numericInput("number_PCR_rxns","Number of PCR reactions:", value = 36, min = 1),
      numericInput("DNA_per_PCR",    "DNA per PCR (ng):",    value = 10, min = 0,  step = 0.1),
      numericInput("max_vol",        "Max volume per sample (µL):", value = 90, min = 1),
      downloadButton("download", "Download processed data")
    ),
    mainPanel(DTOutput("results_table"))
  )
)

server <- function(input, output, session){
  
  process_data <- reactive({
    req(input$sample_files)
    
    dna_concentration <- read_csv(input$sample_files$datapath, show_col_types = FALSE) %>%
      # keep only the original quant OR (if present) its requant
      group_by(dna_extract_tube_id) %>%
      filter(quant_stage == "requant" | !any(quant_stage == "requant")) %>%
      slice_head(n = 1) %>%
      ungroup()
    
    ## 1. pick the concentration column
    conc_cols     <- colnames(dna_concentration) %>% str_subset("^[pn]g_per_ul_mean$")
    dna_conc_var  <- conc_cols[1]   # take the first match
    
    ## 2. convert requested DNA-per-PCR to pg if the file is in pg/µL
    dna_per_pcr <- if (str_detect(dna_conc_var, "^pg")) input$DNA_per_PCR * 1000 else input$DNA_per_PCR
    
    ## 3. do the volume maths
    dna_concentration <- dna_concentration %>%
      mutate(goal_volume_ul = (input$number_PCR_rxns * dna_per_pcr) / .data[[dna_conc_var]],
             transfer_volume = pmin(goal_volume_ul, input$max_vol),
             ul_to_add       = input$number_PCR_rxns * input$ul_per_PCR - transfer_volume,
             !!str_glue("actual_{str_extract(dna_conc_var, '^[pn]g')}_dna_per_pcr") :=
               .data[[dna_conc_var]] * transfer_volume / input$number_PCR_rxns) %>%
      select(-goal_volume_ul)
    
    dna_concentration
  })
  
  output$results_table <- renderDT({
    df <- process_data()
    num <- vapply(df, is.numeric, logical(1))
    df[num] <- lapply(df[num], round, 1)
    datatable(df, options = list(pageLength = 25))
  })
  
  output$download <- downloadHandler(
    filename = function() glue::glue("processed_samples_{Sys.Date()}.csv"),
    content  = function(file) write_csv(process_data(), file)
  )
}

shinyApp(ui, server)
