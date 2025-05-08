library(tidyverse)
library(janitor)
library(readxl)


sample_metadata_file <- '../example_data/rbd_edna-extraction_plate-map_updated.xlsx'
sheet_metadata <- "Plate-map-tidy"
sample_concentration_file <- list.files('../example_data/dna_concentrations', 
                                pattern = '*csv$',
                                full.names = TRUE)

##PCR Setting
ul_per_PCR <- 2 #Z2

##PCR Goals Settings
number_PCR_rxns <- 36 #X2
DNA_per_PCR <- 10 #Y2

## Safety Settings
max_vol <- 90 #U2

#### Functions ####
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
    rename(dna_plate_id = plate_id,
           dna_plate_col = plate_col,
           dna_plate_row = plate_row)
}

#### Data ####
dna_plate_map <- read_plate_map(sample_metadata_file, sheet_metadata)

dna_concentration <- read_csv(sample_concentration_file, 
                         show_col_types = FALSE)  %>%
  #If there is a requant remove the original otherwise keep original
  group_by(dna_extract_tube_id) %>%
  filter(quant_stage == "requant" | !("requant" %in% quant_stage)) %>%
  slice_head(n = 1) %>%
  ungroup()


target_dna_quantity <- number_PCR_rxns * DNA_per_PCR
eluted_volume <- number_PCR_rxns * ul_per_PCR

dna_concentration %>%
  select(dna_extract_tube_id, ng_per_ul_mean) %>%
  mutate(goal_volume_ul = target_dna_quantity / ng_per_ul_mean) %>%
  mutate(transfer_volume = if_else(goal_volume_ul > max_vol,
                                   max_vol,
                                   goal_volume_ul),
         ul_to_add = eluted_volume - transfer_volume,
         actual_ng_dna_per_pcr = (ng_per_ul_mean * transfer_volume) / number_PCR_rxns) %>%
  select(-goal_volume_ul) %>%
  
  ggplot(aes(x = actual_ng_dna_per_pcr)) +
  geom_histogram()
