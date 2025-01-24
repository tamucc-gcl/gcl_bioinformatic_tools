## Perform in Silico digestion of genome using pairs of restriction enzymes
##TODO - Finalize output tables & graphs
##TODO - Run digestions in parallel
##TODO - implement more than just ddRAD


if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  outdir <- args[1]
  digestion_protocol <- args[2] #ddrad
  genome_file <- args[3]
  restriction_enzymes <- args[4] #Named following New England Biolabs Naming
  size_window <- as.integer(args[5])
  
  # .libPaths(.libPaths()[grep(pattern = 'conda', .libPaths())]) #Fix strange bug with sbatch jobs
} else {
  outdir <- '~/Google Drive/TAMUCC-CORE/tmp_dir'
  digestion_protocol <- 'ddrad'
  genome_file <- '~/Google Drive/TAMUCC-CORE/tmp_dir/methylation_rad/GCF_002234675.1_ASM223467v1_genomic.fna'
  restriction_enzymes <- c("EcoRI", 'HindIII', 'BstBI', 'SphI')
  size_window <- c(400, 600)
}

#### Libraries ####
library(DECIPHER) |> suppressMessages() |> suppressWarnings()
library(tidyverse) |> suppressMessages() |> suppressWarnings()

#### Functions ####
digest_to_tibble <- function(dna_string_set){
  tibble(name = names(dna_string_set),
         frag_length = width(dna_string_set)) %>%
    mutate(side = str_extract(name, 'top|bottom'),
           name = str_remove(name, side))
}

ddrad_digestion <- function(genome, enzymes){
  # in Silico digestion to only have sequences with both restriction enzymes not the same one on each side.
  # i.e. ddRAD
  # genome <- the_genome; enzymes <- RESTRICTION_ENZYMES[restriction_enzymes]
  
  digestion <- DigestDNA(enzymes,
                         genome,
                         type = "fragments",
                         strand = "top") %>%
    unlist()
  full_digestion <- digestion
  
  for(cutter in enzymes){
    partial_digest <- DigestDNA(cutter,
                                the_genome,
                                type = "fragments",
                                strand = "top") %>%
      unlist()
    
    digestion <- digestion[!digestion %in% partial_digest]
  }
  
  digestion 
}

# ezrad_digestion

# ograd_digestion

perform_digest <- function(genome, enzymes, digest_type){
  if(str_to_lower(digest_type) == 'ddrad'){
    the_digestion <- ddrad_digestion(genome, enzymes)
  }
  
  digest_to_tibble(the_digestion)
}

#### Prep inputs ####
data(RESTRICTION_ENZYMES)

the_genome <- readDNAStringSet(genome_file)

#### Digestions #### 
processed_digestion <- expand_grid(enzyme1 = restriction_enzymes,
            enzyme2 = restriction_enzymes) %>%
  filter(enzyme1 < enzyme2) %>%
  rowwise(enzyme1, enzyme2) %>%
  reframe(perform_digest(the_genome, 
                         RESTRICTION_ENZYMES[c(enzyme1, enzyme2)], 
                         digestion_protocol))

#### Outputs ####
processed_digestion %>%
  filter(frag_length > 200,
         frag_length < 1000) %>%
  ggplot(aes(x = frag_length)) +
  geom_histogram(binwidth = 25) +
  geom_vline(xintercept = size_window, 
             linetype = 'dashed',
             colour = 'red') +
  facet_grid(enzyme1 ~ enzyme2)

processed_digestion %>%
  filter(frag_length > 200,
         frag_length < 1000) %>%
  mutate(bins = cut_width(frag_length, width = 25, boundary = 0)) %>%
  summarise(n_frag = n(),
            length = sum(frag_length),
            mid_length = median(frag_length),
            .by = c(enzyme1, enzyme2, bins)) %>%
  ggplot(aes(x = mid_length, y = n_frag)) +
  geom_point() +
  facet_grid(enzyme1 ~ enzyme2)
