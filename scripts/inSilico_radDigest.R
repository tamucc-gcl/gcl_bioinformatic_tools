## Perform in Silico digestion of genome using pairs of restriction enzymes
##TODO - Finalize output tables & graphs
##TODO - Run digestions in parallel
##TODO - implement more than just ddRAD
##TODO - add in custom restriction enzymes - i.e. MspI = 'C/CGG'

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
  genome_file <- '~/../Downloads/GCA_034783695.1_ASM3478369v1_genomic.fna'
  restriction_enzymes <- c("NlaIII", 'MluCI')
  size_window <- c(300, 500)
}

#### Libraries ####
library(DECIPHER) |> suppressMessages() |> suppressWarnings()
library(tidyverse) |> suppressMessages() |> suppressWarnings()
library(ggtext) |> suppressMessages() |> suppressWarnings()

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
RESTRICTION_ENZYMES <- c(RESTRICTION_ENZYMES, 'MspI' = 'C/CGG')

the_genome <- readDNAStringSet(genome_file)

# the_genome <- the_genome[width(the_genome) > 10000]
sum(width(the_genome))
sum(width(the_genome) > 10000)

tibble(contig_length = width(the_genome)) %>%
  arrange(-contig_length) %>%
  mutate(cum_length = cumsum(contig_length),
         id = row_number()) %>%
  ggplot(aes(x = id,
             y = cum_length)) +
  geom_line() +
  geom_vline(xintercept = sum(width(the_genome) > 10000),
             linetype = 'dashed') +
  geom_hline(yintercept = sum(width(the_genome[width(the_genome) > 10000])),
             linetype = 'dashed') +
  scale_x_continuous(labels = scales::comma_format()) +
  scale_y_continuous(labels = scales::comma_format(scale = 1/1e6, suffix = 'M')) +
  labs(x = 'Fragment ID',
       y = 'Cumulative Genome Length',
       title = 'Contig length distribution',
       caption = str_c('Reference genome used: ', 
                       str_extract(genome_file, 'GC[FA]_.*\\.[0-9]'))) +
  theme_classic(base_size = 16) +
  theme(panel.background = element_rect(colour = 'black'),
        strip.background = element_blank(),
        plot.title = element_markdown())

#### Digestions #### 
processed_digestion <- expand_grid(enzyme1 = restriction_enzymes,
                                   enzyme2 = restriction_enzymes) %>%
  filter(enzyme1 < enzyme2) %>%
  
  #Temp addition - need to rework to make possible to specify special enzyme
  # filter(enzyme1 == 'EcoRI', enzyme2 == 'MspI') %>%
  rowwise(enzyme1, enzyme2) %>%
  reframe(perform_digest(the_genome, 
                         RESTRICTION_ENZYMES[c(enzyme1, enzyme2)],
                         # c(cut_site1, cut_site2),
                         digestion_protocol))

#### Outputs ####
n_fragments <- processed_digestion %>%
  filter(frag_length >= min(size_window),
         frag_length <= max(size_window)) %>%
  summarise(n = n(),
            .by = c(enzyme1, enzyme2)) %>%
  mutate(position = mean(size_window))

processed_digestion %>%
  filter(frag_length > 200,
         frag_length < 1000) %>%
  ggplot(aes(x = frag_length)) +
  geom_histogram(binwidth = 25) +
  geom_vline(xintercept = size_window, 
             linetype = 'dashed',
             colour = 'red') +
  geom_text(data = n_fragments,
            aes(x = position, 
                y = Inf,
                label = str_c('Number of target\nsize fragments\n',
                              scales::comma(n, accuracy = 1))),
            vjust = 1.5) +
  facet_grid(enzyme1 ~ enzyme2) +
  scale_x_continuous(labels = scales::comma_format()) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = 'Fragment Length',
       y = 'Number of Fragments',
       title = 'ddRAD <i>in silico</i> digestion of <i>Sardinella gibbosa</i> ',
       caption = str_c('Reference genome used: ', 
                       str_extract(genome_file, 'GC[FA]_.*\\.[0-9]'))) +
  theme_classic(base_size = 16) +
  theme(panel.background = element_rect(colour = 'black'),
        strip.background = element_blank(),
        plot.title = element_markdown())

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



