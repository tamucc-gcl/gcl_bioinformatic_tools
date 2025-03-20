##TODO - 

if(!interactive()){
  args <- commandArgs(trailingOnly = TRUE)
  rainbowbridge_dir <- args[1]
  param_yaml <- args[2]
  .libPaths(.libPaths()[grep(pattern = 'conda', .libPaths())]) #Fix strange bug with sbatch jobs
} else {
  rainbowbridge_dir <- '/scratch/group/p.bio240270.000/prj_sheehy-metabarcoding/testing_highPIDENT_DIFF_lowerPIDENT/'
  param_yaml <- 'metabarcode_rainbowbridge_paired.yml'
}

setwd(rainbowbridge_dir)

library(ShortRead) |> suppressMessages() |> suppressWarnings()
library(yaml) |> suppressMessages() |> suppressWarnings()
library(tidyverse) |> suppressMessages() |> suppressWarnings()
library(ggtext) |> suppressMessages() |> suppressWarnings()
library(ggnested) |> suppressMessages() |> suppressWarnings()
library(SummarizedExperiment) |> suppressMessages() |> suppressWarnings()

# lca_table <- read_delim('output/taxonomy/lca/qcov70_pid70_diff1/lca_taxonomy.tsv')
# zotu_table <- read_delim('output/zotus/zotu_table.tsv')
# lulu_table <- read_delim('output/lulu/lulu_zotu_table.tsv')
# 
# read_delim('output/taxonomy/lca/qcov70_pid70_diff1/lca_intermediate.tsv') %>%
#   filter(zotu == 'Zotu1')

rainbowbridge_yaml <- read_yaml(param_yaml)
zotus_final <- read_delim('output/final/zotu_table_final_curated.tsv', 
                          show_col_types = FALSE)

blast_hits <- str_c('output/blast/pid', rainbowbridge_yaml$`percent-identity`, 
                    '_eval', rainbowbridge_yaml$evalue, 
                    '_qcov', rainbowbridge_yaml$qcov, 
                    '_max', rainbowbridge_yaml$`max-query-results`, 
                    '/blast_result_merged.tsv') %>%
  read_delim(col_names = c('zotu', 'seqid', 'taxid', 'blast_species', 'commonname', 'blast_domain',
                           'pident', 'length', 'qlen', 'slen', 'mismatch', 'gapopen', 'gaps',
                           'qstart', 'wend', 'sstart', 'send', 'stitle', 'evalue', 'bitscore',
                           'qcov', 'qcovhsp'), 
             show_col_types = FALSE)

lca_taxonomy <- str_c('output/taxonomy/lca/qcov', rainbowbridge_yaml$`lca-qcov`, 
                      '_pid', rainbowbridge_yaml$`lca-pid`, 
                      '_diff', rainbowbridge_yaml$`lca-diff`, 
                      '/lca_taxonomy.tsv') %>%
  read_delim(show_col_types = FALSE)


final_zotu_sequences <- readDNAStringSet(list.files(path = 'output/zotus', pattern = '_zotus.fasta$', full.names = TRUE))[zotus_final$zotu]

#### Functions ####
get_sample_readcounts <- function(sample_df, sample_dir){
  read_counts <- countFastq(sample_dir, 
                            pattern = 'fastq(.gz)?') %>%
    as_tibble(rownames = 'read_id') %>%
    dplyr::rename(reads = records)
  
  out <- left_join(sample_df,
                   select(read_counts, 
                          read_id, reads, nucleotides),
                   by = c('r1' = 'read_id'))
  
  if(any(colnames(sample_df) == 'r2')){
    out <-  out %>%
      left_join(select(read_counts, 
                       read_id, reads, nucleotides),
                by = c('r2' = 'read_id'),
                suffix = c('.r1', '.r2'))
  }
  out
}


output_sunburst_OLD <- function(zotu_data, sample_names, yvar, outfile){
  ## Make the sunburst plot
  
  #yvar can be either "n_zotu" or "total_reads"
  sun_data <- zotu_data %>%
    group_by(across(domain:species)) %>%
    select(-unique_hits) %>%
    summarise(n_zotu = n_distinct(zotu),
              across(where(is.numeric), sum),
              .groups = 'drop') %>%
    rowwise %>%
    mutate(total_reads = sum(c_across(any_of(sample_names))),
           .keep = 'unused') %>%
    ungroup %>%
    mutate(across(where(is.character), 
                  ~if_else(. == 'LCA_dropped', NA_character_, .)),
           
           across(where(is.character),
                  ~str_replace_na(., replacement = '')),
           
           path = str_c(domain, kingdom, phylum, class, order, family, genus, species, sep = '-') %>%
             str_remove('-+$')) %>%
    select(path, !!sym(yvar))
  
  sundata_tmp <- tempfile(fileext = ".csv")
  write_csv(sun_data, sundata_tmp)
  
  ## Output file
  library(rmarkdown)
  library(pagedown)
  
  vector_Text_RMD <- c('---',
                       str_c('title: "', yvar, '"'),
                       'output: html_document',
                       '---',
                       '```{r setup, include=FALSE}',
                       'knitr::opts_chunk$set(echo = TRUE)',
                       '```',
                       '```{r cars, echo=FALSE}',
                       'library(sunburstR)',
                       str_c('df <- read.csv("', sundata_tmp, '", header = FALSE,stringsAsFactors = FALSE)'),
                       'sun <- sunburst(df, legend = FALSE)',
                       'htmlwidgets::onRender(sun,',
                       "    'function(el, x) {",
                       '    d3.selectAll(".sunburst-legend text").attr("font-size", "10px");',
                       '    d3.select(el).select(".sunburst-togglelegend").property("checked",true); // force show the legend, check legend',
                       '    d3.select(el).select(".sunburst-togglelegend").on("click")(); // simulate click',
                       '    d3.select(el).select(".sunburst-togglelegend").remove() // remove the legend toggle',
                       '    }',
                       "    ')",
                       '```')
  
  zzfil <- tempfile(fileext = ".Rmd")
  writeLines(text = vector_Text_RMD, con = zzfil)
  
  render(input = zzfil,
         output_file = outfile,
         output_dir = rainbowbridge_dir)
}

output_sunburst <- function(zotu_data, sample_names, yvar, outfile){
  ## Make the sunburst plot
  
  #yvar can be either "n_zotu" or "total_reads"
  sun_data <- zotu_data %>%
    group_by(across(domain:species)) %>%
    select(-unique_hits) %>%
    summarise(n_zotu = n_distinct(zotu),
              across(where(is.numeric), sum),
              .groups = 'drop') %>%
    rowwise %>%
    mutate(total_reads = sum(c_across(any_of(sample_names))),
           .keep = 'unused') %>%
    ungroup %>%
    mutate(across(where(is.character), 
                  ~if_else(. == 'LCA_dropped', NA_character_, .)),
           
           across(where(is.character),
                  ~str_replace_na(., replacement = '')),
           
           path = str_c(domain, kingdom, phylum, class, order, family, genus, species, sep = '-') %>%
             str_remove('-+$')) %>%
    select(path, !!sym(yvar))
  
  library(htmlwidgets); library(sunburstR)
  sun <- sunburst(sun_data, legend = FALSE)
  
  # Export the widget as SVG using the `saveWidget` function and PhantomJS
  saveWidget(sun, file = str_c(rainbowbridge_dir, '/', outfile, '.html'), selfcontained = TRUE)
  
  # Convert HTML to SVG
  # rsvg::rsvg_svg(str_c(rainbowbridge_dir, '/', outfile, '.html'), 
  #                str_c(rainbowbridge_dir, '/', outfile, '.svg'))
  # file.remove(str_c(rainbowbridge_dir, '/', outfile, '.html'))
}

ggnested_jds <- function (data, 
                          mapping = aes(), 
                          ..., 
                          legend_labeling = c("sub", "join", "main"), 
                          join_str = " - ", 
                          legend_title = NULL, 
                          main_keys = TRUE, 
                          nested_aes = c("fill", "color"), 
                          gradient_type = c("both", "shades", "tints"), 
                          min_l = 0.05, 
                          max_l = 0.95, 
                          main_palette = NULL, 
                          base_clr = "#008CF0", 
                          NA_option = 'black') 
{
  #Taken from https://github.com/gmteunisse/ggnested modified to allow NA colour options
  aes_args <- names(mapping)
  if (!"main_group" %in% aes_args) {
    stop("Error: provide the main_group in the aesthetic mapping argument. For non-nested data, use the regular ggplot2 function.")
  }
  if (!"sub_group" %in% aes_args) {
    stop("Error: provide a subgroup in the aesthetic mapping argument. For non-nested data, use the regular ggplot2 function.")
  }
  if ("fill" %in% aes_args & "fill" %in% nested_aes) {
    warning("Warning: fill aesthetics will be ignored in the main ggnested function. Please specify non-nested fill in the geom_* layer. Alternatively,\n            remove 'fill' from mapping_aes.")
    mapping$fill <- NULL
  }
  if (("colour" %in% aes_args | "color" %in% aes_args) & ("colour" %in% 
                                                          nested_aes | "color" %in% nested_aes)) {
    warning("Warning: colour aesthetics will be ignored in the main ggnested function. Please specify non-nested colour in the geom_* layer. Alternatively,\n            remove 'colour' from mapping_aes.")
    mapping$colour <- NULL
    mapping$color <- NULL
  }
  group <- quo_name(mapping$main_group)
  subgroup <- quo_name(mapping$sub_group)
  pal <- nested_palette(data, group, subgroup, gradient_type, 
                        min_l, max_l, main_palette, base_clr, join_str) %>%
    mutate(group_subgroup = if_else(is.na(!!sym(subgroup)), 
                                    str_replace(group_subgroup, '- NA$', '- Unknown'),
                                    group_subgroup))
  
  colours <- pal %>%
    filter(!(is.na(!!sym(group)) & is.na(!!sym(subgroup)))) %>%
    rename(sublabel = !!subgroup, label = !!group) %>% 
    as.data.frame()
  
  if (main_keys) {
    colours <- colours %>%
      group_by(label) %>% 
      group_modify(~add_row(.x, .before = 0)) %>% 
      ungroup() %>% 
      mutate(subgroup_colour = ifelse(is.na(subgroup_colour), "#FFFFFF", subgroup_colour), 
             sublabel = case_when(is.na(sublabel) & is.na(group_subgroup) ~ sprintf("**%s**", as.character(label)),
                                  is.na(sublabel) & !is.na(group_subgroup) ~ 'Unknown',
                                  TRUE ~ as.character(sublabel)),
             group_subgroup = ifelse(is.na(group_subgroup), sprintf("**%s**", as.character(label)), group_subgroup)) %>% 
      as.data.frame()
  }
  vals <- colours$subgroup_colour
  names(vals) <- colours$group_subgroup
  df <- left_join(data, pal, by = c(group, subgroup)) %>% 
    arrange(group, subgroup) %>% 
    mutate(group_subgroup = factor(group_subgroup, ordered = T, levels = colours$group_subgroup), 
           `:=`(!!subgroup, factor(!!sym(subgroup), ordered = T)), 
           `:=`(!!group, factor(!!sym(group), ordered = T))) %>% 
    ungroup() %>% 
    arrange(group_subgroup)
  
  if (legend_labeling[1] == "join") {
    labels <- colours$group_subgroup
    leg_title <- sprintf("%s%s%s", group, join_str, subgroup)
  } else if (legend_labeling[1] == "main") {
    labels <- colours$label
    leg_title <- group
  } else if (legend_labeling[1] == "sub") {
    labels <- colours$sublabel
    leg_title <- subgroup
  } else {
    stop("Invalid option for legend_labeling. Pick one of c('join', 'main', 'sub')")
  }
  
  if (!is.null(legend_title)) {
    leg_title <- legend_title
  }
  
  nested_scale <- scale_discrete_manual(..., aesthetics = nested_aes, 
                                        name = leg_title, values = vals, labels = labels, drop = F,
                                        na.value = NA_option)
  if ("fill" %in% nested_aes) {
    mapping$fill <- quo(group_subgroup)
  }
  if ("colour" %in% nested_aes | "color" %in% nested_aes) {
    mapping$colour <- quo(group_subgroup)
  }
  p <- ggplot(df, mapping, ...) + nested_scale
  if (main_keys) {
    p <- p + theme_nested(theme)
  }
  return(p)
}

#### Get Summary Numbers ####
raw_zotu_table <- read_delim('output/zotus/zotu_table.tsv',
                             show_col_types = FALSE)

lulu_table <- read_delim('output/lulu/lulu_zotu_table.tsv',
                         show_col_types = FALSE)

summary_metrics <- tibble(maxhits = rainbowbridge_yaml$`max-query-results`,
                          qcov = rainbowbridge_yaml$qcov,
                          pid = rainbowbridge_yaml$`percent-identity`,
                          eval = rainbowbridge_yaml$evalue,
                          lca_qcov = rainbowbridge_yaml$`lca-qcov`,
                          lca_pid = rainbowbridge_yaml$`lca-pid`,
                          lca_diff = rainbowbridge_yaml$`lca-diff`,
                          n_zotu_raw = n_distinct(raw_zotu_table$`#OTU ID`),
                          n_hits_blast = nrow(blast_hits),
                          nhits_per_rawzotu = n_hits_blast / n_zotu_raw,
                          n_taxa_blast = n_distinct(blast_hits$blast_species),
                          n_zotu_lulu = n_distinct(lulu_table$zotu),
                          n_zotu_lca = n_distinct(lca_taxonomy$zotu),
                          n_zotu_final = n_distinct(zotus_final$zotu),
                          n_taxa_final = select(zotus_final, domain:species) %>%
                            distinct %>%
                            nrow)

summary_metrics <- count(zotus_final, taxid_rank) %>%
  pivot_wider(names_from = 'taxid_rank',
              values_from = 'n',
              names_prefix = 'n_') %>%
  bind_cols(summary_metrics, .)

write_csv(summary_metrics, 'summary_metrics.csv')


#### Get Sample Info ####
if(rainbowbridge_yaml$`demultiplexed-by` == 'index'){
  #If reads were already demultiplexed
  samples <- read_delim(rainbowbridge_yaml$`sample-map`, 
                        show_col_types = FALSE,
                        col_names = c('sample_id', 'r1', 'r2')) %>%
    mutate(sample_id = str_replace_all(sample_id, '-', '_')) %>%
    get_sample_readcounts(rainbowbridge_yaml$reads) %>%
    mutate(in_final = sample_id %in% colnames(zotus_final)) %>%
    rename(reads = reads.r1)
  
} else {
  
  samples <- list.files(path = 'preprocess/split_samples', pattern = 'fastq') %>%
    tibble(r1 = .) %>%
    mutate(sample_id = str_remove(r1, '__split__') %>%
             str_remove('.fastq')) %>%
    full_join(read_delim(rainbowbridge_yaml$metadata, 
                         show_col_types = FALSE),
              by = c('sample_id' = 'sample')) %>%
    select(-sample_name) %>%
    get_sample_readcounts('preprocess/split_samples') %>%
    mutate(in_final = sample_id %in% colnames(zotus_final),
           reads = replace_na(reads, 0L))
  
}



samples %>%
  select(sample_id, 
         reads, 
         in_final) %>%
  mutate(notes = if_else(in_final, NA_character_, 'Removed from Further Analyses'),
         .keep = 'unused') %>%
  write_csv('sample_sequencing_summary.csv', na = '')

sample_filtering_plot <- samples %>%
  ggplot(aes(x = in_final, y = reads)) +
  geom_boxplot() +
  geom_label(data = . %>% 
               summarise(n = n_distinct(sample_id), 
                         IDs = str_c(sample_id, collapse = '; '),
                         .by = in_final),
             aes(y = Inf, 
                 label = if_else(in_final,
                                 as.character(n),
                                 str_c(n, IDs, sep = '\n'))),
             vjust = 1, size = 12) + 
  scale_y_continuous(trans = scales::log10_trans(),
                     labels = scales::comma_format(),
                     limits = c(1, NA)) +
  scale_x_discrete(breaks = c('TRUE', 'FALSE'),
                   labels = c('Samples Retained\nin Analysis', 
                              'Samples Removed\nDue to Lack of Reads'),
                   limits = c('TRUE', 'FALSE')) +
  labs(x = NULL,
       y = 'log<sub>10</sub>(# Reads)') +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        panel.background = element_rect(colour = 'black'),
        legend.position = 'right',
        legend.text = element_markdown(),
        legend.key = element_blank())
ggsave('sample_filtering.png', 
       plot = sample_filtering_plot,
       height = 7,
       width = 7)

reads_per_sample_plot <- samples %>%
  mutate(sample_id = fct_reorder(sample_id, reads)) %>%
  ggplot(aes(y = sample_id, x = reads)) +
  geom_col() + 
  geom_vline(xintercept = median(samples$reads), colour = 'red') +
  scale_x_continuous(labels = scales::comma_format(),
                     trans = scales::log10_trans()) +
  annotation_logticks(sides = 'b') +
  labs(y = NULL,
       x = 'log<sub>10</sub>(# Reads)') +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        panel.background = element_rect(colour = 'black'),
        legend.position = 'right',
        legend.text = element_markdown(),
        legend.key = element_blank())
ggsave('reads_per_sample.png', 
       plot = reads_per_sample_plot,
       height = 15,
       width = 7)


#### Summarize Taxonomic Results ####
output_sunburst(zotus_final, samples$sample_id, 'n_zotu', 'taxonomy_nZOTU')
output_sunburst(zotus_final, samples$sample_id, 'total_reads', 'taxonomy_nReads')

#### Visualize Classification Quality ####
summarized_blast <- select(blast_hits, 
                           zotu, pident, qcov) %>%
  summarise(across(c(pident, qcov), list(mean = mean, median = median, 
                                         min = min, max = max)),
            .by = 'zotu')

summarized_blast_plot <- select(zotus_final, 
                                zotu, taxid_rank) %>%
  left_join(summarized_blast,
            by = 'zotu') %>%
  pivot_longer(cols = -c(zotu, taxid_rank),
               names_to = c('metric', 'measure'),
               names_pattern = '(.*)_(.*)') %>%
  mutate(across(c(measure, taxid_rank), str_to_sentence),
         taxid_rank = factor(taxid_rank, 
                             levels = c('Domain', 'Kingdom', 'Phylum',
                                        'Class', 'Order', 'Family',
                                        'Genus', 'Species'))) %>%
  ggplot(aes(x = taxid_rank, y = value)) +
  stat_summary(fun.data = mean_sdl) +
  facet_grid(measure ~ metric) +
  labs(x = NULL, 
       y = 'Value') +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        panel.background = element_rect(colour = 'black'),
        legend.position = 'right',
        legend.text = element_markdown(),
        legend.key = element_blank(),
        strip.background = element_blank())


summarized_blast_plot <- select(zotus_final, 
                                zotu, taxid_rank) %>%
  left_join(select(blast_hits, 
                   zotu, pident, qcov),
            by = 'zotu') %>%
  
  pivot_longer(cols = c(pident, qcov),
               names_to = 'metric', 
               values_to = 'value') %>%
  summarise(mean_value = mean(value),
            sd_value = sd(value),
            .by = c(taxid_rank, metric)) %>%
  
  mutate(taxid_rank = factor(str_to_title(taxid_rank), 
                             levels = c('Domain', 'Kingdom', 'Phylum',
                                        'Class', 'Order', 'Family',
                                        'Genus', 'Species'))) %>%
  ggplot(aes(y = mean_value, x = taxid_rank)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_value - sd_value,
                    ymax = mean_value + sd_value),
                width = 0.5) +
  facet_grid( ~ metric,
              scales = 'free_y',
              space = 'free_y',
              switch = "both") +
  labs(x = NULL, 
       y = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        panel.background = element_rect(colour = 'black'),
        legend.position = 'right',
        legend.text = element_markdown(),
        legend.key = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16),
        strip.placement = "outside",
        strip.clip = "off")


ggsave('summarized_blast_taxa.png', 
       plot = summarized_blast_plot,
       height = 7,
       width = 10)


#### Summarize Sample Results ####
sample_composition <- zotus_final %>%
  select(-unique_hits) %>%
  mutate(taxid = as.character(taxid),
         across(where(is.character), 
                ~if_else(. == 'LCA_dropped', NA_character_, .))) %>%
  group_by(across(domain:species)) %>%
  summarise(across(where(is.numeric), sum),
            across(c(zotu, taxid), ~unique(.) %>% str_c(collapse = '; ')),
            .groups = 'drop') %>% 
  mutate(lowest_level = case_when(!is.na(species) ~ str_c('s_', species),
                                  !is.na(genus) ~ str_c('g_', genus),
                                  !is.na(family) ~ str_c('f_', family),
                                  !is.na(order) ~ str_c('o_', order),
                                  !is.na(class) ~ str_c('c_', class),
                                  !is.na(phylum) ~ str_c('p_', phylum),
                                  !is.na(kingdom) ~ str_c('k_', kingdom),
                                  !is.na(domain) ~ str_c('d_', domain),
                                  TRUE ~ 'Unknown'),
         .after = species) %>%
  pivot_longer(cols = any_of(samples$sample_id),
               names_to = 'sample_id',
               values_to = 'n_reads') %>%
  filter(n_reads > 0) 



sample_composition_plot <- sample_composition %>%
  mutate(sample_id = fct_reorder(sample_id, n_reads)) %>%
  ggnested_jds(aes(y = sample_id, x = n_reads, 
                   main_group = class, sub_group = lowest_level),
               legend_labeling = 'sub', legend_title = 'Lowest Taxonomic\nClassification',
               main_keys = TRUE, nested_aes = c("fill"), 
               gradient_type = 'both',
               NA_option = 'black') +
  geom_col(position = 'fill') +
  scale_x_continuous(labels = scales::percent_format()) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(y = NULL, 
       x = 'Relative Number of Reads (%)') +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        panel.background = element_rect(colour = 'black'),
        legend.position = 'right',
        legend.text = element_markdown(),
        legend.key = element_blank())

ggsave('sample_composition.png', 
       plot = sample_composition_plot,
       height = 15,
       width = 7)

sample_composition_plot2 <- sample_composition %>%
  mutate(sample_id = fct_reorder(sample_id, n_reads)) %>%
  ggnested_jds(aes(y = sample_id, x = n_reads, 
                   main_group = class, sub_group = lowest_level),
               legend_labeling = 'sub', legend_title = 'Lowest Taxonomic\nClassification',
               main_keys = TRUE, nested_aes = c("fill"), 
               gradient_type = 'both',
               NA_option = 'black') +
  geom_col(position = 'stack') +
  scale_x_continuous(labels = scales::comma_format()) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(y = NULL, 
       x = 'Number of Reads') +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        panel.background = element_rect(colour = 'black'),
        legend.position = 'right',
        legend.text = element_markdown(),
        legend.key = element_blank())

ggsave('sample_composition2.png', 
       plot = sample_composition_plot2,
       height = 15,
       width = 7)


blast_classification_plot <- sample_composition %>%
  select(zotu, lowest_level) %>%
  distinct %>%
  rowwise(lowest_level) %>%
  reframe(zotu = str_split(zotu, ';') %>% unlist %>% str_trim) %>%
  
  left_join(select(blast_hits, 
                   zotu, pident, qcov),
            by = 'zotu') %>%
  left_join(select(zotus_final, zotu, taxid_rank),
            by = 'zotu') %>%
  pivot_longer(cols = c(pident, qcov),
               names_to = 'metric', 
               values_to = 'value') %>%
  summarise(mean_value = mean(value),
            sd_value = sd(value),
            .by = c(taxid_rank, lowest_level, metric)) %>%
  
  mutate(lowest_level = fct_reorder(lowest_level, mean_value),
         taxid_rank = factor(str_to_title(taxid_rank), 
                             levels = c('Domain', 'Kingdom', 'Phylum',
                                        'Class', 'Order', 'Family',
                                        'Genus', 'Species'))) %>%
  ggplot(aes(x = mean_value, y = lowest_level)) +
  geom_col() +
  geom_errorbar(aes(xmin = mean_value - sd_value,
                    xmax = mean_value + sd_value),
                width = 0.5) +
  facet_grid(taxid_rank ~ metric,
             scales = 'free_y',
             space = 'free_y',
             switch = "both") +
  labs(x = NULL, 
       y = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        panel.background = element_rect(colour = 'black'),
        legend.position = 'right',
        legend.text = element_markdown(),
        legend.key = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16),
        strip.placement = "outside",
        strip.clip = "off")

ggsave('blast_classification_plot.png', 
       plot = blast_classification_plot,
       height = 15,
       width = 7)

#### Output zOTUs FASTA ####
zotu_sequence_names <- select(sample_composition, zotu, domain:species, lowest_level) %>%
  distinct %>%
  rowwise(domain:species, lowest_level) %>%
  reframe(zotu = str_split(zotu, '; ') %>%
            unlist) %>%
  distinct %>%
  arrange(parse_number(zotu)) %>%
  mutate(across(domain:species,
                ~str_c(str_sub(cur_column(), 1, 1), ., sep = '_'))) %>%
  select(domain:species, zotu) %>%
  pivot_longer(cols = c(domain:species),
               names_to = 'taxon_level',
               values_to = 'taxa_name',
               values_drop_na = TRUE) %>%
  summarise(taxonomy = str_c(taxa_name, collapse = '; '),
            .by = zotu) %>%
  mutate(seq_name = str_c(zotu, taxonomy, sep = ': ')) %>%
  select(-taxonomy) %>%
  arrange(match(zotu, names(final_zotu_sequences)))

names(final_zotu_sequences) <- zotu_sequence_names$seq_name
Biostrings::writeXStringSet(final_zotu_sequences, 'zotu_sequences.fasta')

#### Create Taxonomy Table (Summarize over zOTUs) ####
zotus_updated <- left_join(zotus_final,
                           summarized_blast,
                           by = c('zotu')) %>%
  relocate(all_of(colnames(summarized_blast)[-1]),
           .after = unique_hits) %>%
  rowwise() %>%
  mutate(across(domain:species, ~ if_else(. == "LCA_dropped" & any(c_across(cur_column():species) != "LCA_dropped"), NA_character_, .))) %>%
  ungroup()
write_csv(zotus_updated, 'zotu_table.csv')

taxa_table <- summarise(zotus_final,
                        across(c(where(is.numeric), -unique_hits), 
                               sum),
                        n_zotu = n_distinct(zotu),
                        .by = c(domain:species)) %>%
  relocate(n_zotu, .after = taxid)
write_csv(taxa_table, 'lowest_taxonomy_table.csv')

#### Output Files with metadata for zOTU and lowest taxonomy tables ####
sample_columns <- colnames(zotus_final)[colnames(zotus_final) %in% samples$sample_id]

c('zotu: A unique identifier for each Zero-radius Operational Taxonomic Unit (ZOTU), representing a unique sequence variant.',
  'taxid: A numeric taxonomic identifier (from NCBI) that corresponds to the organism’s entry in a taxonomy database.',
  'unique_hits: Number of unique BLAST sequence hits.',
  'pident (mean, median, min, max) & qcov (mean, median, min, max): mean/median/min/max values of the BLAST hits passing filtering used to taxonomically classify the zOTU',
  'domain, kingdom, phylum, class, order, family, genus, species: Taxonomic classification of the zOTU.',
  'taxid_rank: The lowest level of taxonomic classification.',
  str_c(str_c(sample_columns[1:3], collapse = ', '), 
        ', ..., ', 
        str_c(sample_columns[(length(sample_columns) - 2):length(sample_columns)], collapse = ', '), 
        ': Number of reads of each zOTU in the specified sample')) %>%
  write_lines('zotu_table_metadata.txt')

c('domain, kingdom, phylum, class, order, family, genus, species: Taxonomic classification of the zOTU.',
  'taxid_rank: The lowest level of taxonomic classification.',
  'n_zotu: Number of zOTUs classified to the specified lowest level of taxonomic classification',
  str_c(str_c(sample_columns[1:3], collapse = ', '), 
        ', ..., ', 
        str_c(sample_columns[(length(sample_columns) - 2):length(sample_columns)], collapse = ', '), 
        ': Number of reads of each lowest taxonomic level in the specified sample')) %>%
  write_lines('lowest_taxonomy_table_metadata.txt')

#### Summarize BLAST -> Taxonomy ####
# tmp <- lca_taxonomy %>%
#   inner_join(filter(zotus_final, Sausage > 0) %>%
#                select(zotu),
#              by = 'zotu') %>%
#   left_join(nest(blast_hits,
#                  blast_hits = -c(zotu)),
#             by = 'zotu') %>%
#   rowwise %>%
#   mutate(n_chicken = sum(str_detect(blast_hits$blast_species, 'Gallus gallus')),
#          blast_hits_filtered = filter(blast_hits, pident >= 70, qcov >= 70) %>%
#            list,
#          n_chicken_post = sum(str_detect(blast_hits_filtered$blast_species, 'Gallus gallus')))
# 
# tmp %>%
#   # filter(n_chicken > 0) %>%
#   ungroup %>%
#   # slice(4) %>%
#   select(zotu, blast_hits_filtered) %>%
#   unnest(blast_hits_filtered) %>%
#   mutate(diff = max(pident) - pident) %>%
#   # filter(blast_species == 'Sepiella maindroni')
#   filter(diff <= 1) %>%
#   count(blast_species)
# 
# #### Recheck Metazoan ####
# zotus_final %>%
#   filter(phylum == 'LCA_dropped') %>%
#   select(zotu) %>%
#   distinct %>%
#   left_join(nest(blast_hits,
#                  blast_hits = -c(zotu)),
#             by = 'zotu') %>%
#   unnest(blast_hits) %>%
#   group_by(zotu, blast_species) %>%
#   summarise(n = n()) %>%
#   pivot_wider(names_from = blast_species, values_from = n)
# 
# 
# zotus_final %>%
#   filter(phylum == 'LCA_dropped') %>%
#   select(zotu) %>%
#   distinct %>%
#   left_join(nest(blast_hits,
#                  blast_hits = -c(zotu)),
#             by = 'zotu') %>%
#   unnest(blast_hits) %>%
#   # group_by(zotu) %>%
#   # filter(qcov == max(qcov)) %>%
#   # filter(blast_species == 'Sepiella maindroni') %>%
#   count(blast_species)
# 
# 
# 
# tmp <- zotus_final %>%
#   filter(phylum == 'LCA_dropped') %>%
#   select(zotu) %>%
#   distinct %>%
#   left_join(nest(blast_hits,
#                  blast_hits = -c(zotu)),
#             by = 'zotu') %>%
#     rowwise %>%
#     mutate(blast_hits_filtered = filter(blast_hits, pident >= 90, qcov >= 90) %>%
#              mutate(diff = max(pident) - pident) %>%
#              arrange(diff) %>%
#              # filter(qcov == max(qcov)) %>%
#              list) %>%
#     select(zotu, blast_hits_filtered) 
# 
# tmp %>%
#   unnest(blast_hits_filtered) %>%
#   # filter(diff < 0.75) %>%
#   # filter(blast_species %in% c('Sepiella maindroni',
#   #                             'Cheylostigmaeus occultatus',
#   #                             'Typhlodromus (Anthoseius) sp. X-BIOTU Ph-89',
#   #                             'Bubalus bubalis',
#   #                             'Equus caballus'))
#   count(zotu, blast_species) %>%
#   count(blast_species)
