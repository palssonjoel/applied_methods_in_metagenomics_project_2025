############################################################
# Title: Metagenomic MapStat Processing
# Author: Joel Pålsson
# Date: 26-10-2025
# Description: Loads .mapstat files, cleans nested data,
#              extracts taxonomy, and computes basepair totals.
############################################################

# 1. Setup
#===========================================================
library(tidyverse)
library(tools)
library(compositions)
library(zCompositions)

taxonomy = c(
  "Domain",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species")

#Enter path to files
path_resFinder = "data/kma_resfinder/"
path_silva = "data/kma_silva/"
path_resFinder_genes = "data/resfinder_genes.csv"

#Enter taxonomic level to plot
tax_level = "Phylum"

# 2. Function definitions
#===========================================================

# ---- 3.1 Load .mapstat files ----
clean = function(nested_data, taxonomy) {
  nested_data <- nested_data %>% 
    mutate(data = map(data, ~ .x %>%                              #Take each tibble in data column and apply functions  
                        separate_wider_delim("refSequence",                            #Separate Template into Taxonomy and Template
                                             names = c("refSequence", "Taxonomy"),
                                             delim = " ",
                                             too_many = "merge") %>% 
                        separate_wider_delim("Taxonomy",                            #Separate Taxonomy into each taxonomic group
                                             names = taxonomy,
                                             delim = ";",
                                             too_few = "align_start",
                                             too_many = "merge")))
  
  #Calculate total basepairs for each sample
  bp_total <- nested_data %>%
    unnest() %>% 
    group_by(sample) %>% 
    filter(Domain == "Bacteria") %>% 
    summarize(bpTotal_sum = sum(bpTotal))
  
  #Merge total basepairs with corresponding sample
  nested_data <- nested_data %>% 
    merge(bp_total)
}

# ---- 3.2 Clean nested Silva data ----
load_mapstat = function(path) {
  files <- dir(path = path,               #Get all files in directory
               recursive = TRUE,
               full.names = TRUE,
               pattern = "\\.mapstat$")
  samples <- basename(file_path_sans_ext(files))  #Extract sample names from file names
  
  #Load file names
  files <- dir(path = path,               #Get all files in directory
               recursive = TRUE,
               full.names = TRUE,
               pattern = "\\.mapstat$")
  samples <- basename(file_path_sans_ext(files))  #Extract sample names from file names
  
  #Create nested table with the data
  data_nested <- tibble(sample = files,                 
                        data = map(files,
                                   read_delim,
                                   delim = "\t", 
                                   skip = 6)) %>% 
    mutate(sample = samples) %>% #Rename file column contents to sample names
    mutate(data = map(data,
                      ~ rename(.x,
                               refSequence = `# refSequence`)))
}

# ---- 3.3 Apply closure ----
closure_table = function(unnested_data) {
  
  #Adjust bpTotals with gene lengths and apply closure
  unnested_data <- unnested_data %>% 
    group_by(sample) %>% 
    mutate(closure = clo(bpTotal / geneLength, total = 100)) %>% 
    group_by(sample, geneClass) %>%
    summarise(classClosure = sum(closure, na.rm = TRUE)) %>%
    ungroup() %>%
    rename(`AR_Class` = geneClass,
           Sample = sample)
  #use cmultRepl() if zeroes are in the data.
}

# 4. Data Loading
# ==========================================================

resFinder_data_nested <- load_mapstat(path_resFinder)

silva_data_nested <- load_mapstat(path_silva) %>% 
  clean(taxonomy)

resFinder_genes <- read_csv(path_resFinder_genes) %>% 
  rename(geneLength = Length,
         geneClass = Class)

# 5. Data Wrangling
# ==========================================================

# ---- 5.1 Merge resFinder data ----

#Add bpTotal to resFinder tibble
resFinder_data_nested <- resFinder_data_nested %>% 
  mutate(bpTotal_sum = silva_data_nested$bpTotal_sum)

# ---- 5.2 Unnest Data ----
resFinder_data_unnested <- resFinder_data_nested %>% 
  unnest()

# ---- 5.3 Merge resFinder data with gene data, and get closures by AMR class ----

#Merge resFinder data with resFinder gene data
resFinder_data_unnested <- merge(resFinder_data_unnested,
                                 resFinder_genes)

#Get closures. Wide data is necessary for alr() later
closures = closure_table(resFinder_data_unnested) %>% 
  pivot_wider(names_from = AR_Class,
              values_from = classClosure)

#ALR transformation
ALR <- closures %>% 
  column_to_rownames("Sample") %>%
  acomp() %>%         # convert to compositional data object
  alr() %>%           # apply Additive Log-Ratio transformation
  as_tibble(rownames = "Sample")
  
CLR <- closures %>% 
  column_to_rownames("Sample") %>%
  acomp() %>%         # convert to compositional data object
  clr() %>%           # apply Additive Log-Ratio transformation
  as_tibble(rownames = "Sample")


#Concept plot
closures %>%
  pivot_longer(
    cols = -Sample,               # pivot all class columns except Sample
    names_to = "AR_Class",        # new column for AMR class names
    values_to = "classClosure"    # new column for closure values
  ) %>% 
  ggplot(aes(x = Sample, y = classClosure, fill = AR_Class)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    x = "Sample",
    y = "Class Closure Value",
    fill = "AMR Class",
    title = "AMR Class Distribution per Sample"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


