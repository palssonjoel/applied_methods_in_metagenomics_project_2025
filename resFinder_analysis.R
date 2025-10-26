library(tidyverse)
library(tools)

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

#Enter taxonomic level to plot
tax_level = "Phylum"

#This function takes a nested tibble and applies data cleaning procedures to each
#tibble in the nested tibble. Assumes that the data in the nested_tibble is found in 
#a column called "data"
clean = function(nested_data) {
  nested_data <- nested_data %>% 
    mutate(data = map(data, ~ .x %>%                              #Take each tibble in data column and apply functions  
                        rename("Template" = `#Template`) %>%                        #Change `#Template` column name
                        separate_wider_delim("Template",                            #Separate Template into Taxonomy and Template
                                             names = c("Template", "Taxonomy"),
                                             delim = " ",
                                             too_many = "merge") %>% 
                        separate_wider_delim("Taxonomy",                            #Separate Taxonomy into each taxonomic group
                                             names =taxonomy,
                                             delim = ";",
                                             too_few = "align_start",
                                             too_many = "merge") %>% 
                        mutate(Template_Identity = as.numeric(Template_Identity),   #Make variables numeric for plotting
                               Template_Coverage = as.numeric(Template_Coverage),
                               Query_Identity = as.numeric(Query_Identity),
                               Query_Coverage = as.numeric(Query_Coverage),
                               Depth = as.numeric(Depth),
                               q_value = as.numeric(q_value),
                               Expected = as.numeric(Expected)) 
    )
    )
}

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

#Load resFinderdata
resFinder_data_nested <- load_mapstat(path_resFinder)

#Load Silva data
silva_data_nested <- load_mapstat(path_silva)

#Unnest data
resFinder_data_unnested <- resFinder_data_nested %>% 
  unnest()

silva_data_unnested <- silva_data_nested %>% 
  unnest()


