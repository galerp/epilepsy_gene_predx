library(tidyverse, quietly = T)


start <- Sys.time()
message(" \n Begin propagation... \n ")

hpo_ancs <- read_csv(input.yaml$hpo_ancestor)
base_hpo <- read_csv(input.yaml$base_hpo)


prop_hpo<- base_hpo %>% 
    left_join(hpo_ancs) %>% 
    dplyr::select(ID, start_year, finish_year, ancs) %>% 
    dplyr::rename(HPO = ancs) %>% 
    separate_rows(HPO, sep = ";") %>% 
    #Remove duplicated HPO terms in each patient
    distinct()


write_csv(tot_base, paste0(input.yaml$file_path,"example_bin_prop.csv"))
