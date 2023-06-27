library(tidyverse, quietly = T)


setwd("/Users/galerp/Documents/manuscripts/cube3/git_repo/Cube3/Files/")

#Set output directory
output_dir = "/Users/galerp/Desktop/"

hpo_ancs <- read_csv("HPO_ancs_rl_2020-10-12_dl_2021-08-03.csv")
base_hpo <- read_csv("example_bin_base.csv")


prop_hpo<- base_hpo %>% 
    left_join(hpo_ancs) %>% 
    select(ID, start_year, finish_year, ancs) %>% 
    rename(HPO = ancs) %>% 
    separate_rows(HPO, sep = ";") %>% 
    #Remove duplicated HPO terms in each patient
    distinct()


write_csv(prop_hpo, paste0(output_dir,"example_bin_prop.csv"))
