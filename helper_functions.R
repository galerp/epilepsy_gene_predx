library(tidyverse, quietly = T)

pat_prop <- function(pat_base){
  pat_table_prop <- pat_base %>% 
    left_join(hpo_ancs %>% select(-definition)) %>% 
    select(famID, Ancestors) %>% 
    rename(HPO = Ancestors) %>% 
    separate_rows(HPO, sep = ";") %>% 
    #Remove duplicated HPO terms in each patient
    distinct()
}
