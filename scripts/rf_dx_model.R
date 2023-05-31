library(tidyverse)
library(Hmisc)
library(hablar)
library(magrittr)
library(mltools)
library(data.table)
library(gbm)
library(caTools)
library(randomForest)
library(caret)
library(ggcorrplot)
library(factoextra)
library(FactoMineR)
library(pROC)
library(ROCR)


#HPO helper files
hpo_def <- read_csv(input.yaml$hpo_tree)
hpo_isa <- read_csv(input.yaml$hpo_isa)

#prop file
prop_hpo_full <- read_csv(input.yaml$prop_hpo)


gene_class <- read_csv(input.yaml$gene_class)


gene_dx <- read_csv(input.yaml$gene_dx)
gene_class <- read_csv(input.yaml$gene_class)


#Remove individuals with genetic dx but no age of dx
gene_nodx_age <- gene_dx %>%
  select(ID, Gene,age_genetic_dx) %>%
  filter(!is.na(Gene),!is.na(ID), Gene !="") %>%
  filter(is.na(age_genetic_dx)) %>% 
  unique() 

#Diagnosed individuals with age of dx
prop_hpo_dx <- prop_hpo_full %>% 
  filter(ID %nin% gene_nodx_age$ID) %>% 
  #Ensure that removing non diagnosed individuals
  filter(!is.na(Gene), Gene !="") %>% 
  #Subtract small number to prevent possible hidden floats 
  mutate(age_genetic_dx = age_genetic_dx-0.00001) %>% 
  #Remove terms after diagnosis
  filter(finish_year<age_genetic_dx)

#Non-diagnosed individuals
prop_hpo_nodx <- prop_hpo_dx %>% 
  filter(ID %nin% gene_nodx_age$ID) %>% 
  filter(ID %nin% prop_hpo_dx$ID) 
  
#Combined prop table
prop_hpo <- prop_hpo_dx %>% 
  rbind(prop_hpo_nodx)

#Compose full genetic diagnoses
mono_gene <- gene_dx %>% 
  filter(ID %in% prop_hpo$ID) %>% 
  filter(!is.na(age_genetic_dx))

all_gene <- mono_gene %>% 
  mutate(Gene = "All")

gene_classes <- mono_gene %>% 
  left_join(gene_class %>% 
              select(Gene, Class)) %>% 
  select(-Gene) %>%
  rename(Gene = Class) %>%
  filter(!is.na(Gene)) %>% 
  select(ID, age_genetic_dx, Gene)

#All diagnoses
full_dx <- mono_gene %>% 
  rbind(all_gene) %>% 
  rbind(gene_classes)



#Features
filt_accord <- read_csv(input.yaml$sig_feats) %>%
  #The filters below was applied prior to file upload
  # filter(pval<0.05) %>% 
  # filter(yes_npats > 1) %>% 
  mutate(t_end = round(t_end, digits=5)) %>% 
  mutate(t_start = round(t_start, digits=5))



rf_run <- function(filt_accord, prop_hpo, genotype, max_age, n_feats){
  
  geno_pat <- full_dx %>% 
    filter(Gene == genotype)
  
  gene_prop_y <- prop_hpo %>% 
    filter(ID %in% geno_pat$ID)
  
  control_prop_y <- prop_hpo %>% 
    filter(ID %nin% prop_hpo$ID)
  
  #Features with filter for genetic dx
  gene_accord <- filt_accord %>% 
    filter(gene == genotype)
  
  #Get top features for this gene and time period
  gene_top_pval <- gene_accord %>% 
    filter(t_end<=max_age) %>% 
    dplyr::group_by(HPO, def) %>% 
    dplyr::summarize(pval = min(pval)) %>%
    left_join(gene_accord %>% select(HPO, pval, t_start,t_end,yes_npats)) %>%
    arrange(pval)
  
  #Prune features to remove redundant terms
  gene_top_pval_prune <- gene_top_pval %>%
    left_join(hpo_isa) %>%
    left_join(gene_top_pval %>% rename(is.a = HPO) %>% select(-def, -pval) %>% mutate(prune = "Yes")) %>%
    distinct() %>%
    filter(!is.na(prune))
  
  gene_top_pval_f <- gene_top_pval %>%
    left_join(gene_top_pval_prune %>% 
                ungroup(HPO) %>%
                select(is.a, pval, t_start, t_end,prune) %>% rename(HPO = is.a))%>%
    distinct() %>%
    filter(is.na(prune)) %>%
    arrange(pval)
  
  gene_top_terms <- gene_top_pval_f[1:n_feats,]
  
  #Need next line as at low ages there are less than 30 significant terms
  gene_top_terms <- gene_top_terms %>% filter(!is.na(HPO))
  
  #This is a quick and dirty way to get features into the appropriate format 
  # for ml functions - an "encoder helper"
  encod_helper <- as.data.frame(matrix(nrow=10000,ncol = 4))
  names(encod_helper) <- c("HPO","def","t_start","t_end")
  
  cur_r <- 1
  for(i in 1:nrow(gene_top_terms)){
    
    i_start <- gene_top_terms$t_start[i]
    i_end <- gene_top_terms$t_end[i]
    
    hp <- gene_top_terms$HPO[i]
    hp_def <- gene_top_terms$def[i]
    
    #Rounding error can cause it to miss the last bin
    segs <- unique(round((c(seq(i_start, i_end, 1/12),i_end)),digits=5))%>% sort()
    segs_2 <- segs[1:(length(segs)-1)]
    
    encod_helper$t_start[cur_r:(cur_r+length(segs_2)-1)] <- segs_2
    encod_helper$t_end[cur_r:(cur_r+length(segs_2)-1)] <- segs[2:length(segs)]
    encod_helper$HPO[cur_r:(cur_r+length(segs_2)-1)] <- hp
    encod_helper$def[cur_r:(cur_r+length(segs_2)-1)] <- hp_def
    cur_r <- cur_r+length(segs_2)
  }
  
  
  encod_helper <- encod_helper %>% 
    filter(!is.na(HPO))

  
  ###########################
  # One-Hot Encode
  ###########################
  
  #Gene Features One Hot Encoding
  gene_long_y <- gene_prop_y %>% 
    rename(t_start = start_year, t_end = finish_year) %>% 
    mutate(t_start = as.factor(round(t_start*12))) %>% 
    mutate(t_end = as.factor(round(t_end*12))) %>% 
    mutate(HPO = as.factor(HPO)) %>% 
    select(-def)%>% 
    dplyr::left_join(encod_helper %>% mutate(top_feats = "Yes") %>% 
                       mutate(t_start = as.factor(round(t_start*12))) %>% 
                       mutate(t_end = as.factor(round(t_end*12))) %>% 
                       mutate(HPO = as.factor(HPO))) %>% 
    filter(top_feats == "Yes") %>% 
    left_join(gene_top20 %>% mutate(feat = paste(HPO,t_start,t_end,sep="_")) %>% 
                select(HPO, def, feat))
  
  
  gene_long_sub <- gene_long_y %>% 
    select(pat_id, feat) %>% 
    mutate(feat = as.factor(feat)) %>% 
    distinct()
  
  #Generate feature that is a sum of all features
  gene_feat_sum <- gene_long_sub %>%
    count(pat_id) %>%
    rename(total_feats = n)
  
  gene_long_sub_tb <- as.data.table(gene_long_sub)
  
  gene_wide <- dcast(gene_long_sub_tb[, list(V1=1, pat_id, feat)], 
                     pat_id ~ feat, fun=sum, value.var="V1", 
                     drop=c(TRUE, FALSE))#%>% 
  
  #Take into account individuals with no features
  pat_miss <- gene_prop_y$pat_id[gene_prop_y$pat_id %nin% gene_long_y$pat_id] %>% unique
  
  fill_miss <- as.data.frame(matrix(nrow=length(pat_miss), ncol = length(gene_wide)))
  names(fill_miss) <- names(gene_wide)
  fill_miss$pat_id <- pat_miss
  fill_miss[is.na(fill_miss)] <- 0
  gene_wide <- gene_wide %>% rbind(fill_miss) %>% 
    left_join(gene_feat_sum)
  
  
  ################
  # CONTROLS One-Hot Encode
  ################
  
  control_long <- control_prop_y %>% 
    rename(t_start = start_year, t_end = finish_year) %>% 
    mutate(t_start = as.factor(round(t_start*12))) %>% 
    mutate(t_end = as.factor(round(t_end*12))) %>% 
    mutate(HPO = as.factor(HPO)) %>% 
    select(-def) %>% 
    dplyr::left_join(encod_helper %>% mutate(top_feats = "Yes") %>% 
                       mutate(t_start = as.factor(round(t_start*12))) %>% 
                       mutate(t_end = as.factor(round(t_end*12))) %>% 
                       mutate(HPO = as.factor(HPO))) %>% 
    filter(top_feats == "Yes") %>% 
    left_join(gene_top20 %>% mutate(feat = paste(HPO,t_start,t_end,sep="_")) %>% 
                select(HPO, def, feat))
  

    control_long_sub <- control_long %>% 
    select(pat_id, feat) %>% 
    mutate(feat = as.factor(feat)) %>% 
    distinct()
    
  #Generate feature that is a sum of all features
  control_feat_sum <- control_long_sub %>%
    count(pat_id) %>%
    rename(total_feats = n)
  
  
  control_long_sub_tb <- as.data.table(control_long_sub)
  
  control_wide <- dcast(control_long_sub_tb[, list(V1=1, pat_id, feat)], 
                        pat_id ~ feat, fun=sum, value.var="V1", 
                        drop=c(TRUE, FALSE))%>% 
    left_join(control_feat_sum)
  
  #Take into account individuals with no features
  pat_miss <- control_prop_y$pat_id[control_prop_y$pat_id %nin% control_long$pat_id] %>% unique
  
  fill_miss <- as.data.frame(matrix(nrow=length(pat_miss), ncol = length(gene_wide)))
  names(fill_miss) <- names(gene_wide)
  fill_miss$pat_id <- pat_miss
  fill_miss[is.na(fill_miss)] <- 0
  control_wide <- control_wide %>% rbind(fill_miss)
  
  ### Combine ALL DATA ####
  tot_wide <- control_wide %>% mutate(gene = 0) %>% 
    rbind(gene_wide %>% mutate(gene = 1))
  tot_wide$gene <- as.factor(tot_wide$gene)
  tot_wide[is.na(tot_wide)] <- 0

  ####################
  # Random Pulls
  ####################
  
  control_pats <- control_wide$pat_id %>% unique
  gene_pats <- gene_wide$pat_id %>% unique
  n_gene <- length(unique(gene_pats))
  
  #Make sure this is True
  # n_gene == length(unique(gene_wide$pat_id))
  
  test_perc <- 0.3
  ntot <- round(n_gene/0.4) #data imbalance 60/40
  
  
  
  test_n <- round(ntot * test_perc)
  train_n <- round(ntot * (1-test_perc))
  
  test_gene_n <- floor(test_perc*n_gene)
  test_control_n <- test_n-test_gene_n
  
  train_gene_n <- n_gene - test_gene_n
  train_control_n <- train_n-train_gene_n
  
  
  #PARAMETERS
  n_its <- 1000 #number of iterations
  
  #Below some lines are commented out that extract feature importance
  #Uncomment to get this data. Can considerably slow function for large datasets
  
  
  #OUTPUT MATRIX
  # Note: it is a matrix for efficiency
  cnames_1 = c("rf_Accuracy","rf_Precision","rf_Recall","rf_AUC","rf_F1","rf_NPV")
  # If you want feature importance, uncomment
  # feat_cnames = paste0("feat_",seq(1,tot_feats+1,1))
  # cnames_all <- c(cnames_1, feat_cnames)
  
  cnames_all <- cnames_1
  xgb_its <- as.data.frame(matrix(nrow = n_its,ncol = length(cnames_all)))
  
  names(xgb_its) <- c(cnames_all)
  

  for(i in 1:n_its){
    
    control_samp = sample(control_pats, ntot - n_gene)
    
    control_train <- sample(control_samp, train_control_n)
    control_test <- control_samp[control_samp %nin% control_train]
    
    gene_train <- sample(gene_pats, train_gene_n)
    gene_test <- gene_pats[gene_pats %nin% gene_train]
    
    train_pats <- c(gene_train, control_train)
    test_pats <- c(gene_test, control_test)
    
    train_set <- tot_wide %>% filter(pat_id %in% train_pats) %>% distinct()
    
    test_set <- tot_wide %>% filter(pat_id %in% test_pats) %>% distinct()
    
    y_train <- (train_set$gene %>% as.numeric) - 1
    y_test <- (test_set$gene %>% as.numeric) - 1
    
    X_train <- train_set %>% 
      select(-gene, -pat_id)
    X_test <- test_set %>% 
      select(-gene,-pat_id)
    
    ###############
    # TEST Model RandomForest
    ###############
    
    rf <- randomForest(X_train, y = as.factor(as.character(y_train)))
    
    
    rf_pred <- predict(rf, X_test, type="prob") # for class probabilities
    rf_pr_test <- prediction(rf_pred[,2], as.character(y_test))
    r_auc_test1 <- performance(rf_pr_test, measure = "auc")@y.values[[1]] 
    rf_AUC <- r_auc_test1
    
    
    rf_pred <- as.data.frame(rf_pred)
    colnames(rf_pred) <- levels(tot_wide$gene)
    
    rf_pred$PredictedClass <- apply(rf_pred, 1, function(y) colnames(rf_pred)[which.max(y)])
    rf_pred$ActualClass <- as.character(y_test)
    
    rf_pred_acc <- rf_pred %>% 
      select("0","1","PredictedClass", "ActualClass") %>% 
      mutate(FP = case_when(((PredictedClass != ActualClass) & 
                               (PredictedClass==1))~1,
                            TRUE~0)) %>% 
      mutate(TP = case_when(((PredictedClass == ActualClass) &
                               (PredictedClass == 1))~1,
                            TRUE~0)) %>% 
      mutate(FN = case_when(((PredictedClass != ActualClass) & 
                               (ActualClass==1))~1,
                            TRUE~0)) %>% 
      mutate(TN = case_when(((PredictedClass == ActualClass) & 
                               (PredictedClass==0))~1,
                            TRUE~0))
    
    rf_Accuracy <- sum(rf_pred$PredictedClass == rf_pred$ActualClass) / nrow(rf_pred)
    rf_Precision <- sum(rf_pred_acc$TP)/
      (sum(rf_pred_acc$TP) + sum(rf_pred_acc$FP))
    rf_Recall <- sum(rf_pred_acc$TP)/
      (sum(rf_pred_acc$TP) + sum(rf_pred_acc$FN))
    rf_F1 <- 2 * (rf_Recall * rf_Precision) / 
      (rf_Recall + rf_Precision)
    rf_NPV <-  sum(rf_pred_acc$TN)/
      (sum(rf_pred_acc$FN) + sum(rf_pred_acc$TN))

    ############
    #For feature importance
    ############
    
    # #Conditional=True, adjusts for correlations between predictors.
    # i_scores <- varImp(rf, conditional=TRUE)
    # i_scores <- i_scores %>% tibble::rownames_to_column("var")
    # i_scores_comp <- paste(i_scores$var,i_scores$Overall,sep=";")
    #For features
    # xgb_its[i,1:length(cnames_all)] <- c(rf_Accuracy,rf_Precision,rf_Recall,rf_AUC, rf_F1,rf_NPV,i_scores_comp)
    
    xgb_its[i,1:length(cnames_all)] <- c(rf_Accuracy,rf_Precision,rf_Recall,rf_AUC, rf_F1,rf_NPV)

    
  }
  xgb_it_tot <- as.data.frame(xgb_its)
  names(xgb_it_tot) <- cnames_all
  xgb_its$Iteration <- seq(1,n_its)
  xgb_it_tot$ngene = length(gene_pats)
  xgb_it_tot$ncontrol = length(control_pats)
  xgb_it_tot$age_max = max_age
  xgb_it_tot$n_feats <- n_feats
  return(xgb_it_tot)
}


#

#Ages to test
start_age = 2/12 + 0.001
end_age = 5.084

tst_ages <- seq((start_age + 1/12), 5.084, 1/12)
tst_ages = tst_ages + 0.001 #Prevents potential issues with hidden floats

#Gene to test
genotype = "SCN1A"

#Features to test
feat_seq = c(5,10,15,20)

for (n_feats in feat_seq){
  print(n_feats)

  rf_it_tot <- rf_run(filt_accord, genotype, start_age,n_feats)
  
  for(max_age in tst_ages){
    print(max_age)
    rf_it_cur <- rf_run(filt_accord, genotype, max_age,n_feats)
    rf_it_tot <- bind_rows(rf_it_tot, rf_it_cur)
  }
  
  
  output_file = paste0(input.yaml$output_dir, genotype, "_rf_model_",
                       n_feats,"feat.csv")
  write_csv(rf_it_tot, output_file)
}


