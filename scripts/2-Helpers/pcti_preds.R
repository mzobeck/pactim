# This script creates a helper function that makes predictions for the PACTIM model.
# The function takes in a data frame, the level of predictions (total, group, or subgroup),
# the disease to predict (optional), and the number of draws to make.
# It also accepts (optional) a value for the percent of the country living in a malaria endemic region (mpr0)
# for a Burkitt lymphoma (BL) and Kaposi sarcoma (KS) adjustment.
# Note: "PCTI" was the initial name for PACTIM. It lives on in function form. 

# Uncomment and run this section to install required packages if not already installed
# packages <- c("tidyverse", "readxl", "dtplyr", "brms", "tidybayes")
# new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
# if(length(new_packages) > 0) install.packages(new_packages)

pcti_preds <- function(
    data = NULL,
    level = NULL,
    disease = NULL,
    numb.draws = 1000
) {
  suppressWarnings({
    
  library(tidyverse)
  library(readxl)
  library(dtplyr)
  library(tidybayes)
  library(brms)
  
  # This code loads the augments.R helper script. This path needs to be changed if its in a different folder. 
  source("scripts/2-Helpers/augments.R")
  

  # Input validation --------------------------------------------------------

  if (!(level %in% c("total", "group", "subgroup"))) {
    stop("'level' must have value 'total', 'group', or 'subgroup'.")
  }
  
  if (length(which(names(data) %in% c("geo_unit", "population", "logpop", "age", "sex"))) != 5) {
    stop("'data' must have unique variables named 'geo_unit', 'age', 'sex', 'population', and 'logpop'.")
  }
  
  if (is.null(data$mpr0)) {
    message("No 'mpr0' variable. No adjustments for BL/KS will be performed.")
  }
    
  if (is.null(disease) & level %in% c("group", "subgroup")) {
      message("No 'disease' variable. Predictions for all diseases will be provided.")
    }
  
    suppressMessages({
      data <- data %>% ungroup
      data.names <- names(data)
      
      # Total Preds -------------------------------------------------------------
      if (level == "total") {
        
        print("Making initial predictions.")
        
        preds.initial <- data %>% 
          #Predicting from one model for every row in "data" data frame  
          add_predicted_draws(readRDS("models/total_model.rds"), ndraws = numb.draws)
        
        if (is.null(data$mpr0)| all(is.na(data$mpr0))) {
          
          print("Making final predictions.")
          
          preds.final <- preds.initial %>% ungroup 
          
        } else {
          
          print("Adjusting for BL and KS.")
          
          #BL Adjustments
          mod.bl <- readRDS("models/top_lvl_cats.rds") %>% 
            filter(disease == "lymphomas and reticuloendothelial neoplasms") %>% 
            .$model %>% 
            .[[1]]
          
          post.bl <- as_draws_df(mod.bl) %>% 
            select(shape)
          
          adj.samp <- post.bl %>% 
            sample_n(numb.draws) %>% 
            mutate(.draw = 1:numb.draws)
          
          bl.adj <- data %>% 
            select(geo_unit, mpr0) %>% 
            filter(!is.na(mpr0)) %>% 
            distinct(geo_unit, mpr0) %>% 
            mutate(data = map(mpr0, ~bl_augs(x = .x))) %>% 
            unnest(data) %>% 
            select(-x)
          
          bl.augments<- data %>% 
            left_join(bl.adj) %>% 
            filter(!is.na(bl_inc)) %>% 
            expand_grid(tibble(.draw = 1:numb.draws)) %>% 
            left_join(adj.samp) %>% 
            rowwise() %>% 
            mutate(scale = bl_inc/shape,
                   lambda = rgamma(1, shape = shape, scale = scale),
                   poisson = rpois(1, lambda),
                   .prediction = poisson/1e6 * population) %>% 
            ungroup %>% 
            select(all_of(c(data.names, ".draw", ".prediction")))
 
          #KS Adjustments
          mod.ks <- readRDS("models/top_lvl_cats.rds") %>% 
            filter(disease == "soft tissue and other extraosseous sarcomas") %>% 
            .$model %>% 
            .[[1]]
          
          post.ks <- as_draws_df(mod.ks)%>% 
            select(shape)
          
          adj.samp <- post.ks %>% 
            sample_n(numb.draws) %>% 
            mutate(.draw = 1:numb.draws)
          
          ks.adj <- data %>% 
            select(geo_unit, mpr0) %>% 
            filter(!is.na(mpr0)) %>% 
            distinct(geo_unit, mpr0) %>% 
            mutate(data = map(mpr0, ~ks_augs(x = .x))) %>% 
            unnest(data) %>% 
            select(-x)
          
          ks.augments<- data %>% 
            left_join(ks.adj) %>% 
            filter(!is.na(ks_inc)) %>% 
            expand_grid(tibble(.draw = 1:numb.draws)) %>% 
            left_join(adj.samp) %>% 
            rowwise() %>% 
            mutate(scale = ks_inc/shape,
                   lambda = rgamma(1, shape = shape, scale = scale),
                   poisson = rpois(1, lambda),
                   .prediction = poisson/1e6 * population) %>% 
            ungroup %>% 
            select(all_of(c(data.names, ".draw", ".prediction")))
   
        
          
          #Make final predictions with bl augments 
          print("Making final predictions.")
          
          preds.final <- preds.initial %>% 
            bind_rows(bl.augments) %>% 
            bind_rows(ks.augments) %>% 
            group_by(across(all_of(c(data.names, ".draw")))) %>% 
            summarise(.prediction = sum(.prediction)) %>% 
            ungroup
        }
      }
      
      # Group Preds -------------------------------------------------------------
      if (level == "group") {
        
        print("Making initial predictions.")
        
        # If no disease provided, produce predictions for all diseases
        if (is.null(disease)) {
          
          diseases <- readRDS("models/top_lvl_cats.rds")$disease
          
          preds.initial <- data %>% 
            expand_grid(disease = diseases) %>% 
            group_by(disease) %>% 
            nest() %>% 
            left_join(readRDS("models/top_lvl_cats.rds") %>% 
                        select(disease, model)) %>% 
            group_by(disease) %>% 
            mutate(preds = map2(model, data, function(model, new_data)
              add_predicted_draws(new_data, model, ndraws = numb.draws))) %>%  
            select(disease, preds) %>% 
            unnest(preds) 
        }
        
        # If one ore more disease groups are provided, predict only those groups 
        if (!is.null(disease)) {
          
          preds.initial <- data %>% 
            expand_grid(disease = disease) %>% 
            group_by(disease) %>% 
            nest() %>% 
            left_join(readRDS("models/top_lvl_cats.rds") %>% 
                        select(disease, model)) %>% 
            group_by(disease) %>% 
            mutate(preds = map2(model, data, function(model, new_data)
              add_predicted_draws(new_data, model, ndraws = numb.draws))) %>%  
            select(disease, preds) %>% 
            unnest(preds) 
        }
        
        # If no mpr0 data or all NA, make final predictions
        if (is.null(data$mpr0) | all(is.na(data$mpr0)) | 
            !is.null(disease) & !any(disease %in% c("lymphomas and reticuloendothelial neoplasms",
                                                 "soft tissue and other extraosseous sarcomas"))) {
          
          print("Making final predictions.")
          
          
          preds.final <- preds.initial %>% ungroup 
          
        } else {
        
          # BL Adjustment
          # If mpr0 provided and BL group is provided in disease, make BL predictions
          if(is.null(disease) | "lymphomas and reticuloendothelial neoplasms" %in% disease) {
            
            print("Adjusting for BL.")
          mod.bl <- readRDS("models/top_lvl_cats.rds") %>% 
            filter(disease == "lymphomas and reticuloendothelial neoplasms") %>% 
            .$model %>% 
            .[[1]]
          
          post.bl <- as_draws_df(mod.bl) %>% 
            select(shape)
          
          adj.samp <- post.bl %>% 
            sample_n(numb.draws) %>% 
            mutate(.draw = 1:numb.draws)
          
          
          bl.adj <- data %>% 
            select(geo_unit, mpr0) %>% 
            filter(!is.na(mpr0)) %>% 
            distinct(geo_unit, mpr0) %>% 
            mutate(data = map(mpr0, ~bl_augs(x = .x))) %>% 
            unnest(data) %>% 
            select(-x)
          
          bl.augments<- data %>% 
            left_join(bl.adj) %>% 
            filter(!is.na(bl_inc)) %>% 
            expand_grid(tibble(.draw = 1:numb.draws)) %>% 
            left_join(adj.samp) %>% 
            rowwise() %>% 
            mutate(scale = bl_inc/shape,
                   lambda = rgamma(1, shape = shape, scale = scale),
                   poisson = rpois(1, lambda),
                   .prediction = poisson/1e6 * population) %>% 
            ungroup %>% 
            select(all_of(c(data.names, ".draw", ".prediction"))) %>% 
            mutate(disease = "lymphomas and reticuloendothelial neoplasms")
          }
          
         
          # KS adjustments
          # If mpr0 provided and KS group is provided in disease, make BL predictions
          if(is.null(disease) | "soft tissue and other extraosseous sarcomas" %in% disease) {
            
          print("Adjusting for KS.")
              
          mod.ks <- readRDS("models/top_lvl_cats.rds") %>% 
            filter(disease == "soft tissue and other extraosseous sarcomas") %>% 
            .$model %>% 
            .[[1]]
          
          post.ks <- as_draws_df(mod.ks)%>% 
            select(shape)
          
          adj.samp <- post.ks %>% 
            sample_n(numb.draws) %>% 
            mutate(.draw = 1:numb.draws)
          
          ks.adj <- data %>% 
            select(geo_unit, mpr0) %>% 
            filter(!is.na(mpr0)) %>% 
            distinct(geo_unit, mpr0) %>% 
            mutate(data = map(mpr0, ~ks_augs(x = .x))) %>% 
            unnest(data) %>% 
            select(-x)
          
          ks.augments<- data %>% 
            left_join(ks.adj) %>% 
            filter(!is.na(ks_inc)) %>% 
            expand_grid(tibble(.draw = 1:numb.draws)) %>% 
            left_join(adj.samp) %>% 
            rowwise() %>% 
            mutate(scale = ks_inc/shape,
                   lambda = rgamma(1, shape = shape, scale = scale),
                   poisson = rpois(1, lambda),
                   .prediction = poisson/1e6 * population) %>% 
            ungroup %>% 
            select(all_of(c(data.names, ".draw", ".prediction"))) %>% 
            mutate(disease = "soft tissue and other extraosseous sarcomas")
          
          }
          
          #Make final predictions with augments 
          print("Making final predictions.")
          
          preds.final <- preds.initial
          
          if(is.null(disease) | "lymphomas and reticuloendothelial neoplasms" %in% disease) {
            
            preds.final <- preds.final %>% 
            bind_rows(bl.augments)
        }
            
          if(is.null(disease) | "soft tissue and other extraosseous sarcomas" %in% disease) {
            
            preds.final <- preds.final %>% 
            bind_rows(ks.augments)
          }
            
          preds.final <- preds.final %>% 
            group_by(across(all_of(c(data.names, "disease", ".draw")))) %>% 
            summarise(.prediction = sum(.prediction)) %>% 
            ungroup
          
        }
      }
      
      
      # Subgroup Preds ----------------------------------------------------------
      if (level == "subgroup") {
        
        print("Making initial predictions.")
        # If no disease provided, produce predictions for all diseases
        if (is.null(disease)) {
          
        diseases <- readRDS("models/baseline_model.rds")$disease
        
        preds.initial <- data %>% 
          #Expand data to make disease specific predictions
          expand_grid(disease = diseases) %>% 
          group_by(disease) %>% 
          nest() %>% 
          left_join(readRDS("models/baseline_model.rds") %>% 
                      select(disease, model)) %>% 
          group_by(disease) %>% 
          mutate(preds = map2(model, data, function(model, new_data)
            add_predicted_draws(new_data, model, ndraws = numb.draws))) %>%  
          select(disease, preds) %>% 
          unnest(preds) 
        
        } 
        
        # If one ore more disease subgroups are provided, predict only those subgroups 
        if (!is.null(disease)) {
          
          preds.initial <- data %>% 
            #Expand data to make disease specific predictions
            expand_grid(disease = disease) %>% 
            group_by(disease) %>% 
            nest() %>% 
            left_join(readRDS("models/baseline_model.rds") %>% 
                        select(disease, model)) %>% 
            group_by(disease) %>% 
            mutate(preds = map2(model, data, function(model, new_data)
              add_predicted_draws(new_data, model, ndraws = numb.draws))) %>%  
            select(disease, preds) %>% 
            unnest(preds) 
          
        } 

        # If no mpr0 data or all NA, make final predictions
        if (is.null(data$mpr0) | all(is.na(data$mpr0)) | 
            !is.null(disease) & !any(disease %in% c("burkitt lymphoma",
                                                    "kaposi sarcoma"))) {
          
          print("Making final predictions.")
          
          preds.final <- preds.initial %>% ungroup 
          
        } else {
          
          if(is.null(disease) | "burkitt lymphoma" %in% disease) {
            
            print("Adjusting for BL.")
          
          # BL Adjustment
          # If mpr0 provided and BL group is provided in disease, make BL predictions
          mod.bl <- readRDS("models/top_lvl_cats.rds") %>% 
            filter(disease == "lymphomas and reticuloendothelial neoplasms") %>% 
            .$model %>% 
            .[[1]]
          
          post.bl <- as_draws_df(mod.bl) %>% 
            select(shape)
          
          
          adj.samp <- post.bl %>% 
            sample_n(numb.draws) %>% 
            mutate(.draw = 1:numb.draws)
          
          
          bl.adj <- data %>% 
            select(geo_unit, mpr0) %>% 
            filter(!is.na(mpr0)) %>% 
            distinct(geo_unit, mpr0) %>% 
            mutate(data = map(mpr0, ~bl_augs(x = .x))) %>% 
            unnest(data) %>% 
            select(-x)
          
          bl.augments<- data %>% 
            left_join(bl.adj) %>% 
            filter(!is.na(bl_inc)) %>% 
            expand_grid(tibble(.draw = 1:numb.draws)) %>% 
            left_join(adj.samp) %>% 
            rowwise() %>% 
            mutate(scale = bl_inc/shape,
                   lambda = rgamma(1, shape = shape, scale = scale),
                   poisson = rpois(1, lambda),
                   .prediction = poisson/1e6 * population) %>% 
            ungroup %>% 
            select(all_of(c(data.names, ".draw", ".prediction"))) %>% 
            mutate(disease = "burkitt lymphoma")
          }
          
          # KS adjustments
          # If mpr0 provided and KS group is provided in disease, make BL predictions
          if(is.null(disease) | "kaposi sarcoma" %in% disease) {
            
            print("Adjusting for KS.")
          
          mod.ks <- readRDS("models/top_lvl_cats.rds") %>% 
            filter(disease == "soft tissue and other extraosseous sarcomas") %>% 
            .$model %>% 
            .[[1]]
          
          post.ks <- as_draws_df(mod.ks)%>% 
            select(shape)
          
          adj.samp <- post.ks %>% 
            sample_n(numb.draws) %>% 
            mutate(.draw = 1:numb.draws)
          
          ks.adj <- data %>% 
            select(geo_unit, mpr0) %>% 
            filter(!is.na(mpr0)) %>% 
            distinct(geo_unit, mpr0) %>% 
            mutate(data = map(mpr0, ~ks_augs(x = .x))) %>% 
            unnest(data) %>% 
            select(-x)
          
          ks.augments<- data %>% 
            left_join(ks.adj) %>% 
            filter(!is.na(ks_inc)) %>% 
            expand_grid(tibble(.draw = 1:numb.draws)) %>% 
            left_join(adj.samp) %>% 
            rowwise() %>% 
            mutate(scale = ks_inc/shape,
                   lambda = rgamma(1, shape = shape, scale = scale),
                   poisson = rpois(1, lambda),
                   .prediction = poisson/1e6 * population) %>% 
            ungroup %>% 
            select(all_of(c(data.names, ".draw", ".prediction"))) %>% 
            mutate(disease = "kaposi sarcoma")
          }
          
          
          #Make final predictions with bl augments 
          print("Making final predictions.")
          
          preds.final <- preds.initial
          
          if(is.null(disease) | "burkitt lymphoma" %in% disease) {
            
            preds.final <- preds.final %>% 
              bind_rows(bl.augments)
          }
          
          if(is.null(disease) | "kaposi sarcoma" %in% disease) {
            
            preds.final <- preds.final %>% 
              bind_rows(ks.augments)
          }
          
          preds.final <- preds.final %>% 
            group_by(across(all_of(c(data.names, "disease", ".draw")))) %>% 
            lazy_dt() %>% 
            summarise(.prediction = sum(.prediction)) %>% 
            ungroup %>% 
            as_tibble
          
        }
      }
      
    })
  })
  
  return(preds.final)
  
}

