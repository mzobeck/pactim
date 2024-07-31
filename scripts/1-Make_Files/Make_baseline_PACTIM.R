# This script will use the data obtained from SEEER to create the PACTIM models. 
# The variable names from SEER must mastch the names in the "Load and clean data" section or else the code
# Will need to be modified manually. 
#
# Three levels of modesl are produced. First, the ICCC-3 subgroup model, then the ICCC-3 group model, and finally the total cancer model.
# Cross-validation (not shown) results suggest that an age by sex interaction is needed for germ cell tumors. 
# A Poisson model is used for diseases with very low incidence to improve the stability of model estimates. 

# Load the packages  ---------------------------------------------------

# Uncomment and run this section to install required packages if not already installed
# packages <- c("tidyverse", "readxl", "janitor", "brms", "loo")
# new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
# if(length(new_packages) > 0) install.packages(new_packages)

# Important note: We use cmdstanr to help fit the the `brms` models. This makes the model fitting process faster.
# To use cmdstanr you need to intall the package with the instructions at https://mc-stan.org/cmdstanr/articles/cmdstanr.html
# If you do not want to intall cmdstanr, erase all of the rows "backend = "cmdstanr" below. 

library(tidyverse)
library(readxl)
library(janitor)
library(brms)
library(loo)

# Factor levels 
age.groups <-c("0_4 years", "5_9 years", "10_14 years", "15_19 years")

# Load and clean data -------------------------------------------------------
data.raw <-read_csv("data/raw/seer/data.csv") %>%
  clean_names() %>%
  mutate(across(where(is.character), tolower))

data <- data.raw %>%
  rename(
    race_eth = contains("race"),
    disease = contains("iccc"),
    year = contains("year_of_diagnosis"),
    age = contains("age_groups_fixed")
  ) %>%
  filter(!str_detect(race_eth, "unknown")) %>%
  mutate(age = factor(age, levels = age.groups)) %>%
  separate(
    disease,
    into = c("cat", "disease"),
    sep = "\\s",
    extra = "merge"
  )


# Make ICCC Subgroup model ------------------------------------------------
data.model <- data %>%
  filter(str_detect(cat, "a|b|c|d|e|f|v(?![a-zA-Z])|not")) %>% 
  filter(!str_detect(disease, "neuroblastoma and other peripheral nervous cell tumors")) %>% 
  mutate(cat = str_remove_all(cat, "a|b|c|d|e|f")) %>%
  filter(sex != "male and female") %>%
  mutate(logpop = log(population)) %>% 
  mutate(disease = ifelse(str_detect(disease, "classified by iccc"), 
                          "not classified by iccc or in situ",
                          disease))

# Pull out germ cell tumors because model comparison demonstrated an interaction term is needed
gonadal_gct <- data.model %>%
  filter(str_detect(disease, "malignant gonadal germ cell")) %>%
  mutate(logpop_c = logpop - mean(logpop))

nogonadal_gct <- data.model %>%
  filter(str_detect(disease, "malignant extracranial")) %>%
  mutate(logpop_c = logpop - mean(logpop))

gct_tib <- bind_rows(gonadal_gct, nogonadal_gct)

everything_else <- data.model %>%
  filter(!(disease %in% gct_tib$disease))

# Poisson model - no interactions  
baseline.poisson.file <- "models/intermediates/model_tib.rds"

if (!file.exists(baseline.poisson.file)) {
  model.tib <- data.model %>%
    group_by(disease) %>%
    nest() %>%
    mutate(model = map(
      data,
      ~ brm(
        formula = count ~ age + sex + offset(logpop),
        family = poisson,
        data = .x,
        save_pars = save_pars(all = TRUE),
        backend = "cmdstanr",
        chains = 4,
        cores = 4,
        iter = 4000
      )
    ))
  
  saveRDS(model.tib, baseline.poisson.file)
}

# Negative Binomial (AKA Gamma Poisson)  

nb.file <- "models/intermediates/model_tib_nb_bespoke.rds"

if (!file.exists(nb.file)) {
  gct.models <- gct_tib %>%
    group_by(disease) %>%
    nest() %>%
    mutate(model = purrr::map(
      data,
      ~ brm(
        formula = count ~  sex * age + offset(logpop),
        family = negbinomial(),
        data = .x,
        save_pars = save_pars(all = TRUE),
        backend = "cmdstanr",
        chains = 4,
        cores = 4,
        iter = 4000
      )
    ))
  
  eveything.models <- everything_else %>%
    group_by(disease) %>%
    nest() %>%
    mutate(model = map(
      data,
      ~ brm(
        formula = count ~  sex + age + offset(logpop),
        family = negbinomial(),
        data = .x,
        save_pars = save_pars(all = TRUE),
        backend = "cmdstanr",
        chains = 4,
        cores = 4,
        iter = 4000
      )
    ))
  
  # save
  model.tib.nb <- bind_rows(gct.models, eveything.models)
  
  saveRDS(model.tib.nb, nb.file)
  
}

# This code finds the diagnoses that have very unstable estimates due to the small incidence.
#> These are the diseases that require the Poisson model

model.tib.nb <- readRDS(nb.file)

excluders <- model.tib.nb %>%
  mutate(feffs = map(model, function(mod)
    abs(fixef(mod)["Intercept", ]))) %>%
  select(disease, feffs) %>%
  unnest_wider(feffs) %>%
  filter(Estimate > 20)


# Final Disease Model Tib 
#> Poisson: unspec malignant renal tumors, unspecified malignant hepatic tumors, kaposi sarcoma, nasopharyngeal carcinomas
#> Negative binomial - Age*Sex: gonadal gct, extragonadal gct
#> Negative binomial - no interaction: Everything not listed above

final_model_tib <- readRDS(nb.file) %>%
  #Get rid of all adjusted models
  filter(!(disease %in% excluders$disease)) %>%
  #bind poisson for excluders
  bind_rows(readRDS(baseline.poisson.file) %>% filter(disease %in% excluders$disease))

saveRDS(final_model_tib, "models/baseline_model.rds")



# Make ICCC Group Model ---------------------------------------------------
data.toplvls <- data %>% 
  filter(!str_detect(cat, "a|b|c|d|e|f")) %>% 
  filter(!str_detect(disease, "classified by iccc")) 

# Make model data - also need interaction term for gcts. 
data.model.top <- data.toplvls %>% 
  filter(sex != "male and female") %>% 
  mutate(logpop = log(population))

gct_tib <- data.model.top %>% 
  filter(str_detect(disease, "germ cell")) 

everything_else <- data.model.top %>% 
  filter(!(disease %in% gct_tib$disease))


# Fit the models
nb.top.file <- "models/top_lvl_cats.rds"

if (!exists(nb.top.file)) {
  top.gct.models <- gct_tib %>% 
    group_by(disease) %>% 
    nest() %>%  
    mutate(model = map(data, 
                       ~ brm(
                         formula = count ~  sex*age + offset(logpop),
                         family = negbinomial(), data = .x, save_pars = save_pars(all = TRUE),
                         backend = "cmdstanr", chains = 4, cores = 4, iter = 4000)
    ))
  
  top.eveything.models <- everything_else %>% 
    group_by(disease) %>% 
    nest() %>%  
    mutate(model = map(data, 
                       ~ brm(
                         formula = count ~  sex + age + offset(logpop),
                         family = negbinomial(), data = .x, save_pars = save_pars(all = TRUE),
                         backend = "cmdstanr", chains = 4, cores = 4, iter = 4000)
    ))
  
  model.top.nb <- bind_rows(top.gct.models, top.eveything.models)
  
  # Save the model
  saveRDS(model.top.nb, nb.top.file)
}

# Total cancers model -----------------------------------------------------

# Get the population data for each age, sex, year, race/ethnicity group 
pops <- data %>% 
  distinct(age, sex, year, race_eth, population) %>% 
  mutate(logpop = log(population)) %>% ungroup %>% 
  filter(sex != "male and female")  


# Make Model Data
data.model.tote <- data %>%
  filter(str_detect(cat, "a|b|c|d|e|f|v(?![a-zA-Z])|not")) %>% 
  filter(!str_detect(disease, "neuroblastoma and other peripheral nervous cell tumors")) %>% 
  mutate(cat = str_remove_all(cat, "a|b|c|d|e|f")) %>%
  filter(sex != "male and female") %>% 
  group_by(sex, age, year, race_eth) %>% 
  summarise(count = sum(count, na.rm = T))  %>% 
  left_join(pops)


# Fit the total model
total.file <- "models/total_model.rds"

if (!file.exists(total.file)) {
  
  total.model <- brm(
    formula = count ~  sex + age + offset(logpop),
    family = negbinomial(),
    data = data.model.tote,
    save_pars = save_pars(all = TRUE),
    backend = "cmdstanr",
    chains = 4,
    cores = 4,
    iter = 4000
    
  )
  
  saveRDS(total.model, total.file)
}



