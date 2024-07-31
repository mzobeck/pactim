# This script will augment the incidence of BL and KS based on the proportion of 
# the population in a country that lives in a Malaria endemic region.
# The methods and adjustment coefficients are drawn from the supplement to 
# Johnston et al. 2021 (https://www.sciencedirect.com/science/article/pii/S1877782119301729)

# Load the packages  ---------------------------------------------------

# Uncomment and run this section to install required packages if not already installed
# packages <- c("dplyr", "tidyr")
# new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
# if(length(new_packages) > 0) install.packages(new_packages)


library(dplyr)
library(tidyr)

# Factor levels
age.groups <- c("0_4 years", "5_9 years", "10_14 years", "15_19 years")
sex.groups <- c("male", "female")

# Make the augmented Burkitt lymphoma inc per million estimates 
bl_augs <- function(x = 0, #this is where the mpr0 value will be passed in the pcti_preds function - 0 is a placeholder
                    age = c("0_4 years", "5_9 years", "10_14 years", "15_19 years"),
                    sex =  c("male", "female")) {
  
  expand_grid(age = age.groups, sex = sex) %>% 
    mutate(x = x) %>% 
    mutate(bl_inc = case_when(
      age == "0_4 years" & sex == "male" ~ 2.7 +2.1*x + 10.5*x^2,
      age == "5_9 years" & sex == "male" ~ 2.8 +18.2*x + 28.3*x^2,
      age == "10_14 years" & sex == "male" ~ 1.7 +10.1*x + 10.3*x^2,
      age == "15_19 years" & sex == "male" ~ 1.7 +10.1*x + 10.3*x^2,
      
      age == "0_4 years" & sex == "female" ~ 1.1 +2.5*x + 10.1*x^2,
      age == "5_9 years" & sex == "female" ~ 1 + 10.6*x + 22.9*x^2,
      age == "10_14 years" & sex == "female" ~ 0.5 +2.4*x + 14.5*x^2,
      age == "15_19 years" & sex == "female" ~ 0.5 +2.4*x + 14.5*x^2
    ))
  
} 

# Make the augmented kaposi sarcoma (ks) inc per million estimates 
ks_augs <- function(x = 0, #this is where the mpr0 value will be passed in the pcti_preds function - 0 is a placeholder
                    age = c("0_4 years", "5_9 years", "10_14 years", "15_19 years"),
                    sex =  c("male", "female")) {
  
  expand_grid(age = age.groups, sex = sex) %>% 
    mutate(x = x) %>% 
    mutate(ks_inc = case_when(
      age == "0_4 years" & sex == "male" ~ 0.3 +46.6*x - 30.8*x^2,
      age == "5_9 years" & sex == "male" ~ 0.5 +93*x - 69.8*x^2,
      age == "10_14 years" & sex == "male" ~ 0.2 +76.3*x - 50.3*x^2,
      age == "15_19 years" & sex == "male" ~ 0.2 +76.3*x - 50.3*x^2,
      
      age == "0_4 years" & sex == "female" ~ 0.2 +24.2*x - 16.5*x^2,
      age == "5_9 years" & sex == "female" ~ 0.4 +47*x - 33.5*x^2,
      age == "10_14 years" & sex == "female" ~ 0.2 +44*x - 28.7*x^2,
      age == "15_19 years" & sex == "female" ~ 0.2 +44*x - 28.7*x^2,
    ))
}
