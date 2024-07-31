# This script provides an example for how to use the helper functions and 
# models in the "model" folder to make predictions for a given population.

# This script will load the pcti_preds() function.
# pcti_preds() also calls the augments.R script to make adjustments for 
# Burkitt lymphoma and Kaposi sarcoma as needed


# Load packages and augments helper functions -----------------------------

# Uncomment and run this section to install required packages if not already installed
# packages <- c("tidyverse", "tidybayes")
# new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
# if(length(new_packages) > 0) install.packages(new_packages)

library(tidyverse)
library(tidybayes)

# This loads the pcti_preds helper function. 
# You will need to change this path if you saved the script in another folder. 
source("scripts/2-Helpers/pcti_preds.R")

# Make example data 
standard_pop <- expand_grid(
  # All predictions using pcti_preds requires one or more regions (labeled geo_units)
  geo_unit = "Hobbiton",
  # Create sex categories 
  sex = c("male", "female"),
  # Create age categories
  age = c("0_4 years", "5_9 years", "10_14 years", "15_19 years"),
  # Specify the population size - 1 million for each age group 
  population = 1e6,
  # The model requires a log(population) value 
  logpop = log(population),
  # If adjustments for KS or BL are desired - this represents the proportion of the pop in a given 
  # country that lives in a malaria endemic region. The variable can take a value between 0 and 1.
  # If no adjustments are required then the mpr0 variable can be absent from the tibble or take a value of NA.
  mpr0 = 0.5
)

# Example 1 - Predicting all cancers --------------------------------------

# pcti_preds takes four arguments:
#  data: data frame as above in standard_pop. 
#           The variables "geo_unit", "age", "sex", "population", and "logpop" must be present in the data frame. 
#  level: "total", "group", or "subgroup" - representing the levels of the ICCC-3
#  disease: a specific disease name for single disease predictions. 
#           If no disease is specified, then predictions for each disease within the group or subgroup will be provided.
#           The text has to exactly match the lower case disease name. If you would like to find a disease name
#           then run the pcti_preds() function without specifying a disease and then inspect the unique
#           values of the `disease` column of the resulting tibble.
#  numb.draws: how many simulations to run. The default is 1000

# The result is a tibble with `numb.draws` number of rows for each
# geo_unit-, age-, sex-specific population for each disease indexed by the `.draw' column. 
# Setting level = "total" gives predictions for all cancers combined. 
pred_tib <- pcti_preds(
  data = standard_pop,
  level = "total",
  numb.draws = 500
) 

print(pred_tib)

# Summarize the median and 95% prediction interval for total cancers. 
# Because there were one million people in each age and sex category in our `standard_pop` data, 
# these results can be interpreted as the incidence per million in that category.  
pred_tib %>% 
  group_by(geo_unit, sex, age) %>% 
  median_qi(.prediction)


# Example 2: Prediction Cancer Groups  ------------------------------------

# The produces is a tibble with `numb.draws` number of rows for each
# geo_unit-, age-, sex-specific population for each disease indexed by the `.draw' column. 
# Setting level = "group" gives predictions for the 12 ICCC-3 cancer groups and the "not classified by iccc or in situ" category 
pred_tib2 <- pcti_preds(
  data = standard_pop,
  level = "group",
  numb.draws = 500
) 

print(pred_tib2)

# The disease groups 
pred_tib2 %>% distinct(disease)

# Now we summarize the median (95% PI) for each age and sex category for each disease group
# Since there is 1 geo_unit, 4 age categories, 2 sexes, and 13 disease groups, there are 1*4*2*13 = 104 rows in the resulting tibble
pred_tib2 %>% 
  group_by(geo_unit, sex, age, disease) %>% 
  median_qi(.prediction)


# Example 3: Prediction Cancer Subgroups  ------------------------------------

# This produces is a tibble with `numb.draws` number of rows for each geo_unit-age-sex-disease combination. 
# There are 49 disease subgroups 
pred_tib3 <- pcti_preds(
  data = standard_pop,
  level = "subgroup",
  numb.draws = 500
) 

print(pred_tib3)

# The disease subgroups
pred_tib3 %>% distinct(disease)

# Now we summarize the median (95% PI) for each age and sex category for each disease group
# Since there is 1 geo_unit, 4 age categories, 2 sexes, and 49 diseases, there are 1*4*2*49 = 392 rows in the resulting tibble
pred_tib3 %>% 
  group_by(geo_unit, sex, age, disease) %>% 
  median_qi(.prediction)


# Example 4: Disease specific predictions  --------------------------------

# This will produce predictions only for kaposi sarcoma and acute myeloid leukemia. 
# The kaposi incidence will be adjusted since we said the Habbiton has an mpr0 of 0.5. 

pred_tib4 <- pcti_preds(
  data = standard_pop,
  level = "subgroup",
  numb.draws = 500,
  disease = c("kaposi sarcoma", "acute myeloid leukemias")
) 

print(pred_tib4)

# The disease subgroups
pred_tib4 %>% distinct(disease)

# Now we summarize the median (95% PI) for each age and sex category for each disease group
# Since there is 1 geo_unit, 4 age categories, 2 sexes, and 49 diseases, there are 1*4*2*49 = 392 rows in the resulting tibble
pred_tib4 %>% 
  group_by(geo_unit, sex, age, disease) %>% 
  median_qi(.prediction)
