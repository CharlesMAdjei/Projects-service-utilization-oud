
#--> Create Random Samples
################################################################################

library(tidyverse)
library(skimr)


# Load the data
original_file_path = "\\\\cdc.gov\\project\\NCIPC_DUIP_HST_Applied_Research\\Buprenorphine_Service_Utilization\\Data\\raw\\top5_services_mixed.csv"
save_path <- "\\\\cdc.gov\\project\\NCIPC_DUIP_HST_Applied_Research\\Buprenorphine_Service_Utilization\\Data\\mixed"


df <- read.csv(original_file_path)

df <- df %>%
  mutate(
    age_cat = cut(
      age_at_init, 
      breaks = c(18, 44, 64, Inf), 
      labels = c("18-44", "45-64", "65+")
    )
  )

# Create a random sample of n = 10,000 without stratification
###############################################################################
set.seed(123)
random_sample_df <- sample_n(df, 10000)
write.csv(random_sample_df, file.path(save_path, "random_sample_df"))




#--> Create 10 different stratified samples (by age, gender, region, and payer) of size 10,000
################################################################################
# Generate 10 random seeds
# random_seeds <- sample(1:10000, 10, replace = FALSE)
random_seed_list <- c(7452, 8016, 7162, 8086, 7269, 9196,  623,  934, 2948, 2146)

# Loop to create 10 samples with different random seeds
for (i in 1:10) {
  set.seed(random_seed_list[i])
  
  sampled_df <- df %>%
    group_by(
      age_cat,
      patient_gender, 
      most_frequent_region,
      most_frequent_known_payer) %>%
    sample_frac(10000 / nrow(df)) %>%
    ungroup()
  
  file_name <- paste0("sample_", i, "_seed_", random_seed_list[i])
  write.csv(sampled_df, file.path(save_path, file_name))
}


