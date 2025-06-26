# --> Libraries 
################################################################################
library(tidyverse)
library(foreach)
library(doParallel)
library(cluster)       # For PAM clustering and Gower distance
library(factoextra)    # For silhouette visualization
library(ggplot2)       # For creating plots
library(tictoc)        # For timing code execution
library(gridExtra)     # For arranging plots
library(Rtsne)         # For t-SNE visualization
library(writexl)       # For exporting results to Excel
library(car)



source("\\\\cdc.gov\\project\\NCIPC_DUIP_HST_Applied_Research\\Buprenorphine_Service_Utilization\\Code\\custom_functions.R")

# Read the non_clustering_variables
################################################################################
data_path0 <- "\\\\cdc.gov\\project\\NCIPC_DUIP_HST_Applied_Research\\Buprenorphine_Service_Utilization\\Data\\raw\\non_clustering_variables_11122024.csv"
non_clustering_vars_df <- read_csv(data_path0, na = c("N/A", ""))



# --> Cluster Experiment Parameters 
################################################################################

# Specify minimum and maximum number of clusters for PAM
min_clusters <- 2
max_clusters <- 20

# Set normalized = 1 if distance metric should be normalized, 0 if not
normalized <- 0

# Define directory paths
path <- "\\\\cdc.gov\\project\\NCIPC_DUIP_HST_Applied_Research\\Buprenorphine_Service_Utilization"
data_type <- "mixed"
sample_nm <- "sample_6_seed_9196"

# --> Input Data 
################################################################################
data_path <- file.path(path, "Data", paste(data_type, sample_nm, sep = "\\"))

# Load the data
sampled_df <- read_csv(data_path)





# --> Data Preparation 
################################################################################

# Convert selected columns to factors
factor_columns <- c(
  "patient_gender", 
  "most_frequent_region", 
  "most_frequent_known_payer",
  "alcohol_and_drug_management_treatment_and_rehabilitation", 
  "consultation_evaluation_and_preventative_care", 
  "laboratory", 
  "physical_occupational_and_speech_therapy_exercises_manipulation_and_other_procedures",
  "psychological_and_psychiatric_evaluation_and_therapy"
)

# Create age groups for stratified sampling and convert selected columns to factors
sampled_df <- sampled_df |> 
  select(!`...1`) |> 
  mutate(across(all_of(factor_columns), as.factor))





# --> Overall Summary 
################################################################################
overall_summary <- sampled_df |> 
  summarise(across(
    where(is.numeric), 
    list(
      median = ~ median(.x, na.rm = TRUE),
      mean = ~ mean(.x, na.rm = TRUE)
    )
  ))





# --> Select Columns for Cluster Analysis
################################################################################
clust_df <- sampled_df |> 
  select(
    -c(patient_id, 
       last_svc_fr_dt, 
       init_date,
       age_at_init,
       age_cat, 
       patient_gender, 
       most_frequent_region,
       most_frequent_known_payer,
       total_other_visits,
       total_psych_mental_health_visits, 
       total_lab_visits
    )
  )






# --> Run the PAM Algorithm and Display Summary Stats
################################################################################

# Compute Gower's Distance
gower_dist <- daisy(clust_df, metric = "gower")
summary(gower_dist)


set.seed(1234)  # For reproducibility

# Perform PAM clustering for k = 11
pam_fit <- pam(gower_dist, diss = TRUE, k = 11)





# --> Cluster Validation Using Silhouette Scores
################################################################################
silhouette_scores <- silhouette(pam_fit$clustering, gower_dist)
silhouette_summary <- summary(silhouette_scores)

# Extract cluster-specific silhouette widths and sizes
cluster_silhouettes <- data.frame(
  Cluster = silhouette_summary$clus.avg.widths,
  Size = silhouette_summary$clus.sizes,
  Avg_Width = silhouette_summary$avg.width
)

print(cluster_silhouettes)  # Display silhouette summary





# --> Display the Medoids Representing Each Cluster
################################################################################
# View the medoids from the sampled_df
medoids_df <- sampled_df[pam_fit$medoids, ]
View(medoids_df)

# Export medoids data to an Excel file
write_xlsx(medoids_df, "pam_medoids_results.xlsx")








# --> Add cluster assignments to the sample_df 
################################################################################
pam_results <- sampled_df |> 
  inner_join(non_clustering_vars_df |> janitor::clean_names(), join_by(patient_id)) |> 
  
  select(
    age_at_init, 
    patient_gender, 
    most_frequent_region, 
    most_frequent_known_payer, 
    provider_specialty_at_init,
    most_frequent_specialty,
    days_prior_init, 
    total_visits, 
    total_procedures, 
    total_outpatient_visit, 
    is_inpatient, 
    total_inpatient_visit, 
    last_inpatient_days_prior_init, 
    is_ed, 
    total_ed_visits, 
    last_ed_days_prior_init, 
    alcohol_and_drug_management_treatment_and_rehabilitation, 
    alcohol_and_drug_management_treatment_and_rehabilitation_days_before_init, 
    consultation_evaluation_and_preventative_care, 
    consultation_evaluation_and_preventative_care_days_before_init, 
    physical_occupational_and_speech_therapy_exercises_manipulation_and_other_procedures, 
    physical_occupational_and_speech_therapy_exercises_manipulation_and_other_procedures_days_before_init, 
    psychological_and_psychiatric_evaluation_and_therapy, 
    psychological_and_psychiatric_evaluation_and_therapy_days_before_init, 
    laboratory, 
    laboratory_days_before_init,
    top_5_alcoho_prc,
    top_5_consult_prc,
    top_5_psycho_prc,
    top_5_physical_prc,
    top_5_lab_prc) |>  
  
  mutate(cluster = pam_fit$clustering,
         is_inpatient = factor(if_else(total_inpatient_visit > 0, "Yes", "No")),
         is_ed = factor(if_else(total_ed_visits > 0, "Yes", "No")))







# --> Results for heatmap (clustering variables only)
################################################################################

# Columns to calculate mean values
numeric_cols <- c(
  "days_prior_init", "total_visits", "total_procedures", "total_outpatient_visit",
  "total_inpatient_visit", "total_ed_visits"
)

categorical_cols <- c(
  "alcohol_and_drug_management_treatment_and_rehabilitation",
  "consultation_evaluation_and_preventative_care",
  "physical_occupational_and_speech_therapy_exercises_manipulation_and_other_procedures",
  "psychological_and_psychiatric_evaluation_and_therapy",
  "laboratory"
)

# Calculate cluster-specific means
cluster_grp_df <- pam_results %>%
  group_by(cluster) %>%
  summarise(
    across(all_of(numeric_cols), ~ mean(.x, na.rm = TRUE)),
    across(all_of(categorical_cols), ~ mean(.x == "Yes", na.rm = TRUE))
  )

# Calculate overall means
overall_means <- pam_results %>%
  summarise(
    across(all_of(numeric_cols), ~ mean(.x, na.rm = TRUE)),
    across(all_of(categorical_cols), ~ mean(.x == "Yes", na.rm = TRUE))
  ) %>%
  mutate(cluster = "overall")

# Combine cluster-specific and overall means for the heat map
clustering_vars_heatmap <- cluster_grp_df %>%
  pivot_longer(cols = -cluster, names_to = "metric", values_to = "Value") %>%
  pivot_wider(names_from = cluster, values_from = Value, names_prefix = "cluster_") %>%
  left_join(
    overall_means %>%
      pivot_longer(cols = -cluster, names_to = "metric", values_to = "Value") %>%
      pivot_wider(names_from = cluster, values_from = Value, names_prefix = "cluster_"),
    by = "metric"
  ) %>%
  relocate(cluster_overall, .after = 1)  

write_xlsx(clustering_vars_heatmap, "\\\\cdc.gov\\project\\NCIPC_DUIP_HST_Applied_Research\\Buprenorphine_Service_Utilization\\Results\\clustering_vars_heatmap.xlsx")







# --> Output for cluster grouping
################################################################################
cluster_grp_tbl_df <- cluster_grp_df |>
  # Reverse the logic for "percentile_number_of_days_between_the_last_service_date_and_the_initiation_date"
  mutate(
    percentile_days_prior_init = case_when(
      days_prior_init <= quantile(days_prior_init, 0.25, na.rm = TRUE) ~ "High",
      
      days_prior_init > quantile(days_prior_init, 0.25, na.rm = TRUE) &
        days_prior_init <= quantile(days_prior_init, 0.75, na.rm = TRUE) ~ "Moderate",
      days_prior_init > quantile(days_prior_init, 0.75, na.rm = TRUE) ~ "Low",
      TRUE ~ NA_character_)) |>
  
  #apply logic to all other columns
  mutate(across(.cols = where(is.numeric) & !contains("days_prior_init") & !contains("cluster"),
                .fns = ~case_when(
                  . <= quantile(., 0.25, na.rm = TRUE) ~ "Low",
                  . > quantile(., 0.25, na.rm = TRUE) & . <= quantile(., 0.75, na.rm = TRUE) ~ "Moderate",
                  . > quantile(., 0.75, na.rm = TRUE) ~ "High",
                  TRUE ~ NA_character_),
                .names = "percentile_{.col}")) |>
  mutate(row_id = row_number())



# Calculate mode for each row's percentile columns
mode_df <- cluster_grp_tbl_df |> 
  pivot_longer(cols = starts_with("percentile_"), names_to = "percentile", values_to = "value") |> 
  group_by(row_id) |> 
  summarize(mode_value = Mode(value), .groups = 'drop')  # Ensure groups are dropped after summarization

# Join back with the original data
cluster_grp_tbl_df <- left_join(cluster_grp_tbl_df, mode_df, by = "row_id")

# Optionally remove the row_id column if no longer needed
cluster_grp_tbl_df <- cluster_grp_tbl_df |> 
  select(cluster, contains("percentile"), mode_value)|> 
  pivot_longer(cols = -cluster, names_to = "variable", values_to = "Value") |> 
  pivot_wider(names_from = cluster, values_from = Value, names_prefix = "cluster_") 

write_xlsx(cluster_grp_tbl_df, "\\\\cdc.gov\\project\\NCIPC_DUIP_HST_Applied_Research\\Buprenorphine_Service_Utilization\\Results\\cluster_grouping_tbl.xlsx")


cluster_grp_assignment_df <- pivot_longer(cluster_grp_tbl_df,
                                          cols = starts_with("cluster"),
                                          names_to = "cluster",
                                          values_to = "cluster_group") |> 
                             filter(variable == "mode_value") |> 
                             separate(cluster, into = c("prefix", "cluster_num"), sep = "_") |> 
                             select(cluster_num, cluster_group) |> 
                             mutate(cluster_num = as.integer(cluster_num)) |> 
                             inner_join(pam_results, by = c("cluster_num" = "cluster")) |> 
                             mutate(cluster_num = as.factor(cluster_num))

write.csv(cluster_grp_assignment_df, "\\\\cdc.gov\\project\\NCIPC_DUIP_HST_Applied_Research\\Buprenorphine_Service_Utilization\\Results\\cluster_grp_assignment_df.csv")



                




# -----> Detail output
################################################################################

#---> Apply the function to each categorical variables
gender_df <- calculate_overall_percent(pam_results, patient_gender, "patient_gender") |> 
  inner_join(
    calculate_cluster_percent(pam_results, cluster, patient_gender, "patient_gender") |> 
      select(-variable) ,
    by="categories"
    )


region_df <- calculate_overall_percent(pam_results, most_frequent_region, "most_frequent_region") |> 
  inner_join(
    calculate_cluster_percent(pam_results, cluster, most_frequent_region, "most_frequent_region") |> 
      select(-variable),
    by="categories"
    )


payer_df <- calculate_overall_percent(pam_results, most_frequent_known_payer, "most_frequent_known_payer") |> 
  inner_join(
    calculate_cluster_percent(pam_results, cluster, most_frequent_known_payer, "most_frequent_known_payer") |> 
      select(-variable),
    by="categories"
    )


specialty_at_init_df <- calculate_overall_percent(pam_results, provider_specialty_at_init, "provider_specialty_at_init") |> 
  inner_join(
    calculate_cluster_percent(pam_results, cluster, provider_specialty_at_init, "provider_specialty_at_init") |> 
      select(-variable),
    by="categories"
    )


frequent_specialty_df <- calculate_overall_percent(pam_results, most_frequent_specialty, "most_frequent_specialty") |> 
  inner_join(
    calculate_cluster_percent(pam_results, cluster, most_frequent_specialty, "most_frequent_specialty") |> 
      select(-variable),
    by="categories"
    )


inpatient_df <- calculate_overall_percent(pam_results, is_inpatient, "is_inpatient") |> 
  inner_join(
    calculate_cluster_percent(pam_results, cluster, is_inpatient, "is_inpatient") |> 
      select(-variable) ,
    by="categories"
    )


ed_df <- calculate_overall_percent(pam_results, is_ed, "is_ed") |> 
  inner_join(
    calculate_cluster_percent(pam_results, cluster, is_ed, "is_ed") |> 
      select(-variable),
    by="categories")


#---> Apply the function to each Numeric variables
# Calculate means
means_df <- calculate_overall_summary(pam_results, mean, "mean") |>
  rename(cluster_overall = Value) |> 
  inner_join(calculate_summary(pam_results, mean, "mean") |> select(-categories), by = "variable") |> 
  relocate(categories, .after = 1)

# Calculate medians
medians_df <- calculate_overall_summary(pam_results, median, "median") |>
  rename(cluster_overall = Value) |> 
  inner_join(calculate_summary(pam_results, median, "median") |> select(-categories), by = "variable") |> 
  relocate(categories, .after = 1)

#---> Combine all summary dataframes into one
detail_summary_df <- bind_rows(
  gender_df, 
  region_df, 
  payer_df, 
  specialty_at_init_df,
  frequent_specialty_df,
  inpatient_df, 
  ed_df,
  medians_df,
  means_df)

write_xlsx(detail_summary_df, "\\\\cdc.gov\\project\\NCIPC_DUIP_HST_Applied_Research\\Buprenorphine_Service_Utilization\\Results\\detail_summary.xlsx")






#---> Top 5 PRC description for the top 5 services
top_5_alcoho_prc_df <- calculate_overall_percent(pam_results, top_5_alcoho_prc, "top_5_alcoho_prc") |> 
  inner_join(
    calculate_cluster_percent(pam_results, cluster, top_5_alcoho_prc, "top_5_alcoho_prc") |> 
      select(-variable) ,
    by="categories"
  )


top_5_consult_prc_df <- calculate_overall_percent(pam_results, top_5_consult_prc, "top_5_consult_prc") |> 
  inner_join(
    calculate_cluster_percent(pam_results, cluster, top_5_consult_prc, "top_5_consult_prc") |> 
      select(-variable),
    by="categories"
  )


top_5_psycho_prc_df <- calculate_overall_percent(pam_results, top_5_psycho_prc, "top_5_psycho_prc") |> 
  inner_join(
    calculate_cluster_percent(pam_results, cluster, top_5_psycho_prc, "top_5_psycho_prc") |> 
      select(-variable),
    by="categories"
  )


top_5_physical_prc_df <- calculate_overall_percent(pam_results, top_5_physical_prc, "top_5_physical_prc") |> 
  inner_join(
    calculate_cluster_percent(pam_results, cluster, top_5_physical_prc, "top_5_physical_prc") |> 
      select(-variable),
    by="categories"
  )


top_5_lab_prc_df <- calculate_overall_percent(pam_results, top_5_lab_prc, "top_5_lab_prc") |> 
  inner_join(
    calculate_cluster_percent(pam_results, cluster, top_5_lab_prc, "top_5_lab_prc") |> 
      select(-variable),
    by="categories"
  )


#---> Combine all summary dataframes into one
top5_prc_des_wihtna_df <- bind_rows(
  top_5_alcoho_prc_df, 
  top_5_consult_prc_df, 
  top_5_psycho_prc_df, 
  top_5_physical_prc_df,
  top_5_lab_prc_df)





#---> Top 5 PRC description for the top 5 services
top_5_alcoho_prc_narm_df <- calculate_overall_percent(pam_results |> filter(!is.na(top_5_alcoho_prc)), top_5_alcoho_prc, "top_5_alcoho_prc") |> 
  inner_join(
    calculate_cluster_percent(pam_results |> filter(!is.na(top_5_alcoho_prc)), cluster, top_5_alcoho_prc, "top_5_alcoho_prc") |> 
      select(-variable) ,
    by="categories"
  )


top_5_consult_prc_narm_df <- calculate_overall_percent(pam_results |> filter(!is.na(top_5_consult_prc)), top_5_consult_prc, "top_5_consult_prc") |> 
  inner_join(
    calculate_cluster_percent(pam_results |> filter(!is.na(top_5_consult_prc)), cluster, top_5_consult_prc, "top_5_consult_prc") |> 
      select(-variable),
    by="categories"
  )


top_5_psycho_prc_narm_df <- calculate_overall_percent(pam_results |> filter(!is.na(top_5_psycho_prc)), top_5_psycho_prc, "top_5_psycho_prc") |> 
  inner_join(
    calculate_cluster_percent(pam_results |> filter(!is.na(top_5_psycho_prc)), cluster, top_5_psycho_prc, "top_5_psycho_prc") |> 
      select(-variable),
    by="categories"
  )


top_5_physical_prc_narm_df <- calculate_overall_percent(pam_results |> filter(!is.na(top_5_physical_prc)), top_5_physical_prc, "top_5_physical_prc") |> 
  inner_join(
    calculate_cluster_percent(pam_results |> filter(!is.na(top_5_physical_prc)), cluster, top_5_physical_prc, "top_5_physical_prc") |> 
      select(-variable),
    by="categories"
  )


top_5_lab_prc_narm_df <- calculate_overall_percent(pam_results |> filter(!is.na(top_5_lab_prc)), top_5_lab_prc, "top_5_lab_prc") |> 
  inner_join(
    calculate_cluster_percent(pam_results |> filter(!is.na(top_5_lab_prc)), cluster, top_5_lab_prc, "top_5_lab_prc") |> 
      select(-variable),
    by="categories"
  )


#---> Combine all summary dataframes into one
top5_prc_des_withoutna_df <- bind_rows(
  top_5_alcoho_prc_narm_df, 
  top_5_consult_prc_narm_df, 
  top_5_psycho_prc_narm_df, 
  top_5_physical_prc_narm_df,
  top_5_lab_prc_narm_df)





# Create a function to reduce code repetition
calculate_and_join <- function(df, column) {
  calculate_overall_percent(df %>%
                              filter(!is.na({{column}})), {{column}}, "{{column}}") %>%
    inner_join(
      calculate_cluster_percent(df %>%
                                  filter(!is.na({{column}})), cluster, {{column}}, "{{column}}"), 
      by = "categories"
    )
}

# Calculate and join for each column
top_5_alcoho_prc_narm_df <- calculate_and_join(pam_results, top_5_alcoho_prc)
top_5_consult_prc_narm_df <- calculate_and_join(pam_results, top_5_consult_prc)
top_5_psycho_prc_narm_df <- calculate_and_join(pam_results, top_5_psycho_prc)
top_5_physical_prc_narm_df <- calculate_and_join(pam_results, top_5_physical_prc)
top_5_lab_prc_narm_df <- calculate_and_join(pam_results, top_5_lab_prc)

# Combine all summary dataframes into one
top5_prc_des_withoutna_df1 <- bind_rows(
  top_5_alcoho_prc_narm_df,
  top_5_consult_prc_narm_df,
  top_5_psycho_prc_narm_df,
  top_5_physical_prc_narm_df,
  top_5_lab_prc_narm_df
)



# Cohort summary
#--> Create Random Samples
################################################################################

library(tidyverse)
library(skimr)


# Load the data
original_file_path = "\\\\cdc.gov\\project\\NCIPC_DUIP_HST_Applied_Research\\Buprenorphine_Service_Utilization\\Data\\raw\\top5_services_mixed.csv"



df <- read.csv(original_file_path)

df <- df %>%
  mutate(
    age_cat = cut(
      age_at_init, 
      breaks = c(18, 44, 64, Inf), 
      labels = c("18-44", "45-64", "65+")
    )
  )

df |> count()

total = nrow(df)
df |> summarise(median_age = median(age_at_init),
                mean_age = mean(age_at_init))
df |> 
  group_by(patient_gender) |> 
  summarise(n = n(), 
            percent = round((n/total)*100,1))

df |> 
  group_by(most_frequent_region) |> 
  summarise(n = n(), 
            percent = round((n/total)*100,1))

df |> 
  group_by(most_frequent_known_payer) |> 
  summarise(n = n(), 
            percent = round((n/total)*100,2))
