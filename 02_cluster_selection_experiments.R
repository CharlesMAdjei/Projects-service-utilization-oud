# --> Libraries 
###############################################################################
library(tidyverse)
library(foreach)
library(doParallel)
library(cluster)       # For PAM clustering and Gower distance
library(factoextra)    # For silhouette visualization
library(ggplot2)       # For saving plots
library(tictoc)
library(gridExtra)


# --> Cluster Experiment Parameters 
###############################################################################

# Specify min and max number of clusters for PAM
min_clusters = 2
max_clusters = 20

# Set normalized = 1 if distance metric should be normalized, 0 if not 
normalized = 0

# Directories
path = "\\\\cdc.gov\\project\\NCIPC_DUIP_HST_Applied_Research\\Buprenorphine_Service_Utilization"
data_type = "mixed"
sample_nm = "sample_10_seed_2146"


# --> Create Directories per above parameters 
###############################################################################
normalized_nm = ifelse(normalized == 0, "Reg", "Norm")

path_full <- paste(path, "Results", data_type, sample_nm, normalized_nm, sep = "\\")
if (dir.exists(path_full)) {
  unlink(path_full, recursive = TRUE)  # Delete the existing directory and its contents
}
dir.create(file.path(path_full), showWarnings = TRUE, recursive = TRUE)  # Create a new directory


# --> Input Data 
###############################################################################
data_path = file.path(path, "Data", paste(data_type, sample_nm, sep = "\\"))
sampled_df <- read_csv(data_path)


# --> Data preparation --------------------------------------------------------
###############################################################################

# --> Convert string columns to factors 
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

# --> Create age groups for stratified sampling
sampled_df <- sampled_df |> 
  select(!`...1`) |> 
  mutate(across(all_of(factor_columns), as.factor))


clust_df <- sampled_df |> 
  select(-c(patient_id, 
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
  ))



# --> Gower distance calculation for mixed data
gower_dist <- daisy(clust_df, metric = "gower")
summary(gower_dist)


# Optional normalization step (not typically used with Gower)
if (normalized == 1) {
  clust_df <- BBmisc::normalize(clust_df, method = "standardize")
}


# --> Parallel Processing Setup -----------------------------------------------
###############################################################################
cores <- detectCores() - 1  # Detect the number of cores and use one less to avoid overloading
cl <- makeCluster(cores)
registerDoParallel(cl)


# Set the seed for reproducibility
set.seed(1234)


# --> Initialize list to store silhouette widths
silhouette_avg_widths <- foreach(k = min_clusters:max_clusters, 
                                 .packages = c('cluster', 'factoextra', 'ggplot2', 'dplyr')) %dopar% {
                                   
                                   # Perform PAM clustering using Gower distance
                                   clust_pam <- pam(gower_dist, k = k, diss = TRUE, cluster.only = FALSE)
                                   
                                   # Cluster validation using silhouette scores
                                   silhouette_scores <- silhouette(clust_pam$clustering, gower_dist)
                                   silhouette_summary <- summary(silhouette_scores)
                                   
                                   # Save cluster solution
                                   cluster_file <- file.path(path_full, paste0("Cluster_solution_", k, ".rdata"))
                                   save(clust_pam, file = cluster_file)
                                   
                                   # Save silhouette average width for this k
                                   avg_sil_width <- silhouette_summary$avg.width
                                   
                                   # Write silhouette summary to CSV
                                   sil_widths <- data.frame(
                                     Cluster = silhouette_summary$clus.avg.widths,
                                     Size = silhouette_summary$clus.sizes,
                                     Avg_Width = silhouette_summary$avg.width
                                   )
                                   
                                   csv_file <- file.path(path_full, paste0("Indices_", k, ".csv"))
                                   write.table(as.character(Sys.time()), csv_file, sep = ",", eol = "\n", append = TRUE, row.names = FALSE, col.names = FALSE)
                                   write.table(sil_widths, csv_file, sep = ",", eol = "\n", append = TRUE, row.names = TRUE, col.names = TRUE)
                                   
                                   # Save silhouette plot
                                   file_name_sil <- file.path(path_full, paste0("Silhouette_plot_", k, "_clusters.png"))
                                   g <- fviz_silhouette(silhouette_scores, print.summary = FALSE)
                                   ggsave(file_name_sil, plot = g)
                                   
                                   # Return average silhouette width for this iteration
                                   return(avg_sil_width)
                                 }

# --> Stop parallel cluster
stopCluster(cl)
registerDoSEQ()


# --> Convert silhouette_avg_widths to a data frame
silhouette_df <- data.frame(
  k = min_clusters:max_clusters,
  avg_width = unlist(silhouette_avg_widths)
)

# --> Plot the overall average silhouette width across all k
ggplot(silhouette_df, aes(x = as.factor(k), y = avg_width)) +
  geom_line(group = 1) +  # No color, using group = 1 to connect points
  geom_point() +
  labs(title = "Average Silhouette Width vs Number of Clusters",
       x = "Number of Clusters (k)",
       y = "Average Silhouette Width") +
  theme(
    panel.background = element_rect(fill = "white"),  # White background
    panel.grid = element_blank(),  # Remove gridlines
    axis.line = element_line(color = "black"),  # Add axis lines
    plot.title = element_text(hjust = 0.5)  # Center the plot title
  )

# Save the overall silhouette width plot
ggsave(file.path(path_full, "Overall_Silhouette_Width.png"))
