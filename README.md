# Health Service Utilization and Patient Profiles Among Individuals Initiating Buprenorphine

This repository contains R scripts and documentation for a research study analyzing healthcare service utilization patterns among patients initiating buprenorphine treatment for opioid use disorder (OUD). The objective of this study is to examine the types and frequency of services used prior to treatment initiation and to identify distinct patient profiles using cluster analysis.

## 📚 Study Overview

The study leverages national all-payer linked prescription and healthcare encounter claims data to:

- Characterize the health service utilization patterns of individuals in the 6 months prior to initiating buprenorphine.
- Identify patient subgroups (clusters) based on these utilization patterns using unsupervised machine learning techniques (Partitioning Around Medoids – PAM).
- Inform intervention strategies and targeted treatment pathways based on these profiles.

## 📁 Repository Structure

### `01_sampling_sensitivity_final_analysis.R`

This script performs the following tasks:

- Creates 10 stratified samples for robust model comparison and reproducibility.
- Conducts sensitivity analyses to assess robustness across different samples and stratifications.
- Executes the final analysis pipeline for generating clean datasets and summary tables for reporting.

### `02_cluster_selection_experiments.R`

This script is used to:

- Run experiments for selecting the optimal number of clusters using multiple validity metrics (e.g., silhouette score, gap statistic).
- Generate diagnostic plots to visualize clustering performance across different `k` values.
- Document the rationale for selecting the final number of clusters used in PAM.

### `03_run_final_PAM_clustering.R`

This script implements the final clustering model:

- Applies the PAM algorithm using the previously selected optimal number of clusters.
- Generates cluster membership assignments and appends them to the final analytic dataset.
- Produces summary statistics and visualizations for each identified cluster to support reporting and manuscript development.
