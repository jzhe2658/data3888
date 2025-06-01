# DATA3888 Capstone Final Submission

## Project Overview

This repository contains the final scripts, workflow chart, Shiny app, model objects, and test data for our DATA3888 capstone project on gene expression-based SLE flare risk prediction.

---

## Directory Structure

- `08_feature_model.qmd`: Final main workflow script for analysis and modeling.
- `sle_app.R`: Shiny app code for interactive risk prediction.
- `Workflow_diagram.jpeg`: Workflow chart illustrating data processing, modeling, and deployment steps.
- `df_test2.rds`: Processed metadata for test/validation samples; used for survival analysis, subgrouping, and plotting.
- `lasso_df3_50_model.rds`: The final trained LASSO Cox model (top 50 genes, rate-of-change features) for reproducible prediction.
- `risk_scores_df3_50.rds`: Precomputed risk scores for all samples, used for risk stratification and survival curve plotting.
- `sample_patients/`: Example sample files for Shiny app demonstration and user testing.

---

## Data Availability

- **Original Data Download Links:**
  - [GSE88884 (NCBI GEO)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88884)
  - [GSE65391 (NCBI GEO)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65391)

- **Processed Data (Google Drive):**  
  [Link to processed datasets](https://drive.google.com/drive/folders/1bzE2tSTAHwzK0T1qL-SuERT0yEYOGAbv)

> **Note:**  
> Due to large file sizes, raw and processed datasets are **not included directly in this GitHub repo**. Please use the above links to download the required data as needed.

---

## How to Use

1. **Reproduce the Analysis**
   - Open `08_feature_model.qmd` with RStudio or Quarto.
   - Follow the code chunks to run the workflow.
   - Ensure all required R packages are installed as indicated in the script.

2. **Test the Shiny App**
   - Open `sle_app.R` in R, run the app, and upload sample data from `sample_patients/` for demonstration.

3. **Model and Result Files**
   - `lasso_df3_50_model.rds` can be loaded directly to reproduce or extend predictions.
   - `risk_scores_df3_50.rds` and `df_test2.rds` are used for risk stratification, survival analysis, and visualization.

---

## Notes

- All scripts, models, and app are based on the **GSE65391** dataset.
- For technical issues, please refer to comments in the scripts or contact the project team.
