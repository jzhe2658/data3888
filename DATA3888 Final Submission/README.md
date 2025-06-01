# DATA3888 Capstone Final Submission

## Project Overview

This repository contains the final scripts, Shiny app, and test data for our DATA3888 capstone project based on gene expression datasets for SLE prediction.

---

## Directory Structure

- `08_feature_model.qmd`: Final main workflow script for analysis and modeling.
- `shiny-app/`: Contains the Shiny app code for interactive risk prediction.
- `sample_patients/`: Example patient data files for demonstration and testing.

---

## Data Availability

- **Original Data Download Links:**
  - [GSE88884 (NCBI GEO)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88884)
  - [GSE65391 (NCBI GEO)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65391)

- **Processed Data (Google Drive):**  
  [Link to processed datasets](https://drive.google.com/drive/folders/1bzE2tSTAHwzK0T1qL-SuERT0yEYOGAbv)

**Note:**  
Due to large file sizes, raw and processed datasets are **not included directly in this GitHub repo**. Please use the above links to download the required data as needed.

---

## How to Use

1. **Reproduce the Analysis**
   - Open `08_feature_model.qmd` with RStudio or Quarto.
   - Follow the code chunks to run the analysis workflow.
   - Ensure required R packages are installed as indicated in the script.

2. **Test the Shiny App**
   - Go to the `shiny-app/` directory and open `sle_app.R` in R.
   - Run the app; upload patient data from `sample_patients/` for demonstration.
   - Make sure to install any necessary dependencies as described in the script.

---

## Notes

- All scripts and app are built for the **GSE65391** dataset, with process documentation included.
- For any issues or questions, please refer to script comments or contact the project team.

---
