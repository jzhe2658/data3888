# DATA3888
This project uses gene expression and clinical data from the GSE88884 study to explore biomarkers and develop a risk prediction model for Systemic Lupus Erythematosus (SLE).

## 📂 Project Structure

- **data/**: Raw and cleaned datasets.
- **scripts/**: R code files for data cleaning, analysis, and modelling.
- **results/**: Output files (plots, tables).
- **shiny_app/**: Optional folder for future deployment of interactive risk calculator.

## Data Access
The large raw datasets are NOT uploaded to GitHub due to size limits. Download datasets from this link:
``` perl
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88884
```
Place them in the data/ folder to run the analysis.

## Branch Workflow

We follow a **branch-per-task workflow** to collaborate efficiently:

| Branch Name             | Purpose                                                   |
|------------------------|-----------------------------------------------------------|
| `main`                 | Stable and merged version of all work                     |
| `data-cleaning`        | Data preprocessing (completed)                             |
| `eda-analysis`        | Exploratory Data Analysis ()                 |
| `differential-analysis`| Differential Expression Analysis ()         |
| `model-training`      | Predictive Model Development ()             |
| `shiny-app`           | Shiny App Deployment (later)                    |


## 👥 How to Contribute

1. Clone the repository:
```bash
git clone https://github.com/jzhe2658/data3888.git
cd data3888
```
2. Checkout the data-cleaning branch
```bash
git checkout data-cleaning
```

3. Create your own feature branch
For example, if you are doing EDA: (please follow the table aboove)
```bash
git checkout eda-analysis
```

4. Work on your scripts
Place your .qmd/.Rmd or .R script in the scripts/ folder.

5. Add and commit your changes
```bash
git add scripts/02_eda_analysis.qmd
git commit -m "Add EDA analysis script"
git push -u origin eda-analysis
```

6. Open a Pull Request
Once your analysis is complete, open a pull request to merge your branch into main.

⚠️ Data Size Notice
Do NOT commit large dataset files (.txt, .rds, .csv) to the repository.
Only analysis scripts and documentation should be committed.
