# 🧬 Data Cleaning - GSE88884 SLE Dataset

This branch contains the **data preprocessing and cleaning steps** for our biomedical project using the GSE88884 dataset.

---

## 📄 Overview

We processed and cleaned the following datasets:

1. **Gene Expression Data**  
   Pre-processed microarray gene expression matrix, filtered and matched to clinical samples.

2. **Phenotypic (Clinical) Data**  
   Metadata including age, sex, race, complement levels, anti-dsDNA antibody status, and SLE Disease Activity Index (SLEDAI) scores.

---

## 📂 Files in this branch

| File                                         | Description                                       |
|----------------------------------------------|---------------------------------------------------|
| `scripts/01_data_cleaning.qmd`              | Quarto file containing the data cleaning workflow |
| `.gitignore`                                | To ignore large data files and temporary files    |
| `README.md`                                 | This file (branch-specific instructions)          |

---

## ✅ Cleaned Data

We created the following cleaned data objects:

| Object               | Description                                             |
|----------------------|---------------------------------------------------------|
| `expr_df_matched`    | Cleaned gene expression matrix (probes × samples)       |
| `pheno_clean`        | Cleaned phenotypic (clinical) data                      |

---

## 📥 Access to Cleaned Data

The cleaned datasets (`.rds` format) are **NOT stored in GitHub** due to size limits.

Please download them from the following link:
➡️ https://drive.google.com/drive/folders/1X8KMvZ5hE9LekCa_Tg7TOucRaRggEeFA?usp=drive_link

Once downloaded, place them in the `data/` folder:
