<h1>DATA3888 Capstone Final Submission</h1>

<hr>

<h2>Project Overview</h2>
<p>
  This repository contains the final scripts, workflow chart, Shiny app, trained model objects, and demonstration/test data for our DATA3888 capstone project on gene expression-based SLE flare risk prediction.
</p>

<hr>

<h2>Directory Structure</h2>
<ul>
  <li><b>SLE_report.qmd / SLE_report.html</b>: Full reproducible report (Quarto/RMarkdown) for analysis and visualization.</li>
  <li><b>Workflow_diagram.jpeg</b>: Workflow chart showing the steps from data processing to deployment.</li>
  <li>
    <b>Shiny app/</b>
    <ul>
      <li><b>sle_app.R</b>: Shiny app for interactive risk prediction and visualization.</li>
      <li><b>df_test2.rds</b>: Processed metadata for test/validation samples (used in app and plotting).</li>
      <li><b>lasso_df3_50_model.rds</b>: Final trained LASSO Cox model (top 50 genes, rate-of-change features).</li>
      <li><b>risk_scores_df3_50.rds</b>: Precomputed risk scores for risk stratification and plotting.</li>
      <li><b>gene_info.csv</b>: Gene symbol annotation for interpretation.</li>
      <li><b>references.bib</b>: References for all code and data sources.</li>
    </ul>
  </li>
  <li><b>sample patients/</b>: Example patient data for app demonstration (CSV format, e.g., <code>flare_patient1.csv</code>).</li>
</ul>

<hr>

<h2>Data Availability</h2>
<ul>
  <li>
    <b>Original Data Download Links:</b>
    <ul>
      <li><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88884">GSE88884 (NCBI GEO)</a></li>
      <li><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65391">GSE65391 (NCBI GEO)</a></li>
    </ul>
  </li>
  <li>
    <b>Processed Data (Google Drive):</b>
    <ul>
      <li><a href="https://drive.google.com/drive/folders/1bzE2tSTAHwzK0T1qL-SuERT0yEYOGAbv">Link to processed datasets</a></li>
    </ul>
  </li>
</ul>
<blockquote>
  <b>Note:</b><br>
  Due to large file sizes, raw and processed datasets are <b>not included directly in this GitHub repo</b>. Please use the above links to download the required data as needed.
</blockquote>

<hr>

<h2>How to Use</h2>
<ol>
  <li>
    <b>Reproduce the Full Analysis</b><br>
    Open <code>SLE_report.qmd</code> in RStudio or Quarto, run all code chunks to reproduce analysis, figures, and results.  
    Ensure all R packages listed in the script are installed.
  </li>
  <li>
    <b>Run the Shiny App</b><br>
    Go to the <code>Shiny app/</code> directory, open <code>sle_app.R</code> in R, and run the app.<br>
    Upload demonstration files from <code>sample patients/</code> (e.g., <code>flare_patient1.csv</code>) to test predictions and risk visualization.
  </li>
  <li>
    <b>Model and Output Files</b><br>
    <code>lasso_df3_50_model.rds</code>: Load to reproduce or extend predictions.<br>
    <code>risk_scores_df3_50.rds</code>, <code>df_test2.rds</code>: Used for risk stratification, survival analysis, and all visualizations in report/app.
  </li>
</ol>

<hr>

<h2>Notes</h2>
<ul>
  <li>All workflows, models, and app were built based on the <b>GSE65391</b> dataset.</li>
  <li>For technical issues, refer to the comments in the scripts or contact the project team.</li>
</ul>
