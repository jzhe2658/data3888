library(shiny)
library(bs4Dash)
library(tidyverse)
library(DT)
library(glmnet)
library(ggplot2)
library(scales)
library(survival)
library(survminer)
library(stringr)
library(readxl)

risk_scores <- readRDS("risk_scores_df3_50.rds")
cox_model <- readRDS("lasso_df3_50_model.rds")

all_genes <- rownames(coef(cox_model, s = "lambda.min"))
all_genes <- all_genes[all_genes != "(Intercept)"]

df_test2 <- readRDS("df_test2.rds")

# Combine immune type and risk group
df_test2$risk_immune_group <- paste(df_test2$immune_type, df_test2$risk_group, sep = " - ")
df_test2$risk_immune_group <- factor(df_test2$risk_immune_group)
# Survival object
surv_obj <- Surv(df_test2$cumulative_time, df_test2$flare_next)
# Fit Kaplan-Meier model
km_fit_combo <- survfit(surv_obj ~ risk_immune_group, data = df_test2)

# cohort medians from df_cox_filtered (reference the Qin et al. (2015) paper in your Shiny app)
neutrophil_median <- 3.71
lymphocyte_median <- 1.37

ui <- dashboardPage(
  title = "Lupus Flare Predictor",
  fullscreen = TRUE,
  dark = FALSE,  # Light theme
  header = dashboardHeader(
    title = dashboardBrand(
      title = "Lupus Flare Risk App",
      color = "primary"
    )
  ),
  sidebar = dashboardSidebar(
    skin = "light",
    status = "primary",
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home"), selected = TRUE),
      menuItem("Predict Flare Risk", tabName = "predict", icon = icon("heartbeat")),
      menuItem("Immune Profile", tabName = "immune", icon = icon("dna"))
    )
  ),
  body = dashboardBody(
    tabItems(
      tabItem(
        tabName = "home",
        jumbotron(
          title = "Welcome to the Lupus Flare Risk App!",
          status = "primary",
          lead = "Upload gene expression data and predict future flare risk using a validated Cox model.",
          btnName = "Learn More",
          href = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7934302/",
          "This tool helps stratify SLE patients by flare risk using gene expression data and a LASSO-regularized Cox model."
        )
      ),
      tabItem(
        tabName = "predict",
        h2("Predict Flare Risk from Gene Expression"),
        fileInput("patient_csv", "Upload CSV", accept = ".csv"),
        helpText("Please upload a single-row CSV with the *rate of change* in gene expression (Visit2 − Visit1)."),
        br(),
        fluidRow(
          box(
            title = "Uploaded Gene Expression",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            DTOutput("patient_table")
          )
        ),
        fluidRow(
          box(
            title = "Risk Summary",
            status = "danger",
            solidHeader = TRUE,
            width = 6,
            uiOutput("risk_summary")
          ),
          box(
            title = "Top Gene Contributions for you",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            plotOutput("gene_plot")
          ),
          box(
            title = "Top Gene Descriptions",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            DTOutput("gene_info_table")
          )
        )
      ),
      tabItem(
        tabName = "immune",
        h2("Immune Type Classification"),
        fluidRow(
          box(
            title = "Input Immune Cell Counts",
            status = "info",
            solidHeader = TRUE,
            width = 6,
            # Pre-fill default values close to the median to hint to the use
            numericInput("neutrophil", "Neutrophil Count", 
                         value = 3.5, min = 0, step = 0.01),
            helpText("Typical range: 1.5 – 8.0 ×10⁹/L"),
            numericInput("lymphocyte", "Lymphocyte Count", 
                         value = 1.2, min = 0, step = 0.01),
            helpText("Typical range: 1.0 – 4.0 ×10⁹/L"),
            actionButton("classify", "Classify Immune Type", class = "btn-primary")
          ),
          box(
            title = "Classification Result",
            status = "success",
            solidHeader = TRUE,
            width = 6,
            uiOutput("immune_result")
          ),
          box(
            title = "Combined Risk & Immune Group",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            uiOutput("user_explanation")
          )
        ),
        box(
          title = "How Your Group Performed Over Time",
          width = 12,
          status = "primary",
          solidHeader = TRUE,
          plotOutput("survival_plot", height = "600px"),
          uiOutput("survival_note")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  patient_df <- reactive({
    req(input$patient_csv)
    read_csv(input$patient_csv$datapath, show_col_types = FALSE) %>% as.data.frame()
  })
  
  output$patient_table <- renderDT({
    req(patient_df())
    
    input_long <- patient_df() %>%
      pivot_longer(
        cols = everything(),
        names_to = "Gene",
        values_to = "Expression"
      ) %>%
      mutate(Expression = round(Expression, 3))
    datatable(input_long, options = list(pageLength = 3))
  })
  
  output$risk_summary <- renderUI({
    req(patient_df())
    input_data <- patient_df()
    logfc <- as.numeric(input_data[1, ])
    names(logfc) <- colnames(input_data)
    print(round(logfc[1:5], 4))
    # Fill in only the 50 selected genes, set others to 0
    full_row <- setNames(rep(0, length(all_genes)), all_genes)
    for (gene in names(logfc)) {
      if (gene %in% names(full_row)) {
        full_row[gene] <- logfc[gene]
      }
    }
    X_patient <- matrix(full_row, nrow = 1)
    colnames(X_patient) <- names(full_row)
    
    log_risk <- predict(cox_model, newx = X_patient, s = "lambda.min", type = "link")[1]
    print(paste("🔢 log_risk:", round(log_risk, 4)))
    risk_values <- as.numeric(risk_scores)
    percentile <- ecdf(risk_values)(log_risk) * 100
    
    print(summary(risk_values))
    
    risk_category <- case_when(
      percentile >= 75 ~ "<span style='color: darkred;'>High Risk</span>",
      percentile >= 25 ~ "<span style='color: #336699;'>Moderate Risk</span>",
      TRUE ~ "<span style='color: darkgreen;'>Low Risk</span>"
    )
    
    recommendation <- if (percentile >= 75) {
      paste0(
        "<br><b style='color:darkred;'>Recommendation:</b><br>",
        "• Schedule more frequent checkups.<br>",
        "• Consider reviewing your current treatment plan for possible early interventions.<br>",
        "• For a more personalized recommendation, go to the <b>Immune Profile</b> tab and input your ",
        "<i>neutrophil</i> and <i>lymphocyte</i> counts to receive immune-stratified long-term advice."
      )
    } else {
      ""
    }
    
    HTML(paste0(
      "<div style='font-size: 18px; line-height: 1.6;'>",
      "<b>Percentile:</b> <span style='font-size: 24px; color: darkorange;'>", round(percentile, 1), "%</span><br>",
      "<small><i>This percentile compares your predicted risk score to all patients in our study. For example, a percentile of 82.7% means your risk score is higher than 82.7% of patients — placing you in the higher-risk group.</i></small><br><br>",
      "<b>Interpretation:</b> ", risk_category,
      recommendation,
      "</div>"
    ))
  })
  
  output$gene_plot <- renderPlot({
    req(patient_df())
    input_data <- patient_df()
    
    input_data <- patient_df()
    logfc <- as.numeric(input_data[1, ])
    names(logfc) <- colnames(input_data)
    # Fill in only the 50 selected genes, set others to 0
    full_row <- setNames(rep(0, length(all_genes)), all_genes)
    for (gene in names(logfc)) {
      if (gene %in% names(full_row)) {
        full_row[gene] <- logfc[gene]
      }
    }
    coef_all <- coef(cox_model, s = "lambda.min")
    coef_vec <- as.vector(coef_all)
    gene_names <- rownames(coef_all)
    valid_genes <- gene_names[coef_vec != 0 & gene_names != "(Intercept)"]
    nonzero_coefs <- coef_vec[coef_vec != 0 & gene_names != "(Intercept)"]
    names(nonzero_coefs) <- valid_genes
    
    contrib_df <- tibble(
      Gene = names(nonzero_coefs),
      Coefficient = as.numeric(nonzero_coefs),
      Expression = as.numeric(full_row[names(nonzero_coefs)]),
      Contribution = Coefficient * Expression
    )
    
    top_contrib <- contrib_df %>%
      slice_max(order_by = abs(Contribution), n = 10) %>%
      arrange(Contribution) %>%
      mutate(
        label = paste0("Expr: ", round(Expression, 2)),
        Contribution = round(Contribution, 2)
      )
    
    ggplot(top_contrib, aes(x = Contribution, y = reorder(Gene, Contribution), fill = Contribution)) +
      geom_col(width = 0.6, color = "black", alpha = 0.95) +
      geom_text(aes(label = label),
                hjust = ifelse(top_contrib$Contribution > 0, -0.1, 1.1),
                size = 3.5, color = "black", fontface = "italic") +
      scale_fill_gradient2(
        low = "#3B4CC0", mid = "gray95", high = "#B22222", midpoint = 0,
        name = "Risk Contribution"
      ) +
      labs(
        title = "Top Gene Contributions to Flare Risk",
        subtitle = "Positive = higher risk, Negative = protective",
        x = "Contribution Score",
        y = NULL
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 13, face = "italic"),
        axis.text.y = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.position = "right",
        panel.grid.major.y = element_blank()
      ) +
      coord_cartesian(xlim = c(min(top_contrib$Contribution) * 1.3,
                               max(top_contrib$Contribution) * 1.3))
  })
  
  output$gene_info_table <- renderDT({
    req(patient_df())
    
    # Load gene info
    gene_info_all <- read_csv("gene_info.csv", show_col_types = FALSE)
    
    # Prepare gene vector from patient
    input_data <- patient_df()
    logfc <- as.numeric(input_data[1, ])
    names(logfc) <- colnames(input_data)
    
    # Fill in only the 50 selected genes, set others to 0
    full_row <- setNames(rep(0, length(all_genes)), all_genes)
    for (gene in names(logfc)) {
      if (gene %in% names(full_row)) {
        full_row[gene] <- logfc[gene]
      }
    }
    
    coef_all <- coef(cox_model, s = "lambda.min")
    coef_vec <- as.vector(coef_all)
    gene_names <- rownames(coef_all)
    valid_genes <- gene_names[coef_vec != 0 & gene_names != "(Intercept)"]
    nonzero_coefs <- coef_vec[coef_vec != 0 & gene_names != "(Intercept)"]
    names(nonzero_coefs) <- valid_genes
    
    # Get top 10 genes by contribution
    top_genes <- tibble(
      Gene = names(nonzero_coefs),
      Contribution = unname(nonzero_coefs * full_row[names(nonzero_coefs)])
    ) %>%
      slice_max(order_by = abs(Contribution), n = 10) %>%
      pull(Gene)
    
    # Merge with description table
    top_gene_info <- gene_info_all %>%
      dplyr::select(Gene, Function, SLE_Relevance) %>%
      arrange(Gene)
    
    # Render interactive searchable table
    datatable(
      as.data.frame(top_gene_info),  # ensure it's clean
      filter = "top",
      rownames = FALSE,
      escape = TRUE,
      options = list(
        pageLength = 5,
        autoWidth = TRUE,
        columnDefs = list(list(width = '40%', targets = c(1, 2)))  # adjust column width
      )
    )
  })
  
  
  # Immune result logic
  output$immune_result <- renderUI({
    req(input$classify)  # Only run when button is clicked
    isolate({
      neutro <- input$neutrophil
      lympho <- input$lymphocyte
      
      if (is.na(neutro) || is.na(lympho)) {
        return(HTML("<span style='color:red;'>Please enter both neutrophil and lymphocyte counts.</span>"))
      }
      
      immune_type <- case_when(
        neutro > neutrophil_median & lympho <= lymphocyte_median ~ "Innate-dominant",
        lympho > lymphocyte_median & neutro <= neutrophil_median ~ "Adaptive-dominant",
        TRUE ~ "Mixed"
      )
      
      HTML(paste0(
        "<div style='font-size: 18px;'>",
        "<b>Immune Type:</b> <span style='color: #2C3E50; font-size: 24px;'>", immune_type, "</span><br>",
        "<b>Explanation:</b> Based on your input values and cohort medians (neutrophils = ", neutrophil_median, 
        ", lymphocytes = ", lymphocyte_median, 
        "), your immune profile is classified as <i>", immune_type, "</i>.",
        "<br><br><i style='font-size: 14px;'>This classification approach is adapted from cohort-based stratification methods used in SLE studies, such as <a href='https://doi.org/10.1016/j.cca.2015.07.019' target='_blank'>Qin et al. (2015)</a>.</i>",
        "</div>"
      ))
    })
  })
  
  # Use both the immune type and risk percentile to describe their immune_risk_group
  output$user_explanation <- renderUI({
    req(input$classify, patient_df())
    
    isolate({
      # Immune type
      neutro <- input$neutrophil
      lympho <- input$lymphocyte
      
      immune_type <- case_when(
        neutro > neutrophil_median & lympho <= lymphocyte_median ~ "Innate-dominant",
        lympho > lymphocyte_median & neutro <= neutrophil_median ~ "Adaptive-dominant",
        TRUE ~ "Mixed"
      )
      
      # Risk percentile
      input_data <- patient_df()
      logfc <- as.numeric(input_data[1, ])
      names(logfc) <- colnames(input_data)
      # Fill in only the 50 selected genes, set others to 0
      full_row <- setNames(rep(0, length(all_genes)), all_genes)
      for (gene in names(logfc)) {
        if (gene %in% names(full_row)) {
          full_row[gene] <- logfc[gene]
        }
      }
      
      X_patient <- matrix(full_row, nrow = 1)
      colnames(X_patient) <- names(full_row)
      log_risk <- predict(cox_model, newx = X_patient, s = "lambda.min", type = "link")[1]
      
      risk_values <- as.numeric(risk_scores)
      percentile <- ecdf(risk_values)(log_risk) * 100
      
      
      risk_group <- case_when(
        percentile >= 75 ~ "High",
        percentile >= 25 ~ "Moderate",
        TRUE ~ "Low"
      )
      
      group <- paste(immune_type, risk_group, sep = " - ")
      
      advice_msg <- case_when(
        immune_type == "Innate-dominant" & risk_group == "High" ~ 
          "Enhanced monitoring: Regular disease activity assessment such as SLEDAI score to detect changes promptly.<br>
Treatment adjustment: Consider immunomodulatory agents like anti-TNF-α biologics to control excessive innate immune response.<br>
Lifestyle management: Avoid known triggers such as infections and UV exposure, and maintain healthy habits.<br>
Reference: Munroe et al. (2024). <a href='https://www.sciencedirect.com/science/article/pii/S0003496724155392' target='_blank'>Link</a>",
        
        immune_type == "Mixed" & risk_group == "High" ~ 
          "Comprehensive treatment: Combine immunosuppressants and biologics to regulate both innate and adaptive immune responses.<br>
Regular follow-up: Closely monitor disease activity and organ function, adjusting treatment strategies as needed.<br>
Patient education: Enhance understanding of the disease to improve adherence.<br>
Medical Reference: Smith et al. (2021). <a href='https://www.sciencedirect.com/science/article/pii/S0896841121000238' target='_blank'>Link</a>",
        
        immune_type == "Adaptive-dominant" & risk_group == "High" ~ 
          "Targeted therapy: Use B-cell targeted treatments such as anti-CD20 monoclonal antibodies to reduce autoantibody production.<br>
Immune monitoring: Regular testing of autoantibody and complement levels to assess disease activity.<br>
Psychological support: Provide counseling to help patients cope with long-term disease stress.<br>
Reference: Smith et al. (2021). <a href='https://www.sciencedirect.com/science/article/pii/S0896841121000238' target='_blank'>Link</a>",
        
        immune_type == "Innate-dominant" & risk_group == "Low" ~ 
          "Maintenance treatment: Continue current treatment plans, avoiding arbitrary discontinuation.<br>
Periodic check-ups: Conduct comprehensive assessments every 3-6 months to ensure disease stability.<br>
Healthy lifestyle: Maintain balanced diet and moderate exercise to enhance overall health.<br>
Reference: Jones et al. (2021). <a href='https://advancesinrheumatology.biomedcentral.com/articles/10.1186/s42358-021-00202-7' target='_blank'>Link</a>",
        
        immune_type == "Mixed" & risk_group == "Low" ~ 
          "Continued observation: Regularly assess immune function to maintain a balanced immune state.<br>
Preventive measures: Vaccinations to prevent infections and minimize flare triggers.<br>
Lifestyle management: Maintain regular routines and a positive mental state to stabilize the disease.<br>
Medical Reference: Jones et al. (2021). <a href='https://advancesinrheumatology.biomedcentral.com/articles/10.1186/s42358-021-00202-7' target='_blank'>Link</a>",
        
        immune_type == "Adaptive-dominant" & risk_group == "Low" ~ 
          "Simplified treatment: Gradually reduce immunosuppressive agents under medical supervision to avoid overtreatment.<br>
Long-term monitoring: Comprehensive assessments every 6-12 months to ensure sustained remission.<br>
Health promotion: Encourage participation in moderate physical activity to boost physical fitness.<br>
Reference: Jones et al. (2021). <a href='https://advancesinrheumatology.biomedcentral.com/articles/10.1186/s42358-021-00202-7' target='_blank'>Link</a>",
        
        TRUE ~ 
          "Your immune-risk profile is uncommon in our dataset. Please consider personalized consultation for further guidance."
      )
      
    
      
      # Output UI
      HTML(paste0(
        "<div style='font-size: 16px;'>",
        "<b>You are classified as:</b> <span style='font-size: 22px; color: #004085;'>", group, "</span><br><br>",
        "This combination is associated with different probabilities of remaining flare-free over time, as shown in the curve above.<br>",
        "The curve represents each subtype-risk combination’s probability of staying flare-free since first visit.<br><br>",
        "<b>Personalized Advice:</b><br>",
        advice_msg,
        "</div>"
      ))
      
    
    })
  })
  
  #  Kaplan–Meier survival chart
  output$survival_plot <- renderPlot({
    ggsurvplot(
      km_fit_combo,
      data = df_test2,
      palette = "Dark2",
      title = "Kaplan-Meier Curves by Risk and Immune Type",
      xlab = "Time (days)", ylab = "Flare-free Survival Probability"
    )$plot
  })
  
  output$survival_note <- renderUI({
    HTML("<div style='font-size: 14px; color: #444; padding-top: 10px;'>
    <b>Note:</b> The survival curve shows how likely each immune-risk group is to remain flare-free over time.
    Steeper drops (e.g., <span style='color: purple; font-weight: bold;'>purple</span> or <span style='color: green; font-weight: bold;'>green</span> lines) indicate higher flare risk.
    Flatter lines (e.g., <span style='color: orange; font-weight: bold;'>orange</span>) suggest better long-term flare-free survival.
  </div>")
  })
    
}

shinyApp(ui, server)
