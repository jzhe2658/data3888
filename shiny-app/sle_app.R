library(shiny)
library(bs4Dash)
library(tidyverse)
library(DT)
library(glmnet)
library(ggplot2)
library(scales)
library(survival)
library(survminer)
library(plotly)

risk_scores <- readRDS("risk_scores.rds")
cox_model <- readRDS("lasso_cox_model.rds")

all_genes <- rownames(coef(cox_model, s = "lambda.min"))
all_genes <- all_genes[all_genes != "(Intercept)"]

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
      color = "primary",
      image = "nerd.webp"
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
            title = "Top Gene Contributions",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            plotOutput("gene_plot")
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
            numericInput("neutrophil", "Neutrophil Count", value = 3.5, min = 0, step = 0.01),
            numericInput("lymphocyte", "Lymphocyte Count", value = 1.2, min = 0, step = 0.01),
            actionButton("classify", "Classify Immune Type", class = "btn-primary")
          ),
          box(
            title = "Classification Result",
            status = "success",
            solidHeader = TRUE,
            width = 6,
            uiOutput("immune_result")
          )
        ),
        fluidRow(
          box(
            title = "How Your Group Performed Over Time",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            plotlyOutput("survival_plot", height = "600px")
          )
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
    datatable(patient_df(), options = list(scrollX = TRUE))
  })
  
  output$risk_summary <- renderUI({
    req(patient_df())
    full_row <- setNames(rep(0, length(all_genes)), all_genes)
    input_data <- patient_df()
    for (gene in colnames(input_data)) {
      if (gene %in% names(full_row)) {
        full_row[gene] <- input_data[[gene]][1]
      }
    }
    X_patient <- matrix(full_row, nrow = 1)
    colnames(X_patient) <- names(full_row)
    
    log_risk <- predict(cox_model, newx = X_patient, s = "lambda.min", type = "link")[1]
    percentile <- ecdf(risk_scores)(log_risk) * 100
    
    risk_category <- case_when(
      percentile >= 75 ~ "<span style='color: darkred;'>High Risk</span>",
      percentile >= 25 ~ "<span style='color: #336699;'>Moderate Risk</span>",
      TRUE ~ "<span style='color: darkgreen;'>Low Risk</span>"
    )
    
    recommendation <- if (percentile >= 75) {
      "<br><b style='color:darkred;'>Recommendation:</b> Schedule more frequent consultations or early intervention."
    } else {
      ""
    }
    
    HTML(paste0(
      "<div style='font-size: 18px; line-height: 1.6;'>",
      "<b>Percentile:</b> <span style='font-size: 24px; color: darkorange;'>", round(percentile, 1), "%</span><br>",
      "<b>Interpretation:</b> ", risk_category,
      recommendation,
      "</div>"
    ))
  })
  
  output$gene_plot <- renderPlot({
    req(patient_df())
    input_data <- patient_df()
    full_row <- setNames(rep(0, length(all_genes)), all_genes)
    for (gene in colnames(input_data)) {
      if (gene %in% names(full_row)) {
        full_row[gene] <- input_data[[gene]][1]
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
      full_row <- setNames(rep(0, length(all_genes)), all_genes)
      for (gene in colnames(input_data)) {
        if (gene %in% names(full_row)) {
          full_row[gene] <- input_data[[gene]][1]
        }
      }
      X_patient <- matrix(full_row, nrow = 1)
      colnames(X_patient) <- names(full_row)
      log_risk <- predict(cox_model, newx = X_patient, s = "lambda.min", type = "link")[1]
      percentile <- ecdf(risk_scores)(log_risk) * 100
      
      risk_group <- case_when(
        percentile >= 75 ~ "High",
        percentile >= 25 ~ "Moderate",
        TRUE ~ "Low"
      )
      
      group <- paste(immune_type, risk_group, sep = " - ")
      
      HTML(paste0(
        "<div style='font-size: 16px;'>",
        "<b>You are classified as:</b> <span style='font-size: 22px; color: #004085;'>", group, "</span><br><br>",
        "This combination is associated with different probabilities of remaining flare-free over time, as shown in the curve above.",
        "<br>The curve represents each subtype-risk combination’s probability of staying flare-free since first visit.",
        "</div>"
      ))
    })
  })
}

shinyApp(ui, server)
