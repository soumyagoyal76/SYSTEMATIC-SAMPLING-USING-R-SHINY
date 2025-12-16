# ==============================================================================
# Project: Interactive Stratified Sampling Designer
# Type: R Shiny Application
# Author: SOUMYA GOYAL
# ==============================================================================
# Load packages
library(promises)
library(shiny)
library(ggplot2)
library(reshape2)
# -----------------------------
# Helper functions
# -----------------------------

# Fix rounding so sum(nh) == n and nh >= 1
fix_sum_to_n <- function(nh, n) {
  nh[nh < 1] <- 1
  diff <- n - sum(nh)
  if (diff == 0) return(nh)
  
  idx <- seq_along(nh)
  for (k in seq_len(abs(diff))) {
    i <- idx[(k - 1) %% length(idx) + 1]
    nh[i] <- nh[i] + sign(diff)
    if (nh[i] < 1) nh[i] <- 1
  }
  nh
}

# Variance of stratified mean with FPC:
# Var(y_st) = sum Wh^2 * (1 - nh/Nh) * Sh^2 / nh
var_strat_mean <- function(df, nh) {
  Wh <- df$Wh
  Nh <- df$Nh
  Sh <- df$Sh
  fh <- nh / Nh
  sum(Wh^2 * (1 - fh) * (Sh^2) / nh)
}

# Find minimum n such that Var <= V_target for a given allocation fractions p
find_required_n <- function(df, p, V_target, n_max = 200000) {
  for (n in 2:n_max) {
    nh <- round(n * p)
    nh <- fix_sum_to_n(nh, n)
    v <- var_strat_mean(df, nh)
    if (v <= V_target) {
      return(list(n_total = n, nh = nh, var = v))
    }
  }
  stop("No solution found up to n_max. Relax error/confidence or increase n_max.")
}

# Allocation fractions
p_proportional <- function(df) {
  p <- df$Wh
  p / sum(p)
}
p_neyman <- function(df) {
  w <- df$Nh * df$Sh
  w / sum(w)
}
p_opt_cost <- function(df) {
  w <- (df$Nh * df$Sh) / sqrt(df$Ch)
  w / sum(w)
}
p_opt_time <- function(df) {
  w <- (df$Nh * df$Sh) / sqrt(df$Th)
  w / sum(w)
}

# -----------------------------
# UI
# -----------------------------
ui <- fluidPage(
  titlePanel("Stratified Sampling Design: Sample Size & Allocation"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Design Parameters"),
      
      sliderInput("error", "Allowed sampling bias / Margin of Error (E):",
                  min = 0.1, max = 5.0, value = 1.5, step = 0.1),
      
      selectInput("conf", "Confidence Level (Z):",
                  choices = c("90% (1.645)" = 1.645,
                              "95% (1.96)"  = 1.96,
                              "99% (2.58)"  = 2.58),
                  selected = 1.96),
      
      hr(),
      h4("Strata Data (Demo)"),
      helpText("You can later replace these values with your real strata."),
      tableOutput("strataPreview"),
      hr(),
      checkboxInput("showDetails", "Show achieved SE and Z*SE", value = TRUE)
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Allocation Plot", plotOutput("allocPlot", height = "420px")),
        tabPanel("Detailed Table", tableOutput("summaryTable"))
      ),
      br(),
      wellPanel(
        h4("Required Total Sample Size (n)"),
        verbatimTextOutput("total_n_output")
      )
    )
  )
)

# -----------------------------
# Server
# -----------------------------
server <- function(input, output, session) {
  
  # Fixed demo strata (you can convert to editable inputs later)
  base_df <- reactive({
    df <- data.frame(
      Stratum = factor(c(1, 2, 3, 4)),
      Nh = c(4000, 3000, 2000, 1000),
      Sh = c(10, 20, 30, 40),
      Ch = c(4, 6, 8, 10),
      Th = c(1, 1.5, 2, 2.5)
    )
    N <- sum(df$Nh)
    df$Wh <- df$Nh / N
    df
  })
  
  output$strataPreview <- renderTable({
    base_df()[, c("Stratum","Nh","Sh","Ch","Th")]
  })
  
  sampling_results <- reactive({
    df <- base_df()
    
    E <- input$error
    Z <- as.numeric(input$conf)
    
    # Target variance based on margin of error and confidence:
    # Z * sqrt(Var) <= E  => Var <= (E/Z)^2
    V_target <- (E / Z)^2
    
    # Proportional
    res_prop <- find_required_n(df, p_proportional(df), V_target)
    df$n_Proportional <- res_prop$nh
    
    # Neyman
    res_neym <- find_required_n(df, p_neyman(df), V_target)
    df$n_Neyman <- res_neym$nh
    
    # Optimum (Cost)
    res_cost <- find_required_n(df, p_opt_cost(df), V_target)
    df$n_Opt_Cost <- res_cost$nh
    
    # Optimum (Time)
    res_time <- find_required_n(df, p_opt_time(df), V_target)
    df$n_Opt_Time <- res_time$nh
    
    list(
      df = df,
      V_target = V_target,
      Z = Z,
      E = E,
      prop = res_prop,
      neym = res_neym,
      cost = res_cost,
      time = res_time
    )
  })
  
  output$total_n_output <- renderText({
    res <- sampling_results()
    txt <- c(
      sprintf("Proportional n = %d", res$prop$n_total),
      sprintf("Neyman n       = %d", res$neym$n_total),
      sprintf("Optimum (Cost) n = %d", res$cost$n_total),
      sprintf("Optimum (Time) n = %d", res$time$n_total)
    )
    
    if (isTRUE(input$showDetails)) {
      txt <- c(
        txt,
        "",
        "Achieved precision check (Z * SE ≤ E):",
        sprintf("Proportional: SE=%.6f, Z*SE=%.6f",
                sqrt(res$prop$var), res$Z*sqrt(res$prop$var)),
        sprintf("Neyman      : SE=%.6f, Z*SE=%.6f",
                sqrt(res$neym$var), res$Z*sqrt(res$neym$var)),
        sprintf("Opt Cost    : SE=%.6f, Z*SE=%.6f",
                sqrt(res$cost$var), res$Z*sqrt(res$cost$var)),
        sprintf("Opt Time    : SE=%.6f, Z*SE=%.6f",
                sqrt(res$time$var), res$Z*sqrt(res$time$var)),
        "",
        sprintf("Target: E=%.2f at Z=%.3f (Var target = %.6f)", res$E, res$Z, res$V_target)
      )
    }
    
    paste(txt, collapse = "\n")
  })
  
  output$summaryTable <- renderTable({
    res <- sampling_results()
    res$df[, c("Stratum","Nh","Sh","Ch","Th","n_Proportional","n_Neyman","n_Opt_Cost","n_Opt_Time")]
  })
  
  output$allocPlot <- renderPlot({
    res <- sampling_results()
    df <- res$df
    
    df_long <- melt(df,
                    id.vars = "Stratum",
                    measure.vars = c("n_Proportional","n_Neyman","n_Opt_Cost","n_Opt_Time"),
                    variable.name = "Method",
                    value.name = "SampleSize"
    )
    
    # nicer labels
    df_long$Method <- factor(df_long$Method,
                             levels = c("n_Proportional","n_Neyman","n_Opt_Cost","n_Opt_Time"),
                             labels = c("Proportional","Neyman","Optimum (Cost)","Optimum (Time)"))
    
    ggplot(df_long, aes(x = Stratum, y = SampleSize, fill = Method)) +
      geom_col(position = "dodge") +
      theme_minimal() +
      labs(
        title = "Stratum-wise Allocation Comparison",
        x = "Stratum",
        y = "Sample size in stratum (nh)"
      )
  })
}

# -----------------------------
# Run App
# -----------------------------
shinyApp(ui = ui, server = server)
