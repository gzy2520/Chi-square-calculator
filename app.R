#' ==============================================================================
#' 遗传学卡方值计算器 - R Shiny 交互式 Web 应用
#' 适用于遗传学实验、孟德尔分离比验证、适合度检验
#' ==============================================================================

library(shiny)
library(bslib)

source("R/core_chisq.R")

ui <- page_navbar(
  title = "🧬 遗传学卡方值计算器 (Chi-Square Calculator)",
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#1E88E5",
    secondary = "#6c757d",
    success = "#2E7D32",
    warning = "#F57C00",
    danger = "#D32F2F"
  ),
  nav_panel(
    title = "🧮 适合度检验计算器",
    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        title = "输入检验数据",

        # Preset selection
        selectInput(
          "preset",
          "遗传学经典比例预设 (可选快捷填入):",
          choices = c(
            "自定义输入" = "custom",
            "孟德尔自交 3:1 (单对等位基因)" = "m_3_1",
            "孟德尔测交 1:1 (等位基因测交)" = "m_1_1",
            "不完全显性 1:2:1" = "m_1_2_1",
            "孟德尔两对自交 9:3:3:1 (双杂合子自交)" = "m_9_3_3_1",
            "双隐性测交 1:1:1:1" = "m_1_1_1_1"
          ),
          selected = "custom"
        ),

        textAreaInput(
          "obs_input",
          "1. 观测值 (Observed Values):",
          value = "315, 108, 101, 32",
          rows = 2,
          placeholder = "用逗号或空格分隔，如: 315, 108, 101, 32"
        ),

        radioButtons(
          "exp_mode",
          "2. 期望值输入模式:",
          choices = c("直接输入期望数值" = "direct", "输入理论分离比例 (自动换算)" = "ratio"),
          selected = "direct"
        ),

        conditionalPanel(
          condition = "input.exp_mode == 'direct'",
          textAreaInput(
            "exp_input",
            "期望值 (Expected Values):",
            value = "312.75, 104.25, 104.25, 34.75",
            rows = 2,
            placeholder = "与观测值等长，如: 312.75, 104.25, 104.25, 34.75"
          )
        ),

        conditionalPanel(
          condition = "input.exp_mode == 'ratio'",
          textInput(
            "ratio_input",
            "理论遗传比例 (Ratio):",
            value = "9:3:3:1",
            placeholder = "如 9:3:3:1 或 3:1 或 1:1"
          )
        ),

        checkboxInput("auto_df", "自动根据表型数计算自由度 (df = k - 1)", value = TRUE),

        conditionalPanel(
          condition = "!input.auto_df",
          numericInput("custom_df", "指定自由度 (Degrees of Freedom):", value = 3, min = 1, step = 1)
        ),

        conditionalPanel(
          condition = "(input.auto_df && input.obs_input.split(/[,\\s]+/).filter(Boolean).length == 2) || (!input.auto_df && input.custom_df == 1)",
          checkboxInput("yates_corr", "启用 Yates 连续性校正 (推荐在 df = 1 时选用)", value = FALSE)
        ),

        hr(),
        actionButton("btn_calc", "立即重新计算", class = "btn-primary w-100", icon = icon("calculator"))
      ),

      # Main Content Area
      tagList(
        # Summary Row: 3 metric cards
        layout_column_wrap(
          width = 1/3,
          card(
            card_header("实际计算卡方值 (χ²_calc)"),
            card_body(
              uiOutput("metric_chisq")
            )
          ),
          card(
            card_header("p = 0.05 临界表值 (χ²_0.05)"),
            card_body(
              uiOutput("metric_crit005")
            )
          ),
          card(
            card_header("p = 0.01 临界表值 (χ²_0.01)"),
            card_body(
              uiOutput("metric_crit001")
            )
          )
        ),

        # Conclusion Banner
        uiOutput("conclusion_banner"),

        # Visuals & Tables
        navset_card_tab(
          nav_panel(
            title = "📈 卡方分布曲线与拒绝域",
            plotOutput("plot_distribution", height = "360px")
          ),
          nav_panel(
            title = "📊 观测值 vs 期望值对比",
            plotOutput("plot_comparison", height = "360px")
          ),
          nav_panel(
            title = "📋 各组分项明细表",
            tableOutput("table_breakdown")
          )
        )
      )
    )
  ),

  nav_panel(
    title = "📖 遗传学卡方速查与说明",
    card(
      card_header("遗传学卡方适合度检验知识速查"),
      card_body(
        markdown("
### 一、适合度检验原理 (Goodness-of-Fit Test)
在遗传学中，卡方适合度检验用于检验杂交试验后代的表型分离比是否符合特定的孟德尔分离假说（如 3:1、1:1、9:3:3:1 等）。

- **原假设 $H_0$**：观测值与理论预期值无显著差异，符合预期的遗传比率。
- **备择假设 $H_1$**：观测值与理论预期值有显著差异，不符合预期的遗传比率。

### 二、计算公式
$$\\chi^2 = \\sum_{i=1}^k \\frac{(O_i - E_i)^2}{E_i}$$
- $O_i$：第 $i$ 类表型的实际观测数 (Observed)
- $E_i$：第 $i$ 类表型的理论期望数 (Expected)
- $k$：表型类别数
- 自由度 $df = k - 1$ (在仅有总数限制的简单适合度检验中)

### 三、显著性判定准则
| 检验结果 | 概率范围 | 判定结论 | 遗传学解释 |
| :--- | :--- | :--- | :--- |
| $\\chi^2 < \\chi^2_{0.05}$ | $P > 0.05$ | **不显著 (ns)** | 接受 $H_0$：观测结果符合理论遗传规律 |
| $\\chi^2_{0.05} \\le \\chi^2 < \\chi^2_{0.01}$ | $0.01 < P \\le 0.05$ | **差异显著 (*)** | 拒绝 $H_0$：实际数据与理论比率存在显著差异 |
| $\\chi^2 \\ge \\chi^2_{0.01}$ | $P \\le 0.01$ | **差异极显著 (**)** | 极显著拒绝 $H_0$：严重偏离理论比率 (可能存在致死、连锁等) |

### 四、Yates 连续性校正说明
当自由度 $df = 1$ (即只有两类表型) 且样本量较小时，离散的二项分布近似为连续的卡方分布时可能偏高，可采用 Yates 连续性校正：
$$\\chi^2_c = \\sum_{i=1}^2 \\frac{(|O_i - E_i| - 0.5)^2}{E_i}$$
        ")
      )
    )
  )
)

server <- function(input, output, session) {

  # Handle presets
  observeEvent(input$preset, {
    if (input$preset == "m_3_1") {
      updateTextAreaInput(session, "obs_input", value = "305, 95")
      updateRadioButtons(session, "exp_mode", selected = "ratio")
      updateTextInput(session, "ratio_input", value = "3:1")
    } else if (input$preset == "m_1_1") {
      updateTextAreaInput(session, "obs_input", value = "85, 75")
      updateRadioButtons(session, "exp_mode", selected = "ratio")
      updateTextInput(session, "ratio_input", value = "1:1")
    } else if (input$preset == "m_1_2_1") {
      updateTextAreaInput(session, "obs_input", value = "52, 98, 50")
      updateRadioButtons(session, "exp_mode", selected = "ratio")
      updateTextInput(session, "ratio_input", value = "1:2:1")
    } else if (input$preset == "m_9_3_3_1") {
      updateTextAreaInput(session, "obs_input", value = "315, 108, 101, 32")
      updateRadioButtons(session, "exp_mode", selected = "ratio")
      updateTextInput(session, "ratio_input", value = "9:3:3:1")
    } else if (input$preset == "m_1_1_1_1") {
      updateTextAreaInput(session, "obs_input", value = "48, 52, 49, 51")
      updateRadioButtons(session, "exp_mode", selected = "ratio")
      updateTextInput(session, "ratio_input", value = "1:1:1:1")
    }
  })

  # Parse helpers
  parse_vector <- function(txt) {
    if (is.null(txt) || nchar(trimws(txt)) == 0) return(NULL)
    clean <- gsub("[,;，；\t:]+", " ", txt)
    vals <- suppressWarnings(as.numeric(unlist(strsplit(trimws(clean), "\\s+"))))
    vals[!is.na(vals)]
  }

  calc_results <- reactive({
    # Dependency on calculate button or inputs
    input$btn_calc

    isolate({
      obs <- parse_vector(input$obs_input)
      if (is.null(obs) || length(obs) == 0) return(NULL)

      exp_val <- NULL
      ratio_val <- NULL

      if (input$exp_mode == "direct") {
        exp_val <- parse_vector(input$exp_input)
      } else {
        ratio_val <- parse_vector(input$ratio_input)
      }

      df_val <- if (input$auto_df) NULL else input$custom_df
      yates_val <- isTRUE(input$yates_corr)

      tryCatch({
        calculate_genetics_chisq(
          observed = obs,
          expected = exp_val,
          ratio = ratio_val,
          df = df_val,
          yates = yates_val
        )
      }, error = function(e) {
        list(error = e$message)
      })
    })
  })

  # Metric: chisq calc
  output$metric_chisq <- renderUI({
    res <- calc_results()
    if (is.null(res) || !is.null(res$error)) {
      return(h3("—", class = "text-muted"))
    }
    p_str <- if (res$p_value < 1e-4) sprintf("%.2e", res$p_value) else sprintf("%.4f", res$p_value)
    tagList(
      h2(sprintf("%.4f", res$chisq_calc), class = "text-primary fw-bold mb-0"),
      p(class = "text-muted mt-1 mb-0", sprintf("P-value = %s (df = %d)", p_str, res$df))
    )
  })

  # Metric: 0.05
  output$metric_crit005 <- renderUI({
    res <- calc_results()
    if (is.null(res) || !is.null(res$error)) return(h3("—", class = "text-muted"))
    tagList(
      h2(sprintf("%.4f", res$crit_005), class = "text-info fw-bold mb-0"),
      p(class = "text-muted mt-1 mb-0", "α = 0.05 临界值 (5% 界限)")
    )
  })

  # Metric: 0.01
  output$metric_crit001 <- renderUI({
    res <- calc_results()
    if (is.null(res) || !is.null(res$error)) return(h3("—", class = "text-muted"))
    tagList(
      h2(sprintf("%.4f", res$crit_001), class = "text-warning fw-bold mb-0"),
      p(class = "text-muted mt-1 mb-0", "α = 0.01 临界值 (1% 界限)")
    )
  })

  # Conclusion banner
  output$conclusion_banner <- renderUI({
    res <- calc_results()
    if (is.null(res)) return(NULL)
    if (!is.null(res$error)) {
      return(div(class = "alert alert-danger", icon("triangle-exclamation"), " 错误：", res$error))
    }

    alert_class <- switch(res$significance_level,
      "not_significant" = "alert-success",
      "significant"     = "alert-warning",
      "highly_significant" = "alert-danger"
    )

    icon_name <- switch(res$significance_level,
      "not_significant" = "circle-check",
      "significant"     = "triangle-exclamation",
      "highly_significant" = "circle-xmark"
    )

    div(
      class = paste("alert", alert_class, "p-3 my-3 shadow-sm"),
      h4(class = "alert-heading fw-bold mb-2",
         icon(icon_name), " 检验判定：", res$conclusion_short),
      p(class = "mb-1", res$conclusion_detail),
      if (length(res$warnings) > 0) {
        div(class = "mt-2 pt-2 border-top border-secondary-subtle small",
            lapply(res$warnings, function(w) div(icon("circle-info"), " ", w)))
      }
    )
  })

  # Breakdown table
  output$table_breakdown <- renderTable({
    res <- calc_results()
    if (is.null(res) || !is.null(res$error)) return(NULL)
    tb <- res$breakdown
    colnames(tb) <- c("类别 (Category)", "观测数 (O)", "期望数 (E)", "偏差 (O-E)", "偏差平方 (O-E)²", "卡方分项 (χ²_i)", "贡献占比 (%)")
    tb
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # Distribution Plot
  output$plot_distribution <- renderPlot({
    res <- calc_results()
    if (is.null(res) || !is.null(res$error)) return(NULL)

    df <- res$df
    x_max <- max(res$crit_001 * 1.5, res$chisq_calc * 1.25, 15)
    x_vals <- seq(0.01, x_max, length.out = 500)
    y_vals <- dchisq(x_vals, df = df)

    # Plot setup
    par(mar = c(4.5, 4.5, 3, 2), bg = "#FAFAFA")
    plot(x_vals, y_vals, type = "l", lwd = 2.5, col = "#37474F",
         xlab = "卡方统计量 (χ²)", ylab = "概率密度 (Density)",
         main = sprintf("卡方分布概率密度曲线 (df = %d)", df),
         xlim = c(0, x_max), ylim = c(0, max(y_vals) * 1.15))
    grid(col = "#E0E0E0")

    # Shading 0.05 region
    x_05 <- seq(res$crit_005, x_max, length.out = 200)
    y_05 <- dchisq(x_05, df = df)
    polygon(c(res$crit_005, x_05, x_max), c(0, y_05, 0), col = rgb(1, 0.6, 0.2, 0.35), border = NA)

    # Shading 0.01 region
    x_01 <- seq(res$crit_001, x_max, length.out = 200)
    y_01 <- dchisq(x_01, df = df)
    polygon(c(res$crit_001, x_01, x_max), c(0, y_01, 0), col = rgb(0.9, 0.1, 0.1, 0.45), border = NA)

    # Redraw line
    lines(x_vals, y_vals, lwd = 2.5, col = "#263238")

    # Mark critical values
    abline(v = res$crit_005, col = "#E65100", lty = 2, lwd = 2)
    abline(v = res$crit_001, col = "#B71C1C", lty = 2, lwd = 2)

    # Mark calculated chisq
    abline(v = res$chisq_calc, col = "#00796B", lwd = 3.5)

    legend("topright",
           legend = c(
             sprintf("计算卡方值 χ² = %.4f", res$chisq_calc),
             sprintf("p=0.05 表值 (χ²=%.4f)", res$crit_005),
             sprintf("p=0.01 表值 (χ²=%.4f)", res$crit_001),
             "α=0.05 拒绝域 (P ≤ 0.05)",
             "α=0.01 拒绝域 (P ≤ 0.01)"
           ),
           col = c("#00796B", "#E65100", "#B71C1C", rgb(1, 0.6, 0.2, 0.5), rgb(0.9, 0.1, 0.1, 0.6)),
           lty = c(1, 2, 2, NA, NA),
           pch = c(NA, NA, NA, 15, 15),
           lwd = c(3.5, 2, 2, NA, NA),
           bg = "white", box.col = "#B0BEC5")
  })

  # Comparison Plot
  output$plot_comparison <- renderPlot({
    res <- calc_results()
    if (is.null(res) || !is.null(res$error)) return(NULL)

    par(mar = c(4.5, 4.5, 3, 2), bg = "#FAFAFA")
    mat <- rbind(res$observed, res$expected)
    colnames(mat) <- paste0("类别 ", seq_along(res$observed))
    barplot(mat, beside = TRUE, col = c("#1E88E5", "#78909C"),
            ylab = "频数 (Counts)", xlab = "表型类别",
            main = "观测值 (Observed) 与 理论期望值 (Expected) 对比",
            ylim = c(0, max(mat) * 1.2))
    grid(nx = NA, ny = NULL, col = "#E0E0E0")
    legend("topright", legend = c("观测值 (O)", "期望值 (E)"),
           fill = c("#1E88E5", "#78909C"), bg = "white", box.col = "#B0BEC5")
  })
}

shinyApp(ui = ui, server = server)
