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
        width = 390,
        title = "输入检验数据与设置",

        # Preset selection with full genetic ratios
        selectInput(
          "preset",
          "遗传学分离比预设 (或选自由设置):",
          choices = list(
            "自由设置" = c("自由自定义比例 (Custom Ratio)" = "custom"),
            "基础孟德尔分离比" = c(
              "3:1 (单杂合自交)" = "m_3_1",
              "1:1 (等位基因测交)" = "m_1_1",
              "1:2:1 (不完全显性/共显性)" = "m_1_2_1",
              "2:1 (显性纯合致死)" = "m_2_1",
              "9:3:3:1 (双杂合自交自由组合)" = "m_9_3_3_1",
              "1:1:1:1 (双隐性测交)" = "m_1_1_1_1"
            ),
            "两对基因相互作用 (9:3:3:1 变型)" = c(
              "9:7 (互补作用 / 香豌豆花色)" = "m_9_7",
              "9:3:4 (隐性上位 / 家兔毛色)" = "m_9_3_4",
              "12:3:1 (显性上位 / 西葫芦皮色)" = "m_12_3_1",
              "13:3 (显性抑制 / 家鸡羽色)" = "m_13_3",
              "15:1 (重叠基因 / 荠菜果形)" = "m_15_1",
              "9:6:1 (积加作用 / 南瓜果形)" = "m_9_6_1",
              "6:3:2:1 (单对显性致死自交)" = "m_6_3_2_1"
            )
          ),
          selected = "m_9_7"
        ),

        radioButtons(
          "exp_mode",
          "期望值设定模式:",
          choices = c("自由设置理论比例 (Ratio)" = "ratio", "直接输入期望数值 (Counts)" = "direct"),
          selected = "ratio"
        ),

        conditionalPanel(
          condition = "input.exp_mode == 'ratio'",
          div(
            textInput(
              "ratio_input",
              "理论分离比 (支持冒号/逗号/空格分隔):",
              value = "9:7",
              placeholder = "如 9:7, 9:3:4, 13:3, 1:4:6:4:1 等"
            ),
            uiOutput("ratio_hint")
          )
        ),

        conditionalPanel(
          condition = "input.exp_mode == 'direct'",
          textAreaInput(
            "exp_input",
            "期望值 (Expected Values):",
            value = "180, 140",
            rows = 2,
            placeholder = "与观测值等长，如: 180, 140"
          )
        ),

        textAreaInput(
          "obs_input",
          "实际观测值 (Observed Values, O):",
          value = "188, 140",
          rows = 2,
          placeholder = "用逗号或空格分隔，如: 188, 140"
        ),

        checkboxInput("auto_df", "自动计算自由度 (df = 表型数 k - 1)", value = TRUE),

        conditionalPanel(
          condition = "!input.auto_df",
          numericInput("custom_df", "手动指定自由度 (df):", value = 1, min = 1, step = 1)
        ),

        conditionalPanel(
          condition = "(input.auto_df && input.obs_input.split(/[,\\s]+/).filter(Boolean).length == 2) || (!input.auto_df && input.custom_df == 1)",
          checkboxInput("yates_corr", "启用 Yates 连续性校正 (推荐在 df = 1 且两类别时选用)", value = FALSE)
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
    title = "📖 遗传学常见分离比速查表",
    card(
      card_header("遗传学经典基因互作分离比大全 (均源于 9:3:3:1 的演变)"),
      card_body(
        markdown("
### 一、孟德尔经典分离与两对基因上位/互作比率汇总

| 作用类型 | 理论比率 | 自由度 $df$ | 典型生物学实例 | 分子/遗传机理简述 |
| :--- | :--- | :--- | :--- | :--- |
| **单杂合自交** | **3:1** | 1 | 豌豆高茎与矮茎 (Dd × Dd) | 完全显性 |
| **测交分离比** | **1:1** | 1 | 豌豆测交 (Dd × dd) | 等位基因分离 |
| **不完全显性** | **1:2:1** | 2 | 金鱼草花色 (红:粉:白) | 杂合子表型介于纯合子之间 |
| **显性致死** | **2:1** | 1 | 小鼠黄毛基因 ($A^y A \times A^y A$) | 显性纯合致死 ($A^y A^y$ 死亡) |
| **自由组合** | **9:3:3:1** | 3 | 豌豆黄圆:绿圆:黄皱:绿皱 | 两对独立等位基因无互作 |
| **双杂测交** | **1:1:1:1** | 3 | 测交后代表型比例 | 两对等位基因各自独立分离 |
| **互补作用** | **9:7** | 1 | 香豌豆花色 (紫花:白花) | 需两对显性基因同时存在 ($A\\_B\\_$) 才能成色，缺少任一均为白色 ($9 : (3+3+1)$) |
| **隐性上位** | **9:3:4** | 2 | 家兔毛色 (灰:黑:白) | 一对隐性纯合 ($bb$) 对另一对表型有遮盖效应 ($9 : 3 : (3+1)$) |
| **显性上位** | **12:3:1** | 2 | 西葫芦皮色 (白:黄:绿) | 显性基因 ($A\\_$) 遮盖另一对基因的表现 ($(9+3) : 3 : 1$) |
| **显性抑制** | **13:3** | 1 | 家鸡羽色 (白羽:着色羽) | 一个显性基因本身不显色，但能抑制另一显性基因表现 ($(9+3+1) : 3$) |
| **重叠作用** | **15:1** | 1 | 荠菜果形 (三角形:卵圆形) | 只要含任一显性基因即表现同一性状 ($(9+3+3) : 1$) |
| **积加作用** | **9:6:1** | 2 | 南瓜果形 (扁盘:圆球:长圆) | 双显性一种形状，单显性另一种形状，双隐性第三种形状 ($9 : (3+3) : 1$) |
| **单对显性致死**| **6:3:2:1** | 3 | 两对基因其中一对纯合致死 | $(2:1) \\times (3:1) = 6:3:2:1$ |
        ")
      )
    )
  )
)

server <- function(input, output, session) {

  # Handle presets
  observeEvent(input$preset, {
    req(input$preset)
    if (input$preset == "m_9_7") {
      updateTextAreaInput(session, "obs_input", value = "188, 140")
      updateRadioButtons(session, "exp_mode", selected = "ratio")
      updateTextInput(session, "ratio_input", value = "9:7")
    } else if (input$preset == "m_9_3_4") {
      updateTextAreaInput(session, "obs_input", value = "180, 62, 78")
      updateRadioButtons(session, "exp_mode", selected = "ratio")
      updateTextInput(session, "ratio_input", value = "9:3:4")
    } else if (input$preset == "m_12_3_1") {
      updateTextAreaInput(session, "obs_input", value = "241, 58, 21")
      updateRadioButtons(session, "exp_mode", selected = "ratio")
      updateTextInput(session, "ratio_input", value = "12:3:1")
    } else if (input$preset == "m_13_3") {
      updateTextAreaInput(session, "obs_input", value = "262, 58")
      updateRadioButtons(session, "exp_mode", selected = "ratio")
      updateTextInput(session, "ratio_input", value = "13:3")
    } else if (input$preset == "m_15_1") {
      updateTextAreaInput(session, "obs_input", value = "302, 18")
      updateRadioButtons(session, "exp_mode", selected = "ratio")
      updateTextInput(session, "ratio_input", value = "15:1")
    } else if (input$preset == "m_9_6_1") {
      updateTextAreaInput(session, "obs_input", value = "182, 118, 20")
      updateRadioButtons(session, "exp_mode", selected = "ratio")
      updateTextInput(session, "ratio_input", value = "9:6:1")
    } else if (input$preset == "m_6_3_2_1") {
      updateTextAreaInput(session, "obs_input", value = "152, 78, 48, 26")
      updateRadioButtons(session, "exp_mode", selected = "ratio")
      updateTextInput(session, "ratio_input", value = "6:3:2:1")
    } else if (input$preset == "m_3_1") {
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
    } else if (input$preset == "m_2_1") {
      updateTextAreaInput(session, "obs_input", value = "142, 68")
      updateRadioButtons(session, "exp_mode", selected = "ratio")
      updateTextInput(session, "ratio_input", value = "2:1")
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

  output$ratio_hint <- renderUI({
    req(input$ratio_input)
    r <- parse_vector(input$ratio_input)
    if (is.null(r) || length(r) < 2) {
      return(p(class = "text-danger small mt-1", "⚠️ 请输入至少两项比例 (例如 9:7)"))
    }
    pcts <- paste0(round(r / sum(r) * 100, 1), "%", collapse = " : ")
    p(class = "text-info small mt-1",
      sprintf("✓ 识别到 %d 项比例 (%s) | 理论比率: %s",
              length(r), paste(r, collapse = " : "), pcts))
  })

  calc_results <- reactive({
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
      h2(sprintf("%.4f", res$crit_005), class = "text-warning fw-bold mb-0"),
      p(class = "text-muted mt-1 mb-0", "α = 0.05 临界值 (5% 界限)")
    )
  })

  # Metric: 0.01
  output$metric_crit001 <- renderUI({
    res <- calc_results()
    if (is.null(res) || !is.null(res$error)) return(h3("—", class = "text-muted"))
    tagList(
      h2(sprintf("%.4f", res$crit_001), class = "text-danger fw-bold mb-0"),
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

    lines(x_vals, y_vals, lwd = 2.5, col = "#263238")
    abline(v = res$crit_005, col = "#E65100", lty = 2, lwd = 2)
    abline(v = res$crit_001, col = "#B71C1C", lty = 2, lwd = 2)
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
