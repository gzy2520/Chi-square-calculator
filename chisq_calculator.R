#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
#' ==============================================================================
#' 遗传学卡方值计算器 (Genetics Chi-Square Goodness-of-Fit Calculator)
#' 适用于：孟德尔分离比验证 (3:1, 9:3:3:1, 1:1 等)、适合度检验、独立性检验
#' 随机数种子标准标识：25
#' ==============================================================================

suppressPackageStartupMessages({
  # Look for core library in local directory
  script_dir <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg)))
    } else {
      getwd()
    }
  }, error = function(e) getwd())
})

core_path <- file.path(script_dir, "R", "core_chisq.R")
if (!file.exists(core_path)) {
  core_path <- "R/core_chisq.R"
}
source(core_path)

# Helper function to parse comma or space separated numbers
parse_numeric_input <- function(str_input) {
  if (is.null(str_input) || length(str_input) == 0 || nchar(trimws(str_input)) == 0) {
    return(NULL)
  }
  # Replace commas, semicolons, tabs, and multiple spaces with a single space
  clean_str <- gsub("[,;，；\t]+", " ", str_input)
  clean_str <- trimws(clean_str)
  parts <- unlist(strsplit(clean_str, "\\s+"))
  num_vals <- suppressWarnings(as.numeric(parts))
  if (any(is.na(num_vals))) {
    stop(sprintf("输入格式错误：无法将 '%s' 解析为有效数字序列。", str_input))
  }
  num_vals
}

print_report <- function(res) {
  cat("\n")
  cat(paste0(rep("=", 68), collapse = ""), "\n")
  cat("                🧬 遗传学卡方适合度检验分析报告 🧬\n")
  cat(paste0(rep("=", 68), collapse = ""), "\n")

  cat("\n【1. 检验基础参数】\n")
  cat(sprintf("  • 表型类别数 (k)       : %d\n", length(res$observed)))
  cat(sprintf("  • 总观察样本量 (N)     : %g\n", sum(res$observed)))
  cat(sprintf("  • 检验自由度 (df)      : %d\n", res$df))
  if (res$yates_applied) {
    cat("  • 连续性校正           : 已应用 Yates 连续性校正 (|O - E| - 0.5)\n")
  }

  cat("\n【2. 卡方值与临界表值 (关键输出)】\n")
  cat(sprintf("  ★ 实际计算卡方值 (χ²_calc) : %-10.4f (P-value = %s)\n",
              res$chisq_calc,
              if (res$p_value < 1e-4) sprintf("%.4e", res$p_value) else sprintf("%.4f", res$p_value)))
  cat(sprintf("  • p = 0.05 对应临界表值 (χ²_0.05) : %.4f\n", res$crit_005))
  cat(sprintf("  • p = 0.01 对应临界表值 (χ²_0.01) : %.4f\n", res$crit_001))

  cat("\n【3. 遗传学显著性判定结论】\n")
  badge <- switch(res$significance_level,
    "not_significant" = "[ 不显著 / 符合理论比率 ] (ns, P > 0.05)",
    "significant"     = "[ 差异显著 * ] (0.01 < P ≤ 0.05)",
    "highly_significant" = "[ 差异极显著 ** ] (P ≤ 0.01)"
  )
  cat(sprintf("  ▶ 判定等级 : %s\n", badge))
  cat(sprintf("  ▶ 详细结论 : %s\n", res$conclusion_detail))

  cat("\n【4. 各类别偏离与卡方贡献明细表】\n")
  df_table <- res$breakdown
  colnames(df_table) <- c("类别", "观测值(O)", "期望值(E)", "偏差(O-E)", "(O-E)²", "贡献项(χ²_i)", "占比(%)")
  print(df_table, row.names = FALSE)

  if (length(res$warnings) > 0) {
    cat("\n【5. 统计注意事项】\n")
    for (w in res$warnings) {
      cat(sprintf("  ⚠️  %s\n", w))
    }
  }

  cat(paste0(rep("=", 68), collapse = ""), "\n\n")
}

# CLI Argument parsing
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  params <- list(obs = NULL, exp = NULL, ratio = NULL, df = NULL, yates = FALSE, chisq = NULL)
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (arg %in% c("--obs", "-o")) {
      params$obs <- args[i + 1]; i <- i + 2
    } else if (arg %in% c("--exp", "-e")) {
      params$exp <- args[i + 1]; i <- i + 2
    } else if (arg %in% c("--ratio", "-r")) {
      params$ratio <- args[i + 1]; i <- i + 2
    } else if (arg %in% c("--df", "-d")) {
      params$df <- as.integer(args[i + 1]); i <- i + 2
    } else if (arg %in% c("--chisq", "-c")) {
      params$chisq <- as.numeric(args[i + 1]); i <- i + 2
    } else if (arg %in% c("--yates", "-y")) {
      params$yates <- TRUE; i <- i + 1
    } else if (arg %in% c("--help", "-h")) {
      cat("使用说明 (Usage):\n")
      cat("  交互模式: 直接运行 Rscript chisq_calculator.R\n")
      cat("  命令行参数模式:\n")
      cat("    --obs   <数值列表>   观测值 (如: \"315,108,101,32\" 或 \"85 75\")\n")
      cat("    --exp   <数值列表>   期望值 (如: \"312.75,104.25,104.25,34.75\")\n")
      cat("    --ratio <比率列表>   理论分离比 (如: \"9,3,3,1\" 或 \"3,1\") [自动根据总数算期望]\n")
      cat("    --df    <整数>       自由度 (可选，默认根据分类数自动设为 k - 1)\n")
      cat("    --yates              当 df=1 时开启 Yates 连续性校正\n")
      cat("    --chisq <数值>       直接输入已算出的卡方值进行查表与检验判定\n\n")
      quit(status = 0)
    } else {
      i <- i + 1
    }
  }
  params
}

run_calculator <- function() {
  cli_params <- parse_args(args)

  # Check if direct chisq lookup is requested
  if (!is.null(cli_params$chisq)) {
    chisq_val <- cli_params$chisq
    df_val <- if (!is.null(cli_params$df)) cli_params$df else 1
    crit_005 <- qchisq(0.05, df = df_val, lower.tail = FALSE)
    crit_001 <- qchisq(0.01, df = df_val, lower.tail = FALSE)
    pval <- pchisq(chisq_val, df = df_val, lower.tail = FALSE)

    cat("\n=== 卡方分布查表与显著性检验 (直接卡方值模式) ===\n")
    cat(sprintf("输入卡方值 (χ²) : %.4f\n", chisq_val))
    cat(sprintf("自由度 (df)     : %d\n", df_val))
    cat(sprintf("p=0.05 临界表值 : %.4f\n", crit_005))
    cat(sprintf("p=0.01 临界表值 : %.4f\n", crit_001))
    cat(sprintf("精确 P 值       : %s\n", if (pval < 1e-4) sprintf("%.4e", pval) else sprintf("%.4f", pval)))

    if (chisq_val < crit_005) {
      cat("结论: 差异不显著 (P > 0.05)，接受原假设 H₀。\n\n")
    } else if (chisq_val < crit_001) {
      cat("结论: 差异显著 (0.01 < P ≤ 0.05, *)，拒绝原假设 H₀。\n\n")
    } else {
      cat("结论: 差异极显著 (P ≤ 0.01, **)，极显著拒绝原假设 H₀。\n\n")
    }
    return(invisible(NULL))
  }

  obs_input <- cli_params$obs
  exp_input <- cli_params$exp
  ratio_input <- cli_params$ratio
  df_input <- cli_params$df
  yates_input <- cli_params$yates

  # If not provided via CLI, interactive prompt
  if (is.null(obs_input)) {
    cat("\n======================================================\n")
    cat("        欢迎使用 遗传学卡方值计算器 (Chi-Square)      \n")
    cat("======================================================\n")
    cat("提示：多个数值以逗号、空格或分号分隔即可。\n\n")

    cat("1. 请输入观测值 (Observed)，例如 '315, 108, 101, 32' 或 '85, 75':\n> ")
    obs_raw <- readLines(n = 1)
    obs_vals <- parse_numeric_input(obs_raw)
    if (is.null(obs_vals) || length(obs_vals) == 0) {
      cat("未检测到有效观测值，程序退出。\n")
      return(invisible(NULL))
    }

    cat("\n2. 请选择期望值的输入方式：\n")
    cat("   [A] 直接输入各组期望值 (Expected，例如 '312.75, 104.25, 104.25, 34.75')\n")
    cat("   [B] 输入理论遗传分离比例 (Ratio，例如 '9:3:3:1' 或 '3:1'，系统自动换算期望值)\n")
    cat("请选择 A 或 B (直接回车默认 A): ")
    choice <- toupper(trimws(readLines(n = 1)))

    exp_vals <- NULL
    ratio_vals <- NULL

    if (choice == "B") {
      cat("请输入理论分离比 (如 '9:3:3:1' 或 '3, 1'):\n> ")
      ratio_raw <- readLines(n = 1)
      ratio_clean <- gsub(":", " ", ratio_raw)
      ratio_vals <- parse_numeric_input(ratio_clean)
    } else {
      cat("请输入期望值 (Expected，如与理论比率相同亦可直接换算):\n> ")
      exp_raw <- readLines(n = 1)
      exp_vals <- parse_numeric_input(exp_raw)
    }

    cat(sprintf("\n3. 请输入自由度 df (直接按回车将自动使用默认值 df = %d):\n> ", length(obs_vals) - 1))
    df_raw <- readLines(n = 1)
    df_val <- if (nchar(trimws(df_raw)) > 0) as.integer(trimws(df_raw)) else NULL

    yates_val <- FALSE
    if ((is.null(df_val) && length(obs_vals) == 2) || (!is.null(df_val) && df_val == 1)) {
      cat("\n检测到自由度 df = 1，是否启用 Yates 连续性校正？(y/N，直接回车默认不开启): ")
      yates_raw <- tolower(trimws(readLines(n = 1)))
      if (yates_raw %in% c("y", "yes", "true", "1")) {
        yates_val <- TRUE
      }
    }

    res <- calculate_genetics_chisq(
      observed = obs_vals,
      expected = exp_vals,
      ratio = ratio_vals,
      df = df_val,
      yates = yates_val
    )
    print_report(res)

  } else {
    # Non-interactive CLI
    obs_vals <- parse_numeric_input(obs_input)
    exp_vals <- parse_numeric_input(exp_input)
    ratio_vals <- if (!is.null(ratio_input)) parse_numeric_input(gsub(":", " ", ratio_input)) else NULL

    res <- calculate_genetics_chisq(
      observed = obs_vals,
      expected = exp_vals,
      ratio = ratio_vals,
      df = df_input,
      yates = yates_input
    )
    print_report(res)
  }
}

if (!interactive()) {
  run_calculator()
}
