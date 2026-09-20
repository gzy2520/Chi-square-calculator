#' Core Chi-Square Calculation Functions for Genetics
#' Author: Pair-programmed with Antigravity
#' Seed symbol: 25

#' Calculate Genetics Chi-Square Goodness-of-Fit Test
#'
#' @param observed Numeric vector of observed counts
#' @param expected Numeric vector of expected counts (or NULL if ratio is provided)
#' @param ratio Optional numeric vector of theoretical ratios (e.g., c(3, 1) or c(9, 3, 3, 1))
#' @param df Degrees of freedom. If NULL, defaults to length(observed) - 1.
#' @param yates Logical, whether to apply Yates' continuity correction (applicable when df == 1)
#'
#' @return A list containing detailed chi-square test results
calculate_genetics_chisq <- function(observed, expected = NULL, ratio = NULL, df = NULL, yates = FALSE) {
  # Input validation
  if (!is.numeric(observed) || any(is.na(observed)) || any(observed < 0)) {
    stop("观测值 (observed) 必须是非负数值向量。")
  }
  k <- length(observed)
  if (k < 1) {
    stop("观测值向量不能为空。")
  }

  total_obs <- sum(observed)

  # Calculate expected if ratio is supplied
  if (is.null(expected) && !is.null(ratio)) {
    if (!is.numeric(ratio) || length(ratio) != k || any(ratio <= 0)) {
      stop("理论比例 (ratio) 的长度必须与观测值相同且大于0。")
    }
    expected <- total_obs * (ratio / sum(ratio))
  }

  if (is.null(expected)) {
    stop("必须提供期望值 (expected) 或理论比例 (ratio)。")
  }

  if (!is.numeric(expected) || length(expected) != k || any(expected <= 0)) {
    stop("期望值 (expected) 必须是与观测值等长且大于0的数值向量。")
  }

  # Degrees of freedom
  if (is.null(df)) {
    df <- max(1, k - 1)
  } else {
    df <- as.integer(df)
    if (df < 1) {
      stop("自由度 (df) 必须为正整数 (>= 1)。")
    }
  }

  # Yates' correction check
  apply_yates <- yates && (df == 1) && (k == 2)

  # Component-wise calculations
  deviations <- observed - expected
  abs_dev <- abs(deviations)

  if (apply_yates) {
    # Yates correction: (|O - E| - 0.5)^2 / E
    # If |O - E| < 0.5, term is 0
    adj_dev <- pmax(0, abs_dev - 0.5)
    contributions <- (adj_dev^2) / expected
  } else {
    contributions <- (deviations^2) / expected
  }

  chisq_calc <- sum(contributions)

  # Critical values from Chi-square distribution
  # p=0.05 critical value: qchisq(0.05, df, lower.tail = FALSE) or qchisq(0.95, df)
  crit_005 <- qchisq(0.05, df = df, lower.tail = FALSE)
  crit_001 <- qchisq(0.01, df = df, lower.tail = FALSE)

  # Exact p-value
  p_val <- pchisq(chisq_calc, df = df, lower.tail = FALSE)

  # Significance conclusion
  # H0: 观测值与期望值无显著差异，符合遗传理论预期
  # H1: 观测值与期望值有显著差异，不符合遗传理论预期
  if (chisq_calc < crit_005) {
    significance_level <- "not_significant"
    conclusion_short <- "差异不显著 (P > 0.05)"
    conclusion_detail <- sprintf("计算卡方值 (%.4f) < 临界值 χ²(0.05)=%.4f。差异不显著 (P = %.4f > 0.05)。接受原假设 H₀，观测数据符合理论遗传分离规律/预期值。",
                                 chisq_calc, crit_005, p_val)
    sig_stars <- "ns (不显著)"
  } else if (chisq_calc < crit_001) {
    significance_level <- "significant"
    conclusion_short <- "差异显著 (P ≤ 0.05, 显著*)"
    conclusion_detail <- sprintf("计算卡方值 (%.4f) 介于 χ²(0.05)=%.4f 与 χ²(0.01)=%.4f 之间。差异显著 (P = %.4f ≤ 0.05)。拒绝原假设 H₀，观测数据与理论预期存在显著差异（可能存在分离偏离、活力差异等）。",
                                 chisq_calc, crit_005, crit_001, p_val)
    sig_stars <- "* (显著, P ≤ 0.05)"
  } else {
    significance_level <- "highly_significant"
    conclusion_short <- "差异极显著 (P ≤ 0.01, 极显著**)"
    conclusion_detail <- sprintf("计算卡方值 (%.4f) ≥ 临界值 χ²(0.01)=%.4f。差异极显著 (P = %.4e ≤ 0.01)。极显著拒绝原假设 H₀，观测数据严重偏离理论预期比例（可能存在致死基因、强连锁或外部选择压力）。",
                                 chisq_calc, crit_001, p_val)
    sig_stars <- "** (极显著, P ≤ 0.01)"
  }

  # Build breakdown table
  breakdown <- data.frame(
    Category = seq_len(k),
    Observed = observed,
    Expected = round(expected, 4),
    Deviation = round(deviations, 4),
    Dev_Squared = round(deviations^2, 4),
    Chisq_Contribution = round(contributions, 4),
    Contribution_Pct = round(contributions / max(chisq_calc, 1e-12) * 100, 2),
    stringsAsFactors = FALSE
  )

  # Check minimum expected frequency rule of thumb (Cochran's rule)
  warning_msgs <- character(0)
  if (any(expected < 5)) {
    warning_msgs <- c(warning_msgs, "提示：存在期望频数 < 5 的组别。根据卡方检验常规准则，当期望频数过小时，卡方近似可能不准确，建议合并表型或改用 Fisher 精确检验。")
  }
  if (total_obs < 30) {
    warning_msgs <- c(warning_msgs, "提示：总样本量 < 30，样本量较小，卡方检验统计功效可能受限。")
  }

  list(
    observed = observed,
    expected = expected,
    df = df,
    chisq_calc = chisq_calc,
    crit_005 = crit_005,
    crit_001 = crit_001,
    p_value = p_val,
    significance_level = significance_level,
    conclusion_short = conclusion_short,
    conclusion_detail = conclusion_detail,
    sig_stars = sig_stars,
    yates_applied = apply_yates,
    breakdown = breakdown,
    warnings = warning_msgs
  )
}
