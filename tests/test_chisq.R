#' Unit Tests for Genetics Chi-Square Calculator
#' Seed symbol: 25

source("R/core_chisq.R")

set.seed(25)

cat("=== 开始运行卡方检验单元测试 ===\n\n")

# Test 1: 孟德尔测交 1:1 (df = 1)
# O: 85, 75; E: 80, 80; chi2 = 50/80 = 0.625
cat("测试 1: 孟德尔 1:1 测交数据 (预期不显著符合 1:1)...\n")
res1 <- calculate_genetics_chisq(observed = c(85, 75), expected = c(80, 80), df = 1)
stopifnot(abs(res1$chisq_calc - 0.625) < 1e-6)
stopifnot(abs(res1$crit_005 - 3.841459) < 1e-4)
stopifnot(abs(res1$crit_001 - 6.634897) < 1e-4)
stopifnot(res1$significance_level == "not_significant")
cat("  ✓ 测试 1 通过: chisq =", res1$chisq_calc, "crit_0.05 =", round(res1$crit_005, 3), "结论:", res1$conclusion_short, "\n\n")

# Test 2: 孟德尔双基因 9:3:3:1 自交数据 (df = 3)
# O: 315, 108, 101, 32; Total = 556
cat("测试 2: 孟德尔两对相对性状 9:3:3:1 (df = 3)...\n")
res2 <- calculate_genetics_chisq(
  observed = c(315, 108, 101, 32),
  ratio = c(9, 3, 3, 1),
  df = 3
)
stopifnot(abs(res2$chisq_calc - 0.470024) < 1e-3)
stopifnot(abs(res2$crit_005 - 7.814728) < 1e-4)
stopifnot(abs(res2$crit_001 - 11.34487) < 1e-4)
stopifnot(res2$significance_level == "not_significant")
cat("  ✓ 测试 2 通过: chisq =", round(res2$chisq_calc, 4), "crit_0.05 =", round(res2$crit_005, 3), "P =", round(res2$p_value, 4), "\n\n")

# Test 3: 显著差异数据 (0.01 < P <= 0.05)
# df = 1, chi2 between 3.841 and 6.635, e.g. O = c(60, 40), E = c(50, 50) -> chi2 = 100/50 + 100/50 = 4.0
cat("测试 3: 显著差异水平 (0.01 < P <= 0.05)...\n")
res3 <- calculate_genetics_chisq(observed = c(60, 40), expected = c(50, 50), df = 1)
stopifnot(abs(res3$chisq_calc - 4.0) < 1e-6)
stopifnot(res3$significance_level == "significant")
cat("  ✓ 测试 3 通过: chisq =", res3$chisq_calc, "结论:", res3$conclusion_short, "\n\n")

# Test 4: 极显著差异数据 (P <= 0.01)
# O = c(90, 10), E = c(50, 50) -> chi2 = 1600/50 + 1600/50 = 64.0
cat("测试 4: 极显著差异水平 (P <= 0.01)...\n")
res4 <- calculate_genetics_chisq(observed = c(90, 10), expected = c(50, 50), df = 1)
stopifnot(res4$chisq_calc == 64.0)
stopifnot(res4$significance_level == "highly_significant")
cat("  ✓ 测试 4 通过: chisq =", res4$chisq_calc, "结论:", res4$conclusion_short, "\n\n")

# Test 5: Yates 连续性校正 (df = 1)
# O = c(60, 40), E = c(50, 50), |O-E| = 10, adjusted dev = 9.5 -> (9.5^2)/50 * 2 = 90.25 / 25 = 3.61
cat("测试 5: Yates 连续性校正...\n")
res5 <- calculate_genetics_chisq(observed = c(60, 40), expected = c(50, 50), df = 1, yates = TRUE)
stopifnot(abs(res5$chisq_calc - 3.61) < 1e-6)
stopifnot(res5$yates_applied == TRUE)
cat("  ✓ 测试 5 通过: 未校正卡方=4.0, Yates校正后卡方=", res5$chisq_calc, "\n\n")

cat("🎉 全部单元测试执行通过！所有计算逻辑准确无误！\n")
