library(dplyr)
library(ggplot2)

# 対象データフレームをリスト化
survey_list <- list(survey1 = survey1, survey2 = survey2, survey3 = survey3, survey4 = survey4)

# 結果格納用リスト
results <- list()

# 各データフレームに対して処理
for (name in names(survey_list)) {
  df <- survey_list[[name]]
  duration <- df$Duration..in.seconds.
  
  # 箱ひげ図の表示
  boxplot(duration, main = paste(name, "- Duration (in seconds)"), ylab = "Seconds")
  
  # 四分位数とIQRの計算
  Q1 <- quantile(duration, 0.25, na.rm = TRUE)
  Q3 <- quantile(duration, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  # 外れ値の数と割合
  n_total <- sum(!is.na(duration))
  n_outliers <- sum(duration > upper_bound, na.rm = TRUE)
  pct_outliers <- round(100 * n_outliers / n_total, 2)
  
  # 結果を保存
  results[[name]] <- list(
    Q1 = Q1,
    Q3 = Q3,
    IQR = IQR,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    n_total = n_total,
    n_outliers = n_outliers,
    pct_outliers = pct_outliers
  )
}

# 結果の出力
for (name in names(results)) {
  cat("\n=============================\n")
  cat(" ", name, "のDuration分析\n")
  cat("=============================\n")
  print(results[[name]])
}
