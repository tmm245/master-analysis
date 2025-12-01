# 必要なパッケージ
library(dplyr)

# 比較パターン定義
patterns <- data.frame(
  X_Var = c("averageNumerical", "averageNumerical", "averageNonNumerical", "averageNonNumerical",
            "averageNumerical", "averageNumerical", "averageNonNumerical", "averageNonNumerical"),
  Y_Var = c("HBM_average_benefit_3", "HBM_average_benefit_3", "HBM_average_benefit_3", "HBM_average_benefit_3",
            "HBM_average_risk_3", "HBM_average_risk_3", "HBM_average_risk_3", "HBM_average_risk_3"),
  stringsAsFactors = FALSE
)

# 相関と Fisher Z を格納するデータフレーム
results <- data.frame()

for (i in seq_len(nrow(patterns))) {
  x_var <- patterns$X_Var[i]
  y_var <- patterns$Y_Var[i]
  
  for (grp in c("BL", "BH")) {
    df_sub <- merged_data_complete_for_path %>%
      filter(berlin_group == grp) %>%
      select(all_of(c(x_var, y_var))) %>%
      na.omit()
    
    r_val <- cor(df_sub[[x_var]], df_sub[[y_var]], method = "pearson")
    n <- nrow(df_sub)
    z <- 0.5 * log((1 + r_val) / (1 - r_val))  # Fisher変換
    
    results <- rbind(results, data.frame(
      Group = grp,
      X_Var = x_var,
      Y_Var = y_var,
      r = r_val,
      n = n,
      z = z
    ))
  }
}

# グループ間の Fisher Z 検定
results$pair_id <- paste(results$X_Var, results$Y_Var, sep = "_")

test_results <- results %>%
  group_by(pair_id) %>%
  summarise(
    X_Var = first(X_Var),
    Y_Var = first(Y_Var),
    r_BL = r[Group == "BL"] %>% first(),
    r_BH = r[Group == "BH"] %>% first(),
    n_BL = n[Group == "BL"] %>% first(),
    n_BH = n[Group == "BH"] %>% first(),
    z_BL = z[Group == "BL"] %>% first(),
    z_BH = z[Group == "BH"] %>% first(),
    z_diff = (z_BH - z_BL) / sqrt(1 / (n_BL - 3) + 1 / (n_BH - 3)),
    p_value = 2 * (1 - pnorm(abs((z_BH - z_BL) / sqrt(1 / (n_BL - 3) + 1 / (n_BH - 3)))))
  ) %>%
  ungroup()

# ② 効果量の追加
test_results <- test_results %>%
  mutate(
    cohen_q = z_BH - z_BL,
    cohen_q_abs = abs(z_BH - z_BL)
  )

# 表示
print(test_results)



#--------------------------------#プロット
# 相関ペアの定義
cor_pairs <- list(
  list(x = "averageNumerical", y = "HBM_average_benefit_3"),
  list(x = "averageNonNumerical", y = "HBM_average_benefit_3"),
  list(x = "averageNumerical", y = "HBM_average_risk_3"),
  list(x = "averageNonNumerical", y = "HBM_average_risk_3")
)

# 各ペアに対して散布図を作成
for (pair in cor_pairs) {
  x_var <- pair$x
  y_var <- pair$y
  
  p <- ggplot(merged_data_complete_for_path %>% filter(berlin_group %in% c("BL", "BH")),
              aes_string(x = x_var, y = y_var, color = "berlin_group")) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(
      title = paste("Scatterplot of", y_var, "vs", x_var),
      x = x_var,
      y = y_var,
      color = "Group"
    ) +
    theme_minimal()
  
  print(p)
}

