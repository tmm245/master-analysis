# ------------------------------------------------------------------
# ドロップアウト分析 (Attrition Analysis) [最終確定版]
# 対象: 初期回答者 571名 vs 最終解析対象者 220名
# ------------------------------------------------------------------

# 1. 分析用データの作成
attrition_data <- df1 %>%
  mutate(
    # 最終サンプル(merged_data_complete)にIDがあれば "Retained", なければ "Dropped"
    Status = ifelse(id %in% merged_data_complete$id, "Retained", "Dropped"),
    
    # 属性変数（df1にある変数名を使用）
    Age = as.numeric(Q21), 
    Gender = as.factor(Q18),      # 1=男性, 2=女性
    Effectiveness = Effectiveness_1,
    Safeness = Safeness_1
  )

# 人数確認 (Retained=220, Dropped=351 になるはず)
cat("【分析対象人数確認】\n")
print(table(attrition_data$Status))

# 2. 統計検定の実行とp値の出力
cat("\n【統計検定結果 (p値) 】\n")

# 年齢 (t検定)
p_age <- t.test(Age ~ Status, data = attrition_data)$p.value
cat(sprintf("1. 年齢 (Age)            : p = %.3f\n", p_age))

# 性別 (カイ二乗検定)
chisq_res <- chisq.test(table(attrition_data$Gender, attrition_data$Status))
cat(sprintf("2. 性別 (Gender)         : p = %.3f\n", chisq_res$p.value))

# 有効性認知 (t検定)
p_eff <- t.test(Effectiveness ~ Status, data = attrition_data)$p.value
cat(sprintf("3. 有効性 (Effectiveness): p = %.3f\n", p_eff))

# 安全性認知 (t検定)
p_safe <- t.test(Safeness ~ Status, data = attrition_data)$p.value
cat(sprintf("4. 安全性 (Safeness)     : p = %.3f\n", p_safe))