# ファイルの読み込み（Shift-JISエンコード）
mdf <- read.csv("merged2.csv", fileEncoding = "shift-jis")


#VIFの確認
# ――――――――――――――
# 文字ラベルで呼び出せるVIFチェック関数
# ――――――――――――――
get_vif_by_label <- function(fit, data, group_label) {
  # lavaanに登録されているグループラベルを取得（例: c("L","H")）
  labels <- inspect(fit, "group.label")
  
  # labelが見つからなければエラー
  if (!group_label %in% labels) {
    stop("group_label must be one of: ", paste(labels, collapse = ", "))
  }
  
  # グループ番号（1 or 2）を探す
  group_number <- which(labels == group_label)
  
  # 回帰パスだけを抽出
  est <- parameterEstimates(fit)
  regs <- subset(est, op == "~" & group == group_number)
  
  cat("===== VIFチェック：Group", group_label, "(lavaan group", group_number, ") =====\n\n")
  
  for (dv in unique(regs$lhs)) {
    preds <- regs$rhs[regs$lhs == dv]
    if (length(preds) < 2) next  # 変数が2つ以上ないとVIF不可
    
    # lm式を作成
    form <- as.formula(paste(dv, "~", paste(preds, collapse = " + ")))
    cat("VIF for:", deparse(form), "\n")
    
    # 該当グループのデータだけで回帰 → VIF
    df_sub <- subset(data, numAttitude_group == group_label)
    vif_val <- vif(lm(form, data = df_sub))
    print(vif_val)
    cat("\n")
  }
}

# ――――――――――――――
# 使い方（ラベル"L"／"H"をそのまま指定）
# ――――――――――――――
get_vif_by_label(fit, merged_data_complete_for_path, "L")
get_vif_by_label(fit, merged_data_complete_for_path, "H")




#クロンバックa計算
subjective_numeracy_alpha <- df2 %>% select(Q12_1, Q12_2, Q12_3, Q12_4, Q12_5)
numAttitude_alpha <- df4 %>% select(数値態度_1_reversed,	数値態度_2_reversed,	数値態度_3_reversed,	数値態度_4_reversed,	数値態度_5_reversed,数値態度_6, 数値態度_7, 数値態度_8, 数値態度_9, 数値態度_10)
HBM_Acceptance3_alpha <- df3 %>% select(Acceptance_1_reversed, Acceptance_2_reversed)
HBM_Benefit4_alpha <- df4 %>% select(HBM_3, HBM_4)
HBM_Risk4_alpha <- df4 %>% select(HBM_5, HBM_6)
Info_Info_Average <- df3 %>% select("Info_Info_1", "Info_Info_2", "Info_Info_3", "Info_Info_4", "Info_Info_5", "Info_Info_6", "Info_Info_7", "Info_Info_8", "Info_Info_9", "Info_Info_10", "Info_Info_11", "Info_Info_12")
averageNumerical <- df3 %>% select("Info_Info_1", "Info_Info_2", "Info_Info_3", "Info_Info_4", "Info_Info_5", "Info_Info_6", "Info_Info_7", "Info_Info_8", "Info_Info_9", "Info_Info_10", "Info_Info_11", "Info_Info_12", "Info_Cred_1", "Info_Cred_2", "Info_Cred_3", "Info_Cred_4", "Info_Cred_5", "Info_Cred_6")
averageNonNumerical <- df3 %>% select("Info_Info_7", "Info_Info_8", "Info_Info_9", "Info_Info_10", "Info_Info_11", "Info_Info_12", "Info_Cred_1", "Info_Cred_2", "Info_Cred_3", "Info_Cred_4", "Info_Cred_5", "Info_Cred_6", "Info_Cred_7", "Info_Cred_8", "Info_Cred_9", "Info_Cred_10", "Info_Cred_11", "Info_Cred_12", "Info_Cred_13")
HBM_Acceptance4_alpha <- df4 %>% select(HBM_1_reversed, HBM_2_reversed)
HBM_Benefit3_alpha <- df3 %>% select(Benefit_1_1, Benefit_2_1)
HBM_Risk3_alpha <- df3 %>% select(Risk_1_1, Risk_2_1)

# クロンバックのαを計算
alpha_result <- alpha(HBM_Risk3_alpha)
print(alpha_result)

#相関行列dfの作成
merged_df_for_Cor <- Reduce(function(x, y) merge(x, y, by = "id", all = TRUE), list(
  df1[, c("id", "berlin_correct_count")], 
  df2[, c("id", "mean_subjective_numeracy", "berlin_oneItem","mean_cliticalThinking_attitude")], 
  df3[, c("id", 
          "averageNumerical", "averageNonNumerical",
          "averageInformativenessPositiveNumerical", "averageInformativenessNetagtiveNumerical",
          "averageInformativenessPositiveNonNumerical", "averageInformativenessNetagtiveNonNumerical",
          "averageCredibilityPositiveNumerical", "averageCredibilityNetagtiveNumerical",
          "averageCredibilityPositiveNonNumerical", "averageCredibilityNetagtiveNonNumerical",
          "averageCredibility", "averageInformativeness",
          "averageInformativenessNumerical", "averageInformativenessNonNumerical",
          "averageCredibilityNumerical", "averageCredibilityNonNumerical","health_correct_count", "HBM_average_benefit_3", "HBM_average_risk_3")], 
  df4[, c("id", "numAttitude_average", "HBM_average_acceptance_4", "HBM_average_benefit_4", "HBM_average_risk_4", "num_vaccination_4")] 
))

#4回回答者の相関df
merged_data_complete_for_cor <- 
  merged_data_complete[, c("berlin_correct_count","mean_subjective_numeracy", "berlin_oneItem", 
                            "mean_cliticalThinking_attitude", "averageNumerical", "averageNonNumerical","Effectiveness_3","Omicron_3","Safeness_3",
                           "Possibility_of_infec_3","Regret_3","health_correct_count", "HBM_average_benefit_3", "HBM_average_risk_3",
                                 "numAttitude_average", "HBM_average_acceptance_4", "HBM_average_benefit_4", "HBM_average_risk_4", "num_vaccination_4", "num_vaccination_3")]

# `全体` の相関行列
# 数値でない列を削除
numeric_df <- merged_data_complete_for_cor[, sapply(merged_data_complete_for_cor, is.numeric)]
# 相関行列を計算
cor_matrix <- cor(numeric_df, use = "complete.obs")
# 小数点第2位まで丸める
cor_matrix <- round(cor_matrix, 2)
# 上三角（対角成分を含む）を NA にする
cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA
# 相関行列を出力
write.csv(cor_matrix, "merged_df_for_Cor_complete0723.csv", row.names = FALSE, fileEncoding = "CP932")

#compデータの相関
numeric_df <- subset(merged_data_complete_for_cor, mean_subjective_numeracy >= 4.1)
# 相関行列を計算
cor_matrix <- cor(numeric_df, use = "complete.obs")
# 小数点第2位まで丸める
cor_matrix <- round(cor_matrix, 2)
# 上三角（対角成分を含む）を NA にする
cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA
# 相関行列を出力
write.csv(cor_matrix, "merged_df_for_Subjective_comp.csv", row.names = FALSE, fileEncoding = "CP932")



# berlin_correct_count による高低
# berlin_correct_count が 0 or 1 のデータを抽出
merged_df_for_berlinLowNumeracy <- subset(merged_df_for_Cor, berlin_correct_count %in% c(0, 1))
# `merged_df_for_Cor` の数値列のみを抽出（`id` を除外）
numeric_df <- merged_df_for_berlinLowNumeracy[, -1]
# 相関行列を計算
cor_matrix <- cor(numeric_df, use = "complete.obs")
# 小数点第2位まで丸める
cor_matrix <- round(cor_matrix, 2)
# 上三角（対角成分を含む）を NA にする
cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA
# 相関行列を出力
#write.csv(cor_matrix, "merged_df_for_berlinLowNumeracy.csv", row.names = FALSE, fileEncoding = "CP932")

# berlin_correct_count が 2 のデータを抽出
merged_df_for_berlinMidNumeracy <- subset(merged_df_for_Cor, berlin_correct_count %in% c(2))
# `merged_df_for_Cor` の数値列のみを抽出（`id` を除外）
numeric_df <- merged_df_for_berlinMidNumeracy[, -1]
# 相関行列を計算
cor_matrix <- cor(numeric_df, use = "complete.obs")
# 小数点第2位まで丸める
cor_matrix <- round(cor_matrix, 2)
# 上三角（対角成分を含む）を NA にする
cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA
# 相関行列を出力
#write.csv(cor_matrix, "merged_df_for_berlinMidNumeracy.csv", row.names = FALSE, fileEncoding = "CP932")

# berlin_correct_count が 3 or 4 のデータを抽出
merged_df_for_berlinHighNumeracy <- subset(merged_df_for_Cor, berlin_correct_count %in% c(3, 4))
# `merged_df_for_Cor` の数値列のみを抽出（`id` を除外）
numeric_df <- merged_df_for_berlinHighNumeracy[, -1]
# 相関行列を計算
cor_matrix <- cor(numeric_df, use = "complete.obs")
# 小数点第2位まで丸める
cor_matrix <- round(cor_matrix, 2)
# 上三角（対角成分を含む）を NA にする
cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA
# 相関行列を出力
#write.csv(cor_matrix, "merged_df_for_berlinHighNumeracy.csv", row.names = FALSE, fileEncoding = "CP932")


#berlin_oneItemによる高低
# berlin_oneItem が 0 のデータを抽出
merged_df_for_berlinOneItemlowNumeracy <- subset(merged_df_for_Cor, berlin_oneItem == 0)
# `merged_df_for_Cor` の数値列のみを抽出（`id` を除外）
numeric_df <- merged_df_for_berlinOneItemlowNumeracy[, -1]
# 相関行列を計算
cor_matrix <- cor(numeric_df, use = "complete.obs")
# 小数点第2位まで丸める
cor_matrix <- round(cor_matrix, 2)
# 上三角（対角成分を含む）を NA にする
cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA
# 相関行列を出力
#write.csv(cor_matrix, "merged_df_for_berlinOneItemlowNumeracy.csv", row.names = FALSE, fileEncoding = "CP932")

# berlin_oneItem が 1 のデータを抽出
merged_df_for_berlinOneItemHighNumeracy <- subset(merged_df_for_Cor, berlin_oneItem == 1)
# `merged_df_for_Cor` の数値列のみを抽出（`id` を除外）
numeric_df <- merged_df_for_berlinOneItemHighNumeracy[, -1]
# 相関行列を計算
cor_matrix <- cor(numeric_df, use = "complete.obs")
# 小数点第2位まで丸める
cor_matrix <- round(cor_matrix, 2)
# 上三角（対角成分を含む）を NA にする
cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA
# 相関行列を出力
#write.csv(cor_matrix, "merged_df_for_berlinOneItemHighNumeracy.csv", row.names = FALSE, fileEncoding = "CP932")

# subjective_numeracy が 1～2.9 のデータを抽出
merged_df_for_subjectiveLownumeracy <- subset(merged_df_for_Cor, mean_subjective_numeracy <= 4.0)
# `merged_df_for_Cor` の数値列のみを抽出（`id` を除外）
numeric_df <- merged_df_for_subjectiveLownumeracy[, -1]
# 相関行列を計算
cor_matrix <- cor(numeric_df, use = "complete.obs")
# 小数点第2位まで丸める
cor_matrix <- round(cor_matrix, 2)
# 上三角（対角成分を含む）を NA にする
cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA
# 相関行列を出力
#write.csv(cor_matrix, "merged_df_for_subjectiveLownumeracy.csv", row.names = FALSE, fileEncoding = "CP932")

# subjective_numeracy が 3～6 のデータを抽出
merged_df_for_subjectiveHighnumeracy <- subset(merged_df_for_Cor, mean_subjective_numeracy >= 4.1)
# `merged_df_for_Cor` の数値列のみを抽出（`id` を除外）
numeric_df <- merged_df_for_subjectiveHighnumeracy[, -1]
# 相関行列を計算
cor_matrix <- cor(numeric_df, use = "complete.obs")
# 小数点第2位まで丸める
cor_matrix <- round(cor_matrix, 2)
# 上三角（対角成分を含む）を NA にする
cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA
# 相関行列を出力
#write.csv(cor_matrix, "merged_df_for_subjectiveHighnumeracy.csv", row.names = FALSE, fileEncoding = "CP932")

# 数値態度 が 1～2.4 のデータを抽出
merged_df_for_lowNumAttitude <- subset(merged_df_for_Cor, numAttitude_average <= 3.4)
# `merged_df_for_Cor` の数値列のみを抽出（`id` を除外）
numeric_df <- merged_df_for_lowNumAttitude[, -1]
# 相関行列を計算
cor_matrix <- cor(numeric_df, use = "complete.obs")
# 小数点第2位まで丸める
cor_matrix <- round(cor_matrix, 2)
# 上三角（対角成分を含む）を NA にする
cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA
# 相関行列を出力
#write.csv(cor_matrix, "merged_df_for_lowNumAttitude.csv", row.names = FALSE, fileEncoding = "CP932")

# 数値態度 が 2.5～5 のデータを抽出
merged_df_for_highNumAttitude <- subset(merged_df_for_Cor, numAttitude_average >= 3.5)
# `merged_df_for_Cor` の数値列のみを抽出（`id` を除外）
numeric_df <- merged_df_for_highNumAttitude[, -1]
# 相関行列を計算
cor_matrix <- cor(numeric_df, use = "complete.obs")
# 小数点第2位まで丸める
cor_matrix <- round(cor_matrix, 2)
# 上三角（対角成分を含む）を NA にする
cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA
# 相関行列を出力
#write.csv(cor_matrix, "merged_df_for_highNumAttitude.csv", row.names = FALSE, fileEncoding = "CP932")

#相関
# 必要な列を抽出
cor_data <- merged_data_complete_for_path %>%
  select(
    numAttitude_average,
    health_correct_count,
    berlin_correct_count,
    mean_subjective_numeracy,
    averageNonNumerical,
    averageNumerical
  )
# 相関行列を作成（欠損値を除外）
cor_matrix <- cor(cor_data, use = "pairwise.complete.obs", method = "pearson")
cor_matrix_rounded <- round(cor_matrix, 2)
# 相関行列を表示
print(cor_matrix_rounded)



