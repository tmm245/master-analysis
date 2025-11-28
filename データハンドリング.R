#testt
# 必要なパッケージの読み込み
library(dplyr)
library(tidyverse)
library(purrr)
library(psych)
library(ggplot2)
library(lavaan)
library(corrplot)
library(car)
# install.packages("rcompanion") # 初回のみ
library(rcompanion)

#-------------------------------------------------------------------------------------------------#
#前処理

# ファイルの読み込み（Shift-JISエンコード）
# データの読み込み & 重複IDの削除
survey1 <- read.csv("survey1.csv", fileEncoding = "shift-jis") 

survey2 <- read.csv("survey2.csv", fileEncoding = "shift-jis") 

survey3 <- read.csv("survey3.csv", fileEncoding = "shift-jis") 

survey4 <- read.csv("survey4.csv", fileEncoding = "shift-jis") 

survey1demo <- read.csv("survey1DemoData.csv", fileEncoding = "UTF-8") 


#-------------------------------------------------------------------------------------------------#
#1回目データ処理

# 1回目の解答必須項目
target_columns1 <- c(
  "Q34", "Q57_1", "Q47", "Q16_1", "Q39_1", "Q40_1", "Q41_1", "Q42_1", 
  "Q53", "Q56", "Q55", "Q54", "Q21", "Q18", "Q22", "Q43_1", 
  "Q28_1", "Q28_2", "Q28_3", "Q28_4", "Q28_5", "Q28_6", "Q28_7", "Q28_8", "Q28_9", "Q28_10", 
  "Q51_1", "Q59_1", "Q52_1", "Q60_1", 
  "Q15_1", "Q15_2", "Q15_3", "Q15_4", "Q15_5", "Q15_6", "Q15_7", "Q15_8", "Q15_9", "Q15_10", 
  "Q58_1", "Q58_2", "Q58_3", "Q58_4", "Q58_5", "Q58_6", "Q58_7", "Q58_8", "Q58_9", 
  "Q61_1", "Q50", "Q46_1", "Q29"
)

# Q28の列（値の一様性をチェックするための対象）
Q28_columns <- c("Q28_1", "Q28_2", "Q28_3", "Q28_4", "Q28_5", 
                 "Q28_6", "Q28_7", "Q28_8", "Q28_9", "Q28_10")

# 1回目用のフィルター関数
filter_data_forOne <- function(df) {
  df <- df %>%
    # 不適切な回答を削除
    filter(!(DistributionChannel == "preview" | #anonymous以外の削除
               Duration..in.seconds. >= 1769 | #回答時間での足切り
               Duration..in.seconds. <= 180 |  #回答時間での足切り
               Q15_5 != 1 |#アテンションチェック項目の確認
               Finished != 1 | #回答完了していないもの除外
               Q34 !=2) #基礎疾患ありor答えたくないの除外
    ) %>%
    # id列がNAまたは空文字の行を削除
    filter(!is.na(id) & id != "") %>%
    # 必須項目のいずれかが NA の場合は削除
    filter(!if_any(all_of(target_columns1), is.na)) %>%
    # Q28_1 ～ Q28_10 のストレートライン
    rowwise() %>%
    filter(length(unique(na.omit(c_across(all_of(Q28_columns))))) > 1) %>%
    ungroup()
  
  return(df)
}

# 1回目ベルリンニューメラシーの正答数計算とグループ分類
berlin_correct_answers <- function(df) {
  # 各設問の正答をリスト化
  correct_answers <- c(Q53 = 3, Q56 = 2, Q55 = 1, Q54 = 3)
  
  # 行ごとに各列の正誤を比較し、正答数を計算
  df$berlin_correct_count <- rowSums(sweep(df[, names(correct_answers)], 2, correct_answers, FUN = "=="), na.rm = TRUE)
  
  # 正答数に基づいてグループを分類(2群)
  df$berlin_group <- ifelse(df$berlin_correct_count >= 2, "BH", "BL")
  
  # 正答数に基づいてグループを分類(3群群)
  #df <- df %>%
  #mutate(berlin_group = case_when(
  #berlin_correct_count %in% c(0, 1) ~ "BL",
  #berlin_correct_count == 2        ~ "BM",
  #berlin_correct_count %in% c(3, 4) ~ "BH",
  #TRUE                              ~ NA_character_
  #))
  
  
  return(df)
}

#-------------------------------------------------------------------------------------------------#
#2回目データ処理

# 2回目の解答必須項目
target_columns2 <- c(
  "QID1", "Q5_1", "Q6_1", "Q7_1", "Q8_1", "Q9", "COVID_pos", "COVID_pos_others", 
  "COVID_concern_oth", "COVID_gen_harm", "COVID_comply_percent", "COVID_comply_pers", 
  "CS_downloaded", "gov_approval", "COVID_info_source", "Q60_1", "Q60_2", "Q60_3", 
  "Q60_4", "Q60_5", "Q60_6", "Q60_7", "Q60_8", "Q60_9", "Q60_10", "Q60_11", 
  "Q60_12", "Q150", "Q63_1", "Q63_2", "Q63_3", "Q63_4", "Q63_5", "Q65_1", 
  "Q65_2", "Q65_3", "Q65_4", "Q65_5", "Q65_6", "Q65_7", "Q65_8", "Q65_9", 
  "Q65_10", "Q65_11", "Q65_12", "Q65_13", "Q65_14", "Q65_15", "Q65_16", "Q12_1", 
  "Q12_2", "Q12_3", "Q12_4", "Q12_5", "Q13"
)

# 2回目用のフィルター関数
filter_data_forTwo <- function(df) {
  df <- df %>%
    # 不適切な回答を削除
    filter(!(DistributionChannel == "preview" | 
               Finished != 1 | #回答完了していないもの除外
               Duration..in.seconds. >= 1702 |  
               Duration..in.seconds. <= 180)
    ) %>%
    # id列がNAまたは空文字の行を削除
    filter(!is.na(id) & id != "") %>%
    # 必須項目のいずれかが NA の場合は削除
    filter(!if_any(all_of(target_columns2), is.na))
  
  return(df)
}

#2回目主観的ニューメラシー計算
mean_subjective_numeracy <- function(df) {
  # 主観的ニューメラシーの平均を計算
  df$mean_subjective_numeracy <- rowMeans(df[, c("Q12_1", "Q12_2", "Q12_3", "Q12_4", "Q12_5")], na.rm = TRUE)
  
  # mean_subjective_numeracy に基づいて subjectiveNumeracy_group を作成
  df$subjectiveNumeracy_group <- ifelse(df$mean_subjective_numeracy <= 4.0, "L", "H")
  
  return(df)
}


#2回目批判的思考態度の計算
mean_cliticalThinking_attitude <- function(df) {
  df$mean_cliticalThinking_attitude <- rowMeans(df[, c("Q65_1", "Q65_2", "Q65_3", "Q65_4", "Q65_5", "Q65_6", "Q65_7", "Q65_8", "Q65_9", "Q65_10", "Q65_11", "Q65_12", "Q65_13", "Q65_14", "Q65_15", "Q65_16")], na.rm = TRUE)
  return(df)
}

#2回目のベルリンワンアイテムの正答計算
berlin_oneItem <- function(df) {
  df$berlin_oneItem <- ifelse(df$Q13 == 2, 1, 0)
  return(df)
}

#-------------------------------------------------------------------------------------------------#
#3回目データ処理

# 3回目の解答必須項目
target_columns3 <- c(
  "Q1", "Q2", "Q3", "Q4_1", "Q5_1", "Q6_1", "Q7_1", "Q8_1", 
  "Attitude_1_1", "Attitude_2_1", "Benefit_1_1", "Benefit_2_1", 
  "Risk_1_1", "Risk_2_1", "Q9", "Q10", "Q11", "Q12", "Q13", 
  "Info_Info_1", "Info_Info_2", "Info_Info_3", "Info_Info_4", "Info_Info_5", 
  "Info_Info_6", "Info_Info_7", "Info_Info_8", "Info_Info_9", "Info_Info_10", 
  "Info_Info_11", "Info_Info_12", "Info_Cred_1", "Info_Cred_2", "Info_Cred_3", 
  "Info_Cred_4", "Info_Cred_5", "Info_Cred_6", "Info_Cred_7", "Info_Cred_8", 
  "Info_Cred_9", "Info_Cred_10", "Info_Cred_11", "Info_Cred_12", "Info_Cred_13", 
  "H_Num1", "H_Num2", "H_Num3", "H_Num4", "H_Num5"
)

# 3回目用のフィルター関数
filter_data_forThree <- function(df) {
  df <- df %>%
    # 不適切な回答を削除
    filter(!(DistributionChannel == "preview" | #anonymous以外の削除
               Finished != 1 | #回答完了していないもの除外
               Duration..in.seconds. >= 1111 | #回答時間での足切り
               Duration..in.seconds. <= 180 |  #回答時間での足切り
               Info_Cred_7 != 1) #アテンションチェック項目の確認
    ) %>%
    # id列がNAまたは空文字の行を削除
    filter(!is.na(id) & id != "") %>%
    # 必須項目のいずれかが NA の場合は削除
    filter(!if_any(all_of(target_columns3), is.na))
  
  return(df)
}

#3回目のヘルスニューメラシー計算
health_correct_answers <- function(df) {
  # 各設問の正答をリスト化
  h_correct_answers <- c(H_Num1 = 4, H_Num2 = 2, H_Num3 = 3, H_Num4 = 2, H_Num5 = 2)
  
  # 行ごとに各列の正誤を比較し、正答数を計算
  df$health_correct_count <- rowSums(sweep(df[, names(h_correct_answers)], 2, h_correct_answers, FUN = "=="))
  
  return(df)
}

#3回目変数平均値計算

#数値情報整理
averageNonNumericalVars <- c(
  "Info_Info_7", "Info_Info_8", "Info_Info_9",
  "Info_Info_10", "Info_Info_11", "Info_Info_12",
  "Info_Cred_8", "Info_Cred_9", "Info_Cred_10",
  "Info_Cred_11", "Info_Cred_12", "Info_Cred_13"
)
averageNumericalVars <- c(
  "Info_Info_1", "Info_Info_2", "Info_Info_3",
  "Info_Info_4", "Info_Info_5", "Info_Info_6",
  "Info_Cred_1", "Info_Cred_2", "Info_Cred_3",
  "Info_Cred_4", "Info_Cred_5", "Info_Cred_6"
)

calculate_infoUse <- function(df) {
  df$averageInformativenessPositiveNumerical <- rowMeans(df[, c("Info_Info_1", "Info_Info_2", "Info_Info_3")], na.rm = TRUE)
  df$averageInformativenessNetagtiveNumerical <- rowMeans(df[, c("Info_Info_4", "Info_Info_5", "Info_Info_6")], na.rm = TRUE)
  df$averageInformativenessPositiveNonNumerical <- rowMeans(df[, c("Info_Info_7", "Info_Info_8", "Info_Info_9")], na.rm = TRUE)
  df$averageInformativenessNetagtiveNonNumerical <- rowMeans(df[, c("Info_Info_10", "Info_Info_11", "Info_Info_12")], na.rm = TRUE)
  
  df$averageCredibilityPositiveNumerical <- rowMeans(df[, c("Info_Cred_1", "Info_Cred_2", "Info_Cred_3")], na.rm = TRUE)
  df$averageCredibilityNetagtiveNumerical <- rowMeans(df[, c("Info_Cred_4", "Info_Cred_5", "Info_Cred_6")], na.rm = TRUE)
  df$averageCredibilityPositiveNonNumerical <- rowMeans(df[, c("Info_Cred_8", "Info_Cred_9", "Info_Cred_10")], na.rm = TRUE)
  df$averageCredibilityNetagtiveNonNumerical <- rowMeans(df[, c("Info_Cred_11", "Info_Cred_12", "Info_Cred_13")], na.rm = TRUE)
  
  df$averageCredibility <- rowMeans(df[, c("Info_Cred_1", "Info_Cred_2", "Info_Cred_3", "Info_Cred_4", "Info_Cred_5", "Info_Cred_6", "Info_Cred_7", "Info_Cred_8", "Info_Cred_9", "Info_Cred_10", "Info_Cred_11", "Info_Cred_12", "Info_Cred_13")], na.rm = TRUE)
  df$averageInformativeness <- rowMeans(df[, c("Info_Info_1", "Info_Info_2", "Info_Info_3", "Info_Info_4", "Info_Info_5", "Info_Info_6", "Info_Info_7", "Info_Info_8", "Info_Info_9", "Info_Info_10", "Info_Info_11", "Info_Info_12")], na.rm = TRUE)
  
  df$averageInformativenessNumerical <- rowMeans(df[, c("Info_Info_1", "Info_Info_2", "Info_Info_3", "Info_Info_4", "Info_Info_5", "Info_Info_6")], na.rm = TRUE)
  df$averageInformativenessNonNumerical <- rowMeans(df[, c("Info_Info_7", "Info_Info_8", "Info_Info_9", "Info_Info_10", "Info_Info_11", "Info_Info_12")], na.rm = TRUE)
  
  df$averageCredibilityNumerical <- rowMeans(df[, c("Info_Cred_1", "Info_Cred_2", "Info_Cred_3", "Info_Cred_4", "Info_Cred_5", "Info_Cred_6")], na.rm = TRUE)
  df$averageCredibilityNonNumerical <- rowMeans(df[, c("Info_Cred_8", "Info_Cred_9", "Info_Cred_10", "Info_Cred_11", "Info_Cred_12", "Info_Cred_13")], na.rm = TRUE)
  
  df$averageNumerical  <- rowMeans(df[, averageNumericalVars], na.rm = TRUE)
  df$averageNonNumerical <- rowMeans(df[, averageNonNumericalVars], na.rm = TRUE)
  
  df$HBM_average_benefit_3 <- rowMeans(df[, c("Benefit_1_1", "Benefit_2_1")], na.rm = TRUE)
  df$HBM_average_risk_3 <- rowMeans(df[, c("Risk_1_1", "Risk_2_1")], na.rm = TRUE)
  
  # --- 逆転項目の計算（指定フォーマットに合わせる） ---
  df$Acceptance_1_reversed <- 6 - as.numeric(df$Attitude_1_1)
  df$Acceptance_2_reversed <- 6 - as.numeric(df$Attitude_2_1)
  # --- 逆転後の平均値計算 ---
  df$HBM_average_acceptance_3 <- rowMeans(
    df[, c("Acceptance_1_reversed", "Acceptance_2_reversed")],
    na.rm = TRUE
  )
  return(df)
}

#-------------------------------------------------------------------------------------------------#
#4回目データ処理

# 4回目の解答必須項目
target_columns4 <- c(
  "Q1", "Effectiveness_1", "Omicron_1", "Safeness_1", "Possibility.of.infec_1", 
  "Regret_1", "Q9", "Q12", "Q10", "Q2", "HBM_1", "HBM_2", "HBM_3", "HBM_4", "HBM_5", "HBM_6", 
  "数値態度_1", "数値態度_2", "数値態度_3", "数値態度_4", "数値態度_5", 
  "数値態度_6", "数値態度_7", "数値態度_8", "数値態度_9", "数値態度_10"
)

# 4回目用のフィルター関数
filter_data_forFour <- function(df) {
  df <- df %>%
    # 不適切な回答を削除
    filter(!(DistributionChannel == "preview" | #anonymous以外の削除
               Finished != 1 | #回答完了していないもの除外
               Duration..in.seconds. >= 1232 |   #回答時間での足切り
               Duration..in.seconds. <= 180) #回答時間での足切り
    ) %>%
    # id列がNAまたは空文字の行を削除
    filter(!is.na(id) & id != "") %>%
    # 必須項目のいずれかが NA の場合は削除
    filter(!if_any(all_of(target_columns4), is.na))
  
  return(df)
}

#4回目の数値態度操作
reverse_items <- function(df) {
  # 逆転項目の列 (数値態度_1 ～ 数値態度_5)
  reverse_columns <- c("数値態度_1", "数値態度_2", "数値態度_3", "数値態度_4", "数値態度_5")
  
  # そのままの値を使う列 (数値態度_6 ～ 数値態度_10)
  original_columns <- c("数値態度_6", "数値態度_7", "数値態度_8", "数値態度_9", "数値態度_10")
  
  # 逆転項目の計算
  df <- df %>%
    mutate(across(all_of(reverse_columns), ~ 6 - ., .names = "{.col}_reversed"))
  
  # 逆転後の値を集める列 (数値態度_1_reversed ～ 数値態度_5_reversed)
  reversed_columns <- paste0(reverse_columns, "_reversed")
  
  # 逆転後の値 + 数値態度_6 ～ 数値態度_10 の平均を計算
  df <- df %>%
    rowwise() %>%
    mutate(数値態度_平均 = mean(c_across(all_of(c(reversed_columns, original_columns))), na.rm = TRUE)) %>%
    ungroup() %>%
    # 数値態度_平均に基づいて numAttitude_group を作成
    mutate(numAttitude_group = ifelse(数値態度_平均 <= 3.4, "L", "H"))
  
  return(df)
}


#4回目HBM計算
create_HBM_scores <- function(df) {
  df %>%
    mutate(
      # 逆転項目の計算
      HBM_1_reversed = 6 - as.numeric(HBM_1),
      HBM_2_reversed = 6 - as.numeric(HBM_2)
    ) %>%
    mutate(
      # 平均値の計算（all_of() を使用）
      HBM_average_acceptance_4 = rowMeans(select(., all_of(c("HBM_1_reversed", "HBM_2_reversed"))), na.rm = TRUE),
      HBM_average_benefit_4 = rowMeans(select(., all_of(c("HBM_3", "HBM_4"))), na.rm = TRUE),
      HBM_average_risk_4 = rowMeans(select(., all_of(c("HBM_5", "HBM_6"))), na.rm = TRUE)
    )
}


#-------------------------------------------------------------------------------------------------#
#共通処理


# 重複削除処理 Duration..in.seconds. が長い方を残す関数
keep_longest_duration_by_id <- function(df) {
  df %>%
    group_by(id) %>%
    slice_max(order_by = Duration..in.seconds., with_ties = FALSE) %>%
    ungroup()
}



# ハンドリング関数を適用
df1 <- berlin_correct_answers(filter_data_forOne(survey1))
df2 <- mean_cliticalThinking_attitude(berlin_oneItem(mean_subjective_numeracy(filter_data_forTwo(survey2))))
df3 <- calculate_infoUse(health_correct_answers(filter_data_forThree(survey3)))
df4 <- create_HBM_scores(reverse_items(filter_data_forFour(survey4)))

write.csv(df3, "重複チェック前.csv", row.names = FALSE, fileEncoding = "CP932")

#重複IDリスト
# リストでまとめて適用
list(df1 = df1, df2 = df2, df3 = df3, df4 = df4) %>%
  lapply(function(df) {
    df %>%
      count(id) %>%
      filter(n > 1)
  })


#重複削除前行数
cat("行数一覧:\n")
cat("df1:", nrow(df1), "\n")
cat("df2:", nrow(df2), "\n")
cat("df3:", nrow(df3), "\n")
cat("df4:", nrow(df4), "\n")


# 重複削除関数を適用
df1 <- keep_longest_duration_by_id(df1)
df2 <- keep_longest_duration_by_id(df2)
df3 <- keep_longest_duration_by_id(df3)
df4 <- keep_longest_duration_by_id(df4)

write.csv(df3, "重複チェック後.csv", row.names = FALSE, fileEncoding = "CP932")

#重複削除後行数
cat("行数一覧:\n")
cat("df1:", nrow(df1), "\n")
cat("df2:", nrow(df2), "\n")
cat("df3:", nrow(df3), "\n")
cat("df4:", nrow(df4), "\n")


#リネーム
df1 <- df1 %>%
  rename(Effectiveness_1 = Q40_1, 
         Safeness_1 = Q42_1)
df2 <- df2 %>%
  rename(Effectiveness_2 = Q5_1, 
         Safeness_2 = Q6_1,
         num_vaccination_2 = QID1)
df4 <- df4 %>%
  rename(numAttitude_average = 数値態度_平均, 
         num_vaccination_4 = Q1,
         Effectiveness_4 = Effectiveness_1,
         Omicron_4 = Omicron_1,
         Safeness_4 = Safeness_1,
         POS_4  = Possibility.of.infec_1,
         Regret_4  = Regret_1  )
df3 <- df3 %>%
  rename(num_vaccination_3 = Q1,
         Effectiveness_3 = Q4_1,
         Omicron_3 = Q5_1,
         Safeness_3 = Q6_1,
         Possibility_of_infec_3 = Q7_1,
         Regret_3 =Q8_1)


# # 列名にサフィックスを追加
# colnames(df1) <- paste0(colnames(survey1), "_S1")
# colnames(df2) <- paste0(colnames(survey2), "_S2")
# colnames(df3) <- paste0(colnames(survey3), "_S3")
# colnames(df4) <- paste0(colnames(survey4), "_S4")

#結合1
# 4つのデータフレームを id をキーに結合
# merged_df <- Reduce(function(x, y) merge(x, y, by = "id", all = TRUE), list(df1, df2, df3, df4))
# merged_df を "merged.csv" として保存 (Shift-JIS エンコーディング)
# write.csv(merged_df, "merged1.csv", row.names = FALSE, fileEncoding = "CP932")


#結合2
# 1. semi_join用のフィルタ処理（id列はそのまま）
filtered_dfs <- list(df2, df3, df4, survey1demo) %>% 
  map(~ semi_join(.x, df1, by = "id"))

# 2. left_joinでdf1に順に統合（id列を共有キーに）
merged_data <- reduce(filtered_dfs, left_join, by = "id", .init = df1)

# 3. 各元dfのid列を保持するため、結合元のid列を付け直す
# （各dfにもとのidが入っていたと仮定してid列をsuffix付きで残す）
merged_data <- merged_data %>%
  mutate(
    id_df2 = ifelse(id %in% df2$id, id, NA),
    id_df3 = ifelse(id %in% df3$id, id, NA),
    id_df4 = ifelse(id %in% df4$id, id, NA),
    id_demo = ifelse(id %in% survey1demo$id, id, NA)
  )

#メタデータ用意
# 必要な列だけ抽出
berlin_group_df <- df1 %>% select(id, berlin_group)
berlin_oneItem_df <- df2 %>% select(id, berlin_oneItem)
subjectiveNumeracy_df <- df2 %>% select(id, subjectiveNumeracy_group)
numAttitude_df <- df4 %>% select(id, numAttitude_group)

# 結合していく
meta_data <- survey1demo %>%
  left_join(berlin_group_df, by = "id") %>%
  left_join(berlin_oneItem_df, by = "id") %>%
  left_join(subjectiveNumeracy_df, by = "id") %>%
  left_join(numAttitude_df, by = "id")

#-------------------------------------------------------------------------------------------------#
#df結合

# データの準備（例：id列が共通と仮定）
dfs <- list(df1 = df1, df2 = df2, df3 = df3, df4 = df4)

# 名前付きの全組み合わせ（2つ以上）
combinations <- unlist(lapply(2:4, function(n) combn(names(dfs), n, simplify = FALSE)), recursive = FALSE)

# 結果格納用リスト
merged_sets <- list()

# 各組み合わせについて結合処理
for (combo in combinations) {
  name <- paste(combo, collapse = "_")
  
  # full_joinを順に適用（reduce）
  merged <- reduce(dfs[combo], full_join, by = "id")
  
  # meta_data を必ず left_join
  merged_with_demo <- left_join(merged, meta_data, by = "id")
  
  # 結果をリストに格納
  merged_sets[[name]] <- merged_with_demo
}

# 例：1つ確認
# View(merged_sets[["df1_df3_df4"]])

#-------------------------------------------------------------------------------------------------#
#マージ処理

# merged_df を "merged.csv" として保存 (Shift-JIS エンコーディング)
write.csv(merged_data, "merged2.csv", row.names = FALSE, fileEncoding = "CP932")

#4回すべて回答したサンプルのdf
merged_data_complete <- subset(merged_data, 
                               !is.na(Q34) & !is.na(num_vaccination_2) & !is.na(num_vaccination_3) & !is.na(num_vaccination_4))

#Subjective Numeracyのグループ分け
# データフレームのコピー
merged_data_complete_for_path <- merged_data_complete
# 条件に基づいて新しい列を追加
merged_data_complete_for_path$subjectiveNumeracy_group <- ifelse(merged_data_complete_for_path$mean_subjective_numeracy <= 4.0, "SL", "SH")
# 条件に基づいて新しい列を追加erged_data_complete_for_path$numAttitude_average <= 3.4, "NL", "NH")