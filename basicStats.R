############################################################
# 修士論文：基礎統計 & 相関行列 出力コード
# （psychのみ使用 + 年齢・性別・学歴 + Acceptance + t1/t2/t3/t4変数追加）
#
# 前提：
#   - データフレーム: merged_data_complete が環境に存在している
#   - 主な列:
#       id,
#       age,
#       gender (1=男性, 2=女性),
#       education (1=中学,2=高校,3=高専,4=専門/専修,
#                  5=短大,6=大学,7=大学院,8=その他),
#       Safeness_1, Safeness_2, Safeness_3, Safeness_4,
#       Effectiveness_1, Effectiveness_2, Effectiveness_3, Effectiveness_4,
#       HBM_average_benefit_3, HBM_average_benefit_4,
#       HBM_average_risk_3,    HBM_average_risk_4,
#       HBM_average_acceptance_3, HBM_average_acceptance_4,
#       num_vaccination_2, num_vaccination_3, num_vaccination_4,
#       mean_cliticalThinking_attitude or mean_criticalThinking_attitude,
#       berlin_correct_count, mean_subjective_numeracy,
#       health_correct_count, numAttitude_average
#
#   - LMM と同じ complete case サンプル（long_cc）に基づいて
#     基礎統計と相関を出す。
############################################################

## パッケージ読み込み --------------------------------------------------------
library(psych)  # 基礎統計＆相関

## 0. 前提チェック -----------------------------------------------------------

if (!exists("merged_data_complete")) {
  stop("オブジェクト 'merged_data_complete' が環境にありません。先に読み込んでください。")
}

if (!("id" %in% names(merged_data_complete))) {
  stop("'id' 列が merged_data_complete に存在しません。")
}

## 1. CT列名のゆれを吸収 ----------------------------------------------------

ct_candidates <- c("mean_cliticalThinking_attitude",
                   "mean_criticalThinking_attitude")

ct_name <- ct_candidates[ct_candidates %in% names(merged_data_complete)][1]

if (is.na(ct_name) || is.null(ct_name)) {
  stop("CT列 (mean_cliticalThinking_attitude / mean_criticalThinking_attitude) が見つかりません。")
}

cat("CT列として使用する変数:", ct_name, "\n")

## 2. LMMで使う long データを base R で再構成 -------------------------------
# mk_wave: t3 / t4 のデータを LMM 用の形式に変換する関数
# （Benefit / Risk 両方を含む wide形式）

mk_wave <- function(df, wave, ct_name) {
  if (wave == "t3") {
    out <- data.frame(
      id                = df$id,
      wave              = factor("t3", levels = c("t3", "t4")),
      CT                = df[[ct_name]],
      Effectiveness_lag = df$Effectiveness_2,
      Safeness_lag      = df$Safeness_2,
      Benefit           = df$HBM_average_benefit_3,
      Risk              = df$HBM_average_risk_3
      # Acceptance は LMM では使っていないのでここでは入れない
    )
  } else { # wave == "t4"
    out <- data.frame(
      id                = df$id,
      wave              = factor("t4", levels = c("t3", "t4")),
      CT                = df[[ct_name]],
      Effectiveness_lag = df$Effectiveness_3,
      Safeness_lag      = df$Safeness_3,
      Benefit           = df$HBM_average_benefit_4,
      Risk              = df$HBM_average_risk_4
    )
  }
  out
}

# t3 / t4 の wide形式データを作成
wave_t3 <- mk_wave(merged_data_complete, "t3", ct_name)
wave_t4 <- mk_wave(merged_data_complete, "t4", ct_name)

# Benefit 用の long
long_ben_t3 <- data.frame(
  id                = wave_t3$id,
  wave              = wave_t3$wave,
  CT                = wave_t3$CT,
  Effectiveness_lag = wave_t3$Effectiveness_lag,
  Safeness_lag      = wave_t3$Safeness_lag,
  outcome           = factor("Benefit", levels = c("Benefit", "Risk")),
  y                 = wave_t3$Benefit
)

long_ben_t4 <- data.frame(
  id                = wave_t4$id,
  wave              = wave_t4$wave,
  CT                = wave_t4$CT,
  Effectiveness_lag = wave_t4$Effectiveness_lag,
  Safeness_lag      = wave_t4$Safeness_lag,
  outcome           = factor("Benefit", levels = c("Benefit", "Risk")),
  y                 = wave_t4$Benefit
)

# Risk 用の long
long_risk_t3 <- data.frame(
  id                = wave_t3$id,
  wave              = wave_t3$wave,
  CT                = wave_t3$CT,
  Effectiveness_lag = wave_t3$Effectiveness_lag,
  Safeness_lag      = wave_t3$Safeness_lag,
  outcome           = factor("Risk", levels = c("Benefit", "Risk")),
  y                 = wave_t3$Risk
)

long_risk_t4 <- data.frame(
  id                = wave_t4$id,
  wave              = wave_t4$wave,
  CT                = wave_t4$CT,
  Effectiveness_lag = wave_t4$Effectiveness_lag,
  Safeness_lag      = wave_t4$Safeness_lag,
  outcome           = factor("Risk", levels = c("Benefit", "Risk")),
  y                 = wave_t4$Risk
)

# 全部まとめて long データに（outcome: Benefit / Risk）
long_all <- rbind(long_ben_t3, long_ben_t4,
                  long_risk_t3, long_risk_t4)

# 必要な変数がすべてNAでない行だけを採用（LMMと同じ complete case 原則）
complete_rows <- !is.na(long_all$id) &
  !is.na(long_all$CT) &
  !is.na(long_all$Effectiveness_lag) &
  !is.na(long_all$Safeness_lag) &
  !is.na(long_all$y)

long_cc <- long_all[complete_rows, ]

cat("LMM用 long データ（complete case）の行数:", nrow(long_cc), "\n")

# ここで使われている id を抽出
analysis_ids <- unique(long_cc$id)

cat("LMM解析に使われた id 数:", length(analysis_ids), "\n")

## 3. 基礎統計用のデータセットを LMMサンプルに絞る -----------------------

# merged_data_complete を analysis_ids でフィルタ
merged_data_analysis <- merged_data_complete[merged_data_complete$id %in% analysis_ids, ]

cat("基礎統計対象データ（LMMサンプル絞り込み後）の次元:",
    paste(dim(merged_data_analysis), collapse = " x "), "\n")

## 4. 基礎統計に含める変数リスト -------------------------------------------
# ★ 年齢・性別・学歴 + Acceptance + t1/t2/t3/t4 の指標を含める

var_names <- c(
  # ---- デモグラフィック ----
  "age",
  "gender",     # 1=男性, 2=女性
  "education",  # 1=中学,2=高校,3=高専,4=専門/専修,5=短大,6=大学,7=大学院,8=その他
  
  # ---- ワクチン認知（t1〜t4）----
  "Safeness_1", "Safeness_2", "Safeness_3", "Safeness_4",
  "Effectiveness_1", "Effectiveness_2", "Effectiveness_3", "Effectiveness_4",
  "HBM_average_benefit_3",    "HBM_average_benefit_4",
  "HBM_average_risk_3",       "HBM_average_risk_4",
  "HBM_average_acceptance_3", "HBM_average_acceptance_4",
  
  # ---- 接種回数（t2/t3/t4）----
  "num_vaccination_2",
  "num_vaccination_3",
  "num_vaccination_4",
  
  # ---- TIC（個人特性） ----
  ct_name,
  "berlin_correct_count",
  "mean_subjective_numeracy",
  "health_correct_count",
  "numAttitude_average"
)

# 実際に存在する列だけに絞る
var_names <- var_names[var_names %in% names(merged_data_analysis)]

if (length(var_names) == 0) {
  stop("指定した変数が merged_data_analysis に存在しません。")
}

cat("基礎統計を出す変数:\n")
print(var_names)

## 5. desc_dat を作成し、CT列名を 'CT' に統一 ------------------------------

# base R で列抽出
desc_dat <- merged_data_analysis[, var_names, drop = FALSE]

# CT列の名称を 'CT' に変更
name_vec <- names(desc_dat)
name_vec[name_vec == ct_name] <- "CT"
names(desc_dat) <- name_vec

# 念のため行数と列数を表示
cat("desc_dat の次元:", paste(dim(desc_dat), collapse = " x "), "\n")

## 6. psych::describe による基礎統計 ---------------------------------------

desc_table <- psych::describe(desc_dat)

cat("\n=== 基礎統計 (psych::describe) ===\n")
print(desc_table)

# CSV出力（Table 1 用の下地）
write.csv(
  desc_table,
  file = "basic_descriptives_table1_psych_LMMsample.csv",
  row.names = TRUE
)

cat("基礎統計（psych::describe）の結果を 'basic_descriptives_table1_psych_LMMsample.csv' に出力しました。\n")

## 7. psych::corr.test による相関 & p値（補正なし） ------------------------

cor_test <- psych::corr.test(
  desc_dat,
  use    = "pairwise",
  method = "pearson",
  adjust = "none"   # 報告には補正なしp値を使用
)

cor_r <- cor_test$r  # 相関係数
cor_p <- cor_test$p  # その r に対応する「素の p値」

cat("\n=== 相関係数 (r) ===\n")
print(cor_r)

cat("\n=== 相関の p値 (unadjusted) ===\n")
print(cor_p)

# CSV出力
write.csv(
  cor_r,
  file = "correlations_r_psych_LMMsample.csv",
  row.names = TRUE
)

write.csv(
  cor_p,
  file = "correlations_p_psych_LMMsample.csv",
  row.names = TRUE
)

cat("相関係数 (r) を 'correlations_r_psych_LMMsample.csv' に、p値を 'correlations_p_psych_LMMsample.csv' に出力しました。\n")

############################################################
# 作成されるファイル：
#  - basic_descriptives_table1_psych_LMMsample.csv
#      → age, gender, education, t1〜t4指標, Acceptance, TIC等を含む
#         N, Mean, SD, Min, Max などの基礎統計。
#  - correlations_r_psych_LMMsample.csv
#      → 上記変数すべての Pearson 相関係数 (r)
#  - correlations_p_psych_LMMsample.csv
#      → 各 r に対応する unadjusted p値
#
# いずれも LMM で使った complete case サンプルに限定。
# gender / education は数値として相関に入ります
# （Table の注釈で「gender: 1=male, 2=female」
#   「education: 1=junior high ... 8=other」などと説明）。
############################################################

#-----------------------------------------------------------------------#
library(dplyr)
library(tidyr)

# 対象変数名
vacc_vars <- paste0("num_vaccination_", 2:4)

# 2～4回目すべてをロング型にして集計
tab_vacc_2to4 <- merged_data_complete %>%
  # num_vaccination_2,3,4 だけ取り出し
  dplyr::select(dplyr::all_of(vacc_vars)) %>%
  # ロング型に変換（wave 列と num_vaccination 列を作る）
  pivot_longer(
    cols = everything(),
    names_to = "wave",
    values_to = "num_vaccination"
  ) %>%
  # 欠損は除外
  filter(!is.na(num_vaccination)) %>%
  # wave（2〜4）ごと×接種回数ごとに件数を集計
  group_by(wave, num_vaccination) %>%
  summarise(n = n(), .groups = "drop") %>%
  # waveごとに割合（％）を計算
  group_by(wave) %>%
  mutate(
    percent = n / sum(n) * 100
  ) %>%
  ungroup() %>%
  arrange(wave, num_vaccination)

# 結果を表示
tab_vacc_2to4

#-----------------------------------------------------------------------#
