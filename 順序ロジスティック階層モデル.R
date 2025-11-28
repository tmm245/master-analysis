############################################################
# 順序ロジスティック階層モデル（num_vaccination の縦断GLMM）
#
# 前提：
#   - merged_data_complete が環境に存在していること
#   - 主な列：
#       id
#       num_vaccination_1 ～ num_vaccination_4
#       Effectiveness_1 ～ Effectiveness_4
#       Safeness_1      ～ Safeness_4
#       HBM_average_benefit_3, HBM_average_benefit_4（あれば）
#       HBM_average_risk_3,    HBM_average_risk_4（あれば）
#       HBM_average_acceptance_3, HBM_average_acceptance_4（あれば）
#       mean_cliticalThinking_attitude or mean_criticalThinking_attitude
#       berlin_correct_count, mean_subjective_numeracy,
#       health_correct_count, numAttitude_average
#       age, gender, education
#
# やること：
#   A) num_vaccination_t を Eff/Safe_lag + CT で予測する
#      階層順序ロジスティックモデル（clmm, (1|id)）
#   B) t4 時点だけを使って num_vaccination_4 を
#      Benefit_4 / Risk_4 / Acceptance_4 と CT で説明する順序ロジスティック（clm）
############################################################

## パッケージ ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ordinal)  # clmm, clm
library(rlang)

## 0. 前提チェック -----------------------------------------------------------

if (!exists("merged_data_complete")) {
  stop("オブジェクト 'merged_data_complete' が環境にありません。先に読み込んでください。")
}

if (!("id" %in% names(merged_data_complete))) {
  stop("'id' 列が merged_data_complete に存在しません。")
}

## 1. CT列名のゆれを処理 ----------------------------------------------------

ct_candidates <- c("mean_cliticalThinking_attitude",
                   "mean_criticalThinking_attitude")

ct_name <- ct_candidates[ct_candidates %in% names(merged_data_complete)][1]

if (is.na(ct_name) || is.null(ct_name)) {
  stop("CT列 (mean_cliticalThinking_attitude / mean_criticalThinking_attitude) が見つかりません。")
}

cat("CT列として使用する変数:", ct_name, "\n")

## 2. ロング形式データの作成 ------------------------------------------------
# num_vaccination_* / Effectiveness_* / Safeness_* / HBM_*_* をロングにする

# タイムインバリアント（個人特性・属性）
base_vars <- c(
  "id",
  ct_name,
  "berlin_correct_count",
  "mean_subjective_numeracy",
  "health_correct_count",
  "numAttitude_average",
  "age",
  "gender",
  "education"
)
base_vars <- base_vars[base_vars %in% names(merged_data_complete)]

# 時点変化する変数のパターン
timevarying_pattern <- paste0(
  "^(num_vaccination",
  "|Effectiveness",
  "|Safeness",
  "|HBM_average_benefit",
  "|HBM_average_risk",
  "|HBM_average_acceptance)",
  "_\\d+$"
)

# ロング化
vacc_long <- merged_data_complete %>%
  dplyr::select(all_of(base_vars), dplyr::matches(timevarying_pattern)) %>%
  dplyr::rename(CT = !!sym(ct_name)) %>%
  tidyr::pivot_longer(
    cols = dplyr::matches(timevarying_pattern),
    names_to = c(".value", "wave"),
    names_pattern = "^(num_vaccination|Effectiveness|Safeness|HBM_average_benefit|HBM_average_risk|HBM_average_acceptance)_(\\d+)$"
  ) %>%
  dplyr::mutate(
    wave = as.integer(wave)
  ) %>%
  dplyr::arrange(id, wave)

cat("ロング化後の次元:", paste(dim(vacc_long), collapse = " x "), "\n")
cat("wave:", paste(sort(unique(vacc_long$wave)), collapse = ", "), "\n")

## 3. ラグ変数の作成（Effectiveness, Safeness） -----------------------------

vacc_long <- vacc_long %>%
  dplyr::group_by(id) %>%
  dplyr::arrange(wave, .by_group = TRUE) %>%
  dplyr::mutate(
    Effectiveness_lag = dplyr::lag(Effectiveness, 1),
    Safeness_lag      = dplyr::lag(Safeness, 1)
  ) %>%
  dplyr::ungroup()

## 4A. 階層順序ロジスティック（Eff/Safe_lag + CT + wave） -----------------

# wave >= 2 で、必要な変数が欠損でない行のみ
model_data_core <- vacc_long %>%
  dplyr::filter(
    wave >= 2,
    !is.na(num_vaccination),
    !is.na(Effectiveness_lag),
    !is.na(Safeness_lag),
    !is.na(CT)
  ) %>%
  dplyr::mutate(
    num_vaccination = factor(num_vaccination, ordered = TRUE),
    wave_factor     = factor(wave, levels = sort(unique(wave))),
    CT_c      = as.numeric(scale(CT, center = TRUE, scale = FALSE)),
    Eff_lag_c = as.numeric(scale(Effectiveness_lag, center = TRUE, scale = FALSE)),
    Safe_lag_c= as.numeric(scale(Safeness_lag,      center = TRUE, scale = FALSE))
  )

cat("model_data_core 次元:", paste(dim(model_data_core), collapse = " x "), "\n")
cat("wave_factor 水準:", paste(levels(model_data_core$wave_factor), collapse = ", "), "\n")

# --- モデル A1: 交互作用なし（コアの安定モデル） ---
m_ord_core_simple <- clmm(
  num_vaccination ~ Eff_lag_c + Safe_lag_c + wave_factor + (1 | id),
  data  = model_data_core,
  link  = "logit",
  Hess  = TRUE,
  nAGQ  = 5
)

cat("\n==== モデル A1: num_vacc_t ~ Eff_lag + Safe_lag + wave + (1|id) ====\n")
print(summary(m_ord_core_simple))

# --- モデル A2: Eff_lag × CT の交互作用を入れたモデル ---
m_ord_core_int <- clmm(
  num_vaccination ~ Eff_lag_c * CT_c + Safe_lag_c + wave_factor + (1 | id),
  data  = model_data_core,
  link  = "logit",
  Hess  = TRUE,
  nAGQ  = 5
)

cat("\n==== モデル A2: num_vacc_t ~ Eff_lag * CT + Safe_lag + wave + (1|id) ====\n")
print(summary(m_ord_core_int))

## 4B. t4時点だけの順序ロジスティック（Benefit/Risk/Acceptance × CT） ----

# HBM 系は t3,t4 しかない前提なので、t4だけで最終回数を説明するモデル

model_data_t4 <- vacc_long %>%
  dplyr::filter(
    wave == 4,
    !is.na(num_vaccination),
    !is.na(HBM_average_benefit),
    !is.na(HBM_average_risk),
    !is.na(HBM_average_acceptance),
    !is.na(CT)
  ) %>%
  dplyr::mutate(
    num_vaccination = factor(num_vaccination, ordered = TRUE),
    CT_c         = as.numeric(scale(CT, center = TRUE, scale = FALSE)),
    Benefit_4_c  = as.numeric(scale(HBM_average_benefit,    center = TRUE, scale = FALSE)),
    Risk_4_c     = as.numeric(scale(HBM_average_risk,       center = TRUE, scale = FALSE)),
    Accept_4_c   = as.numeric(scale(HBM_average_acceptance, center = TRUE, scale = FALSE))
  )

cat("\nmodel_data_t4 次元:", paste(dim(model_data_t4), collapse = " x "), "\n")

# --- モデル B: num_vaccination_4 ~ Benefit_4 * CT + Risk_4 + Accept_4 ---
m_t4 <- clm(
  num_vaccination ~ Benefit_4_c * CT_c + Risk_4_c + Accept_4_c,
  data = model_data_t4,
  link = "logit",
  Hess = TRUE
)

cat("\n==== モデル B: t4 の順序ロジスティック ====\n")
print(summary(m_t4))

############################################################
# ここまでで得られるもの：
#   - m_ord_core_simple:
#       num_vacc_t ~ Eff_lag_c + Safe_lag_c + wave_factor + (1|id)
#   - m_ord_core_int:
#       num_vacc_t ~ Eff_lag_c * CT_c + Safe_lag_c + wave_factor + (1|id)
#   - m_t4:
#       num_vaccination_4 ~ Benefit_4_c * CT_c + Risk_4_c + Accept_4_c
#
# A系が「縦断×階層での接種回数モデル」、
# B系が「最終時点の接種回数と態度の関係」を見るクロスセクションモデル。
############################################################

############################################################
# Acceptance（HBM_average_acceptance_3 / 4）を目的変数にした LMM
#
# 目的：
#   t3 / t4 の Acceptance を、
#   1時点前の Effectiveness / Safeness と CT で予測する
#
#   Acceptance_t ~ Effectiveness_(t-1) * CT + Safeness_(t-1) + wave + (1|id)
############################################################

library(dplyr)
library(lme4)
library(lmerTest)

## 0. CT列名のゆれ処理 ---------------------------------------

ct_candidates <- c("mean_cliticalThinking_attitude",
                   "mean_criticalThinking_attitude")
ct_name <- ct_candidates[ct_candidates %in% names(merged_data_complete)][1]

if (is.na(ct_name)) {
  stop("CT列 (mean_cliticalThinking_attitude / mean_criticalThinking_attitude) が見つかりません。")
}
cat("CT列として使用する変数:", ct_name, "\n")

## 1. t3 / t4 の Acceptance 用の long データを作成 ----------------

mk_wave_acc <- function(df, wave, ct_name){
  if (wave == "t3") {
    dplyr::transmute(
      df,
      id,
      wave = factor("t3", levels = c("t3","t4")),
      CT   = .data[[ct_name]],
      Effectiveness_lag = Effectiveness_2,
      Safeness_lag      = Safeness_2,
      Acceptance        = HBM_average_acceptance_3
    )
  } else { # wave == "t4"
    dplyr::transmute(
      df,
      id,
      wave = factor("t4", levels = c("t3","t4")),
      CT   = .data[[ct_name]],
      Effectiveness_lag = Effectiveness_3,
      Safeness_lag      = Safeness_3,
      Acceptance        = HBM_average_acceptance_4
    )
  }
}

acc_t3 <- mk_wave_acc(merged_data_complete, "t3", ct_name)
acc_t4 <- mk_wave_acc(merged_data_complete, "t4", ct_name)

acc_long <- dplyr::bind_rows(acc_t3, acc_t4) %>%
  # 欠損のある行を落とす
  dplyr::filter(
    !is.na(id),
    !is.na(CT),
    !is.na(Effectiveness_lag),
    !is.na(Safeness_lag),
    !is.na(Acceptance)
  ) %>%
  dplyr::mutate(
    # 中心化（LMMの解釈をしやすくするため）
    CT_c   = as.numeric(scale(CT, center = TRUE, scale = FALSE)),
    Eff_c  = as.numeric(scale(Effectiveness_lag, center = TRUE, scale = FALSE)),
    Safe_c = as.numeric(scale(Safeness_lag,      center = TRUE, scale = FALSE))
  )

cat("acc_long の次元:", paste(dim(acc_long), collapse = " x "), "\n")
table(acc_long$wave)

## 2. Acceptance の LMM を推定 ---------------------------------

# コアモデル：Eff/Safe の主効果のみ
m_acc_core <- lmer(
  Acceptance ~ Eff_c + Safe_c + wave + (1 | id),
  data = acc_long
)
summary(m_acc_core)

# 交互作用（Eff × CT）を入れたモデル
m_acc_int <- lmer(
  Acceptance ~ Eff_c * CT_c + Safe_c + wave + (1 | id),
  data = acc_long
)
summary(m_acc_int)

# wave_factor の水準確認（今は "3","4" になっているはず）
levels(model_data_core$wave_factor)
# [1] "3" "4"

# # 時点間の差（Eff/Safeの傾きがt3とt4で違うか）を検討するモデル
# m_ord_timeint <- clmm(
#   num_vaccination ~ 
#     Eff_lag_c * wave_factor +   # Eff × wave（傾きの時点差）
#     Safe_lag_c * wave_factor +  # Safe × wave
#     (1 | id),
#   data = model_data_core,
#   link = "logit",
#   Hess = TRUE,
#   nAGQ = 5
# )
# 
# summary(m_ord_timeint)

############################################################
# 4C. Benefit / Risk だけを説明変数にした
#     階層順序ロジスティックモデル（t3 + t4）
#
# 目的：
#   num_vaccination_t3/t4 を
#   同時点の Benefit_t, Risk_t だけで予測するモデル：
#
#   num_vaccination ~ Ben_c_only + Risk_c_only + wave_factor + (1 | id)
#
#   → 態度だけで行動をどこまで説明できるかを見る
############################################################

library(dplyr)
library(ordinal)  # clmm

# long_data がすでにある前提（t3/t4のロング）
# 念のため整形＆欠損除去
ord_BR_data <- long_data %>%
  dplyr::filter(
    !is.na(num_vaccination),
    !is.na(Benefit),
    !is.na(Risk),
    !is.na(id),
    !is.na(wave)
  ) %>%
  dplyr::mutate(
    # num_vaccination を順序カテゴリとして扱う
    num_vaccination = factor(num_vaccination, ordered = TRUE),
    # Benefit/Risk を平均中心化
    Ben_c_only  = as.numeric(scale(Benefit, center = TRUE, scale = FALSE)),
    Risk_c_only = as.numeric(scale(Risk,    center = TRUE, scale = FALSE)),
    # wave を因子（t3, t4）に
    wave_factor = factor(wave, levels = c("t3", "t4"))
  )

cat("ord_BR_data 次元:", paste(dim(ord_BR_data), collapse = " x "), "\n")
print(table(ord_BR_data$wave_factor, useNA = "ifany"))

# Benefit / Risk だけを説明変数とする階層順序ロジスティックモデル
m_ord_BR <- clmm(
  num_vaccination ~ Ben_c_only + Risk_c_only + wave_factor + (1 | id),
  data  = ord_BR_data,
  link  = "logit",
  Hess  = TRUE,
  nAGQ  = 5
)

cat("\n==== モデル BR: num_vacc_t ~ Benefit + Risk + wave + (1|id) ====\n")
summary(m_ord_BR)

############################################################
# ord_BR_data を使って、
# Benefit / Risk と 接種回数(num_vaccination) の関係を可視化する ggplot コード
############################################################

library(dplyr)
library(ggplot2)

# 1. 可視化用の数値版 接種回数を作る --------------------------------------
ord_BR_plot <- ord_BR_data %>%
  mutate(
    num_vacc_num = as.numeric(num_vaccination),  # 1,2,3,4,5...
    wave_factor = factor(wave_factor, levels = c("t3","t4"))
  )

# 2. Benefit の水準ごとに「平均接種回数」を見るラインプロット -----------
#   → Benefit が上がると num_vacc_num がどれくらい単調に増えるかを見る

# Benefit をビン分け（例: 5分位）して、その中での平均を計算
ben_summary <- ord_BR_plot %>%
  group_by(wave_factor) %>%
  mutate(Ben_bin = ntile(Benefit, 5)) %>%  # 5分位に分ける
  group_by(wave_factor, Ben_bin) %>%
  summarise(
    Ben_mean       = mean(Benefit, na.rm = TRUE),
    num_vacc_mean  = mean(num_vacc_num, na.rm = TRUE),
    num_vacc_sd    = sd(num_vacc_num, na.rm = TRUE),
    n              = n(),
    num_vacc_se    = num_vacc_sd / sqrt(n),
    .groups = "drop"
  )

# プロット：Benefit（横軸） vs 平均接種回数（縦軸）
p_benefit <- ggplot(ben_summary,
                    aes(x = Ben_mean, y = num_vacc_mean,
                        color = wave_factor, group = wave_factor)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = num_vacc_mean - num_vacc_se,
                    ymax = num_vacc_mean + num_vacc_se),
                width = 0) +
  labs(
    x = "Benefit (HBM, 平均値のビン)",
    y = "平均接種回数（カテゴリの数値化）",
    color = "Wave",
    title = "Benefit の水準と接種回数の関係（t3/t4）"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme_bw()

p_benefit


# 3. Risk 版（Risk が高いほど接種回数が下がるか） -----------------------

risk_summary <- ord_BR_plot %>%
  group_by(wave_factor) %>%
  mutate(Risk_bin = ntile(Risk, 5)) %>%      # 5分位に分ける
  group_by(wave_factor, Risk_bin) %>%
  summarise(
    Risk_mean      = mean(Risk, na.rm = TRUE),
    num_vacc_mean  = mean(num_vacc_num, na.rm = TRUE),
    num_vacc_sd    = sd(num_vacc_num, na.rm = TRUE),
    n              = n(),
    num_vacc_se    = num_vacc_sd / sqrt(n),
    .groups = "drop"
  )

p_risk <- ggplot(risk_summary,
                 aes(x = Risk_mean, y = num_vacc_mean,
                     color = wave_factor, group = wave_factor)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = num_vacc_mean - num_vacc_se,
                    ymax = num_vacc_mean + num_vacc_se),
                width = 0) +
  labs(
    x = "Risk (HBM, 平均値のビン)",
    y = "平均接種回数（カテゴリの数値化）",
    color = "Wave",
    title = "Risk の水準と接種回数の関係（t3/t4）"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme_bw()

p_risk


# 4. 2次元（Benefit × Risk）の中で接種カテゴリがどう分かれているか ------
#    → ほぼ決め打ちで分かれているイメージを見る用

# 元データの整理（さっきと同じ）
ord_BR_plot <- ord_BR_data %>%
  mutate(
    num_vacc_num = as.numeric(num_vaccination),  # 1,2,3,4,5...
    wave_factor  = factor(wave_factor, levels = c("t3", "t4"))
  )

# ---- ここがポイント：Benefit / Risk のビン数を「ユニーク値数以下」にする ----

library(dplyr)
library(ggplot2)

# まずは ord_BR_plot を作っておく（すでにあればこのブロックは不要）
ord_BR_plot <- ord_BR_data %>%
  mutate(
    num_vacc_num = as.numeric(num_vaccination),
    wave_factor  = factor(wave_factor, levels = c("t3", "t4"))
  )

# ---- Benefit / Risk をビン分けして多数派カテゴリを計算 ----

# ビン数は最大でも5に抑える
n_ben  <- min(5, dplyr::n_distinct(ord_BR_plot$Benefit))
n_risk <- min(5, dplyr::n_distinct(ord_BR_plot$Risk))

cat("Benefit bins:", n_ben, " / Risk bins:", n_risk, "\n")

region_df <- ord_BR_plot %>%
  mutate(
    Ben_bin  = ggplot2::cut_number(Benefit, n = n_ben),
    Risk_bin = ggplot2::cut_number(Risk,    n = n_risk)
  ) %>%
  group_by(wave_factor, Ben_bin, Risk_bin) %>%
  summarise(
    Ben_mid  = mean(Benefit, na.rm = TRUE),
    Risk_mid = mean(Risk,    na.rm = TRUE),
    majority_vacc = names(which.max(table(num_vaccination))),
    n = n(),
    .groups = "drop"
  )

head(region_df)

# ---- 大きな四角い点で塗るバージョンのプロット ----

p_br_regions_big <- ggplot(region_df,
                           aes(x = Ben_mid, y = Risk_mid,
                               fill = majority_vacc)) +
  # 四角い点（shape 22）を大きめに描く
  geom_point(shape = 22, size = 10, colour = "grey20") +
  facet_wrap(~ wave_factor) +
  scale_fill_brewer(
    palette = "Set1",
    name = "多数派の接種回数カテゴリ"
  ) +
  labs(
    x = "Benefit（ビン中心）",
    y = "Risk（ビン中心）",
    title = "Benefit × Risk 空間における多数派の接種回数カテゴリ（t3/t4）"
  ) +
  coord_fixed() +
  theme_bw()

p_br_regions_big

############################################################

# 1. まずは、ランダム効果なしの ordinal logit をフィット（全体の決定境界用）
#    ※wave も入れておく（t3/t4でシフトが違う可能性）
m_br_clm <- clm(
  num_vaccination ~ Ben_c_only + Risk_c_only + wave_factor,
  data = ord_BR_data,
  link = "logit"
)

summary(m_br_clm)

# 2. Benefit × Risk のグリッドを作成 -------------------------------------

# 描画範囲（少し余裕を持たせるなら extendrange を使ってもOK）
ben_range  <- range(ord_BR_data$Benefit, na.rm = TRUE)
risk_range <- range(ord_BR_data$Risk,    na.rm = TRUE)

ben_seq  <- seq(ben_range[1],  ben_range[2],  length.out = 80)
risk_seq <- seq(risk_range[1], risk_range[2], length.out = 80)

# 中心化に使った平均
ben_mean  <- mean(ord_BR_data$Benefit, na.rm = TRUE)
risk_mean <- mean(ord_BR_data$Risk,    na.rm = TRUE)

grid <- expand.grid(
  Benefit    = ben_seq,
  Risk       = risk_seq,
  wave_factor = levels(ord_BR_data$wave_factor)
)

grid <- grid %>%
  mutate(
    Ben_c_only  = Benefit - ben_mean,
    Risk_c_only = Risk    - risk_mean
  )

# 3. 各グリッド点で、カテゴリごとの予測確率を計算 ------------------------

prob_mat <- predict(m_br_clm, newdata = grid, type = "prob")

# predict() の戻りは "list" のこともあるので対処
if (is.list(prob_mat)) {
  prob_mat <- prob_mat$fit
}

# 列名はカテゴリラベル（num_vaccination の levels）になっているはず
cat("列名（カテゴリ）:", colnames(prob_mat), "\n")

# 一番確率が高いカテゴリを決定境界として採用
max_col <- max.col(prob_mat, ties.method = "first")
grid$pred_class <- factor(
  colnames(prob_mat)[max_col],
  levels = levels(ord_BR_data$num_vaccination)
)

# 4. 実データも可視化用に準備 --------------------------------------------

ord_BR_plot <- ord_BR_data %>%
  mutate(
    wave_factor = factor(wave_factor, levels = c("t3","t4"))
  )

# 5. 図：背景に「予測カテゴリ領域」、上に実データの散布図 ---------------

p_decision <- ggplot() +
  # 背景：ロジスティック回帰が決めた「カテゴリ領域」
  geom_raster(
    data = grid,
    aes(x = Benefit, y = Risk, fill = pred_class),
    alpha = 0.35
  ) +
  # ざっくり境界線（カテゴリを数値化して等高線として引く）
  geom_contour(
    data = grid,
    aes(x = Benefit, y = Risk, z = as.numeric(pred_class)),
    colour = "black",
    linewidth = 0.3
  ) +
  # 実データの点
  geom_point(
    data = ord_BR_plot,
    aes(x = Benefit, y = Risk, color = num_vaccination),
    size = 1.3,
    alpha = 0.8
  ) +
  facet_wrap(~ wave_factor) +
  scale_fill_brewer(
    palette = "Pastel1",
    name = "モデル予測カテゴリ\n(num_vaccination)"
  ) +
  scale_color_brewer(
    palette = "Set1",
    name = "観測カテゴリ\n(num_vaccination)"
  ) +
  labs(
    x = "Benefit",
    y = "Risk",
    title = "Benefit × Risk に基づく接種回数カテゴリの決定境界\n（ordinal logit の予測）"
  ) +
  coord_fixed() +
  theme_bw()

p_decision

############################################################
# t3 の Benefit / Risk が t4 の接種回数カテゴリを予測する
# クロスラグ ordinal logit ＋ 決定境界プロット
############################################################


## 1. t3 の Benefit/Risk と t4 の num_vaccination を1行にまとめる -----

# t3 側（predictor）
t3_BR <- long_data %>%
  dplyr::filter(wave == "t3") %>%
  dplyr::select(
    id,
    Benefit3 = Benefit,
    Risk3    = Risk
  )

# t4 側（outcome）
t4_vacc <- long_data %>%
  dplyr::filter(wave == "t4") %>%
  dplyr::select(
    id,
    num_vacc_t4 = num_vaccination
  )

# id で結合（t3の態度 + t4 の接種）
crosslag_dat <- t3_BR %>%
  dplyr::inner_join(t4_vacc, by = "id") %>%
  dplyr::filter(
    !is.na(Benefit3),
    !is.na(Risk3),
    !is.na(num_vacc_t4)
  ) %>%
  dplyr::mutate(
    # 目的変数を順序カテゴリに
    num_vacc_t4 = factor(num_vacc_t4, ordered = TRUE)
  )

cat("crosslag_dat 次元:", paste(dim(crosslag_dat), collapse = " x "), "\n")
table(crosslag_dat$num_vacc_t4)

# 中心化（t3 の Benefit / Risk）
ben3_mean  <- mean(crosslag_dat$Benefit3, na.rm = TRUE)
risk3_mean <- mean(crosslag_dat$Risk3,    na.rm = TRUE)

crosslag_dat <- crosslag_dat %>%
  dplyr::mutate(
    Ben3_c  = Benefit3 - ben3_mean,
    Risk3_c = Risk3    - risk3_mean
  )

## 2. クロスラグ ordinal logit モデル ------------------------------

m_cross_BR <- clm(
  num_vacc_t4 ~ Ben3_c + Risk3_c,
  data = crosslag_dat,
  link = "logit"
)

summary(m_cross_BR)

## 3. Benefit3 × Risk3 のグリッドを作って予測 ----------------------

ben_seq  <- seq(min(crosslag_dat$Benefit3), max(crosslag_dat$Benefit3), length.out = 80)
risk_seq <- seq(min(crosslag_dat$Risk3),    max(crosslag_dat$Risk3),    length.out = 80)

grid_cross <- expand.grid(
  Benefit3 = ben_seq,
  Risk3    = risk_seq
) %>%
  dplyr::mutate(
    Ben3_c  = Benefit3 - ben3_mean,
    Risk3_c = Risk3    - risk3_mean
  )

# 各グリッド点でカテゴリ別予測確率を計算
prob_mat_cross <- predict(
  m_cross_BR,
  newdata = grid_cross,
  type = "prob"
)

if (is.list(prob_mat_cross)) {
  prob_mat_cross <- prob_mat_cross$fit
}

cat("列名（カテゴリ）:", colnames(prob_mat_cross), "\n")

# 最尤カテゴリを予測クラスとして付与
max_col_cross <- max.col(prob_mat_cross, ties.method = "first")

grid_cross$pred_class <- factor(
  colnames(prob_mat_cross)[max_col_cross],
  levels = levels(crosslag_dat$num_vacc_t4)
)

## 4. プロット：t3 の Benefit×Risk 空間で、t4 接種カテゴリの決定境界 -----

p_cross_decision <- ggplot() +
  # 背景：モデルが予測するカテゴリ領域
  geom_raster(
    data = grid_cross,
    aes(x = Benefit3, y = Risk3, fill = pred_class),
    alpha = 0.35
  ) +
  # カテゴリが切り替わる境界線
  geom_contour(
    data = grid_cross,
    aes(x = Benefit3, y = Risk3, z = as.numeric(pred_class)),
    colour = "black",
    linewidth = 0.3
  ) +
  # 実データ：t3 Benefit/Risk と t4 接種カテゴリ
  geom_point(
    data = crosslag_dat,
    aes(x = Benefit3, y = Risk3, color = num_vacc_t4),
    size = 1.5,
    alpha = 0.8
  ) +
  scale_fill_brewer(
    palette = "Pastel1",
    name = "モデル予測カテゴリ\n(num_vacc_t4)"
  ) +
  scale_color_brewer(
    palette = "Set1",
    name = "観測カテゴリ\n(num_vacc_t4)"
  ) +
  labs(
    x = "Benefit at t3",
    y = "Risk at t3",
    title = "t3 Benefit × Risk に基づく t4 接種回数カテゴリの決定境界"
  ) +
  coord_fixed() +
  theme_bw()

p_cross_decision
############################################################

## 3. グリッド（Benefit3 × Risk3 × CT_level） ------------------------

ben_seq  <- seq(min(crosslag_ct$Benefit3), max(crosslag_ct$Benefit3), length.out = 80)
risk_seq <- seq(min(crosslag_ct$Risk3),    max(crosslag_ct$Risk3),    length.out = 80)

CT_levels <- c("Low", "Mean", "High")
CT_c_values <- c(-CT_sd, 0, CT_sd)   # CT_c = -1SD, 0, +1SD

grid_list <- list()

for (i in seq_along(CT_levels)) {
  this_level <- CT_levels[i]
  this_CT_c  <- CT_c_values[i]
  
  grid_i <- expand.grid(
    Benefit3 = ben_seq,
    Risk3    = risk_seq
  ) %>%
    mutate(
      Ben3_c = Benefit3 - ben3_mean,
      Risk3_c = Risk3   - risk3_mean,
      CT_c    = this_CT_c,
      CT_level = this_level
    )
  
  grid_list[[i]] <- grid_i
}

grid_ct <- bind_rows(grid_list)

## 4. 予測確率 → 予測カテゴリ -------------------------

prob_mat_ct <- predict(
  m_cross_BR_CT,
  newdata = grid_ct,
  type = "prob"
)

if (is.list(prob_mat_ct)) {
  prob_mat_ct <- prob_mat_ct$fit
}

max_col_ct <- max.col(prob_mat_ct, ties.method = "first")

grid_ct$pred_class <- factor(
  colnames(prob_mat_ct)[max_col_ct],
  levels = levels(crosslag_ct$num_vacc_t4)
)

## 5. プロット：CTごとの決定境界 -------------------------

p_cross_CT <- ggplot() +
  # 背景（モデル予測カテゴリ）
  geom_raster(
    data = grid_ct,
    aes(x = Benefit3, y = Risk3, fill = pred_class),
    alpha = 0.35
  ) +
  # カテゴリ境界（1.5, 2.5, 3.5 付近）
  geom_contour(
    data = grid_ct,
    aes(x = Benefit3, y = Risk3, z = as.numeric(pred_class)),
    breaks = c(1.5, 2.5, 3.5),
    colour = "black",
    linewidth = 0.4
  ) +
  # 実データ（色は t4 接種カテゴリ）
  geom_point(
    data = crosslag_ct,
    aes(x = Benefit3, y = Risk3, color = num_vacc_t4),
    size = 1.3,
    alpha = 0.8
  ) +
  facet_wrap(~ CT_level) +
  scale_fill_brewer(
    palette = "Pastel1",
    name = "モデル予測カテゴリ\n(num_vacc_t4)"
  ) +
  scale_color_brewer(
    palette = "Set1",
    name = "観測カテゴリ\n(num_vacc_t4)"
  ) +
  labs(
    x = "Benefit at t3",
    y = "Risk at t3",
    title = "t3 Benefit × Risk に基づく t4 接種回数カテゴリの決定境界\nCT水準ごとの比較"
  ) +
  coord_fixed() +
  theme_bw()

p_cross_CT

#----------------0 vs 2+ の二値ロジスティック回帰コード--------------

library(dplyr)

## 1. データ整形：0 vs 2+ ----
logit_BR_t3t4 <- merged_data_complete %>%
  transmute(
    id,
    Benefit_t3 = HBM_average_benefit_3,
    Risk_t3    = HBM_average_risk_3,
    vacc4      = num_vaccination_4
  ) %>%
  mutate(
    # 0 = 全く打っていない, 2plus = 2回以上（2,3,4,5）
    vacc4_bin = case_when(
      vacc4 == 0        ~ "0",
      vacc4 >= 2        ~ "2plus",
      TRUE              ~ NA_character_
    ),
    vacc4_bin = factor(vacc4_bin, levels = c("0", "2plus")),  # 参照カテゴリは「0」
    
    # 予測子は grand-mean センタリング
    Benefit_t3_c = scale(Benefit_t3, center = TRUE, scale = FALSE)[, 1],
    Risk_t3_c    = scale(Risk_t3,    center = TRUE, scale = FALSE)[, 1]
  ) %>%
  tidyr::drop_na(vacc4_bin, Benefit_t3_c, Risk_t3_c)

## 2. ロジスティック回帰（0 vs 2+） ----
m_logit_BR_t3t4 <- glm(
  vacc4_bin ~ Benefit_t3_c + Risk_t3_c,
  data   = logit_BR_t3t4,
  family = binomial(link = "logit")
)

summary(m_logit_BR_t3t4)

## 3. オッズ比と95%CI ----
coefs <- summary(m_logit_BR_t3t4)$coefficients

OR  <- exp(coef(m_logit_BR_t3t4))
CI  <- exp(confint(m_logit_BR_t3t4))  # プロファイルCI

OR_table <- cbind(OR = OR, CI)
OR_table

#----------------0 vs 2+ の二値ロジスティック回帰コード--------------
## 必要パッケージ ----
library(dplyr)
library(tidyr)
library(MASS)   # polr()

## 1. t3 Benefit / Risk と t4 接種回数を4カテゴリ(0 / 2 / 3 / 4+)に整形 ----
ord_BR_t3t4_A <- merged_data_complete %>%
  transmute(
    id,
    Benefit_t3 = HBM_average_benefit_3,
    Risk_t3    = HBM_average_risk_3,
    vacc4      = num_vaccination_4
  ) %>%
  mutate(
    # 0 / 2 / 3 / 4+ の4カテゴリ
    vacc4_catA = case_when(
      vacc4 == 0 ~ "0",
      vacc4 == 2 ~ "2",
      vacc4 == 3 ~ "3",
      vacc4 >= 4 ~ "4+",
      TRUE       ~ NA_character_
    ),
    vacc4_catA = factor(
      vacc4_catA,
      levels  = c("0", "2", "3", "4+"),
      ordered = TRUE
    ),
    # 予測子を grand-mean センタリング
    Benefit_t3_c = scale(Benefit_t3, center = TRUE, scale = FALSE)[, 1],
    Risk_t3_c    = scale(Risk_t3,    center = TRUE, scale = FALSE)[, 1]
  ) %>%
  drop_na(vacc4_catA, Benefit_t3_c, Risk_t3_c)

# カテゴリ数の確認（任意）
table(ord_BR_t3t4_A$vacc4_catA)

## 2. 順序ロジスティック回帰（A系 0/2/3/4+） ----
m_ord_BR_t3t4_A <- polr(
  vacc4_catA ~ Benefit_t3_c + Risk_t3_c,
  data = ord_BR_t3t4_A,
  Hess = TRUE
)

summary(m_ord_BR_t3t4_A)

## 3. オッズ比と95%CI、p値の算出 ----

# 係数と t 値
(ct_A <- coef(summary(m_ord_BR_t3t4_A)))

# p値（正規近似）
p_values_A <- 2 * pnorm(abs(ct_A[, "t value"]), lower.tail = FALSE)

# 予測子（Benefit / Risk）のみの p を見たい場合
cbind(ct_A[1:2, ], p = p_values_A[1:2])

# オッズ比と95%CI
OR_A <- exp(coef(m_ord_BR_t3t4_A))
CI_A <- exp(confint(m_ord_BR_t3t4_A))   # プロファイルCI

OR_table_A <- cbind(OR = OR_A, CI_A)
OR_table_A

#-------------------------------------------------------------------#
library(dplyr)
library(MASS)   # polr

## 1-1. 4カテゴリ A 系データ（もうあるならここは不要） ----
ord_BR_t3t4_A <- merged_data_complete %>%
  transmute(
    id,
    Benefit_t3 = HBM_average_benefit_3,
    Risk_t3    = HBM_average_risk_3,
    vacc4      = num_vaccination_4
  ) %>%
  mutate(
    vacc4_catA = case_when(
      vacc4 == 0        ~ "0",
      vacc4 == 2        ~ "2",
      vacc4 == 3        ~ "3",
      vacc4 >= 4        ~ "4plus",
      TRUE              ~ NA_character_
    ),
    vacc4_catA = factor(vacc4_catA,
                        ordered = TRUE,
                        levels = c("0","2","3","4plus")),
    Benefit_t3_c = scale(Benefit_t3, center = TRUE, scale = FALSE)[,1],
    Risk_t3_c    = scale(Risk_t3,    center = TRUE, scale = FALSE)[,1]
  ) %>%
  tidyr::drop_na(vacc4_catA, Benefit_t3_c, Risk_t3_c)

m_ord_BR_t3t4_A <- polr(
  vacc4_catA ~ Benefit_t3_c + Risk_t3_c,
  data = ord_BR_t3t4_A,
  Hess = TRUE
)

## 1-2. 3カテゴリ版（0,2,3plus）を作る ----
ord_BR_t3t4_3cat <- merged_data_complete %>%
  transmute(
    id,
    Benefit_t3 = HBM_average_benefit_3,
    Risk_t3    = HBM_average_risk_3,
    vacc4      = num_vaccination_4
  ) %>%
  mutate(
    vacc4_cat3 = case_when(
      vacc4 == 0        ~ "0",
      vacc4 == 2        ~ "2",
      vacc4 >= 3        ~ "3plus",
      TRUE              ~ NA_character_
    ),
    vacc4_cat3 = factor(vacc4_cat3,
                        ordered = TRUE,
                        levels = c("0","2","3plus")),
    Benefit_t3_c = scale(Benefit_t3, center = TRUE, scale = FALSE)[,1],
    Risk_t3_c    = scale(Risk_t3,    center = TRUE, scale = FALSE)[,1]
  ) %>%
  tidyr::drop_na(vacc4_cat3, Benefit_t3_c, Risk_t3_c)

m_ord_BR_t3t4_3cat <- polr(
  vacc4_cat3 ~ Benefit_t3_c + Risk_t3_c,
  data = ord_BR_t3t4_3cat,
  Hess = TRUE
)

## 1-3. 情報量基準の比較 ----
AIC(m_ord_BR_t3t4_A, m_ord_BR_t3t4_3cat)
BIC(m_ord_BR_t3t4_A, m_ord_BR_t3t4_3cat)

#-------------------------------------------------------------------#

############################################################
# 3カテゴリ版（0,2,3plus）の順序ロジスティック
# 目的：
#   vacc4_cat3（0,2,3plus）を従属変数にして、
#   Benefit_t3, Risk_t3, CT, Berlin の効果と
#   Benefit×CT / Benefit×Berlin の交互作用を検討する
#
# カテゴリ定義：
#   0    : num_vaccination_4 == 0
#   2    : num_vaccination_4 == 2
#   3plus: num_vaccination_4 >= 3
############################################################

library(dplyr)
library(tidyr)
library(MASS)  # polr()

## 0. 前提チェック -----------------------------------------------------------

if (!exists("merged_data_complete")) {
  stop("merged_data_complete が環境にありません。先に読み込んでください。")
}

## 1. CT 列名のゆれ処理 ------------------------------------------------------

ct_candidates <- c("mean_cliticalThinking_attitude",
                   "mean_criticalThinking_attitude")
ct_name <- ct_candidates[ct_candidates %in% names(merged_data_complete)][1]

if (is.na(ct_name) || is.null(ct_name)) {
  stop("CT列 (mean_cliticalThinking_attitude / mean_criticalThinking_attitude) が見つかりません。")
}

cat("CT列として使用する変数:", ct_name, "\n")

## 2. 3カテゴリ＋CT＋Berlin を含むデータフレーム ----------------------------

ord_BR_t3t4_3cat_ext <- merged_data_complete %>%
  transmute(
    id,
    Benefit_t3 = HBM_average_benefit_3,
    Risk_t3    = HBM_average_risk_3,
    vacc4      = num_vaccination_4,
    CT_raw     = .data[[ct_name]],
    berlin     = berlin_correct_count
  ) %>%
  mutate(
    # 3カテゴリ：0, 2, 3plus
    vacc4_cat3 = case_when(
      vacc4 == 0        ~ "0",
      vacc4 == 2        ~ "2",
      vacc4 >= 3        ~ "3plus",
      TRUE              ~ NA_character_
    ),
    vacc4_cat3 = factor(
      vacc4_cat3,
      ordered = TRUE,
      levels  = c("0", "2", "3plus")
    ),
    # 中心化（平均のみ引く）
    Benefit_t3_c = scale(Benefit_t3, center = TRUE, scale = FALSE)[, 1],
    Risk_t3_c    = scale(Risk_t3,    center = TRUE, scale = FALSE)[, 1],
    CT_c         = scale(CT_raw,     center = TRUE, scale = FALSE)[, 1],
    berlin_c     = scale(berlin,     center = TRUE, scale = FALSE)[, 1]
  ) %>%
  tidyr::drop_na(vacc4_cat3, Benefit_t3_c, Risk_t3_c, CT_c, berlin_c)

cat("ord_BR_t3t4_3cat_ext 次元:",
    paste(dim(ord_BR_t3t4_3cat_ext), collapse = " x "), "\n")

## 3. ベースライン（主効果のみ）モデル --------------------------------------

m_ord_3cat_main <- polr(
  vacc4_cat3 ~ Benefit_t3_c + Risk_t3_c + CT_c + berlin_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

summary(m_ord_3cat_main)

## 4. Benefit × CT の交互作用モデル -----------------------------------------

m_ord_3cat_CT <- polr(
  vacc4_cat3 ~ Benefit_t3_c * CT_c + Risk_t3_c + berlin_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

summary(m_ord_3cat_CT)

## 5. Benefit × Berlin の交互作用モデル -------------------------------------

m_ord_3cat_berlin <- polr(
  vacc4_cat3 ~ Benefit_t3_c * berlin_c + Risk_t3_c + CT_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

summary(m_ord_3cat_berlin)

## 6. Benefit × CT ＆ Benefit × Berlin 両方を含むモデル ----------------------

m_ord_3cat_CT_berlin <- polr(
  vacc4_cat3 ~ Benefit_t3_c * CT_c + Benefit_t3_c * berlin_c + Risk_t3_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

summary(m_ord_3cat_CT_berlin)

## 7. AIC / BIC によるモデル比較 --------------------------------------------

AIC(m_ord_3cat_main, m_ord_3cat_CT, m_ord_3cat_berlin, m_ord_3cat_CT_berlin)
BIC(m_ord_3cat_main, m_ord_3cat_CT, m_ord_3cat_berlin, m_ord_3cat_CT_berlin)

## 8. 近似的な LR 検定（ネストしたモデル同士のみ） ---------------------------

cat("\n--- LR test: main vs CT interaction ---\n")
print(anova(m_ord_3cat_main, m_ord_3cat_CT))

cat("\n--- LR test: main vs Berlin interaction ---\n")
print(anova(m_ord_3cat_main, m_ord_3cat_berlin))

cat("\n--- LR test: main vs CT+Berlin interaction ---\n")
print(anova(m_ord_3cat_main, m_ord_3cat_CT_berlin))

############################################################
# これで「当てはまりの良かった 3カテゴリ版」をベースに、
#   - 主効果のみ
#   - Benefit×CT
#   - Benefit×Berlin
#   - その両方
# を一気に比較できます。
############################################################

#---------------------------------------------------------------------------#

############################################################
# 3カテゴリ版（0, 2, 3plus）の順序ロジスティック回帰
#  - 従属変数: vacc4_cat3 (0, 2, 3plus)
#  - 予測子: Benefit_t3, Risk_t3, CT, Berlin（すべて平均中心化）
#
# 前提:
#   merged_data_complete:
#     - id
#     - HBM_average_benefit_3
#     - HBM_average_risk_3
#     - num_vaccination_4
#     - mean_cliticalThinking_attitude または mean_criticalThinking_attitude
#     - berlin_correct_count
############################################################

library(dplyr)
library(tidyr)
library(MASS)  # polr()

## 0. 前提チェック -----------------------------------------------------------

if (!exists("merged_data_complete")) {
  stop("merged_data_complete が環境にありません。先に読み込んでください。")
}

## 1. CT 列名のゆれ処理 ------------------------------------------------------

ct_candidates <- c("mean_cliticalThinking_attitude",
                   "mean_criticalThinking_attitude")
ct_name <- ct_candidates[ct_candidates %in% names(merged_data_complete)][1]

if (is.na(ct_name) || is.null(ct_name)) {
  stop("CT列 (mean_cliticalThinking_attitude / mean_criticalThinking_attitude) が見つかりません。")
}

cat("CT列として使用する変数:", ct_name, "\n")

## 2. 3カテゴリ＋CT＋Berlin を含むデータフレーム ----------------------------

ord_BR_t3t4_3cat_ext <- merged_data_complete %>%
  transmute(
    id,
    Benefit_t3 = HBM_average_benefit_3,
    Risk_t3    = HBM_average_risk_3,
    vacc4      = num_vaccination_4,
    CT_raw     = .data[[ct_name]],
    berlin     = berlin_correct_count
  ) %>%
  mutate(
    # 3カテゴリ: 0, 2, 3plus
    vacc4_cat3 = case_when(
      vacc4 == 0        ~ "0",
      vacc4 == 2        ~ "2",
      vacc4 >= 3        ~ "3plus",
      TRUE              ~ NA_character_
    ),
    vacc4_cat3 = factor(
      vacc4_cat3,
      ordered = TRUE,
      levels  = c("0", "2", "3plus")
    ),
    # 平均中心化（scale = FALSE）
    Benefit_t3_c = scale(Benefit_t3, center = TRUE, scale = FALSE)[, 1],
    Risk_t3_c    = scale(Risk_t3,    center = TRUE, scale = FALSE)[, 1],
    CT_c         = scale(CT_raw,     center = TRUE, scale = FALSE)[, 1],
    berlin_c     = scale(berlin,     center = TRUE, scale = FALSE)[, 1]
  ) %>%
  tidyr::drop_na(vacc4_cat3, Benefit_t3_c, Risk_t3_c, CT_c, berlin_c)

cat("ord_BR_t3t4_3cat_ext 次元:",
    paste(dim(ord_BR_t3t4_3cat_ext), collapse = " x "), "\n")

## 3. ベースライン順序ロジスティックモデル（主効果のみ） -------------------

m_ord_3cat_main <- polr(
  vacc4_cat3 ~ Benefit_t3_c + Risk_t3_c + CT_c + berlin_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

summary(m_ord_3cat_main)

## 4. オッズ比（exp(係数)) と 95% CI も見たい場合 ---------------------------

coefs <- coef(summary(m_ord_3cat_main))

# 係数の logit → OR 変換
OR   <- exp(coefs[, "Value"])
SE   <- coefs[, "Std. Error"]
z    <- qnorm(0.975)

# 95% CI（log-odds を OR に戻す）
lower <- exp(coefs[, "Value"] - z * SE)
upper <- exp(coefs[, "Value"] + z * SE)

OR_table <- data.frame(
  term  = rownames(coefs),
  OR    = OR,
  lower = lower,
  upper = upper
)

print(OR_table)

############################################################
# ここまでで：
#  - 0,2,3plus の 3カテゴリ順序ロジスティック
#  - Benefit_t3_c / Risk_t3_c / CT_c / berlin_c の主効果
# が推定される。
#
# 交互作用を入れたいときは、たとえば：
#   vacc4_cat3 ~ Benefit_t3_c * CT_c + Risk_t3_c + berlin_c
# などに置き換えれば OK（先にやった m_ord_3cat_CT 等）。
############################################################

#------------------------------------------------------------------------#

############################################################
# 4カテゴリ (0,2,3,4plus) vs 3カテゴリ (0,2,3plus) の比較
# - 予測子は Benefit_t3_c, Risk_t3_c のみ（CT/berlin は入れない）
# - 係数とオッズ比(OR)を横並びの表にまとめる
############################################################

library(dplyr)
library(tidyr)
library(MASS)  # polr()

## 0. 前提チェック -----------------------------------------------------------

if (!exists("merged_data_complete")) {
  stop("merged_data_complete が環境にありません。先に読み込んでください。")
}

## 1. 4カテゴリ A 系データ（0,2,3,4plus） -----------------------------------

ord_BR_t3t4_A_simple <- merged_data_complete %>%
  transmute(
    id,
    Benefit_t3 = HBM_average_benefit_3,
    Risk_t3    = HBM_average_risk_3,
    vacc4      = num_vaccination_4
  ) %>%
  mutate(
    vacc4_catA = case_when(
      vacc4 == 0        ~ "0",
      vacc4 == 2        ~ "2",
      vacc4 == 3        ~ "3",
      vacc4 >= 4        ~ "4plus",
      TRUE              ~ NA_character_
    ),
    vacc4_catA = factor(
      vacc4_catA,
      ordered = TRUE,
      levels  = c("0","2","3","4plus")
    ),
    Benefit_t3_c = scale(Benefit_t3, center = TRUE, scale = FALSE)[, 1],
    Risk_t3_c    = scale(Risk_t3,    center = TRUE, scale = FALSE)[, 1]
  ) %>%
  drop_na(vacc4_catA, Benefit_t3_c, Risk_t3_c)

cat("4カテゴリデータ次元:",
    paste(dim(ord_BR_t3t4_A_simple), collapse = " x "), "\n")

m_4cat <- polr(
  vacc4_catA ~ Benefit_t3_c + Risk_t3_c,
  data = ord_BR_t3t4_A_simple,
  Hess = TRUE
)

## 2. 3カテゴリデータ（0,2,3plus） ------------------------------------------

ord_BR_t3t4_3cat_simple <- merged_data_complete %>%
  transmute(
    id,
    Benefit_t3 = HBM_average_benefit_3,
    Risk_t3    = HBM_average_risk_3,
    vacc4      = num_vaccination_4
  ) %>%
  mutate(
    vacc4_cat3 = case_when(
      vacc4 == 0        ~ "0",
      vacc4 == 2        ~ "2",
      vacc4 >= 3        ~ "3plus",
      TRUE              ~ NA_character_
    ),
    vacc4_cat3 = factor(
      vacc4_cat3,
      ordered = TRUE,
      levels  = c("0","2","3plus")
    ),
    Benefit_t3_c = scale(Benefit_t3, center = TRUE, scale = FALSE)[, 1],
    Risk_t3_c    = scale(Risk_t3,    center = TRUE, scale = FALSE)[, 1]
  ) %>%
  drop_na(vacc4_cat3, Benefit_t3_c, Risk_t3_c)

cat("3カテゴリデータ次元:",
    paste(dim(ord_BR_t3t4_3cat_simple), collapse = " x "), "\n")

m_3cat <- polr(
  vacc4_cat3 ~ Benefit_t3_c + Risk_t3_c,
  data = ord_BR_t3t4_3cat_simple,
  Hess = TRUE
)

## 3. 係数・OR・95%CI を取り出す関数 -----------------------------------------

extract_or_table <- function(fit, model_name) {
  coefs <- coef(summary(fit))
  z     <- qnorm(0.975)
  
  # 係数行（Benefit, Risk）のみ抽出（cutpoint は除く）
  coefs_slopes <- coefs[rownames(coefs) %in% c("Benefit_t3_c","Risk_t3_c"), , drop = FALSE]
  
  OR    <- exp(coefs_slopes[, "Value"])
  SE    <- coefs_slopes[, "Std. Error"]
  lower <- exp(coefs_slopes[, "Value"] - z * SE)
  upper <- exp(coefs_slopes[, "Value"] + z * SE)
  
  tibble(
    model = model_name,
    term  = rownames(coefs_slopes),
    beta  = coefs_slopes[, "Value"],
    SE    = SE,
    OR    = OR,
    CI_low  = lower,
    CI_high = upper
  )
}

table_4cat <- extract_or_table(m_4cat,  "4cat (0,2,3,4plus)")
table_3cat <- extract_or_table(m_3cat,  "3cat (0,2,3plus)")

## 4. 横並び比較テーブル -----------------------------------------------------

coef_compare <- bind_rows(table_4cat, table_3cat) %>%
  arrange(term, model)

print(coef_compare)

#---------------------------------------------------------------------------#

############################################################
# 3カテゴリ版（0, 2, 3plus）の順序ロジスティック回帰
#  - 従属変数: vacc4_cat3 (0, 2, 3plus)
#  - 予測子: Benefit_t3, Risk_t3, CT, Berlin（すべて平均中心化）
#  - モデル一覧：
#     m_ord_3cat_main         : 主効果のみ
#     m_ord_3cat_BxCT         : Benefit × CT
#     m_ord_3cat_BxBerlin     : Benefit × Berlin
#     m_ord_3cat_BxCT_BxBer   : Benefit × CT & Benefit × Berlin
#     m_ord_3cat_RxCT         : Risk × CT
#     m_ord_3cat_RxBerlin     : Risk × Berlin
#     m_ord_3cat_RxCT_RxBer   : Risk × CT & Risk × Berlin
#     m_ord_3cat_full         : Benefit/Risk × CT/Berlin 全部入り
############################################################

library(dplyr)
library(tidyr)
library(MASS)  # polr()

## 0. 前提チェック -----------------------------------------------------------

if (!exists("merged_data_complete")) {
  stop("merged_data_complete が環境にありません。先に読み込んでください。")
}

## 1. CT 列名のゆれ処理 ------------------------------------------------------

ct_candidates <- c("mean_cliticalThinking_attitude",
                   "mean_criticalThinking_attitude")
ct_name <- ct_candidates[ct_candidates %in% names(merged_data_complete)][1]

if (is.na(ct_name) || is.null(ct_name)) {
  stop("CT列 (mean_cliticalThinking_attitude / mean_criticalThinking_attitude) が見つかりません。")
}

cat("CT列として使用する変数:", ct_name, "\n")

## 2. 3カテゴリ＋CT＋Berlin を含むデータフレーム ----------------------------

ord_BR_t3t4_3cat_ext <- merged_data_complete %>%
  transmute(
    id,
    Benefit_t3 = HBM_average_benefit_3,
    Risk_t3    = HBM_average_risk_3,
    vacc4      = num_vaccination_4,
    CT_raw     = .data[[ct_name]],
    berlin     = berlin_correct_count
  ) %>%
  mutate(
    # 3カテゴリ: 0, 2, 3plus
    vacc4_cat3 = case_when(
      vacc4 == 0        ~ "0",
      vacc4 == 2        ~ "2",
      vacc4 >= 3        ~ "3plus",
      TRUE              ~ NA_character_
    ),
    vacc4_cat3 = factor(
      vacc4_cat3,
      ordered = TRUE,
      levels  = c("0", "2", "3plus")
    ),
    # 平均中心化（scale = FALSE）
    Benefit_t3_c = scale(Benefit_t3, center = TRUE, scale = FALSE)[, 1],
    Risk_t3_c    = scale(Risk_t3,    center = TRUE, scale = FALSE)[, 1],
    CT_c         = scale(CT_raw,     center = TRUE, scale = FALSE)[, 1],
    berlin_c     = scale(berlin,     center = TRUE, scale = FALSE)[, 1]
  ) %>%
  tidyr::drop_na(vacc4_cat3, Benefit_t3_c, Risk_t3_c, CT_c, berlin_c)

cat("ord_BR_t3t4_3cat_ext 次元:",
    paste(dim(ord_BR_t3t4_3cat_ext), collapse = " x "), "\n")

## 3. ベースライン（主効果のみ）モデル --------------------------------------

m_ord_3cat_main <- polr(
  vacc4_cat3 ~ Benefit_t3_c + Risk_t3_c + CT_c + berlin_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

summary(m_ord_3cat_main)

## 4. Benefit 側の交互作用モデル -------------------------------------------

# Benefit × CT
m_ord_3cat_BxCT <- polr(
  vacc4_cat3 ~ Benefit_t3_c * CT_c + Risk_t3_c + berlin_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

summary(m_ord_3cat_BxCT)

# Benefit × Berlin
m_ord_3cat_BxBerlin <- polr(
  vacc4_cat3 ~ Benefit_t3_c * berlin_c + Risk_t3_c + CT_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

summary(m_ord_3cat_BxBerlin)

# Benefit × CT & Benefit × Berlin
m_ord_3cat_BxCT_BxBer <- polr(
  vacc4_cat3 ~ Benefit_t3_c * CT_c + Benefit_t3_c * berlin_c + Risk_t3_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

summary(m_ord_3cat_BxCT_BxBer)

## 5. Risk 側の交互作用モデル ----------------------------------------------

# Risk × CT
m_ord_3cat_RxCT <- polr(
  vacc4_cat3 ~ Benefit_t3_c + Risk_t3_c * CT_c + berlin_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

summary(m_ord_3cat_RxCT)

# Risk × Berlin
m_ord_3cat_RxBerlin <- polr(
  vacc4_cat3 ~ Benefit_t3_c + Risk_t3_c * berlin_c + CT_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

summary(m_ord_3cat_RxBerlin)

# Risk × CT & Risk × Berlin
m_ord_3cat_RxCT_RxBer <- polr(
  vacc4_cat3 ~ Benefit_t3_c + Risk_t3_c * CT_c + Risk_t3_c * berlin_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

summary(m_ord_3cat_RxCT_RxBer)

## 6. Benefit/Risk 両側で CT/Berlin と交互作用するフルモデル ---------------

m_ord_3cat_full <- polr(
  vacc4_cat3 ~ 
    Benefit_t3_c * CT_c +
    Benefit_t3_c * berlin_c +
    Risk_t3_c    * CT_c +
    Risk_t3_c    * berlin_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

summary(m_ord_3cat_full)

## 7. AIC / BIC によるモデル比較 --------------------------------------------

aic_table <- AIC(
  m_ord_3cat_main,
  m_ord_3cat_BxCT,
  m_ord_3cat_BxBerlin,
  m_ord_3cat_BxCT_BxBer,
  m_ord_3cat_RxCT,
  m_ord_3cat_RxBerlin,
  m_ord_3cat_RxCT_RxBer,
  m_ord_3cat_full
)

bic_table <- BIC(
  m_ord_3cat_main,
  m_ord_3cat_BxCT,
  m_ord_3cat_BxBerlin,
  m_ord_3cat_BxCT_BxBer,
  m_ord_3cat_RxCT,
  m_ord_3cat_RxBerlin,
  m_ord_3cat_RxCT_RxBer,
  m_ord_3cat_full
)

cat("\n### AIC 比較 ###\n")
print(aic_table)

cat("\n### BIC 比較 ###\n")
print(bic_table)

## 8. 近似的な LR 検定（主に main vs 各交互作用モデル） ---------------------

cat("\n--- LR test: main vs Benefit×CT ---\n")
print(anova(m_ord_3cat_main, m_ord_3cat_BxCT))

cat("\n--- LR test: main vs Benefit×Berlin ---\n")
print(anova(m_ord_3cat_main, m_ord_3cat_BxBerlin))

cat("\n--- LR test: main vs Benefit×CT & Benefit×Berlin ---\n")
print(anova(m_ord_3cat_main, m_ord_3cat_BxCT_BxBer))

cat("\n--- LR test: main vs Risk×CT ---\n")
print(anova(m_ord_3cat_main, m_ord_3cat_RxCT))

cat("\n--- LR test: main vs Risk×Berlin ---\n")
print(anova(m_ord_3cat_main, m_ord_3cat_RxBerlin))

cat("\n--- LR test: main vs Risk×CT & Risk×Berlin ---\n")
print(anova(m_ord_3cat_main, m_ord_3cat_RxCT_RxBer))

cat("\n--- LR test: main vs full (Benefit/Risk × CT/Berlin 全部) ---\n")
print(anova(m_ord_3cat_main, m_ord_3cat_full))

############################################################
# ここまでで：
#  - Benefit 側のモデレーション（CT, Berlin）
#  - Risk 側のモデレーション（CT, Berlin）
#  - その組み合わせ
# を 3カテゴリモデルで一気に比較できるようになります。
############################################################
