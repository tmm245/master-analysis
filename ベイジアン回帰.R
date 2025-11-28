############################################################
# brms 多変量マルチレベルモデル
#  - Benefit, Risk, Acceptance: LMM（ガウス）
#  - num_vaccination: 累積接種回数の順序ロジスティック(GLMM)
#
# 目的：
#   Eff/Safe_(t-1) → Benefit/Risk/Acceptance_t（態度）
#   Eff/Safe_(t-1), 態度_t → num_vaccination_t（行動）
#   を「同じランダム切片 (1|id)」で同時推定する
############################################################

## パッケージ ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(brms)
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

## 2. t3 / t4 の long データを作成 -----------------------------------------
# t3: (Eff_2, Safe_2) → Benefit_3, Risk_3, Acceptance_3, num_vaccination_3
# t4: (Eff_3, Safe_3) → Benefit_4, Risk_4, Acceptance_4, num_vaccination_4

mk_wave_long <- function(df, wave, ct_name) {
  if (wave == "t3") {
    dplyr::transmute(
      df,
      id,
      wave = factor("t3", levels = c("t3", "t4")),
      CT   = .data[[ct_name]],
      Effectiveness_lag = Effectiveness_2,
      Safeness_lag      = Safeness_2,
      Benefit           = HBM_average_benefit_3,
      Risk              = HBM_average_risk_3,
      Acceptance        = HBM_average_acceptance_3,
      num_vaccination   = num_vaccination_3
    )
  } else { # wave == "t4"
    dplyr::transmute(
      df,
      id,
      wave = factor("t4", levels = c("t3", "t4")),
      CT   = .data[[ct_name]],
      Effectiveness_lag = Effectiveness_3,
      Safeness_lag      = Safeness_3,
      Benefit           = HBM_average_benefit_4,
      Risk              = HBM_average_risk_4,
      Acceptance        = HBM_average_acceptance_4,
      num_vaccination   = num_vaccination_4
    )
  }
}

long_t3 <- mk_wave_long(merged_data_complete, "t3", ct_name)
long_t4 <- mk_wave_long(merged_data_complete, "t4", ct_name)

long_data <- dplyr::bind_rows(long_t3, long_t4) %>%
  # 必要な列が欠損している行を落とす
  dplyr::filter(
    !is.na(id),
    !is.na(CT),
    !is.na(Effectiveness_lag),
    !is.na(Safeness_lag),
    !is.na(Benefit),
    !is.na(Risk),
    !is.na(Acceptance),
    !is.na(num_vaccination)
  ) %>%
  dplyr::mutate(
    # 中心化（scale=FALSE で平均のみ引く）
    CT_c   = as.numeric(scale(CT, center = TRUE, scale = FALSE)),
    Eff_c  = as.numeric(scale(Effectiveness_lag, center = TRUE, scale = FALSE)),
    Safe_c = as.numeric(scale(Safeness_lag,      center = TRUE, scale = FALSE)),
    Ben_c  = as.numeric(scale(Benefit,           center = TRUE, scale = FALSE)),
    Risk_c = as.numeric(scale(Risk,              center = TRUE, scale = FALSE)),
    Acc_c  = as.numeric(scale(Acceptance,        center = TRUE, scale = FALSE)),
    # 接種回数は順序カテゴリとして扱う
    num_vaccination = factor(num_vaccination, ordered = TRUE)
  )

cat("long_data 次元:", paste(dim(long_data), collapse = " x "), "\n")
print(table(long_data$wave))

## 3. brms 用の bf 式を定義 -----------------------------------------------

# 態度3つはガウスLMM
bf_ben <- bf(
  Benefit ~ Eff_c * CT_c + Safe_c + wave + (1 | id)
)

bf_risk <- bf(
  Risk ~ Eff_c * CT_c + Safe_c + wave + (1 | id)
)

bf_acc <- bf(
  Acceptance ~ Eff_c * CT_c + Safe_c + wave + (1 | id)
)

# 行動（接種回数）は ordinal GLMM
# ★ここで Acceptance（Acc_c）を予測子から外した版に変更★
bf_vacc <- bf(
  num_vaccination ~ 
    Eff_c + Safe_c +         # t-1 の効果・安全性
    Ben_c + Risk_c +         # 認知（Benefit / Risk）
    # Acc_c は入れない
    CT_c +                   # CTの主効果（必要なら残す）
    wave +
    (1 | id),
  family = cumulative("logit")
)

## 4. モデルを同時推定 ------------------------------------------------------

options(mc.cores = parallel::detectCores())

fit_mv <- brm(
  bf_ben + bf_risk + bf_acc + bf_vacc + set_rescor(FALSE),
  data   = long_data,
  chains = 4,
  cores  = 4,
  iter   = 4000,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

## 5. 結果の確認 ------------------------------------------------------------

summary(fit_mv)

# 各サブモデルごとに見るなら：
summary(fit_mv, resp = "Benefit")
summary(fit_mv, resp = "Risk")
summary(fit_mv, resp = "Acceptance")
summary(fit_mv, resp = "num_vaccination")

# 事後分布からパス係数や間接効果を取り出したり、
# conditional_effects() で Eff/CT の単純傾きを図示したりも可能
############################################################

####################ベイジアン媒介#########################
############################################################
# brms 多変量モデル fit_mv から
# Eff/Safe → Benefit/Risk → num_vaccination の媒介効果を算出
#
# - 間接効果（Eff/Safe → Benefit → vacc, Risk → vacc）
# - 合計の間接効果（Benefit + Risk）
# - 直接効果（Eff/Safe → vacc）
# - 総効果（直接＋間接）
# - 間接/総効果の比率
############################################################

library(brms)
library(posterior)
library(dplyr)

## 1. 事後サンプルの取得 ------------------------------------
# fit_mv は既に推定済みの brmsfit オブジェクトとする
# （必要なら load() などで読み込んでおく）

draws <- as_draws_df(fit_mv)

# 係数名を一度確認（必要に応じて走らせてください）
# print(names(draws)[grepl("^b_", names(draws))])

## 2. パス係数の取り出し ------------------------------------
# a パス：Eff/Safe → Benefit/Risk
a_eff_ben   <- draws$b_Benefit_Eff_c
a_safe_ben  <- draws$b_Benefit_Safe_c
a_eff_risk  <- draws$b_Risk_Eff_c
a_safe_risk <- draws$b_Risk_Safe_c

# b パス：Benefit/Risk → num_vaccination
#   ※ brms が num_vaccination を内部的に "numvaccination" として扱うので注意
b_ben  <- draws$b_numvaccination_Ben_c
b_risk <- draws$b_numvaccination_Risk_c

# c' パス：Eff/Safe → num_vaccination（直接効果）
c_eff  <- draws$b_numvaccination_Eff_c
c_safe <- draws$b_numvaccination_Safe_c

## 3. 間接効果の計算 ----------------------------------------
# Effectiveness の媒介経路
ind_eff_ben   <- a_eff_ben  * b_ben              # Eff → Benefit → vacc
ind_eff_risk  <- a_eff_risk * b_risk            # Eff → Risk    → vacc
ind_eff_total <- ind_eff_ben + ind_eff_risk     # 合計間接効果

# Safeness の媒介経路
ind_safe_ben   <- a_safe_ben  * b_ben           # Safe → Benefit → vacc
ind_safe_risk  <- a_safe_risk * b_risk         # Safe → Risk    → vacc
ind_safe_total <- ind_safe_ben + ind_safe_risk # 合計間接効果

# 総効果（Total = 直接 + 間接）
total_eff  <- c_eff  + ind_eff_total
total_safe <- c_safe + ind_safe_total

## 4. 要約関数 ------------------------------------------------
summarise_vec <- function(x) {
  qs <- quantile(x, probs = c(.025, .5, .975), na.rm = TRUE)
  tibble(
    mean   = mean(x, na.rm = TRUE),
    q2.5   = qs[[1]],
    median = qs[[2]],
    q97.5  = qs[[3]],
    pr_gt0 = mean(x > 0, na.rm = TRUE)  # P(効果 > 0 | data)
  )
}

## 5. 結果をテーブル化 --------------------------------------
mediation_results <- bind_rows(
  summarise_vec(ind_eff_ben)   |> mutate(effect = "Eff: indirect via Benefit"),
  summarise_vec(ind_eff_risk)  |> mutate(effect = "Eff: indirect via Risk"),
  summarise_vec(ind_eff_total) |> mutate(effect = "Eff: total indirect (Ben+Risk)"),
  summarise_vec(c_eff)         |> mutate(effect = "Eff: direct (c')"),
  summarise_vec(total_eff)     |> mutate(effect = "Eff: total (direct + indirect)"),
  
  summarise_vec(ind_safe_ben)   |> mutate(effect = "Safe: indirect via Benefit"),
  summarise_vec(ind_safe_risk)  |> mutate(effect = "Safe: indirect via Risk"),
  summarise_vec(ind_safe_total) |> mutate(effect = "Safe: total indirect (Ben+Risk)"),
  summarise_vec(c_safe)         |> mutate(effect = "Safe: direct (c')"),
  summarise_vec(total_safe)     |> mutate(effect = "Safe: total (direct + indirect)")
) |>
  relocate(effect)

print(mediation_results)

## 6. 「間接 / 総効果」の割合（参考） --------------------------
ratio_results <- bind_rows(
  summarise_vec(ind_eff_total / total_eff)   |> mutate(effect = "Eff: (indirect / total)"),
  summarise_vec(ind_safe_total / total_safe) |> mutate(effect = "Safe: (indirect / total)")
) |>
  relocate(effect)

print(ratio_results)


############################################################
# ※ 注意：
# - パラメータ名が異なる場合は names(draws) で確認して適宜書き換えてください。
# - num_vaccination は累積ロジットの線形予測子上の効果（ロジットスケール）です。
# - 結果の解釈は「ロジット上の媒介効果」であることを踏まえて行ってください。
############################################################

#################ベイジアン媒介###############################

############################################################
# brms 多変量モデル fit_mv から
# Eff/Safe → Benefit/Risk → num_vaccination の媒介効果を算出
#
# - 間接効果（Eff/Safe → Benefit → vacc, Risk → vacc）
# - 合計の間接効果（Benefit + Risk）
# - 直接効果（Eff/Safe → vacc）
# - 総効果（直接＋間接）
# - 間接/総効果の比率
############################################################

library(brms)
library(posterior)
library(dplyr)

## 1. 事後サンプルの取得 ------------------------------------
# fit_mv は既に推定済みの brmsfit オブジェクトとする
# （必要なら load() などで読み込んでおく）

draws <- as_draws_df(fit_mv)

# 係数名を一度確認（必要に応じて走らせてください）
# print(names(draws)[grepl("^b_", names(draws))])

## 2. パス係数の取り出し ------------------------------------
# a パス：Eff/Safe → Benefit/Risk
a_eff_ben   <- draws$b_Benefit_Eff_c
a_safe_ben  <- draws$b_Benefit_Safe_c
a_eff_risk  <- draws$b_Risk_Eff_c
a_safe_risk <- draws$b_Risk_Safe_c

# b パス：Benefit/Risk → num_vaccination
b_ben  <- draws$b_num_vaccination_Ben_c
b_risk <- draws$b_num_vaccination_Risk_c

# c' パス：Eff/Safe → num_vaccination（直接効果）
c_eff  <- draws$b_num_vaccination_Eff_c
c_safe <- draws$b_num_vaccination_Safe_c

## 3. 間接効果の計算 ----------------------------------------
# Effectiveness の媒介経路
ind_eff_ben   <- a_eff_ben  * b_ben              # Eff → Benefit → vacc
ind_eff_risk  <- a_eff_risk * b_risk            # Eff → Risk    → vacc
ind_eff_total <- ind_eff_ben + ind_eff_risk     # 合計間接効果

# Safeness の媒介経路
ind_safe_ben   <- a_safe_ben  * b_ben           # Safe → Benefit → vacc
ind_safe_risk  <- a_safe_risk * b_risk         # Safe → Risk    → vacc
ind_safe_total <- ind_safe_ben + ind_safe_risk # 合計間接効果

# 総効果（Total = 直接 + 間接）
total_eff  <- c_eff  + ind_eff_total
total_safe <- c_safe + ind_safe_total

## 4. 要約関数 ------------------------------------------------
summarise_vec <- function(x) {
  qs <- quantile(x, probs = c(.025, .5, .975), na.rm = TRUE)
  tibble(
    mean   = mean(x, na.rm = TRUE),
    q2.5   = qs[[1]],
    median = qs[[2]],
    q97.5  = qs[[3]],
    pr_gt0 = mean(x > 0, na.rm = TRUE)  # P(効果 > 0 | data)
  )
}

## 5. 結果をテーブル化 --------------------------------------
mediation_results <- bind_rows(
  summarise_vec(ind_eff_ben)   |> mutate(effect = "Eff: indirect via Benefit"),
  summarise_vec(ind_eff_risk)  |> mutate(effect = "Eff: indirect via Risk"),
  summarise_vec(ind_eff_total) |> mutate(effect = "Eff: total indirect (Ben+Risk)"),
  summarise_vec(c_eff)         |> mutate(effect = "Eff: direct (c')"),
  summarise_vec(total_eff)     |> mutate(effect = "Eff: total (direct + indirect)"),
  
  summarise_vec(ind_safe_ben)   |> mutate(effect = "Safe: indirect via Benefit"),
  summarise_vec(ind_safe_risk)  |> mutate(effect = "Safe: indirect via Risk"),
  summarise_vec(ind_safe_total) |> mutate(effect = "Safe: total indirect (Ben+Risk)"),
  summarise_vec(c_safe)         |> mutate(effect = "Safe: direct (c')"),
  summarise_vec(total_safe)     |> mutate(effect = "Safe: total (direct + indirect)")
) |>
  relocate(effect)

print(mediation_results)

## 6. 「間接 / 総効果」の割合（参考） --------------------------
ratio_results <- bind_rows(
  summarise_vec(ind_eff_total / total_eff)   |> mutate(effect = "Eff: (indirect / total)"),
  summarise_vec(ind_safe_total / total_safe) |> mutate(effect = "Safe: (indirect / total)")
) |>
  relocate(effect)

print(ratio_results)

############################################################
# ※ 注意：
# - パラメータ名が異なる場合は names(draws) で確認して適宜書き換えてください。
# - num_vaccination は累積ロジットの線形予測子上の効果（ロジットスケール）です。
# - 結果の解釈は「ロジット上の媒介効果」であることを踏まえて行ってください。
############################################################

############################################################
# CT による Eff → Benefit/Risk の交互作用（moderated mediation）
# - モデル: fit_mv （brms 多変量モデル）
# - データ: long_data（CT_c, Eff_c, Safe_c などを含む）
#
# やること：
#  1. 事後サンプルから係数（a, a×CT, b, c'）を抽出
#  2. CT = low / mean / high ごとに Eff の「条件付き間接効果」を計算
#  3. その要約（平均・95%CrI・Pr(>0)）をテーブル化
############################################################

library(brms)
library(posterior)
library(dplyr)

## 0. 前提チェック -----------------------------------------------------------
if (!exists("fit_mv")) stop("fit_mv が環境にありません。先に brms モデルを推定してください。")
if (!exists("long_data")) stop("long_data が環境にありません。")

## 1. 事後サンプルを取得 ----------------------------------------------------
draws <- as_draws_df(fit_mv)

# 必要なら、パラメータ名を一度確認（デバッグ用）
# names(draws)[grepl("Eff_c",  names(draws))]
# names(draws)[grepl("Ben_c",  names(draws))]
# names(draws)[grepl("Risk_c", names(draws))]
# names(draws)[grepl("CT_c",   names(draws))]
# names(draws)[grepl("numvaccination", names(draws))]

## 2. パス係数の取り出し -----------------------------------------------------
# a パス：Eff/Safe → Benefit/Risk（平均的な傾き）
a_eff_ben   <- draws[, "b_Benefit_Eff_c" ] |> unlist()
a_safe_ben  <- draws[, "b_Benefit_Safe_c"] |> unlist()
a_eff_risk  <- draws[, "b_Risk_Eff_c"    ] |> unlist()
a_safe_risk <- draws[, "b_Risk_Safe_c"   ] |> unlist()

# a×CT：Eff × CT → Benefit/Risk の交互作用（Eff の傾きが CT でどう変わるか）
a_eff_ben_CT  <- draws[, "b_Benefit_Eff_c:CT_c"] |> unlist()
a_eff_risk_CT <- draws[, "b_Risk_Eff_c:CT_c"   ] |> unlist()

# b パス：Benefit/Risk → num_vaccination
#   ※ brms が応答名 num_vaccination を "numvaccination" として扱う点に注意
b_ben  <- draws[, "b_numvaccination_Ben_c" ] |> unlist()
b_risk <- draws[, "b_numvaccination_Risk_c"] |> unlist()

# c' パス：Eff/Safe → num_vaccination（直接効果）
c_eff  <- draws[, "b_numvaccination_Eff_c" ] |> unlist()
c_safe <- draws[, "b_numvaccination_Safe_c"] |> unlist()

## 3. CT 水準（low / mean / high）を定義 -----------------------------------
# CT_c は「平均のみ引いた」センタリングなので、CT のSDをそのまま使用
ct_sd <- sd(long_data$CT_c, na.rm = TRUE)

ct_levels <- c(
  low  = -ct_sd,   # CT が低い（約 -1SD）
  mean =  0,       # 平均的な CT
  high =  ct_sd    # CT が高い（約 +1SD）
)

print(ct_levels)

## 4. 要約用関数 ------------------------------------------------------------
summarise_vec <- function(x) {
  qs <- quantile(x, probs = c(.025, .5, .975), na.rm = TRUE)
  tibble(
    mean   = mean(x, na.rm = TRUE),
    q2.5   = qs[[1]],
    median = qs[[2]],
    q97.5  = qs[[3]],
    pr_gt0 = mean(x > 0, na.rm = TRUE)  # P(効果 > 0 | data)
  )
}

## 5. CT ごとの「条件付き間接効果」を計算 ----------------------------------
# Eff のみを対象にした moderated mediation（Eff → Benefit/Risk → vacc の経路）
ct_med_list <- list()

for (nm in names(ct_levels)) {
  ct_val <- ct_levels[[nm]]
  
  # --- a パス（Eff→Benefit/Risk）の CT 条件付き傾き ----------------------
  # a_eff_ben(CT)  = a_eff_ben  + a_eff_ben_CT  * CT
  # a_eff_risk(CT) = a_eff_risk + a_eff_risk_CT * CT
  a_eff_ben_ct  <- a_eff_ben  + a_eff_ben_CT  * ct_val
  a_eff_risk_ct <- a_eff_risk + a_eff_risk_CT * ct_val
  
  # --- 間接効果（Eff のみ） -----------------------------------------------
  # Eff → Benefit → num_vaccination
  ind_eff_ben_ct   <- a_eff_ben_ct  * b_ben
  # Eff → Risk → num_vaccination
  ind_eff_risk_ct  <- a_eff_risk_ct * b_risk
  # 合計間接効果
  ind_eff_total_ct <- ind_eff_ben_ct + ind_eff_risk_ct
  
  # --- 総効果（Total = 直接 c' + 間接） ----------------------------------
  total_eff_ct <- c_eff + ind_eff_total_ct
  
  # --- 要約テーブル -------------------------------------------------------
  ct_med_list[[nm]] <- bind_rows(
    summarise_vec(a_eff_ben_ct)      |> mutate(path = "a_Eff→Benefit (conditional)"),
    summarise_vec(a_eff_risk_ct)     |> mutate(path = "a_Eff→Risk (conditional)"),
    
    summarise_vec(ind_eff_ben_ct)    |> mutate(path = "Eff: indirect via Benefit"),
    summarise_vec(ind_eff_risk_ct)   |> mutate(path = "Eff: indirect via Risk"),
    summarise_vec(ind_eff_total_ct)  |> mutate(path = "Eff: total indirect (Ben+Risk)"),
    
    summarise_vec(c_eff)             |> mutate(path = "Eff: direct (c')"),
    summarise_vec(total_eff_ct)      |> mutate(path = "Eff: total (direct + indirect)")
  ) |>
    mutate(
      CT_level = nm,
      CT_value = ct_val
    ) |>
    relocate(CT_level, CT_value, path)
}

ct_moderated_results <- bind_rows(ct_med_list)

## 6. 間接/総効果の割合（CTごと） ------------------------------------------
ct_ratio_list <- list()

for (nm in names(ct_levels)) {
  ct_val <- ct_levels[[nm]]
  
  a_eff_ben_ct  <- a_eff_ben  + a_eff_ben_CT  * ct_val
  a_eff_risk_ct <- a_eff_risk + a_eff_risk_CT * ct_val
  
  ind_eff_ben_ct   <- a_eff_ben_ct  * b_ben
  ind_eff_risk_ct  <- a_eff_risk_ct * b_risk
  ind_eff_total_ct <- ind_eff_ben_ct + ind_eff_risk_ct
  
  total_eff_ct <- c_eff + ind_eff_total_ct
  
  ratio_ct <- ind_eff_total_ct / total_eff_ct
  
  ct_ratio_list[[nm]] <- summarise_vec(ratio_ct) |>
    mutate(
      CT_level = nm,
      CT_value = ct_val,
      path     = "Eff: (indirect / total)"
    ) |>
    relocate(CT_level, CT_value, path)
}

ct_ratio_results <- bind_rows(ct_ratio_list)

## 7. 結果の表示 ------------------------------------------------------------
print(ct_moderated_results)
print(ct_ratio_results)

############################################################
# ct_moderated_results:
#  - CT_level (low/mean/high)
#  - path:
#      * "a_Eff→Benefit (conditional)"
#      * "a_Eff→Risk (conditional)"
#      * "Eff: indirect via Benefit"
#      * "Eff: indirect via Risk"
#      * "Eff: total indirect (Ben+Risk)"
#      * "Eff: direct (c')"
#      * "Eff: total (direct + indirect)"
#  - 各 path ごとに mean, 95%CrI, Pr(>0)
#
# ct_ratio_results:
#  - CTごとの「Eff の間接効果 / 総効果」の割合
############################################################

