############################################################
# コア分析コード（前処理：gender, education を因子化してから分析）
############################################################

## パッケージ ----------------------------------------------------------
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(MASS)   # polr

## 0. 前提チェック ------------------------------------------------------

if (!exists("merged_data_complete")) {
  stop("オブジェクト 'merged_data_complete' が環境にありません。先に読み込んでください。")
}

if (!("id" %in% names(merged_data_complete))) {
  stop("'id' 列が merged_data_complete に存在しません。")
}

## 共変量として使う変数が存在するかチェック
covars <- c("age", "gender", "education")
missing_covars <- covars[ !(covars %in% names(merged_data_complete)) ]
if (length(missing_covars) > 0) {
  stop(paste0("以下の共変量が merged_data_complete に存在しません: ",
              paste(missing_covars, collapse = ", ")))
}

## 前処理：gender, education を因子化 (factor) に  --------------------
merged_data_complete <- merged_data_complete %>%
  # 必要に応じてラベルは調整
  mutate(
    gender = factor(gender),
    education = factor(education)
  )

## 1. CT 列名のゆれ処理 ------------------------------------------------
ct_candidates <- c("mean_cliticalThinking_attitude",
                   "mean_criticalThinking_attitude")
ct_name <- ct_candidates[ct_candidates %in% names(merged_data_complete)][1]

if (is.na(ct_name) || is.null(ct_name)) {
  stop("CT列 (mean_cliticalThinking_attitude / mean_criticalThinking_attitude) が見つかりません。")
}
cat("CT列として使用する変数:", ct_name, "\n")

## 2. ニュメラシー（Numeracy）列の決定 -------------------------------
num_candidates <- c("berlin_correct_count",
                    "mean_subjective_numeracy",
                    "health_correct_count",
                    "numAttitude_average")
num_name <- num_candidates[num_candidates %in% names(merged_data_complete)][1]

if (is.na(num_name) || is.null(num_name)) {
  stop("ニュメラシー関連の列 (berlin_correct_count 等) が見つかりません。")
}
cat("ニュメラシー列として使用する変数:", num_name, "\n")

#######################################################################
# (1) t-1 Effectiveness / Safeness → t Benefit / Risk の LMM
#######################################################################

## 1-1. t3 / t4 用の long データを作成 -------------------------------
mk_wave_BR <- function(df, wave, ct_name, num_name, covar_names = c("age","gender","education")) {
  if (wave == "t3") {
    dplyr::transmute(
      df,
      id,
      wave  = factor("t3", levels = c("t3","t4")),
      CT    = .data[[ct_name]],
      NUM   = .data[[num_name]],
      Effectiveness_lag = Effectiveness_2,
      Safeness_lag      = Safeness_2,
      Benefit           = HBM_average_benefit_3,
      Risk              = HBM_average_risk_3,
      age        = .data[[covar_names[1]]],
      gender     = .data[[covar_names[2]]],
      education  = .data[[covar_names[3]]]
    )
  } else { # wave == "t4"
    dplyr::transmute(
      df,
      id,
      wave  = factor("t4", levels = c("t3","t4")),
      CT    = .data[[ct_name]],
      NUM   = .data[[num_name]],
      Effectiveness_lag = Effectiveness_3,
      Safeness_lag      = Safeness_3,
      Benefit           = HBM_average_benefit_4,
      Risk              = HBM_average_risk_4,
      age        = .data[[covar_names[1]]],
      gender     = .data[[covar_names[2]]],
      education  = .data[[covar_names[3]]]
    )
  }
}

br_t3 <- mk_wave_BR(merged_data_complete, "t3", ct_name, num_name)
br_t4 <- mk_wave_BR(merged_data_complete, "t4", ct_name, num_name)

br_long <- bind_rows(br_t3, br_t4) %>%
  # 欠損を落とす（主要変数 + 共変量も NA の場合落とす）
  filter(
    !is.na(id),
    !is.na(CT),
    !is.na(NUM),
    !is.na(Effectiveness_lag),
    !is.na(Safeness_lag),
    !is.na(Benefit),
    !is.na(Risk),
    !is.na(age),
    !is.na(gender),
    !is.na(education)
  )

cat("br_long の次元:", paste(dim(br_long), collapse = " x "), "\n")
print(table(br_long$wave))

## 1-2. Benefit / Risk を1本のロングにして LMM 用データに -------------
br_long2 <- br_long %>%
  pivot_longer(
    cols = c(Benefit, Risk),
    names_to  = "outcome",
    values_to = "y"
  ) %>%
  filter(!is.na(y)) %>%
  mutate(
    # センタリング（時間不変の CT / NUM とラグ変数）
    CT_c   = as.numeric(scale(CT, center = TRUE, scale = FALSE)),
    NUM_c  = as.numeric(scale(NUM, center = TRUE, scale = FALSE)),
    Eff_c  = as.numeric(scale(Effectiveness_lag, center = TRUE, scale = FALSE)),
    Safe_c = as.numeric(scale(Safeness_lag,      center = TRUE, scale = FALSE))
    # age, gender, education は factor / numeric なのでそのまま
  )

cat("br_long2 の次元:", paste(dim(br_long2), collapse = " x "), "\n")
print(table(br_long2$outcome, br_long2$wave))

## 1-3. Benefit 用 LMM （共変量を追加） ----------------------------------
br_ben <- br_long2 %>% filter(outcome == "Benefit")

m_ben_cov <- lmer(
  y ~ Eff_c * (CT_c + NUM_c) +
    Safe_c * (CT_c + NUM_c) +
    wave +
    age + gender + education +
    (1 | id),
  data = br_ben
)

cat("\n==== LMM (共変量付き): t-1 Eff/Safe → t Benefit ====\n")
print(summary(m_ben_cov))

## 1-4. Risk 用 LMM （共変量を追加） -------------------------------------
br_risk <- br_long2 %>% filter(outcome == "Risk")

m_risk_cov <- lmer(
  y ~ Eff_c * (CT_c + NUM_c) +
    Safe_c * (CT_c + NUM_c) +
    wave +
    age + gender + education +
    (1 | id),
  data = br_risk
)

cat("\n==== LMM (共変量付き): t-1 Eff/Safe → t Risk ====\n")
print(summary(m_risk_cov))

#######################################################################
# (2) t-1 Benefit / Risk → t4 接種回数(3カテゴリ) の順序ロジスティック
#     共変量（age, gender, education）も含める
#######################################################################

## 2-1. t3 Benefit/Risk + t4 接種回数 + CT + NUM + 共変量 を 1 行にまとめる ------
cross_dat <- merged_data_complete %>%
  transmute(
    id,
    Benefit_t3 = HBM_average_benefit_3,
    Risk_t3    = HBM_average_risk_3,
    CT         = .data[[ct_name]],
    NUM        = .data[[num_name]],
    vacc4      = num_vaccination_4,
    age        = age,
    gender     = gender,
    education  = education
  ) %>%
  filter(
    !is.na(id),
    !is.na(Benefit_t3),
    !is.na(Risk_t3),
    !is.na(CT),
    !is.na(NUM),
    !is.na(vacc4),
    !is.na(age),
    !is.na(gender),
    !is.na(education)
  )

cat("\ncross_dat 次元:", paste(dim(cross_dat), collapse = " x "), "\n")
print(table(cross_dat$vacc4))

## 2-2. 接種回数を 3カテゴリ (0, 2, 3plus) にまとめる ------------------
cross_dat <- cross_dat %>%
  mutate(
    vacc4_cat3 = case_when(
      vacc4 == 0 ~ "0",
      vacc4 == 2 ~ "2",
      vacc4 >= 3 ~ "3plus",
      TRUE       ~ NA_character_
    )
  ) %>%
  filter(!is.na(vacc4_cat3)) %>%
  mutate(
    vacc4_cat3 = factor(vacc4_cat3,
                        levels = c("0","2","3plus"),
                        ordered = TRUE),
    # センタリング (Benefit, Risk, CT, NUM)
    Benefit_t3_c = as.numeric(scale(Benefit_t3, center = TRUE, scale = FALSE)),
    Risk_t3_c    = as.numeric(scale(Risk_t3,    center = TRUE, scale = FALSE)),
    CT_c         = as.numeric(scale(CT,         center = TRUE, scale = FALSE)),
    NUM_c        = as.numeric(scale(NUM,        center = TRUE, scale = FALSE))
    # age, gender, education は factor / numeric のまま
  )

cat("cross_dat（3カテゴリ）の次元:", paste(dim(cross_dat), collapse = " x "), "\n")
print(table(cross_dat$vacc4_cat3))

## 2-3. 順序ロジスティック（polr） -------------------------------------
m_vacc3_cov <- polr(
  vacc4_cat3 ~ Benefit_t3_c * (CT_c + NUM_c) +
    Risk_t3_c    * (CT_c + NUM_c) +
    age + gender + education,
  data   = cross_dat,
  method = "logistic",
  Hess   = TRUE
)

cat("\n==== 順序ロジスティック (共変量付き): t3 Benefit/Risk → t4 接種3カテゴリ ====\n")
print(summary(m_vacc3_cov))
