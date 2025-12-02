############################################################
# コア分析コード（前処理：gender を因子化、education を 3 水準カテゴリにまとめてから分析）
############################################################

## パッケージ ----------------------------------------------------------
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(MASS)   # polr

## 0. 前提チェック ------------------------------------------------------
if (!exists("merged_data_complete")) {
  stop("オブジェクト 'merged_data_complete' が環境に存在しません。先に読み込んでください。")
}
if (!("id" %in% names(merged_data_complete))) {
  stop("'id' 列が merged_data_complete に存在しません。")
}
if (!all(c("age","gender","education") %in% names(merged_data_complete))) {
  stop("age, gender, education のいずれかが merged_data_complete に存在しません。")
}

## 前処理：gender を因子化, education を 3 水準カテゴリに再コーディング ----------------
merged_data_complete <- merged_data_complete %>%
  mutate(
    gender = factor(gender),
    # 以下は education の元のコード (1〜7 など) に応じて低/中/高に分類
    education3 = case_when(
      education %in% c(1, 2)             ~ "low",    # 例: 元コード 1,2 → 低学歴
      education %in% c(3, 4, 5)          ~ "middle", # 3–5 → 中学歴
      education %in% c(6, 7)             ~ "high",   # 6–7 → 高学歴
      TRUE                               ~ NA_character_
    ),
    education3 = factor(education3, levels = c("low","middle","high"))
  )

## （必要なら元の education は残したままにするか drop するか選ぶ）

## 1. CT 列名のゆれ処理 ------------------------------------------------
ct_candidates <- c("mean_cliticalThinking_attitude",
                   "mean_criticalThinking_attitude")
ct_name <- ct_candidates[ct_candidates %in% names(merged_data_complete)][1]
if (is.na(ct_name) || is.null(ct_name)) {
  stop("CT列 が見つかりません。")
}
cat("CT列として使用する変数:", ct_name, "\n")

## 2. ニュメラシー（Numeracy）列の決定 -------------------------------
num_candidates <- c("berlin_correct_count",
                    "mean_subjective_numeracy",
                    "health_correct_count",
                    "numAttitude_average")
num_name <- num_candidates[num_candidates %in% names(merged_data_complete)][1]
if (is.na(num_name) || is.null(num_name)) {
  stop("ニュメラシー関連の列 が見つかりません。")
}
cat("ニュメラシー列として使用する変数:", num_name, "\n")

#######################################################################
# (1) t-1 Effectiveness / Safeness → t Benefit / Risk の LMM
#######################################################################

mk_wave_BR <- function(df, wave, ct_name, num_name, covars = c("age","gender","education3")) {
  if (wave == "t3") {
    df %>% transmute(
      id,
      wave = factor("t3", levels = c("t3","t4")),
      CT    = .data[[ct_name]],
      NUM   = .data[[num_name]],
      Effectiveness_lag = Effectiveness_2,
      Safeness_lag      = Safeness_2,
      Benefit           = HBM_average_benefit_3,
      Risk              = HBM_average_risk_3,
      age        = .data[[covars[1]]],
      gender     = .data[[covars[2]]],
      education3 = .data[[covars[3]]]
    )
  } else {
    df %>% transmute(
      id,
      wave = factor("t4", levels = c("t3","t4")),
      CT    = .data[[ct_name]],
      NUM   = .data[[num_name]],
      Effectiveness_lag = Effectiveness_3,
      Safeness_lag      = Safeness_3,
      Benefit           = HBM_average_benefit_4,
      Risk              = HBM_average_risk_4,
      age        = .data[[covars[1]]],
      gender     = .data[[covars[2]]],
      education3 = .data[[covars[3]]]
    )
  }
}

br_t3 <- mk_wave_BR(merged_data_complete, "t3", ct_name, num_name)
br_t4 <- mk_wave_BR(merged_data_complete, "t4", ct_name, num_name)

br_long <- bind_rows(br_t3, br_t4) %>%
  filter(
    !is.na(id),
    !is.na(CT), !is.na(NUM),
    !is.na(Effectiveness_lag), !is.na(Safeness_lag),
    !is.na(Benefit), !is.na(Risk),
    !is.na(age), !is.na(gender), !is.na(education3)
  )

br_long2 <- br_long %>%
  pivot_longer(cols = c(Benefit, Risk),
               names_to = "outcome",
               values_to = "y") %>%
  filter(!is.na(y)) %>%
  mutate(
    CT_c   = as.numeric(scale(CT, center = TRUE, scale = FALSE)),
    NUM_c  = as.numeric(scale(NUM, center = TRUE, scale = FALSE)),
    Eff_c  = as.numeric(scale(Effectiveness_lag, center = TRUE, scale = FALSE)),
    Safe_c = as.numeric(scale(Safeness_lag,      center = TRUE, scale = FALSE))
  )

## Benefit 用 LMM ------------------------------------------------
br_ben <- br_long2 %>% filter(outcome == "Benefit")
m_ben_cov3 <- lmer(
  y ~ Eff_c * (CT_c + NUM_c) +
    Safe_c * (CT_c + NUM_c) +
    wave +
    age + gender + education3 +
    (1 | id),
  data = br_ben
)
print(summary(m_ben_cov3))

## Risk 用 LMM ------------------------------------------------
br_risk <- br_long2 %>% filter(outcome == "Risk")
m_risk_cov3 <- lmer(
  y ~ Eff_c * (CT_c + NUM_c) +
    Safe_c * (CT_c + NUM_c) +
    wave +
    age + gender + education3 +
    (1 | id),
  data = br_risk
)
print(summary(m_risk_cov3))

#######################################################################
# (2) t-1 Benefit / Risk → t4 接種回数(3カテゴリ) の順序ロジスティック
#     共変量 age, gender, education3 を含める
#######################################################################

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
    education3 = education3
  ) %>%
  filter(
    !is.na(id),
    !is.na(Benefit_t3), !is.na(Risk_t3),
    !is.na(CT), !is.na(NUM),
    !is.na(vacc4),
    !is.na(age), !is.na(gender), !is.na(education3)
  )

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
    Benefit_t3_c = as.numeric(scale(Benefit_t3, center = TRUE, scale = FALSE)),
    Risk_t3_c    = as.numeric(scale(Risk_t3,    center = TRUE, scale = FALSE)),
    CT_c         = as.numeric(scale(CT,         center = TRUE, scale = FALSE)),
    NUM_c        = as.numeric(scale(NUM,        center = TRUE, scale = FALSE))
  )

m_vacc3_cov3 <- polr(
  vacc4_cat3 ~ Benefit_t3_c * (CT_c + NUM_c) +
    Risk_t3_c    * (CT_c + NUM_c) +
    age + gender + education3,
  data   = cross_dat,
  method = "logistic",
  Hess   = TRUE
)
print(summary(m_vacc3_cov3))
