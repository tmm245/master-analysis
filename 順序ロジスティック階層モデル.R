############################################################
# 順序ロジスティック階層モデル ＋ 各種モデル結果を
# すべて一つの CSV にまとめて保存する版
#
# 出力ファイル:
#   all_model_results.csv
############################################################

## パッケージ ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ordinal)   # clmm, clm
library(MASS)      # polr
library(lme4)
library(lmerTest)
library(ggplot2)
library(rlang)

############################################################
# 0. 前提チェック ＆ 共通ヘルパー
############################################################

if (!exists("merged_data_complete")) {
  stop("オブジェクト 'merged_data_complete' が環境にありません。先に読み込んでください。")
}

if (!("id" %in% names(merged_data_complete))) {
  stop("'id' 列が merged_data_complete に存在しません。")
}

## CT列名のゆれ処理（共通） -------------------------
ct_candidates <- c("mean_cliticalThinking_attitude",
                   "mean_criticalThinking_attitude")

ct_name <- ct_candidates[ct_candidates %in% names(merged_data_complete)][1]

if (is.na(ct_name) || is.null(ct_name)) {
  stop("CT列 (mean_cliticalThinking_attitude / mean_criticalThinking_attitude) が見つかりません。")
}

cat("CT列として使用する変数:", ct_name, "\n")

############################################################
# 結果保存用ヘルパー
############################################################

model_results <- list()

add_result <- function(fit,
                       model_id,
                       outcome,
                       predictors = NA_character_,
                       notes = NA_character_) {
  model_class <- class(fit)[1]
  n_obs <- tryCatch(nobs(fit), error = function(e) NA_integer_)
  aic   <- tryCatch(AIC(fit),  error = function(e) NA_real_)
  bic   <- tryCatch(BIC(fit),  error = function(e) NA_real_)
  
  # ordinal / glm(binomial) 系（logitスケール → OR 計算）
  if (inherits(fit, "clmm") ||
      inherits(fit, "clm")  ||
      inherits(fit, "polr") ||
      (inherits(fit, "glm") && family(fit)$family == "binomial")) {
    
    cf <- coef(summary(fit))
    colnames_cf <- colnames(cf)
    
    est_col <- colnames_cf[1]
    se_col  <- colnames_cf[2]
    stat_col <- if ("z value" %in% colnames_cf) {
      "z value"
    } else if ("t value" %in% colnames_cf) {
      "t value"
    } else {
      NA_character_
    }
    
    # p値の推定
    if ("Pr(>|z|)" %in% colnames_cf) {
      p_vals <- cf[, "Pr(>|z|)"]
    } else if ("Pr(>|t|)" %in% colnames_cf) {
      p_vals <- cf[, "Pr(>|t|)"]
    } else if (!is.na(stat_col)) {
      p_vals <- 2 * pnorm(abs(cf[, stat_col]), lower.tail = FALSE)
    } else {
      p_vals <- NA_real_
    }
    
    # ordinal の cutpoint (0|1 など) は除外
    if (inherits(fit, "clmm") || inherits(fit, "clm") || inherits(fit, "polr")) {
      slope_idx <- !grepl("\\|", rownames(cf))
    } else {
      slope_idx <- rep(TRUE, nrow(cf))
    }
    
    cf_use <- cf[slope_idx, , drop = FALSE]
    p_use  <- p_vals[slope_idx]
    
    est  <- cf_use[, est_col]
    se   <- cf_use[, se_col]
    stat <- if (!is.na(stat_col)) cf_use[, stat_col] else NA_real_
    
    z <- qnorm(0.975)
    OR      <- exp(est)
    CI_low  <- exp(est - z * se)
    CI_high <- exp(est + z * se)
    
    res <- data.frame(
      model_id   = model_id,
      model_class= model_class,
      outcome    = outcome,
      predictors = predictors,
      term       = rownames(cf_use),
      estimate   = as.numeric(est),
      std_error  = as.numeric(se),
      statistic  = as.numeric(stat),
      p_value    = as.numeric(p_use),
      OR         = as.numeric(OR),
      CI_low     = as.numeric(CI_low),
      CI_high    = as.numeric(CI_high),
      nobs       = n_obs,
      AIC        = aic,
      BIC        = bic,
      notes      = notes,
      stringsAsFactors = FALSE
    )
    
  } else if (inherits(fit, "lmerMod")) {
    # LMM（線形）: ORはNA
    cf <- summary(fit)$coefficients
    est <- cf[, "Estimate"]
    se  <- cf[, "Std. Error"]
    tval<- cf[, "t value"]
    pval<- cf[, "Pr(>|t|)"]
    
    res <- data.frame(
      model_id   = model_id,
      model_class= model_class,
      outcome    = outcome,
      predictors = predictors,
      term       = rownames(cf),
      estimate   = as.numeric(est),
      std_error  = as.numeric(se),
      statistic  = as.numeric(tval),
      p_value    = as.numeric(pval),
      OR         = NA_real_,
      CI_low     = NA_real_,
      CI_high    = NA_real_,
      nobs       = n_obs,
      AIC        = aic,
      BIC        = bic,
      notes      = notes,
      stringsAsFactors = FALSE
    )
    
  } else {
    # それ以外: とりあえず係数だけ
    cf <- coef(summary(fit))
    est <- cf[, 1]
    se  <- cf[, 2]
    stat<- if (ncol(cf) >= 3) cf[, 3] else NA_real_
    pval<- if (ncol(cf) >= 4) cf[, 4] else NA_real_
    
    res <- data.frame(
      model_id   = model_id,
      model_class= model_class,
      outcome    = outcome,
      predictors = predictors,
      term       = rownames(cf),
      estimate   = as.numeric(est),
      std_error  = as.numeric(se),
      statistic  = as.numeric(stat),
      p_value    = as.numeric(pval),
      OR         = NA_real_,
      CI_low     = NA_real_,
      CI_high    = NA_real_,
      nobs       = n_obs,
      AIC        = aic,
      BIC        = bic,
      notes      = notes,
      stringsAsFactors = FALSE
    )
  }
  
  model_results[[length(model_results) + 1]] <<- res
}

############################################################
# 1. 縦断GLMM: num_vaccination を Eff/Safe_lag(+CT)で予測（A系）
############################################################

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

# ラグ変数の作成
vacc_long <- vacc_long %>%
  dplyr::group_by(id) %>%
  dplyr::arrange(wave, .by_group = TRUE) %>%
  dplyr::mutate(
    Effectiveness_lag = dplyr::lag(Effectiveness, 1),
    Safeness_lag      = dplyr::lag(Safeness, 1)
  ) %>%
  dplyr::ungroup()

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

# --- モデル A1: 交互作用なし ---
m_ord_core_simple <- clmm(
  num_vaccination ~ Eff_lag_c + Safe_lag_c + wave_factor + (1 | id),
  data  = model_data_core,
  link  = "logit",
  Hess  = TRUE,
  nAGQ  = 5
)

cat("\n==== モデル A1: num_vacc_t ~ Eff_lag + Safe_lag + wave + (1|id) ====\n")
print(summary(m_ord_core_simple))

add_result(
  m_ord_core_simple,
  model_id   = "A1_ord_core_simple",
  outcome    = "num_vaccination",
  predictors = "Eff_lag_c + Safe_lag_c + wave_factor",
  notes      = "clmm; wave >= 2"
)

# --- モデル A2: Eff_lag × CT の交互作用 ---
m_ord_core_int <- clmm(
  num_vaccination ~ Eff_lag_c * CT_c + Safe_lag_c + wave_factor + (1 | id),
  data  = model_data_core,
  link  = "logit",
  Hess  = TRUE,
  nAGQ  = 5
)

cat("\n==== モデル A2: num_vacc_t ~ Eff_lag * CT + Safe_lag + wave + (1|id) ====\n")
print(summary(m_ord_core_int))

add_result(
  m_ord_core_int,
  model_id   = "A2_ord_core_int",
  outcome    = "num_vaccination",
  predictors = "Eff_lag_c * CT_c + Safe_lag_c + wave_factor",
  notes      = "clmm; wave >= 2"
)

############################################################
# 2. t4 時点だけ: num_vaccination_4 ~ Benefit_4/Risk_4/Acceptance_4 × CT（B系）
############################################################

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

m_t4 <- clm(
  num_vaccination ~ Benefit_4_c * CT_c + Risk_4_c + Accept_4_c,
  data = model_data_t4,
  link = "logit",
  Hess = TRUE
)

cat("\n==== モデル B: t4 の順序ロジスティック ====\n")
print(summary(m_t4))

add_result(
  m_t4,
  model_id   = "B_t4_BRA_CT",
  outcome    = "num_vaccination_4",
  predictors = "Benefit_4_c * CT_c + Risk_4_c + Accept_4_c",
  notes      = "clm; cross-sectional t4"
)

############################################################
# 3. Acceptance（HBM_average_acceptance_3 / 4）を目的変数にした LMM
############################################################

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
  dplyr::filter(
    !is.na(id),
    !is.na(CT),
    !is.na(Effectiveness_lag),
    !is.na(Safeness_lag),
    !is.na(Acceptance)
  ) %>%
  dplyr::mutate(
    CT_c   = as.numeric(scale(CT, center = TRUE, scale = FALSE)),
    Eff_c  = as.numeric(scale(Effectiveness_lag, center = TRUE, scale = FALSE)),
    Safe_c = as.numeric(scale(Safeness_lag,      center = TRUE, scale = FALSE))
  )

cat("acc_long の次元:", paste(dim(acc_long), collapse = " x "), "\n")
print(table(acc_long$wave))

# コアモデル：Eff/Safe の主効果のみ
m_acc_core <- lmer(
  Acceptance ~ Eff_c + Safe_c + wave + (1 | id),
  data = acc_long
)
cat("\n==== m_acc_core: Acceptance ~ Eff_c + Safe_c + wave + (1|id) ====\n")
print(summary(m_acc_core))

add_result(
  m_acc_core,
  model_id   = "C1_acc_core",
  outcome    = "Acceptance",
  predictors = "Eff_c + Safe_c + wave",
  notes      = "lmer; t3,t4"
)

# 交互作用（Eff × CT）
m_acc_int <- lmer(
  Acceptance ~ Eff_c * CT_c + Safe_c + wave + (1 | id),
  data = acc_long
)
cat("\n==== m_acc_int: Acceptance ~ Eff_c * CT_c + Safe_c + wave + (1|id) ====\n")
print(summary(m_acc_int))

add_result(
  m_acc_int,
  model_id   = "C2_acc_EffxCT",
  outcome    = "Acceptance",
  predictors = "Eff_c * CT_c + Safe_c + wave",
  notes      = "lmer; t3,t4"
)

############################################################
# 4C. Benefit / Risk だけを説明変数にした
#     階層順序ロジスティックモデル（t3 + t4）
############################################################

## 4C-1. t3/t4 用の long_data を作成 -----------------------------

mk_wave_BR <- function(df, wave) {
  if (wave == "t3") {
    dplyr::transmute(
      df,
      id,
      wave            = factor("t3", levels = c("t3","t4")),
      num_vaccination = num_vaccination_3,
      Benefit         = HBM_average_benefit_3,
      Risk            = HBM_average_risk_3
    )
  } else {  # wave == "t4"
    dplyr::transmute(
      df,
      id,
      wave            = factor("t4", levels = c("t3","t4")),
      num_vaccination = num_vaccination_4,
      Benefit         = HBM_average_benefit_4,
      Risk            = HBM_average_risk_4
    )
  }
}

long_data <- dplyr::bind_rows(
  mk_wave_BR(merged_data_complete, "t3"),
  mk_wave_BR(merged_data_complete, "t4")
) %>%
  dplyr::filter(
    !is.na(id),
    !is.na(num_vaccination),
    !is.na(Benefit),
    !is.na(Risk)
  )

cat("long_data 次元:", paste(dim(long_data), collapse = " x "), "\n")
print(table(long_data$wave))

## 4C-2. ord_BR_data の作成 -----------------------------

ord_BR_data <- long_data %>%
  dplyr::filter(
    !is.na(num_vaccination),
    !is.na(Benefit),
    !is.na(Risk),
    !is.na(id),
    !is.na(wave)
  ) %>%
  dplyr::mutate(
    num_vaccination = factor(num_vaccination, ordered = TRUE),
    Ben_c_only  = as.numeric(scale(Benefit, center = TRUE, scale = FALSE)),
    Risk_c_only = as.numeric(scale(Risk,    center = TRUE, scale = FALSE)),
    wave_factor = factor(wave, levels = c("t3", "t4"))
  )

cat("ord_BR_data 次元:", paste(dim(ord_BR_data), collapse = " x "), "\n")
print(table(ord_BR_data$wave_factor, useNA = "ifany"))

# 階層順序ロジスティック
m_ord_BR <- clmm(
  num_vaccination ~ Ben_c_only + Risk_c_only + wave_factor + (1 | id),
  data  = ord_BR_data,
  link  = "logit",
  Hess  = TRUE,
  nAGQ  = 5
)

cat("\n==== モデル BR: num_vacc_t ~ Benefit + Risk + wave + (1|id) ====\n")
print(summary(m_ord_BR))

add_result(
  m_ord_BR,
  model_id   = "D1_ord_BR_clmm",
  outcome    = "num_vaccination",
  predictors = "Ben_c_only + Risk_c_only + wave_factor",
  notes      = "clmm; t3,t4"
)

# ランダム効果なしの ordinal logit（決定境界用）
m_br_clm <- clm(
  num_vaccination ~ Ben_c_only + Risk_c_only + wave_factor,
  data = ord_BR_data,
  link = "logit"
)

cat("\n==== m_br_clm: num_vacc ~ Benefit + Risk + wave (fixed only) ====\n")
print(summary(m_br_clm))

add_result(
  m_br_clm,
  model_id   = "D2_BR_clm_fixed",
  outcome    = "num_vaccination",
  predictors = "Ben_c_only + Risk_c_only + wave_factor",
  notes      = "clm; t3,t4"
)

############################################################
# 4C-3. 可視化（p_benefit, p_risk, p_br_regions_big, p_decision）
#   ※結果CSVには関係ないが元のコードを保持
############################################################

# 可視化用 numeric 接種回数
ord_BR_plot <- ord_BR_data %>%
  mutate(
    num_vacc_num = as.numeric(num_vaccination),
    wave_factor  = factor(wave_factor, levels = c("t3","t4"))
  )

# Benefit ビンごとの平均
ben_summary <- ord_BR_plot %>%
  group_by(wave_factor) %>%
  mutate(Ben_bin = ntile(Benefit, 5)) %>%
  group_by(wave_factor, Ben_bin) %>%
  summarise(
    Ben_mean       = mean(Benefit, na.rm = TRUE),
    num_vacc_mean  = mean(num_vacc_num, na.rm = TRUE),
    num_vacc_sd    = sd(num_vacc_num,   na.rm = TRUE),
    n              = n(),
    num_vacc_se    = num_vacc_sd / sqrt(n),
    .groups = "drop"
  )

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
  theme_bw()

# Risk ビンごとの平均
risk_summary <- ord_BR_plot %>%
  group_by(wave_factor) %>%
  mutate(Risk_bin = ntile(Risk, 5)) %>%
  group_by(wave_factor, Risk_bin) %>%
  summarise(
    Risk_mean      = mean(Risk, na.rm = TRUE),
    num_vacc_mean  = mean(num_vacc_num, na.rm = TRUE),
    num_vacc_sd    = sd(num_vacc_num,   na.rm = TRUE),
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
  theme_bw()

# Benefit×Risk 空間の多数派カテゴリ
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

p_br_regions_big <- ggplot(region_df,
                           aes(x = Ben_mid, y = Risk_mid,
                               fill = majority_vacc)) +
  geom_point(shape = 22, size = 10, colour = "grey20") +
  facet_wrap(~ wave_factor) +
  labs(
    x = "Benefit（ビン中心）",
    y = "Risk（ビン中心）",
    title = "Benefit × Risk 空間における多数派の接種回数カテゴリ（t3/t4）"
  ) +
  coord_fixed() +
  theme_bw()

# 決定境界（fixed model m_br_clm を使用）
ben_range  <- range(ord_BR_data$Benefit, na.rm = TRUE)
risk_range <- range(ord_BR_data$Risk,    na.rm = TRUE)

ben_seq  <- seq(ben_range[1],  ben_range[2],  length.out = 80)
risk_seq <- seq(risk_range[1], risk_range[2], length.out = 80)

ben_mean  <- mean(ord_BR_data$Benefit, na.rm = TRUE)
risk_mean <- mean(ord_BR_data$Risk,    na.rm = TRUE)

grid <- expand.grid(
  Benefit    = ben_seq,
  Risk       = risk_seq,
  wave_factor = levels(ord_BR_data$wave_factor)
) %>%
  mutate(
    Ben_c_only  = Benefit - ben_mean,
    Risk_c_only = Risk    - risk_mean
  )

prob_mat <- predict(m_br_clm, newdata = grid, type = "prob")
if (is.list(prob_mat)) prob_mat <- prob_mat$fit

max_col <- max.col(prob_mat, ties.method = "first")
grid$pred_class <- factor(
  colnames(prob_mat)[max_col],
  levels = levels(ord_BR_data$num_vaccination)
)

ord_BR_plot2 <- ord_BR_data %>%
  mutate(
    wave_factor = factor(wave_factor, levels = c("t3","t4"))
  )

p_decision <- ggplot() +
  geom_raster(
    data = grid,
    aes(x = Benefit, y = Risk, fill = pred_class),
    alpha = 0.35
  ) +
  geom_contour(
    data = grid,
    aes(x = Benefit, y = Risk, z = as.numeric(pred_class)),
    colour = "black",
    linewidth = 0.3
  ) +
  geom_point(
    data = ord_BR_plot2,
    aes(x = Benefit, y = Risk, color = num_vaccination),
    size = 1.3,
    alpha = 0.8
  ) +
  facet_wrap(~ wave_factor) +
  coord_fixed() +
  theme_bw()

############################################################
# 5. t3 の Benefit / Risk が t4 の接種回数カテゴリを予測する
#    クロスラグ ordinal logit
############################################################

t3_BR <- long_data %>%
  dplyr::filter(wave == "t3") %>%
  dplyr::select(
    id,
    Benefit3 = Benefit,
    Risk3    = Risk
  )

t4_vacc <- long_data %>%
  dplyr::filter(wave == "t4") %>%
  dplyr::select(
    id,
    num_vacc_t4 = num_vaccination
  )

crosslag_dat <- t3_BR %>%
  dplyr::inner_join(t4_vacc, by = "id") %>%
  dplyr::filter(
    !is.na(Benefit3),
    !is.na(Risk3),
    !is.na(num_vacc_t4)
  ) %>%
  dplyr::mutate(
    num_vacc_t4 = factor(num_vacc_t4, ordered = TRUE)
  )

cat("crosslag_dat 次元:", paste(dim(crosslag_dat), collapse = " x "), "\n")
print(table(crosslag_dat$num_vacc_t4))

ben3_mean  <- mean(crosslag_dat$Benefit3, na.rm = TRUE)
risk3_mean <- mean(crosslag_dat$Risk3,    na.rm = TRUE)

crosslag_dat <- crosslag_dat %>%
  dplyr::mutate(
    Ben3_c  = Benefit3 - ben3_mean,
    Risk3_c = Risk3    - risk3_mean
  )

m_cross_BR <- clm(
  num_vacc_t4 ~ Ben3_c + Risk3_c,
  data = crosslag_dat,
  link = "logit"
)

cat("\n==== m_cross_BR: num_vacc_t4 ~ Ben3_c + Risk3_c ====\n")
print(summary(m_cross_BR))

add_result(
  m_cross_BR,
  model_id   = "E1_cross_BR_t3_to_t4",
  outcome    = "num_vaccination_4",
  predictors = "Ben3_c + Risk3_c",
  notes      = "clm; cross-lag t3 attitudes -> t4 doses"
)

# 決定境界用グリッド
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

prob_mat_cross <- predict(
  m_cross_BR,
  newdata = grid_cross,
  type = "prob"
)
if (is.list(prob_mat_cross)) prob_mat_cross <- prob_mat_cross$fit

max_col_cross <- max.col(prob_mat_cross, ties.method = "first")

grid_cross$pred_class <- factor(
  colnames(prob_mat_cross)[max_col_cross],
  levels = levels(crosslag_dat$num_vacc_t4)
)

p_cross_decision <- ggplot() +
  geom_raster(
    data = grid_cross,
    aes(x = Benefit3, y = Risk3, fill = pred_class),
    alpha = 0.35
  ) +
  geom_contour(
    data = grid_cross,
    aes(x = Benefit3, y = Risk3, z = as.numeric(pred_class)),
    colour = "black",
    linewidth = 0.3
  ) +
  geom_point(
    data = crosslag_dat,
    aes(x = Benefit3, y = Risk3, color = num_vacc_t4),
    size = 1.5,
    alpha = 0.8
  ) +
  coord_fixed() +
  theme_bw()

############################################################
# 6. 0 vs 2+ の二値ロジスティック回帰（t3 Benefit/Risk → t4 接種）
############################################################

logit_BR_t3t4 <- merged_data_complete %>%
  transmute(
    id,
    Benefit_t3 = HBM_average_benefit_3,
    Risk_t3    = HBM_average_risk_3,
    vacc4      = num_vaccination_4
  ) %>%
  mutate(
    vacc4_bin = case_when(
      vacc4 == 0        ~ "0",
      vacc4 >= 2        ~ "2plus",
      TRUE              ~ NA_character_
    ),
    vacc4_bin = factor(vacc4_bin, levels = c("0", "2plus")),
    Benefit_t3_c = scale(Benefit_t3, center = TRUE, scale = FALSE)[, 1],
    Risk_t3_c    = scale(Risk_t3,    center = TRUE, scale = FALSE)[, 1]
  ) %>%
  tidyr::drop_na(vacc4_bin, Benefit_t3_c, Risk_t3_c)

m_logit_BR_t3t4 <- glm(
  vacc4_bin ~ Benefit_t3_c + Risk_t3_c,
  data   = logit_BR_t3t4,
  family = binomial(link = "logit")
)

cat("\n==== m_logit_BR_t3t4: vacc4_bin(0 vs 2+) ~ Benefit_t3_c + Risk_t3_c ====\n")
print(summary(m_logit_BR_t3t4))

add_result(
  m_logit_BR_t3t4,
  model_id   = "F1_logit_0vs2plus",
  outcome    = "vacc4_bin (0 vs 2plus)",
  predictors = "Benefit_t3_c + Risk_t3_c",
  notes      = "glm binomial; sensitivity 0 vs 2+"
)

############################################################
# 7. 4カテゴリ (0/2/3/4plus) vs 3カテゴリ (0/2/3plus) 比較
############################################################

## 7-1. 4カテゴリ A系 ----------------------------------

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
    vacc4_catA = factor(
      vacc4_catA,
      ordered = TRUE,
      levels  = c("0","2","3","4plus")
    ),
    Benefit_t3_c = scale(Benefit_t3, center = TRUE, scale = FALSE)[, 1],
    Risk_t3_c    = scale(Risk_t3,    center = TRUE, scale = FALSE)[, 1]
  ) %>%
  tidyr::drop_na(vacc4_catA, Benefit_t3_c, Risk_t3_c)

cat("4カテゴリ (0,2,3,4plus) データ次元:",
    paste(dim(ord_BR_t3t4_A), collapse = " x "), "\n")

m_ord_BR_t3t4_A <- polr(
  vacc4_catA ~ Benefit_t3_c + Risk_t3_c,
  data = ord_BR_t3t4_A,
  Hess = TRUE
)

cat("\n==== m_ord_BR_t3t4_A: 4cat (0,2,3,4plus) ====\n")
print(summary(m_ord_BR_t3t4_A))

add_result(
  m_ord_BR_t3t4_A,
  model_id   = "G1_4cat_0_2_3_4plus",
  outcome    = "vacc4_catA",
  predictors = "Benefit_t3_c + Risk_t3_c",
  notes      = "polr; 4 categories 0/2/3/4+"
)

## 7-2. 3カテゴリ (0/2/3plus) -------------------------

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
    vacc4_cat3 = factor(
      vacc4_cat3,
      ordered = TRUE,
      levels  = c("0","2","3plus")
    ),
    Benefit_t3_c = scale(Benefit_t3, center = TRUE, scale = FALSE)[, 1],
    Risk_t3_c    = scale(Risk_t3,    center = TRUE, scale = FALSE)[, 1]
  ) %>%
  tidyr::drop_na(vacc4_cat3, Benefit_t3_c, Risk_t3_c)

cat("3カテゴリ (0,2,3plus) データ次元:",
    paste(dim(ord_BR_t3t4_3cat), collapse = " x "), "\n")

m_ord_BR_t3t4_3cat <- polr(
  vacc4_cat3 ~ Benefit_t3_c + Risk_t3_c,
  data = ord_BR_t3t4_3cat,
  Hess = TRUE
)

cat("\n==== m_ord_BR_t3t4_3cat: 3cat (0,2,3plus) ====\n")
print(summary(m_ord_BR_t3t4_3cat))

add_result(
  m_ord_BR_t3t4_3cat,
  model_id   = "G2_3cat_0_2_3plus",
  outcome    = "vacc4_cat3",
  predictors = "Benefit_t3_c + Risk_t3_c",
  notes      = "polr; 3 categories 0/2/3+"
)

############################################################
# 8. 3カテゴリ版（0,2,3plus）の順序ロジスティック
#   Benefit/Risk + CT + Berlin（主効果 & 各種交互作用）
############################################################

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
    Benefit_t3_c = scale(Benefit_t3, center = TRUE, scale = FALSE)[, 1],
    Risk_t3_c    = scale(Risk_t3,    center = TRUE, scale = FALSE)[, 1],
    CT_c         = scale(CT_raw,     center = TRUE, scale = FALSE)[, 1],
    berlin_c     = scale(berlin,     center = TRUE, scale = FALSE)[, 1]
  ) %>%
  tidyr::drop_na(vacc4_cat3, Benefit_t3_c, Risk_t3_c, CT_c, berlin_c)

cat("ord_BR_t3t4_3cat_ext 次元:",
    paste(dim(ord_BR_t3t4_3cat_ext), collapse = " x "), "\n")

## 8-1. ベースライン（主効果のみ） --------------------

m_ord_3cat_main <- polr(
  vacc4_cat3 ~ Benefit_t3_c + Risk_t3_c + CT_c + berlin_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

cat("\n==== m_ord_3cat_main: main effects ====\n")
print(summary(m_ord_3cat_main))

add_result(
  m_ord_3cat_main,
  model_id   = "H1_3cat_main",
  outcome    = "vacc4_cat3",
  predictors = "Benefit_t3_c + Risk_t3_c + CT_c + berlin_c",
  notes      = "polr; main effects only"
)

## 8-2. Benefit 側の交互作用 ----------------------------

# Benefit × CT
m_ord_3cat_BxCT <- polr(
  vacc4_cat3 ~ Benefit_t3_c * CT_c + Risk_t3_c + berlin_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

cat("\n==== m_ord_3cat_BxCT: Benefit × CT ====\n")
print(summary(m_ord_3cat_BxCT))

add_result(
  m_ord_3cat_BxCT,
  model_id   = "H2_3cat_BxCT",
  outcome    = "vacc4_cat3",
  predictors = "Benefit_t3_c * CT_c + Risk_t3_c + berlin_c",
  notes      = "polr; Benefit x CT"
)

# Benefit × Berlin
m_ord_3cat_BxBerlin <- polr(
  vacc4_cat3 ~ Benefit_t3_c * berlin_c + Risk_t3_c + CT_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

cat("\n==== m_ord_3cat_BxBerlin: Benefit × Berlin ====\n")
print(summary(m_ord_3cat_BxBerlin))

add_result(
  m_ord_3cat_BxBerlin,
  model_id   = "H3_3cat_BxBerlin",
  outcome    = "vacc4_cat3",
  predictors = "Benefit_t3_c * berlin_c + Risk_t3_c + CT_c",
  notes      = "polr; Benefit x Berlin"
)

# Benefit × CT & Benefit × Berlin
m_ord_3cat_BxCT_BxBer <- polr(
  vacc4_cat3 ~ Benefit_t3_c * CT_c + Benefit_t3_c * berlin_c + Risk_t3_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

cat("\n==== m_ord_3cat_BxCT_BxBer: Benefit × (CT & Berlin) ====\n")
print(summary(m_ord_3cat_BxCT_BxBer))

add_result(
  m_ord_3cat_BxCT_BxBer,
  model_id   = "H4_3cat_BxCT_BxBerlin",
  outcome    = "vacc4_cat3",
  predictors = "Benefit_t3_c * CT_c + Benefit_t3_c * berlin_c + Risk_t3_c",
  notes      = "polr; Benefit x CT & Berlin"
)

## 8-3. Risk 側の交互作用 -------------------------------

# Risk × CT
m_ord_3cat_RxCT <- polr(
  vacc4_cat3 ~ Benefit_t3_c + Risk_t3_c * CT_c + berlin_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

cat("\n==== m_ord_3cat_RxCT: Risk × CT ====\n")
print(summary(m_ord_3cat_RxCT))

add_result(
  m_ord_3cat_RxCT,
  model_id   = "H5_3cat_RxCT",
  outcome    = "vacc4_cat3",
  predictors = "Benefit_t3_c + Risk_t3_c * CT_c + berlin_c",
  notes      = "polr; Risk x CT"
)

# Risk × Berlin
m_ord_3cat_RxBerlin <- polr(
  vacc4_cat3 ~ Benefit_t3_c + Risk_t3_c * berlin_c + CT_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

cat("\n==== m_ord_3cat_RxBerlin: Risk × Berlin ====\n")
print(summary(m_ord_3cat_RxBerlin))

add_result(
  m_ord_3cat_RxBerlin,
  model_id   = "H6_3cat_RxBerlin",
  outcome    = "vacc4_cat3",
  predictors = "Benefit_t3_c + Risk_t3_c * berlin_c + CT_c",
  notes      = "polr; Risk x Berlin"
)

# Risk × CT & Risk × Berlin
m_ord_3cat_RxCT_RxBer <- polr(
  vacc4_cat3 ~ Benefit_t3_c + Risk_t3_c * CT_c + Risk_t3_c * berlin_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

cat("\n==== m_ord_3cat_RxCT_RxBer: Risk × (CT & Berlin) ====\n")
print(summary(m_ord_3cat_RxCT_RxBer))

add_result(
  m_ord_3cat_RxCT_RxBer,
  model_id   = "H7_3cat_RxCT_RxBerlin",
  outcome    = "vacc4_cat3",
  predictors = "Benefit_t3_c + Risk_t3_c * CT_c + Risk_t3_c * berlin_c",
  notes      = "polr; Risk x CT & Berlin"
)

## 8-4. Benefit/Risk 両側で CT/Berlin と交互作用するフルモデル ----------

m_ord_3cat_full <- polr(
  vacc4_cat3 ~ 
    Benefit_t3_c * CT_c +
    Benefit_t3_c * berlin_c +
    Risk_t3_c    * CT_c +
    Risk_t3_c    * berlin_c,
  data = ord_BR_t3t4_3cat_ext,
  Hess = TRUE
)

cat("\n==== m_ord_3cat_full: full interactions ====\n")
print(summary(m_ord_3cat_full))

add_result(
  m_ord_3cat_full,
  model_id   = "H8_3cat_full",
  outcome    = "vacc4_cat3",
  predictors = "Benefit_t3_c * (CT_c + berlin_c) + Risk_t3_c * (CT_c + berlin_c)",
  notes      = "polr; full interaction model"
)

############################################################
# 9. すべてのモデル結果を1つのCSVに書き出し
############################################################

model_results_df <- dplyr::bind_rows(model_results)

cat("\n=== 全モデルをまとめた結果の行数:", nrow(model_results_df), "===\n")

# ファイル名は必要に応じて変更
write.csv(model_results_df,
          file = "all_model_results.csv",
          row.names = FALSE)

cat("結果を 'all_model_results.csv' として書き出しました。\n")
