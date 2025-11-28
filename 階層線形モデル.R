# =========================
# 必要パッケージ
# =========================
required_pkgs <- c(
  "dplyr","tidyr","stringr","lme4","lmerTest",
  "broom.mixed","purrr","tibble"
)
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# R2用（performance または MuMIn を用意）
has_performance <- requireNamespace("performance", quietly = TRUE)
has_mumin       <- requireNamespace("MuMIn",       quietly = TRUE)
if (!has_performance && !has_mumin) {
  install.packages("performance", dependencies = TRUE)
  has_performance <- TRUE
}

get_r2 <- function(fit) {
  if (has_performance) {
    r2 <- performance::r2(fit)
    if (inherits(r2, "data.frame")) {
      c(marginal = r2$R2_marginal[1], conditional = r2$R2_conditional[1])
    } else {
      c(marginal = r2$R2_marginal, conditional = r2$R2_conditional)
    }
  } else if (has_mumin) {
    r2 <- MuMIn::r.squaredGLMM(fit)
    c(marginal = as.numeric(r2[1]), conditional = as.numeric(r2[2]))
  } else {
    c(marginal = NA_real_, conditional = NA_real_)
  }
}

# =========================
# ユーザー設定
# =========================
df <- merged_data_complete
id_col        <- "id"
include_time_fixed <- TRUE
random_intercept_only <- TRUE

# ★ Numeracy 候補（名前=列名）に新指標を追加 ★
numeracy_list <- list(
  mean_subjective_numeracy      = "mean_subjective_numeracy",
  berlin_correct_count          = "berlin_correct_count",
  health_correct_count          = "health_correct_count",
  numAttitude_average           = "numAttitude_average",
  berlin_oneItem                = "berlin_oneItem",                # 0/1
  mean_cliticalThinking_attitude = "mean_cliticalThinking_attitude" # 連続
)

# 追加共変量（必要なら列名をベクターで）
covariates <- NULL
include_acceptance <- FALSE

# =========================
# データ整形：ワイド→ロング & ラグ作成（t=3,4のみ）
# =========================
saf_long <- df |>
  dplyr::select(dplyr::all_of(c(id_col, paste0("Safeness_", 1:4)))) |>
  tidyr::pivot_longer(dplyr::starts_with("Safeness_"),
                      names_to = "time", values_to = "Safeness") |>
  dplyr::mutate(time = as.integer(stringr::str_extract(time, "\\d+")))

eff_long <- df |>
  dplyr::select(dplyr::all_of(c(id_col, paste0("Effectiveness_", 1:4)))) |>
  tidyr::pivot_longer(dplyr::starts_with("Effectiveness_"),
                      names_to = "time", values_to = "Effectiveness") |>
  dplyr::mutate(time = as.integer(stringr::str_extract(time, "\\d+")))

risk_long <- df |>
  dplyr::select(dplyr::all_of(c(id_col, "HBM_average_risk_3", "HBM_average_risk_4"))) |>
  tidyr::pivot_longer(dplyr::starts_with("HBM_average_risk_"),
                      names_to = "time", values_to = "Risk") |>
  dplyr::mutate(time = as.integer(stringr::str_extract(time, "\\d+")))

benefit_long <- df |>
  dplyr::select(dplyr::all_of(c(id_col, "HBM_average_benefit_3", "HBM_average_benefit_4"))) |>
  tidyr::pivot_longer(dplyr::starts_with("HBM_average_benefit_"),
                      names_to = "time", values_to = "Benefit") |>
  dplyr::mutate(time = as.integer(stringr::str_extract(time, "\\d+")))

acc_long <- NULL
if (include_acceptance) {
  acc_long <- df |>
    dplyr::select(dplyr::all_of(c(id_col, "HBM_average_acceptance_3", "HBM_average_acceptance_4"))) |>
    tidyr::pivot_longer(dplyr::starts_with("HBM_average_acceptance_"),
                        names_to = "time", values_to = "Acceptance") |>
    dplyr::mutate(time = as.integer(stringr::str_extract(time, "\\d+")))
}

base_long <- saf_long |>
  dplyr::left_join(eff_long, by = c(id_col, "time")) |>
  dplyr::arrange(.data[[id_col]], time)

lagged <- base_long |>
  dplyr::group_by(.data[[id_col]]) |>
  dplyr::arrange(time, .by_group = TRUE) |>
  dplyr::mutate(
    Safeness_lag      = dplyr::lag(Safeness, 1),
    Effectiveness_lag = dplyr::lag(Effectiveness, 1)
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(time %in% c(3,4)) |>
  dplyr::left_join(risk_long,    by = c(id_col, "time")) |>
  dplyr::left_join(benefit_long, by = c(id_col, "time"))

if (include_acceptance) {
  lagged <- dplyr::left_join(lagged, acc_long, by = c(id_col, "time"))
}

if (!is.null(covariates) && length(covariates)>0) {
  lagged <- dplyr::left_join(lagged, df |> dplyr::select(dplyr::all_of(c(id_col, covariates))), by = id_col)
}

base_needed <- c("Safeness_lag", "Effectiveness_lag", "Risk", "Benefit")
if (include_acceptance) base_needed <- c(base_needed, "Acceptance")
if (!is.null(covariates)) base_needed <- c(base_needed, covariates)

if (include_time_fixed) lagged$time_f <- factor(lagged$time)

rand_part <- if (random_intercept_only) {
  paste0(" + (1 | ", id_col, ")")
} else {
  paste0(" + (1 + time | ", id_col, ")")
}

# =========================
# 1つのNumeracy指標で Risk/Benefit を当てて要約
# =========================
fit_one_numeracy <- function(nm_name, nm_col) {
  dat <- lagged |>
    dplyr::left_join(df |> dplyr::select(dplyr::all_of(c(id_col, nm_col))), by = id_col)
  
  needed <- c(base_needed, nm_col)
  dat <- dat |> dplyr::filter(dplyr::if_all(dplyr::all_of(needed), ~ !is.na(.)))
  
  vals <- unique(na.omit(dat[[nm_col]]))
  is_categorical <- nm_name == "berlin_oneItem" || (length(vals) <= 2 && all(vals %in% c(0,1)))
  
  if (is_categorical) {
    dat$Numeracy_cat <- factor(dat[[nm_col]])
    rhs_core <- "(Safeness_lag * Numeracy_cat) + (Effectiveness_lag * Numeracy_cat)"
  } else {
    dat$Numeracy_c <- scale(dat[[nm_col]], center = TRUE, scale = FALSE)[,1]
    rhs_core <- "(Safeness_lag * Numeracy_c) + (Effectiveness_lag * Numeracy_c)"
  }
  rhs_time <- if (include_time_fixed) " + time_f" else ""
  rhs_cov  <- if (!is.null(covariates) && length(covariates)>0) paste0(" + ", paste(covariates, collapse = " + ")) else ""
  rhs_acc  <- if (include_acceptance) " + Acceptance" else ""
  
  frm_risk    <- as.formula(paste0("Risk ~ ",    rhs_core, rhs_time, rhs_acc, rhs_cov, rand_part))
  frm_benefit <- as.formula(paste0("Benefit ~ ", rhs_core, rhs_time, rhs_acc, rhs_cov, rand_part))
  
  fit_risk    <- lmer(frm_risk,    data = dat, REML = TRUE)
  fit_benefit <- lmer(frm_benefit, data = dat, REML = TRUE)
  
  tidy_r <- broom.mixed::tidy(fit_risk, effects = "fixed", conf.int = FALSE)
  tidy_b <- broom.mixed::tidy(fit_benefit, effects = "fixed", conf.int = FALSE)
  
  grab_inter <- function(tidy_df, var_left, var_right) {
    r1 <- subset(tidy_df, term == paste0(var_left, ":", var_right))
    r2 <- subset(tidy_df, term == paste0(var_right, ":", var_left))
    if (nrow(r1) == 1) r1 else if (nrow(r2) == 1) r2 else {
      data.frame(estimate=NA_real_, p.value=NA_real_, statistic=NA_real_)
    }
  }
  
  if (is_categorical) {
    pick_max_t <- function(tidy_df, pattern) {
      sub <- tidy_df[grep(pattern, tidy_df$term), , drop=FALSE]
      if (nrow(sub) == 0) c(estimate=NA_real_, p.value=NA_real_) else {
        j <- which.max(abs(sub$statistic))
        c(estimate=sub$estimate[j], p.value=sub$p.value[j])
      }
    }
    s_r <- pick_max_t(tidy_r, "^Safeness_lag:Numeracy_cat")
    e_r <- pick_max_t(tidy_r, "^Effectiveness_lag:Numeracy_cat")
    s_b <- pick_max_t(tidy_b, "^Safeness_lag:Numeracy_cat")
    e_b <- pick_max_t(tidy_b, "^Effectiveness_lag:Numeracy_cat")
  } else {
    ir_s <- grab_inter(tidy_r, "Safeness_lag",      "Numeracy_c")
    ir_e <- grab_inter(tidy_r, "Effectiveness_lag", "Numeracy_c")
    ib_s <- grab_inter(tidy_b, "Safeness_lag",      "Numeracy_c")
    ib_e <- grab_inter(tidy_b, "Effectiveness_lag", "Numeracy_c")
    s_r <- c(estimate = ir_s$estimate[1], p.value = ir_s$p.value[1])
    e_r <- c(estimate = ir_e$estimate[1], p.value = ir_e$p.value[1])
    s_b <- c(estimate = ib_s$estimate[1], p.value = ib_s$p.value[1])
    e_b <- c(estimate = ib_e$estimate[1], p.value = ib_e$p.value[1])
  }
  
  aic_r <- AIC(fit_risk);    bic_r <- BIC(fit_risk);    ll_r <- as.numeric(logLik(fit_risk))
  aic_b <- AIC(fit_benefit); bic_b <- BIC(fit_benefit); ll_b <- as.numeric(logLik(fit_benefit))
  
  r2_r <- get_r2(fit_risk);  r2_b <- get_r2(fit_benefit)
  
  tibble::tibble(
    numeracy_name = nm_name,
    numeracy_col  = nm_col,
    numeracy_type = ifelse(is_categorical, "categorical(0/1)", "continuous(centered)"),
    n_obs_risk    = nobs(fit_risk),
    n_obs_benefit = nobs(fit_benefit),
    
    Risk_int_Safe_est = as.numeric(s_r["estimate"]), Risk_int_Safe_p = as.numeric(s_r["p.value"]),
    Risk_int_Eff_est  = as.numeric(e_r["estimate"]), Risk_int_Eff_p  = as.numeric(e_r["p.value"]),
    Risk_AIC = aic_r, Risk_BIC = bic_r, Risk_logLik = ll_r,
    Risk_R2_marginal = r2_r["marginal"], Risk_R2_conditional = r2_r["conditional"],
    
    Benefit_int_Safe_est = as.numeric(s_b["estimate"]), Benefit_int_Safe_p = as.numeric(s_b["p.value"]),
    Benefit_int_Eff_est  = as.numeric(e_b["estimate"]), Benefit_int_Eff_p  = as.numeric(e_b["p.value"]),
    Benefit_AIC = aic_b, Benefit_BIC = bic_b, Benefit_logLik = ll_b,
    Benefit_R2_marginal = r2_b["marginal"], Benefit_R2_conditional = r2_b["conditional"]
  )
}

# =========================
# 全Numeracyで一括フィット＆要約表
# =========================
summary_tbl <- purrr::imap_dfr(numeracy_list, ~ fit_one_numeracy(.y, .x))

# AICランキング
rank_risk    <- summary_tbl |> dplyr::arrange(Risk_AIC)
rank_benefit <- summary_tbl |> dplyr::arrange(Benefit_AIC)

# 交互作用の最小p
rank_risk_p <- summary_tbl |>
  dplyr::mutate(best_p = pmin(dplyr::coalesce(Risk_int_Safe_p, 1),
                              dplyr::coalesce(Risk_int_Eff_p,  1))) |>
  dplyr::arrange(best_p)

rank_benefit_p <- summary_tbl |>
  dplyr::mutate(best_p = pmin(dplyr::coalesce(Benefit_int_Safe_p, 1),
                              dplyr::coalesce(Benefit_int_Eff_p,  1))) |>
  dplyr::arrange(best_p)

# ΔAIC & Akaike重み
add_aic_weights <- function(tab, aic_col) {
  a <- tab[[aic_col]]
  dAIC <- a - min(a, na.rm=TRUE)
  w <- exp(-0.5 * dAIC) / sum(exp(-0.5 * dAIC), na.rm=TRUE)
  dplyr::bind_cols(tab, tibble::tibble(delta_AIC = dAIC, akaike_weight = w))
}
rank_risk_w    <- add_aic_weights(rank_risk,    "Risk_AIC")
rank_benefit_w <- add_aic_weights(rank_benefit, "Benefit_AIC")

# 主効果のみ（交互作用なし）の一括サマリー
library(broom.mixed)
summ_main <- purrr::imap_dfr(numeracy_list, function(nm_col, nm_name){
  dat <- lagged |>
    dplyr::left_join(df |> dplyr::select(dplyr::all_of(c(id_col, nm_col))), by = "id") |>
    dplyr::filter(dplyr::if_all(dplyr::all_of(c("Safeness_lag","Effectiveness_lag","Risk","Benefit",nm_col)), ~ !is.na(.)))
  
  vals <- unique(na.omit(dat[[nm_col]]))
  is_cat <- nm_name=="berlin_oneItem" || (length(vals) <= 2 && all(vals %in% c(0,1)))
  
  if (is_cat) {
    dat$Numeracy <- factor(dat[[nm_col]])
  } else {
    dat$Numeracy <- scale(dat[[nm_col]], center = TRUE, scale = FALSE)[,1]
  }
  
  fr <- lmer(Risk    ~ Safeness_lag + Effectiveness_lag + Numeracy + time_f + (1|id), data=dat, REML=TRUE)
  fb <- lmer(Benefit ~ Safeness_lag + Effectiveness_lag + Numeracy + time_f + (1|id), data=dat, REML=TRUE)
  
  tr <- tidy(fr, effects="fixed")
  tb <- tidy(fb, effects="fixed")
  
  tibble::tibble(
    numeracy_name    = nm_name,
    numeracy_col     = nm_col,
    numeracy_type    = ifelse(is_cat, "categorical(0/1)", "continuous(centered)"),
    Risk_num_est     = tr$estimate[tr$term=="Numeracy"],
    Risk_num_p       = tr$p.value[tr$term=="Numeracy"],
    Benefit_num_est  = tb$estimate[tb$term=="Numeracy"],
    Benefit_num_p    = tb$p.value[tb$term=="Numeracy"],
    Risk_AIC_main    = AIC(fr),
    Benefit_AIC_main = AIC(fb)
  )
})

# =========================
# CSV 出力
# =========================
utils::write.csv(summary_tbl,      "numeracy_model_summary.csv",                   row.names = FALSE)
utils::write.csv(rank_risk_w,      "numeracy_rank_by_AIC_risk_with_weights.csv",   row.names = FALSE)
utils::write.csv(rank_benefit_w,   "numeracy_rank_by_AIC_benefit_with_weights.csv",row.names = FALSE)
utils::write.csv(rank_risk_p,      "numeracy_rank_by_p_risk.csv",                  row.names = FALSE)
utils::write.csv(rank_benefit_p,   "numeracy_rank_by_p_benefit.csv",               row.names = FALSE)
utils::write.csv(summ_main,        "numeracy_main_effects_summary.csv",            row.names = FALSE)

# =========================
# 画面表示（サマリー）
# =========================
cat("\n=== Summary (first rows) ===\n"); print(utils::head(summary_tbl, 6))
cat("\n=== Rank by AIC: Risk (with ΔAIC & weights) ===\n"); print(rank_risk_w[, c("numeracy_name","Risk_AIC","delta_AIC","akaike_weight")])
cat("\n=== Rank by AIC: Benefit (with ΔAIC & weights) ===\n"); print(rank_benefit_w[, c("numeracy_name","Benefit_AIC","delta_AIC","akaike_weight")])
cat("\n=== Rank by min p (Risk interactions) ===\n"); print(rank_risk_p[, c("numeracy_name","Risk_int_Safe_p","Risk_int_Eff_p","best_p")])
cat("\n=== Rank by min p (Benefit interactions) ===\n"); print(rank_benefit_p[, c("numeracy_name","Benefit_int_Safe_p","Benefit_int_Eff_p","best_p")])
cat("\n=== Main-effect only summary ===\n"); print(summ_main)
