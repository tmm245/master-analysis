# ============================================================
# Longitudinal LMM/HLM 健全性チェック一括スクリプト
# 前提: データフレーム merged_data_complete が環境に存在
# 必須列:
#   id,
#   Safeness_1:Safeness_4, Effectiveness_1:Effectiveness_4,
#   HBM_average_risk_3:4, HBM_average_benefit_3:4
#   (交互作用例示用) mean_cliticalThinking_attitude
# 出力: ./_checks/ 以下にCSV/PNGを保存
# ============================================================

# ---- 0) パッケージ ----
need <- c("dplyr","tidyr","stringr","lme4","lmerTest","performance","DHARMa","misty","ggplot2","broom.mixed","readr")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(need, library, character.only = TRUE))

dir.create("_checks", showWarnings = FALSE)

# ---- 1) 基本メタ & ID一貫性 ----
df <- merged_data_complete
stopifnot(all(c("id") %in% names(df)))
message("# rows: ", nrow(df), " | # unique id: ", dplyr::n_distinct(df$id))
dup_id <- df %>% count(id) %>% filter(n > 1)
if (nrow(dup_id)) {
  warning("重複IDが存在します。_checks/duplicate_ids.csv を確認してください。")
  readr::write_csv(dup_id, "_checks/duplicate_ids.csv")
}

# ---- 2) ワイド→ロング & ラグ作成（t=3/4）----
id_col <- "id"

saf_long <- df |>
  select(all_of(c(id_col, paste0("Safeness_", 1:4)))) |>
  pivot_longer(starts_with("Safeness_"), names_to = "time", values_to = "Safeness") |>
  mutate(time = as.integer(str_extract(time, "\\d+")))

eff_long <- df |>
  select(all_of(c(id_col, paste0("Effectiveness_", 1:4)))) |>
  pivot_longer(starts_with("Effectiveness_"), names_to = "time", values_to = "Effectiveness") |>
  mutate(time = as.integer(str_extract(time, "\\d+")))

risk_long <- df |>
  select(all_of(c(id_col, "HBM_average_risk_3","HBM_average_risk_4"))) |>
  pivot_longer(starts_with("HBM_average_risk_"), names_to = "time", values_to = "Risk") |>
  mutate(time = as.integer(str_extract(time, "\\d+")))

benefit_long <- df |>
  select(all_of(c(id_col, "HBM_average_benefit_3","HBM_average_benefit_4"))) |>
  pivot_longer(starts_with("HBM_average_benefit_"), names_to = "time", values_to = "Benefit") |>
  mutate(time = as.integer(str_extract(time, "\\d+")))

base_long <- saf_long %>% left_join(eff_long, by = c(id_col,"time")) %>% arrange(.data[[id_col]], time)

lagged <- base_long %>%
  group_by(.data[[id_col]]) %>%
  arrange(time, .by_group = TRUE) %>%
  mutate(
    Safeness_lag      = dplyr::lag(Safeness, 1),
    Effectiveness_lag = dplyr::lag(Effectiveness, 1)
  ) %>%
  ungroup() %>%
  filter(time %in% c(3,4)) %>%
  left_join(risk_long,    by = c(id_col,"time")) %>%
  left_join(benefit_long, by = c(id_col,"time")) %>%
  mutate(time_f = factor(time))

# ---- 3) 時点別Nと欠測パターン ----
# t3/t4 の観測状況
tab_obs <- lagged %>%
  summarize(
    n_id            = n_distinct(.data[[id_col]]),
    n_obs           = n(),
    n_t3            = sum(time == 3),
    n_t4            = sum(time == 4),
    n_complete_both = sum(table(.data[[id_col]], time) %>% apply(1, function(x) all(c(3,4) %in% which(x>0))))
  )
readr::write_csv(tab_obs, "_checks/longitudinal_counts_summary.csv")

# id × time の有無テーブル
id_time <- lagged %>% count(.data[[id_col]], time) %>% tidyr::pivot_wider(names_from = time, values_from = n, values_fill = 0)
readr::write_csv(id_time, "_checks/id_time_presence.csv")

# 欠測パターン（主要列のみ）
vars_check <- c("Safeness_lag","Effectiveness_lag","Risk","Benefit")
lagged$miss_pattern <- apply(is.na(lagged[,vars_check]), 1, function(r) paste0(as.integer(r), collapse=""))
miss_table <- lagged %>% count(miss_pattern, time) %>% arrange(desc(n))
readr::write_csv(miss_table, "_checks/missing_patterns_by_time.csv")

# ---- 4) MCAR（Littleのテスト） ※misty::na.test を使用 ----
# timeごとに主要変数のMCARをざっくり確認（数値のみ）
mcar_list <- list()
for (tt in c(3,4)) {
  dat_tt <- lagged %>% filter(time == tt) %>% select(all_of(vars_check))
  # 変数数 >=2 & 欠測ありの時のみ実行
  if (sum(is.na(as.matrix(dat_tt))) > 0 && ncol(dat_tt) >= 2) {
    mcar <- tryCatch(misty::na.test(dat_tt, method = "little"), error = function(e) NULL)
    mcar_list[[as.character(tt)]] <- mcar
  }
}
sink("_checks/mcar_little_results.txt"); print(mcar_list); sink()
# 参考: LittleのMCARテスト（misty::na.test）。:contentReference[oaicite:1]{index=1}

# ---- 5) LMMのベースモデル（主効果のみ） & 適合指標 ----
# Risk
fit_risk_main <- lmer(Risk ~ Safeness_lag + Effectiveness_lag + time_f + (1|id), data = lagged, REML = TRUE)
# Benefit
fit_ben_main  <- lmer(Benefit ~ Safeness_lag + Effectiveness_lag + time_f + (1|id), data = lagged, REML = TRUE)

# AIC/BIC/対数尤度, R2(Nakagawa), ICC
summ_models <- tibble::tibble(
  model = c("Risk_main","Benefit_main"),
  AIC   = c(AIC(fit_risk_main), AIC(fit_ben_main)),
  BIC   = c(BIC(fit_risk_main), BIC(fit_ben_main)),
  logLik= c(as.numeric(logLik(fit_risk_main)), as.numeric(logLik(fit_ben_main)))
)

r2_risk <- performance::r2_nakagawa(fit_risk_main) # marginal/conditional R²（Nakagawa）:contentReference[oaicite:2]{index=2}
r2_ben  <- performance::r2_nakagawa(fit_ben_main)
icc_r   <- performance::icc(fit_risk_main)         # ICC（混合効果の分散比）:contentReference[oaicite:3]{index=3}
icc_b   <- performance::icc(fit_ben_main)

sink("_checks/model_fit_summary.txt")
cat("\n== Base LMM (Risk) ==\n"); print(summary(fit_risk_main)); print(r2_risk); print(icc_r)
cat("\n== Base LMM (Benefit) ==\n"); print(summary(fit_ben_main));  print(r2_ben);  print(icc_b)
sink()
readr::write_csv(summ_models, "_checks/model_fit_table.csv")

# ---- 6) 多重共線性（VIF相当）/ 収束 / シンギュラリティ ----
sink("_checks/diagnostics_collinearity_singularity.txt")

cat("\n== Collinearity (Risk) ==\n")
print(performance::check_collinearity(fit_risk_main))

cat("\n== Collinearity (Benefit) ==\n")
print(performance::check_collinearity(fit_ben_main))

cat("\n== Convergence ==\n")
print(performance::check_convergence(fit_risk_main))
print(performance::check_convergence(fit_ben_main))

cat("\n== Singularity ==\n")
# 正：check_singularity() を使う（旧: check_singular は存在しない）
print(performance::check_singularity(fit_risk_main))
print(performance::check_singularity(fit_ben_main))

# 参考：lme4ネイティブの簡易判定（補助用）
cat("\n== lme4::isSingular (fallback) ==\n")
print(lme4::isSingular(fit_risk_main, tol = 1e-04))
print(lme4::isSingular(fit_ben_main,  tol = 1e-04))

sink()


# ---- 7) 残差診断（DHARMa & performance::check_model）----
# DHARMa（ガウスLMMでも利用可）:contentReference[oaicite:5]{index=5}
res_risk <- DHARMa::simulateResiduals(fittedModel = fit_risk_main, n = 1000)
png("_checks/DHARMa_residuals_Risk.png", width = 900, height = 600); plot(res_risk); dev.off()
res_ben  <- DHARMa::simulateResiduals(fittedModel = fit_ben_main, n = 1000)
png("_checks/DHARMa_residuals_Benefit.png", width = 900, height = 600); plot(res_ben); dev.off()

# check_model 可視診断（残差正規性・等分散・線形性・ランダム効果正規性など）:contentReference[oaicite:6]{index=6}
png("_checks/check_model_Risk.png", width = 1200, height = 900);  performance::check_model(fit_risk_main);  dev.off()
png("_checks/check_model_Benefit.png", width = 1200, height = 900); performance::check_model(fit_ben_main); dev.off()

# ---- 8) 交互作用の要否: LRT（例: Benefit ~ Effectiveness_lag * CT）----
ct_col <- "mean_cliticalThinking_attitude"  # 必要に応じて列名調整
if (ct_col %in% names(df)) {
  dat_ct <- lagged %>%
    left_join(df %>% select(id, all_of(ct_col)), by = "id") %>%
    filter(if_all(all_of(c("Safeness_lag","Effectiveness_lag","Benefit",ct_col)), ~ !is.na(.))) %>%
    mutate(CT_c = scale(.data[[ct_col]], center = TRUE, scale = TRUE)[,1])
  
  m_main <- lmer(Benefit ~ Safeness_lag + Effectiveness_lag + CT_c + time_f + (1|id),
                 data = dat_ct, REML = FALSE)
  m_int  <- lmer(Benefit ~ Safeness_lag + Effectiveness_lag*CT_c + time_f + (1|id),
                 data = dat_ct, REML = FALSE)
  
  lrt <- anova(m_main, m_int) # 交互作用の寄与（尤度比）を検定
  sink("_checks/LRT_Benefit_CT_interaction.txt")
  cat("== LRT: add Effectiveness_lag × CT_c ==\n"); print(lrt)
  cat("\nAIC main / int: ", AIC(m_main), " / ", AIC(m_int), "\n")
  cat("\nR2 (Nakagawa):\n"); print(performance::r2_nakagawa(m_main)); print(performance::r2_nakagawa(m_int))
  sink()
}

# ---- 9) サンプルサイズの把握（検出力は別途 powerlmm 等で）----
n_subjects <- dplyr::n_distinct(lagged$id)
n_obs      <- nrow(lagged)
obs_per_id <- lagged %>% count(id) %>% summarize(mean_obs = mean(n), sd_obs = sd(n), min_obs = min(n), max_obs = max(n))
sink("_checks/sample_size_summary.txt")
cat("N subjects: ", n_subjects, "\nN observations: ", n_obs, "\n"); print(obs_per_id)
cat("\n時点別 N:\n"); print(table(lagged$time))
sink()

message("✅ 完了: 出力は ./_checks/ に保存されました。")
