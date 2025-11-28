
# 分析データ
data <- merged_data_complete_for_path

# モデル情報
models_info <- list(
  list(dep = "HBM_average_acceptance_3", ind = "mean_subjective_numeracy", grp = "berlin_oneItem"),
  list(dep = "HBM_average_acceptance_3", ind = "mean_subjective_numeracy", grp = "berlin_group"),
  list(dep = "HBM_average_acceptance_3", ind = "mean_subjective_numeracy", grp = "numAttitude_group"),
  list(dep = "HBM_average_acceptance_3", ind = "numAttitude_average", grp = "subjectiveNumeracy_group"),
  list(dep = "HBM_average_acceptance_3", ind = "numAttitude_average", grp = "berlin_oneItem"),
  list(dep = "HBM_average_acceptance_3", ind = "numAttitude_average", grp = "berlin_group"),
  list(dep = "HBM_average_acceptance_4", ind = "mean_subjective_numeracy", grp = "berlin_oneItem"),
  list(dep = "HBM_average_acceptance_4", ind = "mean_subjective_numeracy", grp = "berlin_group"),
  list(dep = "HBM_average_acceptance_4", ind = "mean_subjective_numeracy", grp = "numAttitude_group"),
  list(dep = "HBM_average_acceptance_4", ind = "numAttitude_average", grp = "subjectiveNumeracy_group"),
  list(dep = "HBM_average_acceptance_4", ind = "numAttitude_average", grp = "berlin_oneItem"),
  list(dep = "HBM_average_acceptance_4", ind = "numAttitude_average", grp = "berlin_group"),
  list(dep = "num_vaccination_3", ind = "mean_subjective_numeracy", grp = "berlin_oneItem"),
  list(dep = "num_vaccination_3", ind = "mean_subjective_numeracy", grp = "berlin_group"),
  list(dep = "num_vaccination_3", ind = "mean_subjective_numeracy", grp = "numAttitude_group"),
  list(dep = "num_vaccination_3", ind = "numAttitude_average", grp = "subjectiveNumeracy_group"),
  list(dep = "num_vaccination_3", ind = "numAttitude_average", grp = "berlin_oneItem"),
  list(dep = "num_vaccination_3", ind = "numAttitude_average", grp = "berlin_group"),
  list(dep = "num_vaccination_4", ind = "mean_subjective_numeracy", grp = "berlin_oneItem"),
  list(dep = "num_vaccination_4", ind = "mean_subjective_numeracy", grp = "berlin_group"),
  list(dep = "num_vaccination_4", ind = "mean_subjective_numeracy", grp = "numAttitude_group"),
  list(dep = "num_vaccination_4", ind = "numAttitude_average", grp = "subjectiveNumeracy_group"),
  list(dep = "num_vaccination_4", ind = "numAttitude_average", grp = "berlin_oneItem"),
  list(dep = "num_vaccination_4", ind = "numAttitude_average", grp = "berlin_group")
)

# R2の安全な抽出関数（リスト or ベクトル 対応版）
safe_r2 <- function(r2, varname) {
  # マルチグループなら、まず第1グループを取り出す
  if (is.list(r2)) {
    r2 <- r2[[1]]
  }
  if (!is.null(r2) && varname %in% names(r2)) {
    return(round(r2[[varname]], 3))
  } else {
    return(NA_real_)
  }
}

# 結果格納リスト
all_results <- vector("list", length(models_info))

# ループで24モデルを一括実行
for (i in seq_along(models_info)) {
  info <- models_info[[i]]
  dep <- info$dep
  ind <- info$ind
  grp <- info$grp
  
  cat("\n==============================\n")
  cat(paste("モデル", i, ":", dep, "~", ind, "(Group =", grp, ")\n"))
  cat("==============================\n")
  
  # lavaanモデル構文を組み立て
  model_syntax <- paste0("
    averageNumerical       ~ ", ind, "
    averageNonNumerical    ~ ", ind, "
    HBM_average_risk_3     ~ averageNumerical
    HBM_average_risk_3     ~ averageNonNumerical
    HBM_average_benefit_3  ~ averageNumerical
    HBM_average_benefit_3  ~ averageNonNumerical
    ", dep, "               ~ HBM_average_benefit_3
    ", dep, "               ~ HBM_average_risk_3
    averageNumerical ~~ averageNonNumerical
    HBM_average_risk_3 ~~ HBM_average_benefit_3
  ")
  
  # 1. Configuralモデルの推定
  fit_config <- sem(model_syntax, data = data, group = grp)
  
  # 適合度指標の抽出
  fm <- fitMeasures(fit_config,
                    c("chisq", "df", "pvalue",
                      "cfi", "tli", "rmsea", "srmr"))
  
  # R2の抽出
  r2_vals <- inspect(fit_config, "rsquare")
  
  # 全回帰パス等値モデル＆差分検定の準備
  fit_all_equal <- sem(model_syntax, data = data,
                       group = grp, group.equal = "regressions")
  
  # モデル内の回帰パスを抽出
  model_lines <- trimws(strsplit(model_syntax, "\n")[[1]])
  paths_to_test <- grep("~", model_lines, value = TRUE)
  paths_to_test <- paths_to_test[!grepl("~~", paths_to_test)]
  
  # 各パスごとに差分検定
  path_tests <- lapply(paths_to_test, function(path_label) {
    fit_partial <- sem(model_syntax, data = data, group = grp,
                       group.equal = "regressions",
                       group.partial = path_label)
    LRT <- anova(fit_all_equal, fit_partial)
    p <- LRT$`Pr(>Chisq)`[2]
    group_ns <- table(data[[grp]])
    n_info <- paste(names(group_ns), group_ns, sep = "=", collapse = "; ")
    data.frame(
      model_id               = i,
      dependent               = dep,
      independent             = ind,
      group_var               = grp,
      group_n                = n_info,
      fit_chisq               = round(fm["chisq"], 3),
      fit_df                  =     fm["df"],
      fit_p                   = round(fm["pvalue"], 4),
      fit_cfi                 = round(fm["cfi"], 3),
      fit_tli                 = round(fm["tli"], 3),
      fit_rmsea               = round(fm["rmsea"], 3),
      fit_srmr                = round(fm["srmr"], 3),
      r2_averageNumerical    = safe_r2(r2_vals, "averageNumerical"),
      r2_averageNonNumerical = safe_r2(r2_vals, "averageNonNumerical"),
      r2_HBM_risk            = safe_r2(r2_vals, "HBM_average_risk_3"),
      r2_HBM_benefit         = safe_r2(r2_vals, "HBM_average_benefit_3"),
      r2_dep                 = safe_r2(r2_vals, dep),
      path                    = path_label,
      chisq_diff              = round(LRT$Chisq[2] - LRT$Chisq[1], 3),
      df_diff                 = LRT$Df[2] - LRT$Df[1],
      p_value                 = round(p, 4),
      note                    = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        TRUE     ~ ""
      ),
      stringsAsFactors = FALSE
    )
  })
  
  # 結果をリストに格納
  all_results[[i]] <- bind_rows(path_tests)
}

# CSV書き出し
final_results <- bind_rows(all_results)
write.csv(final_results, "24multi_group_sem_24models_results.csv", row.names = FALSE)




# 条件に合うモデルID（CFI, TLI, RMSEA, SRMR, p値 ??? 0.05）を抽出
valid_model_ids <- final_results %>%
  group_by(model_id) %>%
  filter(
    fit_cfi >= 0.95,
    fit_tli >= 0.95,
    fit_rmsea <= 0.10,
    fit_srmr <= 0.08,
    fit_p >= 0.05       # ??? χ2検定のp値が有意でない（モデル適合良好）
  ) %>%
  summarise(.groups = "drop") %>%
  pull(model_id)

# 条件を満たすモデルの結果を抽出
filtered_results <- final_results %>%
  filter(model_id %in% valid_model_ids)

# CSVとして保存
write_csv(filtered_results, "24selected_models_from_final_results.csv")
