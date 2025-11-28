# 分析データ
data <- merged_sets$df1_df3

# モデル情報（simpleModelPatterns.xlsxの内容を手動定義）
base_models_info <- list(
  list(dep = "HBM_average_acceptance_3", grp = "berlin_oneItem"),
  list(dep = "HBM_average_acceptance_3", grp = "berlin_group"),
  list(dep = "HBM_average_acceptance_3", grp = "numAttitude_group"),
  list(dep = "HBM_average_acceptance_3", grp = "subjectiveNumeracy_group"),
  list(dep = "HBM_average_acceptance_4", grp = "berlin_oneItem"),
  list(dep = "HBM_average_acceptance_4", grp = "berlin_group"),
  list(dep = "HBM_average_acceptance_4", grp = "numAttitude_group"),
  list(dep = "HBM_average_acceptance_4", grp = "subjectiveNumeracy_group"),
  list(dep = "num_vaccination_3", grp = "berlin_oneItem"),
  list(dep = "num_vaccination_3", grp = "berlin_group"),
  list(dep = "num_vaccination_3", grp = "numAttitude_group"),
  list(dep = "num_vaccination_3", grp = "subjectiveNumeracy_group"),
  list(dep = "num_vaccination_4", grp = "berlin_oneItem"),
  list(dep = "num_vaccination_4", grp = "berlin_group"),
  list(dep = "num_vaccination_4", grp = "numAttitude_group"),
  list(dep = "num_vaccination_4", grp = "subjectiveNumeracy_group")
)

# 3グループの組み合わせ
# berlin_comparisons <- list(
#   c("BL", "BM"),
#   c("BL", "BH"),
#   c("BM", "BH")
# )

#2グループの組み合わせ
berlin_comparisons <- list(
  c("BL", "BH")
)

# 展開してモデル情報を再構成
models_info <- list()
for (base in base_models_info) {
  if (base$grp == "berlin_group") {
    for (pair in berlin_comparisons) {
      models_info[[length(models_info)+1]] <- list(
        dep = base$dep,
        grp = base$grp,
        levels = pair
      )
    }
  } else {
    models_info[[length(models_info)+1]] <- list(
      dep = base$dep,
      grp = base$grp,
      levels = NULL
    )
  }
}

# R2の安全な抽出関数
safe_r2 <- function(r2, varname) {
  if (is.list(r2)) {
    r2 <- r2[[1]]
  }
  if (!is.null(r2) && varname %in% names(r2)) {
    return(round(r2[[varname]], 3))
  } else {
    return(NA_real_)
  }
}

# ベースモデルテンプレート
base_model_template <- "
  HBM_average_risk_3 ~ averageNumerical
  HBM_average_risk_3 ~ averageNonNumerical
  HBM_average_benefit_3 ~ averageNumerical
  HBM_average_benefit_3 ~ averageNonNumerical
  <dep> ~ HBM_average_benefit_3
  <dep> ~ HBM_average_risk_3
  HBM_average_risk_3 ~~ HBM_average_benefit_3
"

# 結果格納リスト
all_results <- vector("list", length(models_info))

# モデル一括実行ループ
for (i in seq_along(models_info)) {
  info <- models_info[[i]]
  dep <- info$dep
  grp <- info$grp
  subset_data <- data
  
  if (!is.null(info$levels)) {
    subset_data <- subset(data, get(grp) %in% info$levels)
    subset_data[[grp]] <- droplevels(factor(subset_data[[grp]], levels = info$levels))
  }
  
  cat("\n==============================\n")
  cat(paste("モデル", i, ": dependent =", dep, ", group =", grp,
            if (!is.null(info$levels)) paste0(" (", paste(info$levels, collapse = "/"), ")") else "", "\n"))
  cat("==============================\n")
  
  model_syntax <- gsub("<dep>", dep, base_model_template)
  
  fit_config <- sem(model_syntax, data = subset_data, group = grp)
  fm <- fitMeasures(fit_config, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))
  r2_vals <- inspect(fit_config, "rsquare")
  
  fit_all_equal <- sem(model_syntax, data = subset_data, group = grp, group.equal = "regressions")
  model_lines <- trimws(strsplit(model_syntax, "\n")[[1]])
  paths_to_test <- grep("~", model_lines, value = TRUE)
  paths_to_test <- paths_to_test[!grepl("~~", paths_to_test)]
  
  path_tests <- lapply(paths_to_test, function(path_label) {
    fit_partial <- sem(model_syntax, data = subset_data, group = grp,
                       group.equal = "regressions",
                       group.partial = path_label)
    LRT <- anova(fit_all_equal, fit_partial)
    p <- LRT$`Pr(>Chisq)`[2]
    data.frame(
      model_id               = i,
      dependent              = dep,
      group_var              = grp,
      group_levels           = if (is.null(info$levels)) NA else paste(info$levels, collapse = "/"),
      fit_chisq              = round(fm["chisq"], 3),
      fit_df                 = fm["df"],
      fit_p                  = round(fm["pvalue"], 4),
      fit_cfi                = round(fm["cfi"], 3),
      fit_tli                = round(fm["tli"], 3),
      fit_rmsea              = round(fm["rmsea"], 3),
      fit_srmr               = round(fm["srmr"], 3),
      r2_averageNumerical    = safe_r2(r2_vals, "averageNumerical"),
      r2_averageNonNumerical = safe_r2(r2_vals, "averageNonNumerical"),
      r2_HBM_risk            = safe_r2(r2_vals, "HBM_average_risk_3"),
      r2_HBM_benefit         = safe_r2(r2_vals, "HBM_average_benefit_3"),
      r2_dep                 = safe_r2(r2_vals, dep),
      path                   = path_label,
      chisq_diff             = round(LRT$Chisq[2] - LRT$Chisq[1], 3),
      df_diff                = LRT$Df[2] - LRT$Df[1],
      p_value                = round(p, 4),
      note                   = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        TRUE      ~ ""
      ),
      stringsAsFactors = FALSE
    )
  })
  
  all_results[[i]] <- bind_rows(path_tests)
}

# 結果の統合と出力
final_results <- bind_rows(all_results)
write.csv(final_results, "multi_group_sem_simpleModels_results.csv", row.names = FALSE)
cat("\n??? CSVファイル 'multi_group_sem_simpleModels_results.csv' を出力しました。\n")

# 適合度の良いモデルのみ抽出
valid_model_ids <- final_results %>%
  group_by(model_id) %>%
  filter(
    fit_cfi >= 0.95,
    fit_tli >= 0.95,
    fit_rmsea <= 0.10,
    fit_srmr <= 0.08,
    fit_p >= 0.05
  ) %>%
  summarise(.groups = "drop") %>%
  pull(model_id)

filtered_results <- final_results %>%
  filter(model_id %in% valid_model_ids)

write_csv(filtered_results, "selected_simpleModels_from_final_results.csv")
cat("??? 条件を満たすモデルのみを抽出し、'selected_simpleModels_from_final_results.csv' に保存しました。\n")
