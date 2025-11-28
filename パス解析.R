
# ベースモデルテンプレート（<dep>なし）
model <- "
  HBM_average_risk_3 ~ averageNumerical
  HBM_average_risk_3 ~ averageNonNumerical
  HBM_average_benefit_3 ~ averageNumerical
  HBM_average_benefit_3 ~ averageNonNumerical
  HBM_average_risk_3 ~~ HBM_average_benefit_3
"

#数値態度、アクセプタンス3時点目
model <- '
  HBM_average_risk_3 ~ averageNumerical
  HBM_average_risk_3 ~ averageNonNumerical
  HBM_average_benefit_3 ~ averageNumerical
  HBM_average_benefit_3 ~ averageNonNumerical
  num_vaccination_3 ~ HBM_average_benefit_3
  num_vaccination_3 ~ HBM_average_risk_3
      # 共分散（相関）の追加
  HBM_average_risk_3 ~~ HBM_average_benefit_3
'

#数値態度、アクセプタンス3時点目
  model <- '
    averageNumerical ~ numAttitude_average
    averageNonNumerical ~ numAttitude_average
    HBM_average_risk_3 ~ averageNumerical
    HBM_average_risk_3 ~ averageNonNumerical
    HBM_average_benefit_3 ~ averageNumerical
    HBM_average_benefit_3 ~ averageNonNumerical
    num_vaccination_3 ~ HBM_average_benefit_3
    num_vaccination_3 ~ HBM_average_risk_3
        # 共分散（相関）の追加
    averageNumerical ~~ averageNonNumerical
    HBM_average_risk_3 ~~ HBM_average_benefit_3
  '

# リスクベネフィット3時点目、接種回数3時点目
model <- '
  averageNumerical ~ mean_subjective_numeracy
  averageNonNumerical ~ mean_subjective_numeracy
  HBM_average_risk_3 ~ averageNumerical + averageNonNumerical
  HBM_average_benefit_3 ~ averageNumerical + averageNonNumerical
  num_vaccination_3 ~ HBM_average_benefit_3 + HBM_average_risk_3
    # 共分散（相関）の追加
  averageNumerical ~~ averageNonNumerical
  HBM_average_risk_3 ~~ HBM_average_benefit_3
'



#多母集団同時分析パス解析の実行v2


# -----------------------------
# 手順2：母集団ごとの構成モデル（configural model）
# -----------------------------
fit_configural <- sem(model, 
                      data = merged_sets$df1_df3, 
                      group = "berlin_group.y") #subjectiveNumeracy_group #berlin_oneItem #berlin_group #numAttitude_group #berlin_group.y

summary(fit_configural, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

# -----------------------------
# 手順3：構成不変性（測定不変性、参考程度）
# -----------------------------
# ※このモデルは測定モデルではないので group.equal = "loadings" は非推奨（警告出る）
# 代わりに以下のように構造パス制約モデルで比較するのが実践的

# -----------------------------
# 手順4：構造パスの等値制約の差分検定
# -----------------------------

# 4-1. 全ての構造パスを等値制約したモデル（比較基準モデル）
fit_all_equal <- sem(model, 
                     data = merged_sets$df1_df3, 
                     group = "berlin_group.y", #subjectiveNumeracy_group #berlin_oneItem #berlin_group #numAttitude_group
                     group.equal = "regressions")

# 4-2. モデル式から回帰パスだけ抽出（~~は除外）
model_lines <- strsplit(model, "\n")[[1]]
model_lines <- trimws(model_lines)
model_lines <- model_lines[model_lines != "" & !grepl("^#", model_lines)]
paths_to_test <- grep("~", model_lines, value = TRUE)
paths_to_test <- paths_to_test[!grepl("~~", paths_to_test)]

# 4-3. 各パスごとに1本だけ自由にしたモデルと差分検定
results <- lapply(paths_to_test, function(path_label) {
  fit_partial <- sem(model,
                     data = merged_sets$df1_df3,
                     group = "berlin_group.y", #subjectiveNumeracy_group #berlin_oneItem #berlin_group #numAttitude_group
                     group.equal = "regressions",        # すべて等値制約
                     group.partial = path_label)         # 1本だけ除外
  
  anova_result <- anova(fit_all_equal, fit_partial)
  
  data.frame(
    path = path_label,
    chisq_diff = anova_result$Chisq[2] - anova_result$Chisq[1],
    df_diff = anova_result$Df[2] - anova_result$Df[1],
    p_value = anova_result$`Pr(>Chisq)`[2]
  )
})

# 4-4. 検定結果を表示
results_df <- do.call(rbind, results)
print(results_df)




# 多母集団パス解析の実行
fit <- sem(model, data = merged_data_complete_for_path, group="subjectiveNumeracy_group") #subjectiveNumeracy_group #berlin_oneItem #berlin_group #numAttitude_group
# 結果の表示
summary(fit, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
#決定係数
inspect(fit, what="rsquare")


# 多母集団パス係数の差の検定
# 1. 自由モデル(すべてのパス自由)
fit_free <- sem(model, data = merged_data_complete_for_path, group = "subjectiveNumeracy_group")#subjectiveNumeracy_group #berlin_oneItem #berlin_group #numAttitude_group
# 2. 回帰式のパス名のみを抽出
model_lines <- strsplit(model, "\n")[[1]]
# 空行とコメント行を除去
model_lines <- trimws(model_lines)
model_lines <- model_lines[model_lines != "" & !grepl("^#", model_lines)]
# 回帰式のみを抽出（~ を含むが ~~ を含まない行）
paths_to_test <- grep("~", model_lines, value = TRUE)
paths_to_test <- paths_to_test[!grepl("~~", paths_to_test)]
print("抽出された回帰式パス:")
print(paths_to_test)
# 3. 差分検定のループ
results <- lapply(paths_to_test, function(path_label) {
  fit_constrained <- sem(model, data = merged_data_complete_for_path,
                         group = "subjectiveNumeracy_group", #subjectiveNumeracy_group #berlin_oneItem #berlin_group #numAttitude_group
                         group.equal = c("regressions"),
                         group.partial = path_label) # その1本だけ等しくする
  
  anova_result <- anova(fit_free, fit_constrained)
  data.frame(
    path = path_label,
    chisq_diff = anova_result$Chisq[2] - anova_result$Chisq[1],
    df_diff = anova_result$Df[2] - anova_result$Df[1],
    p_value = anova_result$`Pr(>Chisq)`[2]
  )
})
# 4. 結果をまとめて表示
do.call(rbind, results)
parameterEstimates(fit_free)