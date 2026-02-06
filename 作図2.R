  # 作図用パッケージ
  if (!require("ggeffects")) install.packages("ggeffects")
  if (!require("ggplot2")) install.packages("ggplot2")
  if (!require("broom.mixed")) install.packages("broom.mixed") # 係数抽出用
  
  library(ggeffects)
  library(ggplot2)
  library(broom.mixed)
  library(dplyr)
  
  # ==============================================================================
  # 1. 交互作用プロット: Safeness × Wave (plot_safe_wave)
  # ==============================================================================
  # Benefitモデル(m_ben_wave)を使用
  # terms = c("X軸", "グループ分け(色)")
  pred_safe_wave <- ggpredict(m_ben_wave, terms = c("Safe_c", "wave"))
  
  plot_safe_wave <- ggplot(pred_safe_wave, aes(x = x, y = predicted, color = group)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), 
                alpha = 0.15, linetype = 0) +
    geom_line(linewidth = 1) +
    scale_color_brewer(palette = "Set1", name = "Wave") +
    scale_fill_brewer(palette = "Set1", name = "Wave") +
    labs(
      title = "Interaction: Safeness effect by Wave",
      # subtitle を削除しました
      x = "Safeness (Centered)",
      y = "Predicted Benefit"
    ) +
    theme_classic()
  
  print(plot_safe_wave)
  
  # ==============================================================================
  # 2. 交互作用プロット: Effectiveness × Critical Thinking (plot_eff_ct)
  # ==============================================================================
  # terms = c("X軸", "グループ分け[平均±1SD]")
  pred_eff_ct <- ggpredict(m_ben_wave, terms = c("Eff_c", "CT_c [meansd]"))
  
  plot_eff_ct <- ggplot(pred_eff_ct, aes(x = x, y = predicted, color = group)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), 
                alpha = 0.15, linetype = 0) +
    geom_line(linewidth = 1) +
    scale_color_viridis_d(name = "Critical Thinking", labels = c("Low (-1SD)", "Mean", "High (+1SD)")) +
    scale_fill_viridis_d(name = "Critical Thinking", labels = c("Low (-1SD)", "Mean", "High (+1SD)")) +
    labs(
      title = "Interaction: Effectiveness × Critical Thinking",
      # subtitle を削除しました
      x = "Effectiveness (Centered)",
      y = "Predicted Benefit"
    ) +
    theme_classic()
  
  print(plot_eff_ct)
  
  # ==============================================================================
  # 3. 係数プロット (Forest Plot) (plot_forest)
  # ==============================================================================
  # モデルから係数と信頼区間を抽出
  tidy_ben <- tidy(m_ben_wave, conf.int = TRUE) %>%
    filter(effect == "fixed", term != "(Intercept)") # 切片は除外
  
  plot_forest <- ggplot(tidy_ben, aes(x = estimate, y = reorder(term, estimate))) +
    # 基準線（0）
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    # エラーバー（信頼区間）
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    # 点（推定値）
    geom_point(size = 3, color = "darkblue") +
    labs(
      title = "Fixed Effects on Benefit (LMM)",
      # subtitle を削除しました
      x = "Estimate (Coefficient)",
      y = "Predictors"
    ) +
    theme_bw()
  
  print(plot_forest)
  
  # ==============================================================================
  # 4. 順序ロジスティック回帰の予測確率プロット (plot_polr)
  # ==============================================================================
  # 順序ロジスティックモデル(m_vacc3_cov)を使用
  # Benefit_t3_c の変化に応じた、各カテゴリ(0, 2, 3plus)の確率を計算
  pred_polr <- ggpredict(m_vacc3_cov, terms = "Benefit_t3_c [all]")
  
  plot_polr <- ggplot(pred_polr, aes(x = x, y = predicted, color = response.level)) +
    geom_line(linewidth = 1.2) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = response.level), 
                alpha = 0.1, linetype = 0) +
    labs(
      title = "Predicted Probabilities of Vaccination Count",
      # subtitle を削除しました
      x = "Benefit (t3, Centered)",
      y = "Probability",
      color = "Vaccination Count",
      fill = "Vaccination Count"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(plot_polr)