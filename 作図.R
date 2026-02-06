## 必要パッケージ
library(dplyr)
library(tidyr)
library(MASS)     # polr
library(ggplot2)


## 1. A系カテゴリ（0,2,3,4+）のデータフレーム -------------------------

ord_BR_t3t4_A <- merged_data_complete %>%
  transmute(
    id,
    Benefit_t3 = HBM_average_benefit_3,
    Risk_t3    = HBM_average_risk_3,
    vacc4      = num_vaccination_4
  ) %>%
  mutate(
    # A系カテゴリ：0, 2, 3, 4plus（4と5をまとめる）
    vacc4_catA = case_when(
      vacc4 == 0        ~ "0",
      vacc4 == 2        ~ "2",
      vacc4 == 3        ~ "3",
      vacc4 >= 4        ~ "4plus",
      TRUE              ~ NA_character_
    ),
    vacc4_catA = factor(
      vacc4_catA,
      levels  = c("0", "2", "3", "4plus"),
      ordered = TRUE
    ),
    # 予測子は grand-mean センタリング
    Benefit_t3_c = as.numeric(scale(Benefit_t3, center = TRUE, scale = FALSE)),
    Risk_t3_c    = as.numeric(scale(Risk_t3,    center = TRUE, scale = FALSE))
  ) %>%
  tidyr::drop_na(vacc4_catA, Benefit_t3_c, Risk_t3_c)

# 確認
table(ord_BR_t3t4_A$vacc4_catA, useNA = "ifany")

## 2. A系カテゴリの順序ロジスティックモデル ------------------------

m_ord_BR_t3t4_A <- polr(
  formula = vacc4_catA ~ Benefit_t3_c + Risk_t3_c,
  data    = ord_BR_t3t4_A,
  Hess    = TRUE
)

summary(m_ord_BR_t3t4_A)

## 3. 予測用グリッドを作成 ----------------------------------------

ben_range  <- range(ord_BR_t3t4_A$Benefit_t3_c, na.rm = TRUE)
risk_range <- range(ord_BR_t3t4_A$Risk_t3_c,    na.rm = TRUE)

grid_pred <- expand.grid(
  Benefit_t3_c = seq(ben_range[1],  ben_range[2],  length.out = 100),
  Risk_t3_c    = seq(risk_range[1], risk_range[2], length.out = 100)
)

## 4. 各グリッド点でカテゴリ別予測確率 ---------------------------

prob_mat <- predict(
  m_ord_BR_t3t4_A,
  newdata = grid_pred,
  type    = "probs"
)

prob_df <- as.data.frame(prob_mat)

# 各点で「最も確率の高いカテゴリ」をラベル
max_cat <- colnames(prob_df)[max.col(prob_df, ties.method = "first")]

grid_plot <- grid_pred %>%
  bind_cols(prob_df) %>%
  mutate(
    pred_cat = factor(
      max_cat,
      levels = levels(ord_BR_t3t4_A$vacc4_catA)  # c("0","2","3","4plus")
    )
  )

## 5. 決定境界図を描画 --------------------------------------------

ggplot() +
  # 背景：その点で最も確率の高いカテゴリ
  geom_raster(
    data = grid_plot,
    aes(x = Benefit_t3_c, y = Risk_t3_c, fill = pred_cat),
    alpha = 0.7
  ) +
  # 実データ点を重ねる
  geom_point(
    data = ord_BR_t3t4_A,
    aes(x = Benefit_t3_c, y = Risk_t3_c, color = vacc4_catA),
    size  = 2,
    alpha = 0.8
  ) +
  labs(
    x = "Benefit_t3（中心化）",
    y = "Risk_t3（中心化）",
    fill   = "Predicted\ncategory",
    color  = "Observed\ncategory",
    title  = "A系カテゴリ（0,2,3,4+）順序ロジスティック\nBenefit × Risk による決定境界"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

library(dplyr)
library(tidyr)
library(ggplot2)

## 1. Risk を平均に固定 ----------------------------------------
risk_mean <- mean(ord_BR_t3t4_A$Risk_t3_c, na.rm = TRUE)

## 2. Benefit の範囲で予測用データを作成 -----------------------
ben_range <- range(ord_BR_t3t4_A$Benefit_t3_c, na.rm = TRUE)

newdat_ben <- data.frame(
  Benefit_t3_c = seq(ben_range[1], ben_range[2], length.out = 200),
  Risk_t3_c    = risk_mean        # Risk は全点で平均に固定
)

## 3. 各 Benefit でのカテゴリ別予測確率を計算 ------------------
prob_mat_ben <- predict(
  m_ord_BR_t3t4_A,
  newdata = newdat_ben,
  type    = "probs"
)

prob_df_ben <- as.data.frame(prob_mat_ben)

plot_df_ben <- newdat_ben %>%
  bind_cols(prob_df_ben) %>%
  pivot_longer(
    cols      = colnames(prob_df_ben),   # "0","2","3","4plus"
    names_to  = "vacc_cat",
    values_to = "prob"
  ) %>%
  mutate(
    vacc_cat = factor(
      vacc_cat,
      levels = levels(ord_BR_t3t4_A$vacc4_catA)  # c("0","2","3","4plus")
    )
  )

## 4. プロット ---------------------------------------------------
ggplot(plot_df_ben,
       aes(x = Benefit_t3_c, y = prob, color = vacc_cat)) +
  geom_line(size = 1.1) +
  labs(
    x = "Benefit_t3（中心化）",
    y = "カテゴリ別予測確率",
    color = "接種カテゴリ",
    title = "A系カテゴリ（0,2,3,4+）\nRiskを平均に固定したときの Benefit 別予測確率"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
#--------------------------------------------------------------------------------#
library(dplyr)
library(tidyr)
library(ggplot2)

## 1. Risk を平均に固定 ----------------------------------------
risk_mean <- mean(ord_BR_t3t4_A$Risk_t3_c, na.rm = TRUE)

## 2. Benefit の範囲で予測用データを作成 -----------------------
ben_range <- range(ord_BR_t3t4_A$Benefit_t3_c, na.rm = TRUE)

newdat_ben <- data.frame(
  Benefit_t3_c = seq(ben_range[1], ben_range[2], length.out = 200),
  Risk_t3_c    = risk_mean        # Risk は全点で平均に固定
)

## 3. 各 Benefit でのカテゴリ別予測確率を計算 ------------------
prob_mat_ben <- predict(
  m_ord_BR_t3t4_A,
  newdata = newdat_ben,
  type    = "probs"
)

prob_df_ben <- as.data.frame(prob_mat_ben)

plot_df_ben <- newdat_ben %>%
  bind_cols(prob_df_ben) %>%
  pivot_longer(
    cols      = colnames(prob_df_ben),   # "0","2","3","4plus"
    names_to  = "vacc_cat",
    values_to = "prob"
  ) %>%
  mutate(
    vacc_cat = factor(
      vacc_cat,
      levels = levels(ord_BR_t3t4_A$vacc4_catA)  # c("0","2","3","4plus")
    )
  )

## 4. プロット ---------------------------------------------------
ggplot(plot_df_ben,
       aes(x = Benefit_t3_c, y = prob, color = vacc_cat)) +
  geom_line(size = 1.1) +
  labs(
    x = "Benefit_t3（中心化）",
    y = "カテゴリ別予測確率",
    color = "接種カテゴリ",
    title = "A系カテゴリ（0,2,3,4+）\nRiskを平均に固定したときの Benefit 別予測確率"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
#-------------------------------------------------------------------------------#
library(dplyr)
library(tidyr)
library(ggplot2)

## 1. Benefit を平均に固定 -------------------------------------
ben_mean <- mean(ord_BR_t3t4_A$Benefit_t3_c, na.rm = TRUE)

## 2. Risk の範囲で予測用データを作成 -------------------------
risk_range <- range(ord_BR_t3t4_A$Risk_t3_c, na.rm = TRUE)

newdat_risk <- data.frame(
  Benefit_t3_c = ben_mean,                                   # Benefit は平均に固定
  Risk_t3_c    = seq(risk_range[1], risk_range[2], length.out = 200)
)

## 3. 各 Risk でのカテゴリ別予測確率を計算 --------------------
prob_mat_risk <- predict(
  m_ord_BR_t3t4_A,
  newdata = newdat_risk,
  type    = "probs"
)

prob_df_risk <- as.data.frame(prob_mat_risk)

plot_df_risk <- newdat_risk %>%
  bind_cols(prob_df_risk) %>%
  pivot_longer(
    cols      = colnames(prob_df_risk),   # "0","2","3","4plus"
    names_to  = "vacc_cat",
    values_to = "prob"
  ) %>%
  mutate(
    vacc_cat = factor(
      vacc_cat,
      levels = levels(ord_BR_t3t4_A$vacc4_catA)  # c("0","2","3","4plus")
    )
  )

## 4. プロット --------------------------------------------------
ggplot(plot_df_risk,
       aes(x = Risk_t3_c, y = prob, color = vacc_cat)) +
  geom_line(size = 1.1) +
  labs(
    x = "Risk_t3（中心化）",
    y = "カテゴリ別予測確率",
    color = "接種カテゴリ",
    title = "A系カテゴリ（0,2,3,4+）\nBenefitを平均に固定したときの Risk 別予測確率"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
#-------------------------------------------------------------------------------#
library(dplyr)
library(tidyr)
library(ggplot2)

## 1. Benefit を平均に固定 -------------------------------------
ben_mean <- mean(ord_BR_t3t4_A$Benefit_t3_c, na.rm = TRUE)

## 2. Risk の範囲で予測用データを作成 -------------------------
risk_range <- range(ord_BR_t3t4_A$Risk_t3_c, na.rm = TRUE)

newdat_risk <- data.frame(
  Benefit_t3_c = ben_mean,                                   # Benefit は平均に固定
  Risk_t3_c    = seq(risk_range[1], risk_range[2], length.out = 200)
)

## 3. 各 Risk でのカテゴリ別予測確率を計算 --------------------
prob_mat_risk <- predict(
  m_ord_BR_t3t4_A,
  newdata = newdat_risk,
  type    = "probs"
)

prob_df_risk <- as.data.frame(prob_mat_risk)

plot_df_risk <- newdat_risk %>%
  bind_cols(prob_df_risk) %>%
  pivot_longer(
    cols      = colnames(prob_df_risk),   # "0","2","3","4plus"
    names_to  = "vacc_cat",
    values_to = "prob"
  ) %>%
  mutate(
    vacc_cat = factor(
      vacc_cat,
      levels = levels(ord_BR_t3t4_A$vacc4_catA)  # c("0","2","3","4plus")
    )
  )

## 4. プロット --------------------------------------------------
ggplot(plot_df_risk,
       aes(x = Risk_t3_c, y = prob, color = vacc_cat)) +
  geom_line(size = 1.1) +
  labs(
    x = "Risk_t3（中心化）",
    y = "カテゴリ別予測確率",
    color = "接種カテゴリ",
    title = "A系カテゴリ（0,2,3,4+）\nBenefitを平均に固定したときの Risk 別予測確率"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

#-------------------------------------------------------------------#

library(dplyr)
library(tidyr)
library(MASS)    # polr()
library(ggplot2)

# 3カテゴリデータ作成（0,2,3plus）
ord_BR_t3t4_3cat_simple <- merged_data_complete %>%
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

# 順序ロジスティック（主効果のみ）
m_ord_3cat_simple <- polr(
  vacc4_cat3 ~ Benefit_t3_c + Risk_t3_c,
  data = ord_BR_t3t4_3cat_simple,
  Hess = TRUE
)

summary(m_ord_3cat_simple)

## 1. 代表値の取得（Risk の 25, 50, 75%分位）
risk_q <- quantile(
  ord_BR_t3t4_3cat_simple$Risk_t3_c,
  probs = c(.25, .50, .75),
  na.rm = TRUE
)

risk_q
#   25%   50%   75%  くらいの値が出るはず

## 2. Benefit の範囲を決めてシーケンス生成
benefit_range <- range(ord_BR_t3t4_3cat_simple$Benefit_t3_c, na.rm = TRUE)
benefit_seq   <- seq(benefit_range[1], benefit_range[2], length.out = 100)

## 3. 予測用データ（Benefit を動かし、Risk は3水準に固定）
new_ben <- expand.grid(
  Benefit_t3_c = benefit_seq,
  Risk_t3_c    = as.numeric(risk_q)
)

new_ben$Risk_level <- factor(
  rep(names(risk_q), each = length(benefit_seq)),
  levels = c("25%", "50%", "75%"),
  labels = c("Risk: low (25%)", "Risk: mid (50%)", "Risk: high (75%)")
)

## 4. 予測確率の計算
probs_ben <- predict(m_ord_3cat_simple, newdata = new_ben, type = "probs")

probs_ben_df <- cbind(new_ben, probs_ben) |>
  pivot_longer(
    cols      = c("0", "2", "3plus"),
    names_to  = "vacc_cat",
    values_to = "prob"
  )

probs_ben_df$vacc_cat <- factor(
  probs_ben_df$vacc_cat,
  levels = c("0", "2", "3plus"),
  labels = c("0 doses", "2 doses", "3+ doses")
)

## 5. 作図（Benefit を横軸）
p_benefit <- ggplot(
  probs_ben_df,
  aes(x = Benefit_t3_c, y = prob, color = vacc_cat)
) +
  geom_line(size = 1) +
  facet_wrap(~ Risk_level) +
  labs(
    x = "Benefit at t3 (centered)",
    y = "Predicted probability",
    color = "Vaccination\ncategory",
    title = "Predicted probabilities by Benefit (t3)\nRisk fixed at 25/50/75% quantiles",
    subtitle = "Outcome: vaccination category at t4 (0, 2, 3+)"
  ) +
  theme_bw() +
  theme(
    legend.position = "right"
  )

# 図を表示
print(p_benefit)

## 1. 代表値の取得（Benefit の 25, 50, 75%分位）
benefit_q <- quantile(
  ord_BR_t3t4_3cat_simple$Benefit_t3_c,
  probs = c(.25, .50, .75),
  na.rm = TRUE
)

benefit_q

## 2. Risk の範囲を決めてシーケンス生成
risk_range <- range(ord_BR_t3t4_3cat_simple$Risk_t3_c, na.rm = TRUE)
risk_seq   <- seq(risk_range[1], risk_range[2], length.out = 100)

## 3. 予測用データ（Risk を動かし、Benefit は3水準に固定）
new_risk <- expand.grid(
  Risk_t3_c    = risk_seq,
  Benefit_t3_c = as.numeric(benefit_q)
)

new_risk$Benefit_level <- factor(
  rep(names(benefit_q), each = length(risk_seq)),
  levels = c("25%", "50%", "75%"),
  labels = c("Benefit: low (25%)", "Benefit: mid (50%)", "Benefit: high (75%)")
)

## 4. 予測確率の計算
probs_risk <- predict(m_ord_3cat_simple, newdata = new_risk, type = "probs")

probs_risk_df <- cbind(new_risk, probs_risk) |>
  pivot_longer(
    cols      = c("0", "2", "3plus"),
    names_to  = "vacc_cat",
    values_to = "prob"
  )

probs_risk_df$vacc_cat <- factor(
  probs_risk_df$vacc_cat,
  levels = c("0", "2", "3plus"),
  labels = c("0 doses", "2 doses", "3+ doses")
)

## 5. 作図（Risk を横軸）
p_risk <- ggplot(
  probs_risk_df,
  aes(x = Risk_t3_c, y = prob, color = vacc_cat)
) +
  geom_line(size = 1) +
  facet_wrap(~ Benefit_level) +
  labs(
    x = "Risk at t3 (centered)",
    y = "Predicted probability",
    color = "Vaccination\ncategory",
    title = "Predicted probabilities by Risk (t3)\nBenefit fixed at 25/50/75% quantiles",
    subtitle = "Outcome: vaccination category at t4 (0, 2, 3+)"
  ) +
  theme_bw() +
  theme(
    legend.position = "right"
  )

# 図を表示
print(p_risk)


library(dplyr)
library(tidyr)
library(ggplot2)

# ==============================================================================
# 1. データ整形 & 基礎集計
# ==============================================================================

# グループ分け（t4時点の接種回数 0回 vs 1回以上）
plot_data_raw <- merged_data_complete %>%
  mutate(
    Group = case_when(
      num_vaccination_4 == 0 ~ "Unvaccinated (0)",
      num_vaccination_4 >= 1 ~ "Vaccinated (1+)",
      TRUE ~ NA_character_
    )
  ) %>%
  tidyr::drop_na(Group) 

# 検定用データ（dplyr::selectを明示して競合回避）
data_for_test <- plot_data_raw %>%
  dplyr::select(id, Group, Safeness_1:Safeness_4, Effectiveness_1:Effectiveness_4)

# プロット用データ (Long形式)
plot_data_long <- plot_data_raw %>%
  dplyr::transmute(
    id, Group,
    Safeness_t1 = Safeness_1, Safeness_t2 = Safeness_2, 
    Safeness_t3 = Safeness_3, Safeness_t4 = Safeness_4,
    Effectiveness_t1 = Effectiveness_1, Effectiveness_t2 = Effectiveness_2,
    Effectiveness_t3 = Effectiveness_3, Effectiveness_t4 = Effectiveness_4
  ) %>%
  pivot_longer(
    cols      = -c(id, Group),
    names_to  = c("Measure", "TimePoint"),
    names_sep = "_",
    values_to = "Score"
  ) %>%
  mutate(
    TimePoint = factor(TimePoint, levels = c("t1", "t2", "t3", "t4")),
    Group     = factor(Group, levels = c("Unvaccinated (0)", "Vaccinated (1+)"))
  )

# 平均と標準誤差
summary_data <- plot_data_long %>%
  group_by(Group, Measure, TimePoint) %>%
  summarise(
    Mean = mean(Score, na.rm = TRUE),
    SE   = sd(Score, na.rm = TRUE) / sqrt(sum(!is.na(Score))),
    .groups = "drop"
  )

# ==============================================================================
# 2. 検定の自動実行 (対応のあるt検定)
# ==============================================================================

pairs_to_test <- list(c("t1", "t2"), c("t2", "t3"), c("t3", "t4"))
measures <- c("Safeness", "Effectiveness")
groups   <- c("Unvaccinated (0)", "Vaccinated (1+)")

test_results <- data.frame()

for (m in measures) {
  for (g in groups) {
    sub_df <- data_for_test %>% filter(Group == g)
    
    for (pair in pairs_to_test) {
      t_prev <- pair[1]
      t_curr <- pair[2]
      
      idx_prev <- gsub("t", "", t_prev)
      idx_curr <- gsub("t", "", t_curr)
      col_prev <- paste0(m, "_", idx_prev)
      col_curr <- paste0(m, "_", idx_curr)
      
      # 対応のあるt検定
      res <- t.test(sub_df[[col_prev]], sub_df[[col_curr]], paired = TRUE)
      
      p_val <- res$p.value
      stars <- case_when(
        p_val < 0.001 ~ "***",
        p_val < 0.01  ~ "**",
        p_val < 0.05  ~ "*",
        TRUE          ~ "ns"
      )
      
      if (stars != "ns") { 
        test_results <- rbind(test_results, data.frame(
          Measure   = m,
          Group     = g,
          TimeStart = t_prev,
          TimeEnd   = t_curr,
          P_val     = p_val,
          Label     = stars
        ))
      }
    }
  }
}

# ==============================================================================
# 3. アノテーション（ブラケット）座標の設定
# ==============================================================================

time_map <- c("t1" = 1, "t2" = 2, "t3" = 3, "t4" = 4)

# 高さの設定
base_y_unvac <- 6.5   # Unvaccinated用の基準高さ
base_y_vac   <- 7.2   # Vaccinated用の基準高さ
bracket_h    <- 0.2   # ブラケットの足の長さ

annotation_df <- test_results %>%
  mutate(
    x_start = time_map[TimeStart],
    x_end   = time_map[TimeEnd],
    
    y_pos = case_when(
      Group == "Unvaccinated (0)" ~ base_y_unvac,
      Group == "Vaccinated (1+)"  ~ base_y_vac
    ),
    
    Group = factor(Group, levels = c("Unvaccinated (0)", "Vaccinated (1+)"))
  )

# ==============================================================================
# 4. 作図
# ==============================================================================

ggplot() +
  # --- メインの折れ線グラフ ---
  geom_errorbar(data = summary_data, 
                aes(x = TimePoint, ymin = Mean - SE, ymax = Mean + SE, group = Group, color = Group),
                width = 0.1, size = 0.8) +
  
  geom_line(data = summary_data, 
            aes(x = TimePoint, y = Mean, group = Group, color = Group),
            size = 1.2) +
  
  geom_point(data = summary_data, 
             aes(x = TimePoint, y = Mean, group = Group, color = Group),
             size = 4) +
  
  # --- 有意差ブラケット ---
  geom_segment(data = annotation_df,
               aes(x = x_start, xend = x_end, y = y_pos, yend = y_pos, color = Group),
               size = 0.5, show.legend = FALSE) +
  geom_segment(data = annotation_df,
               aes(x = x_start, xend = x_start, y = y_pos, yend = y_pos - bracket_h, color = Group),
               size = 0.5, show.legend = FALSE) +
  geom_segment(data = annotation_df,
               aes(x = x_end, xend = x_end, y = y_pos, yend = y_pos - bracket_h, color = Group),
               size = 0.5, show.legend = FALSE) +
  
  # --- アスタリスク ---
  geom_text(data = annotation_df,
            aes(x = (x_start + x_end) / 2, y = y_pos + 0.1, label = Label, color = Group),
            size = 6, show.legend = FALSE, vjust = 0) +
  
  # --- 設定 ---
  facet_wrap(~ Measure) +
  scale_color_manual(values = c("Unvaccinated (0)" = "#E41A1C", "Vaccinated (1+)" = "#377EB8")) +
  scale_y_continuous(limits = c(1, 8.5), breaks = 1:7) + 
  
  labs(
    y = "Mean Score (1-7)",
    x = "Time Point",
    color = "Group (at t4)",
    title = "Safeness & Effectiveness Trends"
    # subtitle を削除しました
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 14),
    
    # 凡例を右側に配置
    legend.position = "right",
    
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )