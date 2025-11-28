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
