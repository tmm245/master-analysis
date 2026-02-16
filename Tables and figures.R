# 作図用パッケージ
if (!require("ggeffects")) install.packages("ggeffects")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("broom.mixed")) install.packages("broom.mixed") # 係数抽出用

library(ggeffects)
library(ggplot2)
library(broom.mixed)
library(dplyr)


## 必要パッケージ
library(dplyr)
library(tidyr)
library(MASS)     # polr
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
# 4. 作図 (APAスタイル：タイトル中央・凡例右上配置版)
# ==============================================================================

# フォント設定（必要に応じて変更）
apa_base_size <- 12
apa_font_family <- "sans" 

ggplot() +
  # --- メインの折れ線グラフ ---
  geom_errorbar(data = summary_data, 
                aes(x = TimePoint, ymin = Mean - SE, ymax = Mean + SE, 
                    group = Group, color = Group), 
                width = 0.1, size = 0.7) +
  
  geom_line(data = summary_data, 
            aes(x = TimePoint, y = Mean, group = Group, color = Group, linetype = Group),
            size = 1) +
  
  geom_point(data = summary_data, 
             aes(x = TimePoint, y = Mean, group = Group, color = Group, shape = Group),
             size = 3) +
  
  # --- 有意差ブラケット ---
  geom_segment(data = annotation_df,
               aes(x = x_start, xend = x_end, y = y_pos, yend = y_pos, color = Group),
               size = 0.4, show.legend = FALSE) +
  geom_segment(data = annotation_df,
               aes(x = x_start, xend = x_start, y = y_pos, yend = y_pos - bracket_h, color = Group),
               size = 0.4, show.legend = FALSE) +
  geom_segment(data = annotation_df,
               aes(x = x_end, xend = x_end, y = y_pos, yend = y_pos - bracket_h, color = Group),
               size = 0.4, show.legend = FALSE) +
  
  # --- アスタリスク ---
  geom_text(data = annotation_df,
            aes(x = (x_start + x_end) / 2, y = y_pos + 0.1, label = Label, color = Group),
            size = 5, show.legend = FALSE, vjust = 0, family = apa_font_family) +
  
  # --- スケール設定 ---
  scale_color_manual(values = c("Unvaccinated (0)" = "black", "Vaccinated (1+)" = "gray50")) + 
  scale_linetype_manual(values = c("Unvaccinated (0)" = "solid", "Vaccinated (1+)" = "dashed")) +
  scale_shape_manual(values = c("Unvaccinated (0)" = 16, "Vaccinated (1+)" = 17)) +
  
  # Y軸の範囲（凡例が入るスペースを確保するため、上限を少し余裕を持たせるのがコツです）
  scale_y_continuous(limits = c(1, 9.0), breaks = 1:7, expand = c(0, 0)) +
  
  facet_wrap(~ Measure) +
  
  labs(
    y = "Mean Score",
    x = "Time Point",
    color = NULL, linetype = NULL, shape = NULL # 凡例タイトルを削除
  ) +
  
  # --- テーマ設定 ---
  theme_classic(base_size = apa_base_size, base_family = apa_font_family) +
  # --- テーマ設定 (凡例に枠線を追加) ---
  theme_classic(base_size = apa_base_size, base_family = apa_font_family) +
  theme(
    # 軸設定（変更なし）
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    
    # タイトルを中央揃え（t2とt3の間）
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold", hjust = 0.5), 
    
    # 凡例の位置
    legend.position = c(0.975, 0.95),  # 右上（必要に応じて微調整）
    legend.justification = c(1, 1),   # 右上を基準点に
    
    # 【ここが追加点】凡例を枠で囲む設定
    legend.background = element_rect(
      fill = "white",     # 背景を白（グラフ線を隠すため）
      color = "black",    # 枠線の色
      size = 0.5          # 枠線の太さ
    ),
    
    # 枠線と文字の間に余白を入れる（これがないと窮屈になります）
    legend.margin = margin(t = 5, r = 5, b = 5, l = 5),
    
    # 線種がわかるようにキーの幅を確保
    legend.key.width = unit(1.5, "cm")
  )

# ==============================================================================
# 1. 交互作用プロット: Safeness × Wave (APA Style)
# ==============================================================================
library(ggeffects)
library(ggplot2)

# フォント設定
apa_base_size <- 12
apa_font_family <- "sans" # または "serif"

# 予測値の計算 (既にある場合はスキップ可)
# terms = c("Safe_c", "wave") -> x軸がSafe_c, グループ(線)がwave
pred_safe_wave <- ggpredict(m_ben_wave, terms = c("Safe_c", "wave"))

plot_safe_wave <- ggplot(pred_safe_wave, aes(x = x, y = predicted, group = group)) +
  
  # --- 信頼区間 (Ribbon) ---
  # 白黒印刷でも邪魔にならないよう、薄いグレーまたは透明度を高めた色を使用
  # show.legend = FALSE にして、凡例には線だけを表示させるのがスッキリさせるコツです
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), 
              alpha = 0.1, linetype = 0, show.legend = FALSE) +
  
  # --- 回帰直線 ---
  # linetype = group を指定して、実線・点線などを使い分ける
  geom_line(aes(color = group, linetype = group), size = 1) +
  
  # --- スケール設定 (APA推奨: 白黒対応) ---
  # 4時点あると仮定して、線種を割り当てます
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
  
  # 色はグレースケール(黒〜濃いグレー)にするか、視認性の良い色を指定
  # ここでは前の図に合わせて黒とグレーの濃淡にします
  scale_color_grey(start = 0, end = 0.6) + 
  scale_fill_grey(start = 0, end = 0.6) +
  
  # --- ラベル ---
  labs(
    # APAでは図の中にタイトルを書かず、キャプションに書くのが一般的ですが
    # 必要であればここに title = "..." を追加してください
    x = "Perceived Safeness (Centered)",
    y = "Predicted Perceived Benefit",
    color = NULL, linetype = NULL # 凡例タイトルを削除
  ) +
  
  # --- テーマ設定 ---
  theme_classic(base_size = apa_base_size, base_family = apa_font_family) +
  theme(
    # 軸設定（変更なし）
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    
    # 凡例の位置設定
    legend.position = c(0.2, 0.95), 
    legend.justification = c(1, 1),
    
    # 【ここが変更点】凡例に枠線と背景色をつける
    legend.background = element_rect(
      fill = "white",    # 背景を白にしてグラフ線と重ならないようにする
      color = "black",   # 枠線の色（黒）
      size = 0.5         # 枠線の太さ
    ),
    
    # 【推奨】枠線と文字の間に少し余白を入れる
    legend.margin = margin(t = 4, r = 4, b = 4, l = 4),
    
    # テキストサイズ（必要に応じて）
    legend.text = element_text(size = 10),
    legend.key.width = unit(1.5, "cm")
  )

print(plot_safe_wave)
