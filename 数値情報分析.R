#数値情報群別平均値
merged_data_complete_for_path %>%
  group_by(berlin_group) %>%
  summarise(mean_averageNonNumerical = mean(averageNonNumerical, na.rm = TRUE))

#数値情報群別平均値
merged_data_complete_for_path %>%
  group_by(berlin_group) %>%
  summarise(mean_averageNumerical = mean(averageNumerical, na.rm = TRUE))

#--------------------------------#差の検定

#正規性の検定
# グループBHのデータ
group_H <- merged_data_complete_for_path %>%
  filter(berlin_group == "BH") %>%
  pull(averageNonNumerical)

# グループBLのデータ
group_L <- merged_data_complete_for_path %>%
  filter(berlin_group == "BL") %>%
  pull(averageNonNumerical)

# Shapiro-Wilk検定
shapiro.test(group_H)
shapiro.test(group_L)

wilcox.test(averageNonNumerical ~ berlin_group,
            data = merged_data_complete_for_path,
            subset = berlin_group %in% c("BH", "BL"))

# NA除去済みのベクトルを結合
x <- c(group_H, group_L)
g <- c(rep("H", length(group_H)), rep("L", length(group_L)))

# 効果量の計算
wilcoxonR(x = x, g = g)


#--------------------------------#プロット
#averageNonNumerical 
# データフレームの抽出（BHとBLのみ）
plot_data <- merged_data_complete_for_path %>%
  filter(berlin_group %in% c("BH", "BL"))

# グループのラベルを因子にして順序指定（任意）
plot_data$berlin_group <- factor(plot_data$berlin_group, levels = c("BL", "BH"))

# ggplotでバイオリン＋箱ひげ＋平均点をプロット
ggplot(plot_data, aes(x = berlin_group, y = averageNonNumerical, fill = berlin_group)) +
  geom_violin(trim = FALSE, alpha = 0.3) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red") +
  labs(
    title = "Distribution of averageNonNumerical by Berlin Group",
    x = "Berlin Group",
    y = "averageNonNumerical"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

#averageNumerical 
# パッケージ読み込み（まだであれば）
library(ggplot2)
library(dplyr)

# データ抽出（BLとBHのみ）
plot_data_num <- merged_data_complete_for_path %>%
  filter(berlin_group %in% c("BH", "BL"))

# グループラベルの順序指定（任意）
plot_data_num$berlin_group <- factor(plot_data_num$berlin_group, levels = c("BL", "BH"))

# プロット作成
ggplot(plot_data_num, aes(x = berlin_group, y = averageNumerical, fill = berlin_group)) +
  geom_violin(trim = FALSE, alpha = 0.3) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red") +
  labs(
    title = "Distribution of averageNumerical by Berlin Group",
    x = "Berlin Group",
    y = "averageNumerical"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
