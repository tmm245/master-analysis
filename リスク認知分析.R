#--------------------------------#データの準備と正規性の検定
# データ抽出
group_H_risk <- merged_data_complete_for_path %>%
  filter(berlin_group == "BH") %>%
  pull(HBM_average_risk_3)

group_L_risk <- merged_data_complete_for_path %>%
  filter(berlin_group == "BL") %>%
  pull(HBM_average_risk_3)

# Shapiro-Wilk検定
shapiro.test(group_H_risk)
shapiro.test(group_L_risk)

#--------------------------------#群間差の検定（Wilcoxon）
wilcox.test(HBM_average_risk_3 ~ berlin_group,
            data = merged_data_complete_for_path,
            subset = berlin_group %in% c("BH", "BL"))
# NA除去
group_H_risk_clean <- na.omit(group_H_risk)
group_L_risk_clean <- na.omit(group_L_risk)

# ベクトルとグループラベルの作成
x_risk <- c(group_H_risk_clean, group_L_risk_clean)
g_risk <- c(rep("H", length(group_H_risk_clean)), rep("L", length(group_L_risk_clean)))

# 効果量計算
wilcoxonR(x = x_risk, g = g_risk)

#--------------------------------#グループごとの平均値の算出
merged_data_complete_for_path %>%
  filter(berlin_group %in% c("BH", "BL")) %>%
  group_by(berlin_group) %>%
  summarise(mean_risk = mean(HBM_average_risk_3, na.rm = TRUE))
