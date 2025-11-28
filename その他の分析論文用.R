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

#マンホイットニーのU検定
wilcox.test(averageNonNumerical ~ berlin_group,
            data = merged_data_complete_for_path,
            subset = berlin_group %in% c("BH", "BL"))
wilcoxonR(x = group_H, y = group_L)
