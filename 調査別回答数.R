library(dplyr)

answered_summary <- merged_data_complete %>%
  mutate(
    answered_1 = !is.na(Q34),
    answered_2 = !is.na(QID1),
    answered_3 = !is.na(num_vaccination_3),
    answered_4 = !is.na(num_vaccination_4),
    answered_total = rowSums(across(starts_with("answered_")))
  ) %>%
  count(answered_total) %>%
  arrange(desc(answered_total))

print(answered_summary)
