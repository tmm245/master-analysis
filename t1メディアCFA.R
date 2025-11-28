############################################################
# 合成得点（平均）で一気通貫：t1 -> t3/t4 LMM 分析パイプライン
# 前提: merged_data_complete, df_t1_cred, df_t1_use が既に環境に存在
# 依存: tidyverse, lme4, lmerTest, broom.mixed, performance, emmeans
# 出力: 要約テーブル・単純傾きCSV・ログメモ
############################################################

## 0) パッケージ
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(lme4); library(lmerTest)
  library(emmeans)
  library(broom.mixed); library(performance)
})

## 1) ユーティリティ
id_to_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- NA_character_
  x
}
safe_num <- function(x) suppressWarnings(as.numeric(x))
quick_desc <- function(...) {
  d <- list(...)
  as.data.frame(lapply(d, function(v){
    c(Min.=min(v,na.rm=TRUE), `1Q`=quantile(v,.25,na.rm=TRUE),
      Median=median(v,na.rm=TRUE), Mean=mean(v,na.rm=TRUE),
      `3Q`=quantile(v,.75,na.rm=TRUE), Max.=max(v,na.rm=TRUE))
  }))
}

## 2) t1 合成得点（平均）を作成：Cred（Q15_x）/ Use（Q58_x）
# ここは設問対応を固定：GOV=2項目, MASS=2項目, CIT=4項目, EXP=単一
scores_composite_cred <- df_t1_cred %>%
  transmute(
    id           = id_to_chr(id),
    GOV_cred_t1  = rowMeans(across(c(Q15_1, Q15_2)), na.rm = TRUE),
    MASS_cred_t1 = rowMeans(across(c(Q15_3, Q15_4)), na.rm = TRUE),
    CIT_cred_t1  = rowMeans(across(c(Q15_6, Q15_7, Q15_8, Q15_9)), na.rm = TRUE),
    EXP_cred_t1  = Q15_10
  )

scores_composite_use <- df_t1_use %>%
  transmute(
    id           = id_to_chr(id),
    GOV_use_t1   = rowMeans(across(c(Q58_1, Q58_2)), na.rm = TRUE),
    MASS_use_t1  = rowMeans(across(c(Q58_3, Q58_4)), na.rm = TRUE),
    CIT_use_t1   = rowMeans(across(c(Q58_5, Q58_6, Q58_7, Q58_8)), na.rm = TRUE),
    EXP_use_t1   = Q58_9
  )

## 3) 分析用 long データ（t3,t4のみ：ラグは t2→t3, t3→t4）
ct_candidates <- c("mean_cliticalThinking_attitude","mean_criticalThinking_attitude","CT")
ct_name <- ct_candidates[ct_candidates %in% names(merged_data_complete)][1]
if (is.na(ct_name)) stop("CT列が見つからない: ", paste(ct_candidates, collapse=", "))

df_core <- merged_data_complete %>%
  select(id, !!ct_name,
         Effectiveness_1, Effectiveness_2, Effectiveness_3, Effectiveness_4,
         Safeness_1,      Safeness_2,      Safeness_3,      Safeness_4,
         HBM_average_benefit_3, HBM_average_benefit_4,
         HBM_average_risk_3,    HBM_average_risk_4)

mk_wave <- function(df, wave){
  if (wave=="t3"){
    transmute(df,
              id    = id_to_chr(id),
              wave  = factor("t3", levels=c("t3","t4")),
              CT_t2 = .data[[ct_name]],
              Effectiveness_lag = Effectiveness_2,
              Safeness_lag      = Safeness_2,
              Benefit = HBM_average_benefit_3,
              Risk    = HBM_average_risk_3
    )
  } else {
    transmute(df,
              id    = id_to_chr(id),
              wave  = factor("t4", levels=c("t3","t4")),
              CT_t2 = .data[[ct_name]],
              Effectiveness_lag = Effectiveness_3,
              Safeness_lag      = Safeness_3,
              Benefit = HBM_average_benefit_4,
              Risk    = HBM_average_risk_4
    )
  }
}

long <- bind_rows(mk_wave(df_core,"t3"), mk_wave(df_core,"t4"))

## t1 合成得点の結合（id は character）
stopifnot(is.character(scores_composite_cred$id),
          is.character(scores_composite_use$id))

long <- long %>%
  left_join(scores_composite_cred, by="id", multiple="first") %>%
  left_join(scores_composite_use,  by="id", multiple="first")

## 4) センタリング
# within-person centering for lags（推奨：個人差コントロール）
long <- long %>%
  group_by(id) %>%
  mutate(
    Effectiveness_lag_w = Effectiveness_lag - mean(Effectiveness_lag, na.rm=TRUE),
    Safeness_lag_w      = Safeness_lag      - mean(Safeness_lag, na.rm=TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    # grand-mean centering for moderators (t1 composites)
    CT_c           = scale(CT_t2,         center=TRUE, scale=FALSE)[,1],
    GOV_cred_t1_c  = scale(GOV_cred_t1,   center=TRUE, scale=FALSE)[,1],
    MASS_cred_t1_c = scale(MASS_cred_t1,  center=TRUE, scale=FALSE)[,1],
    CIT_cred_t1_c  = scale(CIT_cred_t1,   center=TRUE, scale=FALSE)[,1],
    GOV_use_t1_c   = scale(GOV_use_t1,    center=TRUE, scale=FALSE)[,1],
    MASS_use_t1_c  = scale(MASS_use_t1,   center=TRUE, scale=FALSE)[,1],
    CIT_use_t1_c   = scale(CIT_use_t1,    center=TRUE, scale=FALSE)[,1]
  )

## 5) LMM（合成得点をモデレーターとして使用）
# 5-1) コア主効果（within中心化版）
m_ben_core <- lmer(Benefit ~ Effectiveness_lag_w + Safeness_lag_w + wave + (1|id), data=long)
m_risk_core <- lmer(Risk    ~ Effectiveness_lag_w + Safeness_lag_w + wave + (1|id), data=long)

# 5-2) 既報の Eff×CT（確認）
m_ben_ct  <- lmer(Benefit ~ Effectiveness_lag_w*CT_c + Safeness_lag_w + wave + (1|id), data=long)
m_risk_ct <- lmer(Risk    ~ Effectiveness_lag_w*CT_c + Safeness_lag_w + wave + (1|id), data=long)

# 5-3) 情報信頼（Cred composites）による Eff→Benefit/Risk の調整
m_ben_govCred_comp  <- lmer(Benefit ~ Effectiveness_lag_w*GOV_cred_t1_c  + Safeness_lag_w + wave + (1|id), data=long)
m_ben_massCred_comp <- lmer(Benefit ~ Effectiveness_lag_w*MASS_cred_t1_c + Safeness_lag_w + wave + (1|id), data=long)
m_ben_citCred_comp  <- lmer(Benefit ~ Effectiveness_lag_w*CIT_cred_t1_c  + Safeness_lag_w + wave + (1|id), data=long)

m_risk_govCred_comp  <- lmer(Risk ~ Effectiveness_lag_w*GOV_cred_t1_c  + Safeness_lag_w + wave + (1|id), data=long)
m_risk_massCred_comp <- lmer(Risk ~ Effectiveness_lag_w*MASS_cred_t1_c + Safeness_lag_w + wave + (1|id), data=long)
m_risk_citCred_comp  <- lmer(Risk ~ Effectiveness_lag_w*CIT_cred_t1_c  + Safeness_lag_w + wave + (1|id), data=long)

# 5-4) 情報利用（Use composites）による調整（任意）
m_ben_govUse_comp   <- lmer(Benefit ~ Effectiveness_lag_w*GOV_use_t1_c   + Safeness_lag_w + wave + (1|id), data=long)
m_ben_massUse_comp  <- lmer(Benefit ~ Effectiveness_lag_w*MASS_use_t1_c  + Safeness_lag_w + wave + (1|id), data=long)
m_ben_citUse_comp   <- lmer(Benefit ~ Effectiveness_lag_w*CIT_use_t1_c   + Safeness_lag_w + wave + (1|id), data=long)

m_risk_govUse_comp  <- lmer(Risk ~ Effectiveness_lag_w*GOV_use_t1_c  + Safeness_lag_w + wave + (1|id), data=long)
m_risk_massUse_comp <- lmer(Risk ~ Effectiveness_lag_w*MASS_use_t1_c + Safeness_lag_w + wave + (1|id), data=long)
m_risk_citUse_comp  <- lmer(Risk ~ Effectiveness_lag_w*CIT_use_t1_c  + Safeness_lag_w + wave + (1|id), data=long)

## 6) 要約テーブル（固定効果 + 95%CI + R2）
summ_fixed <- function(mod, name){
  out <- tidy(mod, effects="fixed", conf.int=TRUE)
  rs  <- r2(mod)
  attr(out, "R2") <- rs
  attr(out, "name") <- name
  out
}
mods <- list(
  ben_core = m_ben_core, risk_core = m_risk_core,
  ben_ct = m_ben_ct, risk_ct = m_risk_ct,
  ben_govCred = m_ben_govCred_comp, ben_massCred = m_ben_massCred_comp, ben_citCred = m_ben_citCred_comp,
  risk_govCred = m_risk_govCred_comp, risk_massCred = m_risk_massCred_comp, risk_citCred = m_risk_citCred_comp,
  ben_govUse = m_ben_govUse_comp, ben_massUse = m_ben_massUse_comp, ben_citUse = m_ben_citUse_comp,
  risk_govUse = m_risk_govUse_comp, risk_massUse = m_risk_massUse_comp, risk_citUse = m_risk_citUse_comp
)
tbls <- lapply(names(mods), function(n) summ_fixed(mods[[n]], n))
names(tbls) <- names(mods)

# 代表表示（コンソール）
message("\n== Core models (Benefit/Risk, within-centered) ==")
print(tidy(m_ben_core, effects = "fixed", conf.int = TRUE))
print(r2(m_ben_core))
print(tidy(m_risk_core, effects = "fixed", conf.int = TRUE))
print(r2(m_risk_core))

## 7) 単純傾き（例：GOV_cred_t1_c の -1SD / 0 / +1SD で Eff→Benefit）
em_ben_govCred <- emtrends(
  m_ben_govCred_comp, ~ GOV_cred_t1_c, var = "Effectiveness_lag_w",
  at = list(GOV_cred_t1_c = c(-sd(long$GOV_cred_t1_c,na.rm=TRUE), 0, sd(long$GOV_cred_t1_c,na.rm=TRUE)))
)
simp_gov <- as.data.frame(em_ben_govCred)
write.csv(simp_gov, "simple_slopes_eff_by_govcred_composite.csv", row.names = FALSE, fileEncoding = "CP932")
message("Saved: simple_slopes_eff_by_govcred_composite.csv")

## 8) スケール向き・相関の簡易チェック（t2→t3 の素朴相関）
message("\n== A) 尺度の向きチェック（要約 & 相関：t2→t3） ==")
print(quick_desc(
  Effectiveness_2 = merged_data_complete$Effectiveness_2,
  Safeness_2      = merged_data_complete$Safeness_2,
  Benefit_3       = merged_data_complete$HBM_average_benefit_3,
  Risk_3          = merged_data_complete$HBM_average_risk_3
))
message("Pearson相関（t2→t3）: Effectiveness_2 / Safeness_2 / Benefit_3 / Risk_3")
print(with(merged_data_complete,
           cor(cbind(Effectiveness_2, Safeness_2, HBM_average_benefit_3, HBM_average_risk_3),
               use = "pairwise.complete.obs"))
)

## 9) ログ（採用方針メモ）
cat("
[採用方針（合成得点ベース）]
- t1の情報信頼（GOV/MASS/CIT）は各項目平均で合成。EXPは単一項目のまま。
- LMMの主分析は within-person centered lag（Eff/Safe）＋ grand-mean centered moderators。
- 既報の Eff×CT は確認的に推定。Cred/Use（合成）×Eff の交互作用は本文で、必要に応じて補助解析に回す。

[補足]
- 2項目下位尺度（GOV/MASS）は平均合成が実務的に妥当。2項目相関 r と 2-item α（Spearman–Brown）を結果に併記推奨。
- 重要な交互作用は頑健性確保のため、因子得点版（別途）でも再確認可（付録）。
")

## 10) 参考：2項目相関 & 2-item α を計算・表示（報告用）
sb_alpha <- function(r) 2*r/(1+r)
pair_stats <- function(df, a, b){
  r <- suppressWarnings(cor(df[[a]], df[[b]], use="pairwise.complete.obs"))
  c(pair = paste(a,b,sep="-"), r = r, alpha_2item = sb_alpha(r))
}
message("\n\n== 2項目相関 & 2-item α（報告用） ==")
print(as.data.frame(rbind(
  pair_stats(df_t1_cred, "Q15_1","Q15_2"),
  pair_stats(df_t1_cred, "Q15_3","Q15_4"),
  pair_stats(df_t1_use , "Q58_1","Q58_2"),
  pair_stats(df_t1_use , "Q58_3","Q58_4")
)))
