#
# ============================================================================= #
# ---- ALFA+ LONGITUDINAL ANALYSIS SCRIPT ----
# ============================================================================= #
# Purpose:
#   - Inclusion criteria, cleaning and recoding of participants and variables 
#     of the ALFA+ study to optimize the following steps
#   - Exploratory analysis of the cohort checking for biases
#   - Longitudinal A-IADL-Q stratified analysis with Linear Mixed-effect Models (LME)
#   - Optimization and validation of the -2.2 A-IADL-Q points threshold in the 
#     ALFA+ cohort for binary decliner / no-decliner stratification of participants
# 
# Description:
# This script does data ingestion of the ALFA+ longitudinal dataset (created
# with the previous script "ALFA+_dataset_create_program.R") and performs an
# inferential and longitudinal analysis of biomarkers to investigate wheather 
# A-IADL-Q functional trajectories differ among cognitively unimpaired individuals 
# depending on it's CSF Amyloid-beta and pTau status and if interindividual 
# variability can be explained by several CSF and plasma biomarkers.
#
# The Pipeline includes:
#   - Setup and data load
#   - Inclusion criteria, depuration and recoding of participant's variables
#   - Quality control checks
#   - Exploratory evaluation of the cohort with a demographic table
#   - Selection of optimal baseline LME model with ANOVA and AIC criteria.
#   - Creation and evaluation of Analysis 1: Time x AT_group
#   - Biomarker quality control, normalization and ratio computation
#   - Creation and evaluation of Analysis 2: Time x Biomarker (one model per Bmk)
#   - Creation and evaluation of Analysis 3: Time x Biomarker x AT_group
#   - Visit selection criteria to determine longitudinal MIC (V2-V1 vs. V3-V1)
#   - Correlation and in-depth comparison between decliners and no-decliner
#
# Inputs:
#   - ALFA long dataset (outuput of the "ALFA+_dataset_create_program.R")
#
# Outputs:
#   - Demographic table
#   - Base LME model performance metrics
#   - Plot of estimated A-IADL-Q longitudinal trajectories stratified by AT_group 
#   - Estimates and p-values of the predictors in Model 1
#   - Biomarker p-values (FDR corrected), and 95% CI association with
#     A-IADL-Q functional decline extracted from LME models
#   - Plot of functional trajectories based on biomarker concentration (tertiles)
#   - Forest plot of the 95% CI derived form the Time x Biomarker x AT_Group model
#   - Correlation table between V2-V1 and V3-V1 criteria
#   - Longitudinal trajectories of functional decline stratified by visit selection
#     criteria using a LOESS approach
#   - SCD and SUCOG correlation between decliners and no-decliners
#   - Comparison table between decliners and no-decliners:
#       - Demographics, genetics and frailty
#       - Baseline cognition and function
#       - Baseline biomarkers
#       - Longitudinal cognition
#     
#
# Reproducibility:
#   - Developed in R version 4.5.2
#   - Key Packages:
#     - dplyr (2.5.1)
#     - tidyr (1.3.1)
#     - readr (2.1.6)
#     - lme4 (1.1-38)
#     - lmerTest (3.1-3)
#     - performance (0.15.3)
#     - gtsummary (2.5.0)
#     - gt (1.1.0)
#     - ggplot2 (4.0.1) 
#
# Code: Eloi Jacomet & Federica Anastasi
# Date: May 5th 2026 
# ============================================================================= #


# ============================================================================= #
# 1. SETUP AND DATA LOADING ----
# ============================================================================= #
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(here)
  library(gtsummary)
  library(gt)
  library(flextable)
  library(lme4)
  library(lmerTest)
  library(sjPlot)
  library(performance)
  library(emmeans)
  library(purrr)
  library(ggplot2)
  library(forcats)
  library(broom.mixed)
  library(ggtext)
  library(tidyverse)
  library(scales)
  library(ggthemes)
  library(ggeffects)
})

rm(list = ls())
setwd(here())
wd <- getwd()


## ---- 1.1 Load alfa_long dataset ----
load("data/ALFAplus_long_IADL_EDU_ATN_BIOMK_APOE_SUCOG_PACC_FRAILTY_TMT_Diagnostic20260421_build_20260421.RData")

## ---- 1.2 Check class and format data ----
str(alfa_long)
colSums(is.na(alfa_long))
alfa_long <- alfa_long %>% mutate(visit = as.character(visit))

## ---- 1.3 Ensure date is codified as date ----
if (!inherits(alfa_long$date, "Date")) {
  alfa_long$date <- as.Date(alfa_long$date)
}

# ============================================================================= #
# 2. INCLUSION CRITERIA, CLEANING AND RECODING ----
# ============================================================================= #

## ---- 2.1 Exclude baseline A-T+ ----
alfa_long <- alfa_long %>%
  filter(AT_bl != "A-T+" | is.na(AT_bl)) %>%   
  mutate(AT_bl = droplevels(AT_bl))


## ---- 2.2 Inclusion sets (ID) ----

# Must have IADL at V1 and >=1 follow-up (V2/V3) IADL
iadl_ok <- alfa_long %>%
  group_by(alfa_id) %>%
  summarise(
    has_V1_iadl = any(visit == "V1" & !is.na(iadl_score)),
    has_fu_iadl = any(visit %in% c("V2", "V2_5", "V3", "V3_5") & !is.na(iadl_score)),  # follow-ups here
    .groups = "drop"
  ) %>%
  filter(has_V1_iadl, has_fu_iadl) %>%
  pull(alfa_id)

# Must have baseline AT defined 

at_ok <- alfa_long %>%
  filter(!is.na(AT_bl)) %>%
  pull(alfa_id) %>%
  unique()


## ---- 2.3 Final participants = intersection ----
final_participants <- Reduce(intersect, list(iadl_ok, at_ok))

## ---- 2.4 Filter dataset for final participants ----
alfa_long_final <- alfa_long %>%
  filter(alfa_id %in% final_participants)

## ---- 2.5 QC report ----
cat("\n--- Inclusion/Exclusion summary ---\n")
cat("N with V1 + ≥1 follow-up IADL:", length(iadl_ok), "\n")
cat("N with AT baseline defined (ATN_bl):", length(at_ok), "\n")
cat("N FINAL participants:", length(final_participants), "\n\n")


## ---- 2.6 Build final dataset, drop V2.5 if IADL absent and recode IADL score in range (70-0) ----
alfa_study <- alfa_long %>%
  filter(alfa_id %in% final_participants) %>%
  mutate(iadl_score = iadl_score/10)
  

for (v in c("V2_5", "V3_5")) {
  
  vv <- alfa_study %>%
    filter(visit == v) %>%
    summarise(
      n_rows = n(),
      n_iadl_nonmiss = sum(!is.na(iadl_score)),
      .groups = "drop"
    )
  
  if (nrow(vv) == 1 && vv$n_rows > 0) {
    cat(v, " rows:", vv$n_rows, " | non-missing IADL:", vv$n_iadl_nonmiss, "\n")
    if (vv$n_iadl_nonmiss == 0) {
      message("Dropping ", v, " (IADL fully missing).")
      alfa_study <- alfa_study %>% filter(visit != v)
    }
  }
}

table(alfa_study$visit, !is.na(alfa_study$iadl_score))

## ---- 2.7 Compute time since baseline (years) from V1 date ----

# Check duplicated V1 rows
baseline_ids <- alfa_study %>%
  filter(visit == "V1") %>%
  count(alfa_id) %>%
  filter(n > 1) %>%
  arrange(desc(n))  # confirm that is empty

# baseline date per participant (pick earliest non-NA V1 date)
baseline_dates <- alfa_study %>%
  filter(visit == "V1") %>%
  transmute(alfa_id, v1_date = as.Date(date)) %>%
  filter(!is.na(v1_date)) %>%
  group_by(alfa_id) %>%
  summarise(v1_date = min(v1_date), .groups = "drop")

alfa_study <- alfa_study %>%
  mutate(date = as.Date(date)) %>%
  left_join(baseline_dates, by = "alfa_id") %>%
  mutate(time_years = as.numeric(date - v1_date) / 365.25)

# quick sanity
range(alfa_study$v1_date)


## ---- 2.8 Recode Sex, Visit and ATN status/_bl ----
alfa_study <- alfa_study %>%
  mutate(
    alfa_id = factor(alfa_id),
    visit   = factor(as.character(visit), levels = c("V1","V2","V2_5","V3","V3_5")),
    sex     = factor(sex, levels = c(1, 2), labels = c("Men", "Women")),
    
    # visit-level A/T/N status: 0/1 -> A-/A+ etc.
    A_status = factor(dplyr::case_when(
      A_status == 1 ~ "A+",
      A_status == 0 ~ "A-",
      TRUE ~ NA_character_
    ), levels = c("A-","A+")),
    
    T_status = factor(dplyr::case_when(
      T_status == 1 ~ "T+",
      T_status == 0 ~ "T-",
      TRUE ~ NA_character_
    ), levels = c("T-","T+")),
    
    N_status = factor(dplyr::case_when(
      N_status == 1 ~ "N+",
      N_status == 0 ~ "N-",
      TRUE ~ NA_character_
    ), levels = c("N-","N+")),
    
    # baseline A/T/N: 0/1 -> A-/A+ etc.
    A_bl = factor(dplyr::case_when(
      A_bl == 1 ~ "A+",
      A_bl == 0 ~ "A-",
      TRUE ~ NA_character_
    ), levels = c("A-","A+")),
    
    T_bl = factor(dplyr::case_when(
      T_bl == 1 ~ "T+",
      T_bl == 0 ~ "T-",
      TRUE ~ NA_character_
    ), levels = c("T-","T+")),
    
    N_bl = factor(dplyr::case_when(
      N_bl == 1 ~ "N+",
      N_bl == 0 ~ "N-",
      TRUE ~ NA_character_
    ), levels = c("N-","N+")),
    
    # ensure AT labels are factors with consistent ordering
    AT_status = factor(as.character(AT_status),
                       levels = c("A-T-","A+T-","A-T+","A+T+")),
    
    AT_bl = factor(as.character(AT_bl),
                   levels = c("A-T-","A+T-","A-T+","A+T+")),
    
    # ensure ATN labels are ordered factors
    ATN_status = factor(as.character(ATN_status),
                        levels = c("A-T-N-","A-T+N-","A-T-N+","A-T+N+",
                                   "A+T-N-","A+T+N-","A+T-N+","A+T+N+"),
                        ordered = TRUE),
    
    ATN_bl = factor(as.character(ATN_bl),
                    levels = c("A-T-N-","A-T+N-","A-T-N+","A-T+N+",
                               "A+T-N-","A+T+N-","A+T-N+","A+T+N+"),
                    ordered = TRUE)
  )


## ---- 2.9 Sanity check for LME readiness ----
cat("\nIADL non-missing by visit (final sample):\n")
print(
  alfa_study %>%
    group_by(visit) %>%
    summarise(n_rows = n(),
              n_iadl_nonmiss = sum(!is.na(iadl_score)),
              pct_iadl_nonmiss = round(100 * mean(!is.na(iadl_score)), 1),
              .groups = "drop")
)

cat("\nBaseline participants (distinct IDs): ",
    n_distinct(alfa_study$alfa_id[alfa_study$visit == "V1"]), "\n", sep = "")

cat("\nTime_years range:\n")
print(range(alfa_study$time_years, na.rm = TRUE))

cat("\nTime_years at V1 (should be 0):\n")
print(summary(alfa_study$time_years[alfa_study$visit == "V1"]))

# ============================================================================= #
# 3. EXPLORATORY ANALYSIS (DEMOGRAPHIC TABLE) ----
# ============================================================================= #

## ---- 3.1 Follow-up summary (n_visits & fu_time_years) ----
fu_summary <- alfa_study %>%
  filter(visit %in% c("V1", "V2", "V3")) %>%
  group_by(alfa_id) %>%
  summarise(
    # Baseline date (V1)
    v1_date = min(date[visit == "V1"], na.rm = TRUE),
    
    # Number of visits with non-missing IADL across V1–V3
    n_visits = n_distinct(visit[visit %in% c("V1", "V2", "V3") & !is.na(iadl_score)]),
    
    # Last follow-up date with IADL (V2 or V3)
    fu_last_date = max(date[visit %in% c("V2", "V3") & !is.na(iadl_score)], na.rm = TRUE),
    fu_last_date = if_else(is.infinite(fu_last_date), as.Date(NA), fu_last_date),
    
    # Follow-up time in years (V1 → last FU)
    fu_time_years = if_else(
      !is.na(fu_last_date),
      as.numeric(fu_last_date - v1_date) / 365.25,
      0
    ),
    .groups = "drop"
  ) %>%
  select(alfa_id, n_visits, fu_time_years)

## ---- 3.2 IADL scores at each visit in wide format ----
iadl_wide <- alfa_study %>%
  filter(visit %in% c("V1", "V2", "V3")) %>%
  select(alfa_id, visit, iadl_score) %>%
  tidyr::pivot_wider(
    names_from  = visit,
    values_from = iadl_score,
    names_prefix = "iadl_"
  )

## ---- 3.3 Baseline dataset: 1 row per participant at V1 ----
baseline_tbl <- alfa_study %>%
  filter(visit == "V1") %>%
  distinct(alfa_id, .keep_all = TRUE) %>%
  left_join(fu_summary, by = "alfa_id") %>%
  left_join(iadl_wide,  by = "alfa_id") %>%
  mutate(
    # sex is already text in your dataset
    sex = case_when(
      sex %in% c("Men", "Man", "M", "1") ~ "Men",
      sex %in% c("Women", "Woman", "W", "2") ~ "Women",
      TRUE ~ NA_character_     # anything else becomes NA (no "Unknown")
    ),
    sex = factor(sex, levels = c("Men", "Women")),
    AT_bl = droplevels(AT_bl)
  )

## ---- 3.4 Demographic table by AT_bl (no biomarkers) ----
demo_bl_by_AT <- baseline_tbl %>%
  select(
    AT_bl,
    Age_V1,
    YearsEducation,
    sex,
    iadl_V1,
    iadl_V2,
    iadl_V3,
    n_visits,
    fu_time_years
  ) %>%
  tbl_summary(
    by = AT_bl,
    missing = "no",
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous() ~ 1,
      fu_time_years ~ 2
    ),
    label = list(
      #AT_bl        ~ "AT group (baseline)",
      Age_V1       ~ "Age at baseline (years)",
      YearsEducation ~ "Education (years)",
      sex          ~ "Sex",
      iadl_V1      ~ "A-IADL-Q score at V1 (baseline)",
      iadl_V2      ~ "A-IADL-Q score at V2",
      iadl_V3      ~ "A-IADL-Q score at V3",
      n_visits     ~ "Number of visits with A-IADL-Q",
      fu_time_years ~ "Follow-up time (years)"
    )
  ) %>%
  add_overall() %>%
  add_p(test = list(
    all_continuous()  ~ "kruskal.test",
    all_categorical() ~ "chisq.test"
  )) %>%
  add_n() %>%
  bold_labels()

demo_bl_by_AT

## ---- 3.5 Improve looks ----
tab_gt <-
  demo_bl_by_AT %>%
  as_gt() %>%
  tab_options(
    table.font.names = "Arial",
    table.font.size = 10,
    data_row.padding = px(2),
    heading.align = "left",
    column_labels.font.weight = "bold"
  ) %>%
  opt_row_striping() %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  )

tab_gt

## ---- 3.6 Save demographic table ----
gtsave(tab_gt, "Table1_Demographics.png", vwidth = 1800)

# ============================================================================= #
# 4. INCLUSION CRITERIA, CLEANING AND RECODING ----
# ============================================================================= #

## ---- 4.1 Compact print helper ----
qc_lme <- function(df, label) {
  cat("\n---", label, "---\n")
  cat("Rows:", nrow(df), "| Participants:", dplyr::n_distinct(df$alfa_id), "\n")
  cat("Visits:\n"); print(table(df$visit, useNA = "ifany"))
  if ("AT_bl" %in% names(df)) { cat("AT_bl:\n"); print(table(df$AT_bl, useNA="ifany")) }
  if ("A_bl"  %in% names(df)) { cat("A_bl:\n");  print(table(df$A_bl,  useNA="ifany")) }
  if ("T_bl"  %in% names(df)) { cat("T_bl:\n");  print(table(df$T_bl,  useNA="ifany")) }
  invisible(df)
}

## ---- 4.2 Global LME dataset ----
dt0 <- alfa_study %>%
  filter(
    !is.na(iadl_score),
    !is.na(time_years),
    !is.na(A_bl),
    !is.na(Age_V1),
    !is.na(sex),
    !is.na(YearsEducation)
  ) %>%
  mutate(
    time2 = time_years^2,
    Age_c = Age_V1 - mean(Age_V1, na.rm = TRUE),
    Edu_c = YearsEducation - mean(YearsEducation, na.rm = TRUE)
  )

qc_lme(dt0, "dt0 (global complete-case)")
cat("AT_bl levels:", paste(levels(dt0$AT_bl), collapse = ", "), "\n")

# Optional: quick cross-tab (often useful)
cat("\nA_bl by visit:\n")
print(table(dt0$A_bl, dt0$visit, useNA = "ifany"))

## ---- 4.3 A+ subset LME dataset ----
dt_Apos <- alfa_study %>%
  filter(
    A_bl == "A+",
    !is.na(iadl_score),
    !is.na(time_years),
    !is.na(Age_V1),
    !is.na(sex),
    !is.na(YearsEducation)
  ) %>%
  mutate(
    time2 = time_years^2,
    Age_c = Age_V1 - mean(Age_V1, na.rm = TRUE),
    Edu_c = YearsEducation - mean(YearsEducation, na.rm = TRUE)
  )

qc_lme(dt_Apos, "dt_Apos (A+ complete-case)")
cat("T_bl levels:", paste(levels(dt_Apos$T_bl), collapse = ", "), "\n")

## ---- 4.4 FIXED-EFFECTS SHAPE SELECTION ----

# fit with max likelihood (i.e. REML = FALSE); random intercept only) and 
# compare linear vs quadratic time (nested models)

m1 <- lmer(iadl_score ~ time_years + Age_V1 + sex + YearsEducation + (1|alfa_id),
           data=dt0, REML=FALSE)

m2 <- lmer(iadl_score ~ time_years + Age_V1 + sex + YearsEducation + (1+time_years|alfa_id),
           data=dt0, REML=FALSE)

m2u <- lmer(iadl_score ~ time_years + Age_V1 + sex + YearsEducation + (1+time_years||alfa_id),
            data=dt0, REML=FALSE)

m3 <- lmer(iadl_score ~ time_years + time2 + Age_V1 + sex + YearsEducation + (1|alfa_id),
           data=dt0, REML=FALSE)

m4 <- lmer(iadl_score ~ time_years + time2 + Age_V1 + sex + YearsEducation + (1+time_years||alfa_id),
           data=dt0, REML=FALSE)

anova(m1, m2u)     # random slope?
isSingular(m2)
anova(m1, m3)     # quadratic?
anova(m2u, m4)    # quadratic given slope (uncorrelated for stability)
AIC(m1, m2, m2u, m3, m4)


# Based on the results, linear time + uncorrelated random slope seems the best options
final_structure <- "(1 + time_years || alfa_id)"

## ---- 4.5 Final global model ----
mod_A <- lmer(
  iadl_score ~ time_years * A_bl  + Age_c + sex + Edu_c +
    (1 + time_years || alfa_id),
  data = dt0, REML = TRUE
)

## ---- 4.6 Model evaluation ----
summary(mod_A)
tab_model(mod_A)
plot_model(mod_A, type = "int")

check_model(mod_A)

# ============================================================================= #
# 5. ANALYSIS 1: TIME x AT_GROUP LONGITUDINAL TRAJECTORIES ----
# ============================================================================= #

## ---- 5.1 Model creation and prediction----
mod_AT <- lmer(
  iadl_score ~ time_years * AT_bl + Age_c + sex + Edu_c +
    (1 + time_years || alfa_id),
  data = dt0, REML = TRUE
)
summary(mod_AT)

tab_model(mod_AT)

pp <- ggpredict(mod_AT, terms = c("time_years", "AT_bl"))

## ---- 5.2 Create table with estimates, CI and p-value ----
rng <- range(dt0$time_years, na.rm = TRUE)

sjPlot::plot_model(
  mod_AT,
  type = "int",
  terms = c("time_years", "AT_bl [quart]")
) +
  coord_cartesian(xlim = rng) +
  scale_x_continuous(limits = rng)
check_model(mod_AT)

emm <- emtrends(mod_AT, ~ AT_bl, var = "time_years")
emm

## ---- 5.3 A+ only: T status (baseline) — within A+ baseline + slope effects ----

mod_T_Apos_rs <- lmer(
  iadl_score ~ time_years * T_bl + Age_c + sex + Edu_c +
    (1 | alfa_id),
  data = dt_Apos, REML = TRUE
)
summary(mod_T_Apos_rs)
lme4::isSingular(mod_T_Apos_rs)
plot_model(mod_T_Apos_rs, type = "int")

## ---- 5.4 Plot individual trajectories ----
plot_dt <- dt0 %>%
  filter(visit %in% c("V1", "V2", "V3")) %>%
  mutate(
    AT_bl = droplevels(AT_bl),
    sex   = droplevels(factor(sex))   # ensure factor
  )

# A grid of times to draw smooth fitted lines
time_grid <- seq(
  min(plot_dt$time_years, na.rm = TRUE),
  max(plot_dt$time_years, na.rm = TRUE),
  length.out = 100
)

# Pprediction grid (fixed effects only)
newdata <- expand.grid(
  time_years = time_grid,
  AT_bl      = levels(plot_dt$AT_bl),
  Age_c      = 0,
  Edu_c      = 0,
  sex        = levels(plot_dt$sex)[1]
)

# Enforce factor levels to match plot_dt
newdata$AT_bl <- factor(newdata$AT_bl, levels = levels(plot_dt$AT_bl))
newdata$sex   <- factor(newdata$sex,   levels = levels(plot_dt$sex))

# Population-level predictions (re.form = NA)
newdata$pred <- as.numeric(predict(mod_AT, newdata = newdata, re.form = NA))

# Plot
graphics.off()
cols <- c(
  "A-T-" = "#2A9D8F",
  "A+T-" = "#F28E2B",
  "A+T+" = "#7A5195"
)

gg_pub2 <- ggplot() +
  geom_line(
    data = plot_dt,
    aes(time_years, iadl_score, group = alfa_id, color = AT_bl),
    alpha = 0.25, linewidth = 0.25
  ) +
  geom_point(
    data = plot_dt,
    aes(time_years, iadl_score, color = AT_bl),
    alpha = 0.45, size = 1
  ) +
  geom_ribbon(
    data = pp,
    aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.15, colour = NA,
    show.legend = FALSE
  ) +
  guides(fill = "none") +
  
  geom_line(data = pp, aes(x = x, y = conf.low,  color = group), linewidth = 0.15, alpha = 0.1) +
  geom_line(data = pp, aes(x = x, y = conf.high, color = group), linewidth = 0.15, alpha = 0.1) +
  geom_line(
    data = pp,
    aes(x = x, y = predicted, color = group),
    linewidth = 1.0
  ) +
  scale_x_continuous(limits = c(0, 9), breaks = 0:9) +
  coord_cartesian(ylim = c(40, 70)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  labs(
    x = "Years from baseline",
    y = "A-IADL-Q score"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.1),
    axis.title = element_text(face = "bold")
  )

gg_pub2
## ---- 5.5 Save plot ----
ggsave("A-IADLQ_trajectories_AT_square.pdf", gg_pub2, width = 5.5, height = 5.5)
ggsave("A-IADLQ_trajectories_AT_square.png", gg_pub2, width = 5.5, height = 5.5, dpi = 600)


## ---- 5.6 Table of slope trends by AT group----
tr <- emtrends(mod_AT, specs = "AT_bl", var = "time_years")

# slopes + within-group p-values
slopes <- as.data.frame(summary(tr, infer = c(TRUE, TRUE))) %>%
  dplyr::transmute(
    AT_group = AT_bl,
    slope = time_years.trend,
    lwr = lower.CL,
    upr = upper.CL,
    p_slope = p.value
  )

# differences in slopes vs reference (time×AT interaction contrasts)
pvals_diff <- as.data.frame(
  summary(
    contrast(tr, method = "trt.vs.ctrl", ref = "A-T-"),
    infer = TRUE,
    adjust = "none"   # <-- key
  )
) %>%
  dplyr::transmute(
    AT_group = sub("^\\(([^)]+)\\).*$", "\\1", contrast),
    p_diff_vs_ref = p.value
  )

tab_panel <- slopes %>%
  dplyr::left_join(pvals_diff, by = "AT_group") %>%
  dplyr::mutate(
    CI = sprintf("%.2f to %.2f", lwr, upr),
    p_slope = format.pval(p_slope, digits = 3, eps = 0.001),
    p_diff_vs_ref = ifelse(is.na(p_diff_vs_ref), "Ref",
                           format.pval(p_diff_vs_ref, digits = 3, eps = 0.001))
  ) %>%
  dplyr::select(AT_group, slope, CI, p_slope, p_diff_vs_ref) %>%
  dplyr::mutate(AT_group = factor(AT_group, levels = c("A-T-", "A+T-", "A+T+"))) %>%
  dplyr::arrange(AT_group)

tab_panel

gt_tab <- tab_panel %>%
  gt() %>%
  fmt_number(columns = slope, decimals = 2) %>%
  cols_label(
    AT_group   = "AT group",
    slope      = html("Slope<br>(points/year)"),
    CI         = "95% CI",
    p_diff_vs_ref = html("p (time×AT)<br>vs A−T−")
  ) %>%
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  tab_style(
    style = cell_text(align = "center", v_align = "middle", weight = "bold"),
    locations = cells_column_labels(everything())
  )  %>%
  tab_options(
    table.font.names = "Arial",
    table.font.size = 10,
    data_row.padding = px(2)
  ) %>%
  opt_row_striping()

gt_tab <- gt_tab %>%
  tab_options(column_labels.padding = px(6))

## ---- 5.7 Save tables ----

gtsave(gt_tab, "panel_decline_table.png", vwidth = 1200)
gtsave(gt_tab, "panel_decline_table.pdf")

gt_tab

## ---- 5.8 Plot trajectories by AT group ----

# prediction grid includes BOTH sexes, and facets will split by AT_bl
newdata <- expand.grid(
  time_years = time_grid,
  AT_bl      = levels(plot_dt$AT_bl),
  sex        = levels(plot_dt$sex),
  Age_c      = 0,
  Edu_c      = 0
)

newdata$AT_bl <- factor(newdata$AT_bl, levels = levels(plot_dt$AT_bl))
newdata$sex   <- factor(newdata$sex,   levels = levels(plot_dt$sex))

# Use mod_AT (sex additive). If you fitted mod_AT_sex with interactions, use that instead.
newdata$pred <- as.numeric(predict(mod_AT, newdata = newdata, re.form = NA))

gg_ATfacet <- ggplot(plot_dt, aes(x = time_years, y = iadl_score, group = alfa_id)) +
  scale_x_continuous(limits = c(0, 9), breaks = seq(0, 9, 1)) +
  geom_line(aes(color = sex), alpha = 0.40, linewidth = 0.35) +
  geom_point(aes(color = sex), alpha = 0.4, size = 0.7) +
  geom_line(
    data = newdata,
    aes(x = time_years, y = pred, color = sex, group = sex),
    linewidth = 1.4, alpha = 1
  ) +
  facet_wrap(~ AT_bl, nrow = 1) +
  labs(
    title = "A-IADL-Q trajectories over time (faceted by baseline AT group)",
    x = "Time since baseline (years)",
    y = "A-IADL-Q score",
    color = "Sex"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(gg_ATfacet)

# ============================================================================= #
# 6. ANALYSIS 2: TIME x BIOMARKER LONGITUDINAL TRAJECTORIES ----
# ============================================================================= #

## ---- 6.1 Helpers (data transformation) ----
safe_ratio <- function(num, den) ifelse(!is.na(num) & !is.na(den) & den != 0, num/den, NA_real_)
zscore <- function(x) as.numeric(scale(x))

winsorize <- function(x, p = 0.01) {
  lo <- quantile(x, p, na.rm = TRUE); hi <- quantile(x, 1-p, na.rm = TRUE)
  pmin(pmax(x, lo), hi)
}
iqr_trim <- function(x, k = 1.5) {
  q1 <- quantile(x, 0.25, na.rm = TRUE); q3 <- quantile(x, 0.75, na.rm = TRUE)
  i <- q3 - q1; lo <- q1 - k*i; hi <- q3 + k*i
  ifelse(!is.na(x) & (x < lo | x > hi), NA_real_, x)
}
process_marker <- function(x, mode = c("winsor_z","iqrtrim_z","zonly")) {
  mode <- match.arg(mode)
  if (mode == "winsor_z")  return(zscore(winsorize(x)))
  if (mode == "iqrtrim_z") return(zscore(iqr_trim(x)))
  zscore(x)
}
safe_ratio_cols <- function(df, num_col, den_col, out_name) {
  if (!(num_col %in% names(df)) || !(den_col %in% names(df))) {
    message("Ratio not created (missing cols): ", out_name,
            " | need: ", num_col, " & ", den_col)
    df[[out_name]] <- NA_real_
    return(df)
  }
  df[[out_name]] <- safe_ratio(df[[num_col]], df[[den_col]])
  df
}

# family assignment for raw + ratios (vectorized)
assign_family <- function(x) {
  dplyr::case_when(
    str_detect(x, regex("CSF", ignore_case = TRUE)) ~ "CSF",
    str_detect(x, regex("PLASMA", ignore_case = TRUE)) ~ "PLASMA",
    str_detect(x, "^CSF_R_") ~ "CSF",
    str_detect(x, "^PLASMA_R_") ~ "PLASMA",
    TRUE ~ NA_character_
  )
}

pretty_bio <- function(x) {
  x %>%
    str_remove("^BIO_BL_") %>%
    str_remove("(_w_z|_iqr_z|_z)$") %>%
    str_replace_all("^BIO_", "") %>%
    str_replace_all("_", " ")
}

lmer_ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

## ---- 6.2 Baseline (V1) biomarker dataset: one row per alfa_id ----
baseline_v1 <- alfa_study %>% filter(visit == "V1")

dup_v1 <- baseline_v1 %>% count(alfa_id) %>% filter(n > 1)
if (nrow(dup_v1) > 0) {
  cat("\nWARNING: duplicated V1 rows for some alfa_id:\n")
  print(head(dup_v1, 20))
}

baseline_bio <- baseline_v1 %>%
  distinct(alfa_id, .keep_all = TRUE)

bio_vars <- names(baseline_bio)[str_detect(names(baseline_bio), "^BIO_")]
cat("BIO vars:", length(bio_vars), "\n")

## ---- 6.3 Compute ratios (raw) at baseline ----
baseline_bio <- baseline_bio %>%
  safe_ratio_cols("BIO_Abeta1.42_PLASMA_Elecsys_UGOT_Roche",
                  "BIO_Abeta1.40_PLASMA_Elecsys_UGOT_Roche",
                  "PLASMA_R_ab42_ab40") %>%
  safe_ratio_cols("BIO_pTau181_PLASMA_Elecsys_UGOT_Roche",
                  "BIO_Abeta1.42_PLASMA_Elecsys_UGOT_Roche",
                  "PLASMA_R_pt181_ab42") %>%
  safe_ratio_cols("BIO_pTau217_PLASMA_Elecsys_UGOT_Roche",
                  "BIO_Abeta1.42_PLASMA_Elecsys_UGOT_Roche",
                  "PLASMA_R_pt217_ab42") %>%
  safe_ratio_cols("BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche",
                  "BIO_Abeta1.40_CSF_Elecsys_UGOT_Roche",
                  "CSF_R_ab42_ab40")

PLASMA_ratios <- c("PLASMA_R_ab42_ab40", "PLASMA_R_pt181_ab42", "PLASMA_R_pt217_ab42")
CSF_ratios    <- c("CSF_R_ab42_ab40")

## ---- 6.4 QC baseline availability + histograms (raw markers + ratios) ----
vars_qc <- unique(c(bio_vars, PLASMA_ratios, CSF_ratios))

qc_tbl <- baseline_bio %>%
  summarise(across(any_of(vars_qc), ~ sum(!is.na(.x)), .names = "{.col}")) %>%
  pivot_longer(everything(), names_to = "var", values_to = "n_nonmiss") %>%
  mutate(family = assign_family(var)) %>%
  arrange(family, n_nonmiss)

cat("\nQC non-missing counts (baseline V1):\n")
print(qc_tbl, n = Inf)

p_hist <- baseline_bio %>%
  dplyr::select(any_of(vars_qc)) %>%
  pivot_longer(everything(), names_to = "var", values_to = "value") %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ var, scales = "free", ncol = 3) +
  theme_minimal(base_size = 11) +
  labs(title = "Baseline distributions (raw biomarkers + ratios)", x = NULL, y = "Count")

print(p_hist)

## ---- 6.5 Settings (edit here only) ----
bio_mode <- "winsor_z"   # "zonly" | "winsor_z" | "iqrtrim_z"
min_id   <- 50        # skip model if < min_id participants with non-missing biomarker

# raw biomarkers to exclude
drop_raw <- c(
  "BIO_Abeta1.40_PLASMA_Elecsys_UGOT_Roche",
  "BIO_pTau181_CSF_Elecsys_UGOT_Roche",
  "BIO_tTau_CSF_Elecsys_UGOT_Roche",
  "BIO_Abeta1.40_CSF_Elecsys_UGOT_Roche"
)

# ratios to exclude (if any)
drop_ratios <- c()

# output folder
results_dir <- file.path(getwd(), "results_biomarkers")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

analysis_tag <- paste0("BIO_", bio_mode, "_minID", min_id)

## ---- 6.6 Histogram of biomarkers extracted with all three methods ----
drop_set <- unique(c(drop_raw, drop_ratios))
vars_after <- setdiff(vars_qc, drop_set)

bio_raw_long <- baseline_bio %>%
  dplyr::select(any_of(vars_after)) %>%
  pivot_longer(everything(), names_to = "var", values_to = "value_raw")

plot_mode <- function(mm) {
  dfp <- bio_raw_long %>%
    group_by(var) %>%
    mutate(value = process_marker(value_raw, mode = mm)) %>%
    ungroup() %>%
    filter(!is.na(value))
  
  ggplot(dfp, aes(x = value)) +
    geom_histogram(bins = 30) +
    facet_wrap(~ var, scales = "free", ncol = 3) +
    theme_minimal(base_size = 11) +
    labs(title = paste0("Baseline distributions (AFTER drop) — ", mm),
         x = NULL, y = "Count") +
    theme(strip.text = element_text(size = 8))
}

p1 <- plot_mode("zonly");    print(p1)
p2 <- plot_mode("winsor_z"); print(p2)
p3 <- plot_mode("iqrtrim_z");print(p3)

## ---- 6.7 Apply settings ----
biomarkers_use    <- setdiff(bio_vars, drop_raw)
PLASMA_ratios_use <- setdiff(PLASMA_ratios, drop_ratios)
CSF_ratios_use    <- setdiff(CSF_ratios, drop_ratios)

vars_use <- unique(c(biomarkers_use, PLASMA_ratios_use, CSF_ratios_use))
cat("\nUsing", length(vars_use), "predictors (biomarkers + ratios)\n")

bio_family <- tibble(var = vars_use) %>%
  mutate(family = assign_family(var))

if (any(is.na(bio_family$family))) {
  print(bio_family %>% filter(is.na(family)))
  stop("Some predictors could not be assigned to CSF/PLASMA family.", call.=FALSE)
}

csf_family_raw    <- bio_family %>% filter(family == "CSF")    %>% pull(var)
plasma_family_raw <- bio_family %>% filter(family == "PLASMA") %>% pull(var)

## ---- 6.8 Preprocess baseline predictors (BIO_BL_<var><suffix>) ----
suffix <- switch(bio_mode, winsor_z="_w_z", iqrtrim_z="_iqr_z", zonly="_z")

baseline_bio <- baseline_bio %>%
  mutate(across(
    all_of(vars_use),
    ~ process_marker(.x, mode = bio_mode),
    .names = paste0("BIO_BL_{.col}", suffix)
  ))

csf_family_pre    <- paste0("BIO_BL_", csf_family_raw, suffix)
plasma_family_pre <- paste0("BIO_BL_", plasma_family_raw, suffix)

bio_bl_pre <- baseline_bio %>%
  dplyr::select(alfa_id, any_of(c(csf_family_pre, plasma_family_pre)))

dt_bio <- dt0 %>% left_join(bio_bl_pre, by = "alfa_id")


## ---- 6.9 Time x biomarker screen ----

# A) Manual labels (edit as needed) 
bio_labels_raw <- c(
  "BIO_Abeta1.42_PLASMA_Elecsys_UGOT_Roche" = "Plasma Aβ42",
  "BIO_pTau181_PLASMA_Elecsys_UGOT_Roche"   = "Plasma pTau181",
  "BIO_pTau217_PLASMA_Elecsys_UGOT_Roche"   = "Plasma pTau217 (Roche)",
  "BIO_pTau217_PLASMA_MSD_LILLY_in.house"    = "Plasma pTau217 (Lilly)",
  "BIO_GFAP_PLASMA_Elecsys_UGOT_Roche"      = "Plasma GFAP",
  "BIO_NFL_PLASMA_Elecsys_UGOT_Roche"       = "Plasma NfL",
  "PLASMA_R_pt181_ab42"                     = "Plasma pTau181/Aβ42",
  "PLASMA_R_pt217_ab42"                     = "Plasma pTau217/Aβ42",
  "PLASMA_R_ab42_ab40"                      = "Plasma Aβ42/Aβ40",
  "BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche"     = "CSF Aβ42",
  "CSF_R_ab42_ab40"                         = "CSF Aβ42/Aβ40",
  "BIO_GFAP_CSF_Elecsys_UGOT_Roche"          = "CSF GFAP",
  "BIO_NFL_CSF_Elecsys_UGOT_Roche"           = "CSF NfL",
  "BIO_pTau217_CSF_MSD_LILLY_in.house"       = "CSF pTau217",
  "BIO_pTau231_CSF_ELISA_UGOT_ADx.NeuroSciences" = "CSF pTau231"
)

# build processed-name label map for the CURRENT suffix automatically
bio_labels <- setNames(
  unname(bio_labels_raw),
  paste0("BIO_BL_", names(bio_labels_raw), suffix)
)

pretty_bio <- function(x) dplyr::recode(x, !!!bio_labels, .default = x)

# B) Biomarker list (robust: only keep biomarkers that exist in dt_bio)
all_biomarkers <- unique(c(csf_family_pre, plasma_family_pre))
all_biomarkers <- all_biomarkers[all_biomarkers %in% names(dt_bio)]

# Optional safety: if the vectors were empty, fall back to the label map names present in dt_bio
if (length(all_biomarkers) == 0) {
  all_biomarkers <- intersect(names(bio_labels), names(dt_bio))
}

# C) NA check (participant-level complete-case N for model)
na_check <- purrr::map_dfr(all_biomarkers, function(b) {
  d_model <- dt_bio %>%
    filter(!is.na(time_years), !is.na(iadl_score), !is.na(.data[[b]]))
  tibble(
    biomarker = b,
    biomarker_label = pretty_bio(b),
    family = assign_family(b),
    n_id_complete_for_model = n_distinct(d_model$alfa_id)
  )
})

# D) Fit one model per biomarker; extract time×biomarker interaction
fit_time_bio <- function(dt, bio, min_id = 50, lmer_ctrl = NULL) {
  
  d_model <- dt %>%
    filter(!is.na(time_years), !is.na(iadl_score), !is.na(.data[[bio]]))
  
  n_id <- n_distinct(d_model$alfa_id)
  
  if (n_id < min_id) {
    return(tibble(
      biomarker = bio,
      biomarker_label = pretty_bio(bio),
      family = assign_family(bio),
      n_id = n_id,
      beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_nominal = NA_real_,
      note = "Skipped: too few participants"
    ))
  }
  
  f <- as.formula(paste0(
    "iadl_score ~ time_years * `", bio, "` + Age_c + sex + Edu_c + (1 + time_years || alfa_id)"
  ))
  
  tryCatch({
    m  <- lmerTest::lmer(f, data = d_model, REML = TRUE, control = lmer_ctrl)
    tt <- broom.mixed::tidy(m, effects = "fixed")
    
    # Interaction term can appear as time_years:`bio` (with backticks)
    term_pat <- paste0("^time_years:`", bio, "`$|^`", bio, "`:time_years$")
    
    # robustly find interaction term time_years:bio OR bio:time_years (no backticks assumed)
    row <- tt %>%
      dplyr::filter(term %in% c(paste0("time_years:", bio), paste0(bio, ":time_years"))) %>%
      dplyr::slice(1)
    
    # fallback (if term names contain backticks or other formatting)
    if (nrow(row) == 0) {
      row <- tt %>%
        dplyr::filter(
          stringr::str_detect(term, stringr::fixed("time_years")) &
            stringr::str_detect(term, stringr::fixed(bio)) &
            stringr::str_detect(term, ":")
        ) %>%
        dplyr::slice(1)
    }
    
    tibble(
      biomarker = bio,
      biomarker_label = pretty_bio(bio),
      family = assign_family(bio),
      n_id = n_id,
      beta = row$estimate,
      se = row$std.error,
      ci_low = row$estimate - 1.96 * row$std.error,
      ci_high = row$estimate + 1.96 * row$std.error,
      p_nominal = row$p.value,
      note = NA_character_
    )
  }, error = function(e) {
    tibble(
      biomarker = bio,
      biomarker_label = pretty_bio(bio),
      family = assign_family(bio),
      n_id = n_id,
      beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_nominal = NA_real_,
      note = paste("Model failed:", e$message)
    )
  })
}

tbl_int <- purrr::map_dfr(all_biomarkers, ~ fit_time_bio(dt_bio, .x, min_id = min_id, lmer_ctrl = lmer_ctrl)) %>%
  group_by(family) %>%                                  # FDR separately within family
  mutate(p_fdr = p.adjust(p_nominal, method = "fdr")) %>%
  ungroup()

## ---- 6.10 Build final display table (FDR only) ----
tab_panel <- tbl_int %>%
  mutate(
    beta_ci = if_else(
      is.na(beta),
      NA_character_,
      sprintf("%.3f (%.3f, %.3f)", beta, ci_low, ci_high)
    ),
    p_fdr_fmt = if_else(
      is.na(p_fdr),
      NA_character_,
      format.pval(p_fdr, digits = 3, eps = 0.001)
    )
  ) %>%
  dplyr::select(
    family,
    Biomarker = biomarker_label,
    `N (participants)` = n_id,
    `Time×biomarker β (95% CI)` = beta_ci,
    `p (FDR)` = p_fdr_fmt
  ) %>%
  mutate(
    family = toupper(family),
    family = factor(family, levels = c("PLASMA", "CSF"))
  ) %>%
  arrange(family, Biomarker)

# gt table (two vertical sections: PLASMA / CSF)
gt_tab_bio <- tab_panel %>%
  # keep a numeric p column for styling
  mutate(p_fdr_num = `p (FDR)`) %>%
  gt(groupname_col = "family") %>%
  cols_label(
    Biomarker = "Biomarker",
    `N (participants)` = "N",
    `Time×biomarker β (95% CI)` = html("Time×biomarker<br>β (95% CI)"),
    `p (FDR)` = html("p<br>(FDR)")
  ) %>%
  # hide the numeric helper column
  cols_hide(columns = p_fdr_num) %>%
  
  tab_style(
    style = cell_text(weight = "bold", align = "center", v_align = "middle"),
    locations = cells_column_labels(everything())
  ) %>%
  # bold significant p-values (FDR < 0.05)
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = `p (FDR)`,
      rows = !is.na(p_fdr_num) & p_fdr_num < 0.05
    )
  ) %>%
  
  opt_row_striping()

gt_tab_bio_compact <- gt_tab_bio %>%
  tab_options(
    table.font.size = 12,          # slightly smaller text
    data_row.padding = px(1),     # tighter rows
    row_group.padding = px(1),
    column_labels.padding = px(3)
  )

gt_tab_bio_compact

## ---- 6.11 Save table ----
gtsave(gt_tab_bio, "time_x_biomarker_table.png", vwidth = 1600)
gtsave(gt_tab_bio, "time_x_biomarker_table.pdf")

## ---- 6.12 Plot predicted trajectories by biomarker tertiles ----

# Significant biomarkers from your earlier screen (FDR < 0.05)
sig_bios <- tbl_int %>%
  dplyr::filter(!is.na(p_fdr), p_fdr < 0.05) %>%
  dplyr::pull(biomarker)

tertile_labels <- c("T1 (low)", "T2 (mid)", "T3 (high)")

# Choose colors (edit hex codes): magnitude palette low -> high

tertile_cols <- c(
  "T1 (low)"  = "#FDE725",
  "T2 (mid)"  = "#35B779",
  "T3 (high)" = "#31688E"
)

plot_traj_by_bio_tertiles <- function(dt, bio_pre,
                                      tertile_labels = tertile_labels,
                                      tertile_cols = tertile_cols,
                                      y_lim = c(65, 70),
                                      lmer_ctrl = NULL) {
  
  # complete-case for time/outcome/biomarker
  d <- dt %>%
    dplyr::filter(!is.na(time_years), !is.na(iadl_score), !is.na(.data[[bio_pre]]))
  
  # tertiles defined at PARTICIPANT level from baseline biomarker distribution
  tert <- d %>%
    dplyr::select(alfa_id, !!rlang::sym(bio_pre)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(tertile = dplyr::ntile(.data[[bio_pre]], 3L)) %>%
    dplyr::mutate(tertile = factor(tertile, levels = 1:3, labels = tertile_labels)) %>%
    dplyr::select(alfa_id, tertile)
  
  dB <- d %>%
    dplyr::left_join(tert, by = "alfa_id") %>%
    dplyr::filter(!is.na(tertile)) %>%
    dplyr::mutate(tertile = factor(tertile, levels = tertile_labels))
  
  # model: time × tertile + covariates
  m <- lmerTest::lmer(
    iadl_score ~ time_years * tertile + Age_c + sex + Edu_c + (1 + time_years || alfa_id),
    data = dB, REML = TRUE, control = lmer_ctrl
  )
  
  # predicted trajectories + CI (robust)
  pp <- ggeffects::ggpredict(m, terms = c("time_years [all]", "tertile"))
  
  # Force exactly 3 groups (prevents missing lines due to factor/scale dropping)
  pp$group <- factor(as.character(pp$group), levels = tertile_labels)
  
  # Safety: stop early if palette/levels mismatch
  stopifnot(all(tertile_labels %in% levels(pp$group)))
  stopifnot(all(tertile_labels %in% names(tertile_cols)))
  
  # label position inside the panel (always visible given y_lim)
  lab_x <- min(pp$x, na.rm = TRUE) + 0.02 * diff(range(pp$x, na.rm = TRUE))
  lab_y <- y_lim[2] - 0.3
  
  p <- ggplot(pp, aes(x = x)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
                alpha = 0.12, colour = NA, show.legend = FALSE) +
    geom_line(aes(y = predicted, color = group),
              linewidth = 1.2) +
    annotate("text", x = lab_x, y = lab_y, label = pretty_bio(bio_pre),
             hjust = 0, vjust = 1, fontface = "bold", size = 6) +
    scale_color_manual(values = tertile_cols, limits = tertile_labels, drop = FALSE) +
    scale_fill_manual(values = tertile_cols, limits = tertile_labels, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = 0, add = 0)) +
    coord_cartesian(ylim = y_lim) +
    labs(x = "Years from baseline", y = "Predicted A-IADL-Q", color = NULL) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(linewidth = 0.2),
      axis.title = element_text(face = "bold"),
      plot.margin = margin(4, 4, 4, 4)
    )
  
  p
}

dir.create(file.path(results_dir, "traj_tertiles_sig"), showWarnings = FALSE, recursive = TRUE)

plots_traj <- purrr::map(sig_bios, ~ plot_traj_by_bio_tertiles(
  dt_bio, .x,
  tertile_labels = tertile_labels,
  tertile_cols = tertile_cols,
  y_lim = c(62.5, 70),
  lmer_ctrl = lmer_ctrl
))

purrr::walk2(sig_bios, plots_traj, \(bio, p) {
  fn <- paste0(gsub("[^A-Za-z0-9]+", "_", bio), "_traj_tertiles")
  
  ggsave(file.path(results_dir, "traj_tertiles_sig", paste0(fn, ".pdf")),
         p, width = 5.2, height = 5.2)
  
  ggsave(file.path(results_dir, "traj_tertiles_sig", paste0(fn, ".png")),
         p, width = 5.2, height = 5.2, dpi = 600)
})


## ---- 6.13 Save plots in a specific folder ----
dir.create(file.path(results_dir, "plots_traj_tertiles_sq"), showWarnings = FALSE, recursive = TRUE)


plots_sq <- purrr::map(sig_bios, ~ plot_traj_tertiles_sq(
  dt_bio, .x,
  min_id = min_id,
  lmer_ctrl = lmer_ctrl,
  tertile_labels = tertile_labels,
  tertile_cols = tertile_cols,
  show_individual = FALSE,
  y_lim = c(60, 70)
))

purrr::walk2(sig_bios, plots_sq, \(bio, p) {
  if (!is.null(p)) {
    fn <- paste0(gsub("[^A-Za-z0-9]+", "_", bio), "_tertiles_sq")
    
    ggsave(file.path(results_dir, "plots_traj_tertiles_sq", paste0(fn, ".pdf")),
           p, width = 5.2, height = 5.2)
    
    ggsave(file.path(results_dir, "plots_traj_tertiles_sq", paste0(fn, ".png")),
           p, width = 5.2, height = 5.2, dpi = 600)
  }
})




table(dB$tertile)

# ============================================================================= #
# 7. ANALYSIS 3: TIME X BIOMARKER x AT_GROUP ----
# ============================================================================= #

## ---- 7.1 Fit time x AT_bl x biomarker model ----
fit_one_biomarker_groupmod <- function(bio_pre, data, group_var, ref_level = NULL, min_id = 50) {
  stopifnot(group_var %in% names(data))
  d <- data %>%
    filter(!is.na(.data[[bio_pre]]),
           !is.na(.data[[group_var]]),
           !is.na(time_years),
           !is.na(iadl_score))
  
  n_id <- n_distinct(d$alfa_id)
  if (n_id < min_id) {
    return(list(
      biomarker=bio_pre, group_var=group_var,
      n_id=n_id, n_rows=nrow(d),
      p_3way=NA_real_, slopes=tibble(), note="Skipped: too few participants"
    ))
  }
  
  # ensure factor + reference
  d[[group_var]] <- factor(d[[group_var]])
  if (!is.null(ref_level) && ref_level %in% levels(d[[group_var]])) {
    d[[group_var]] <- relevel(d[[group_var]], ref = ref_level)
  }
  
  base_rhs <- "Age_c + sex + Edu_c + (1 + time_years || alfa_id)"
  
  # FULL: time * group * biomarker
  f_full <- as.formula(paste0(
    "iadl_score ~ time_years * ", group_var, " * `", bio_pre, "` + ", base_rhs
  ))
  
  # REDUCED: remove ONLY 3-way term; keep all lower order terms
  f_red <- as.formula(paste0(
    "iadl_score ~ time_years * ", group_var,
    " + time_years * `", bio_pre, "` + ", group_var, " * `", bio_pre, "` + ",
    base_rhs
  ))
  
  out <- tryCatch({
    # LRT: ML
    m_red  <- lmer(f_red,  data=d, REML=FALSE, control=lmer_ctrl)
    m_full <- lmer(f_full, data=d, REML=FALSE, control=lmer_ctrl)
    p_3way <- anova(m_red, m_full)[["Pr(>Chisq)"]][2]
    
    # refit full with REML for effects/CI
    m_reml <- lmer(f_full, data=d, REML=TRUE, control=lmer_ctrl)
    
    b <- fixef(m_reml)
    V <- as.matrix(vcov(m_reml))
    nms <- names(b)
    bio_pat <- paste0("`?", fixed(bio_pre), "`?")
    
    # time:bio (in reference group)
    t_bio <- nms[str_detect(nms, "^time_years:") & str_detect(nms, bio_pat) & !str_detect(nms, paste0("^time_years:", group_var))]
    if (length(t_bio) != 1) stop("Cannot uniquely find time_years:bio for ", bio_pre, " (", group_var, ")")
    
    # 3-way terms time:groupLevel:bio for each non-ref level
    t_3 <- nms[str_detect(nms, paste0("^time_years:", group_var)) & str_detect(nms, bio_pat)]
    t_3 <- sort(t_3)
    
    g_levels <- levels(model.frame(m_reml)[[group_var]])
    ref <- g_levels[1]
    nonref <- g_levels[-1]
    
    # map 3-way terms to nonref levels (treatment coding)
    if (length(nonref) != length(t_3)) {
      stop("Mismatch: group levels vs 3-way terms for ", bio_pre,
           "\nlevels: ", paste(g_levels, collapse=" | "),
           "\n3way: ", paste(t_3, collapse=" | "))
    }
    map3 <- setNames(t_3, nonref)
    
    lin <- function(terms, w) {
      L <- setNames(rep(0, length(b)), nms)
      L[terms] <- w
      est <- sum(L * b)
      se  <- sqrt(as.numeric(t(L) %*% V %*% L))
      z   <- est / se
      p   <- 2 * pnorm(abs(z), lower.tail = FALSE)
      ci  <- est + c(-1,1) * qnorm(0.975) * se
      tibble(estimate=est, se=se, conf.low=ci[1], conf.high=ci[2], p.value=p)
    }
    
    # within-group slope effects: ref uses time:bio; others add corresponding 3-way
    slopes <- map_dfr(g_levels, function(g) {
      if (g == ref) {
        eff <- lin(t_bio, 1)
      } else {
        eff <- lin(c(t_bio, map3[[g]]), c(1,1))
      }
      tibble(group=g) %>% bind_cols(eff)
    }) %>%
      mutate(group = factor(group, levels = g_levels))
    
    list(
      biomarker=bio_pre, group_var=group_var,
      n_id=n_id, n_rows=nrow(d),
      p_3way=p_3way, slopes=slopes, note=NA_character_
    )
  }, error=function(e) {
    list(
      biomarker=bio_pre, group_var=group_var,
      n_id=n_id, n_rows=nrow(d),
      p_3way=NA_real_, slopes=tibble(),
      note=paste("Model failed:", e$message)
    )
  })
  
  out
}

## ---- 7.2 Moderation by A_bl and by AT_bl ----
all_biomarkers <- unique(c(csf_family_pre, plasma_family_pre))

# Ensure expected references (change if you prefer)
ref_A  <- if ("A-" %in% levels(factor(dt_bio$A_bl)))  "A-"  else NULL
ref_AT <- if ("A-T-" %in% levels(factor(dt_bio$AT_bl))) "A-T-" else NULL

res_A  <- map(all_biomarkers, ~ fit_one_biomarker_groupmod(.x, dt_bio, group_var="A_bl",  ref_level=ref_A,  min_id=min_id))
res_AT <- map(all_biomarkers, ~ fit_one_biomarker_groupmod(.x, dt_bio, group_var="AT_bl", ref_level=ref_AT, min_id=min_id))

make_tables <- function(res_list, group_name) {
  tbl_main <- map_dfr(res_list, \(x) {
    tibble(
      biomarker = x$biomarker,
      family = assign_family(x$biomarker),
      biomarker_label = pretty_bio(x$biomarker),
      n_id = x$n_id, n_rows = x$n_rows,
      p_3way = x$p_3way,
      note = x$note,
      group = group_name
    )
  }) %>%
    group_by(family) %>%
    mutate(p_3way_fdr = p.adjust(p_3way, method="fdr")) %>%
    ungroup()
  
  tbl_slopes <- map_dfr(res_list, \(x) {
    if (nrow(x$slopes) == 0) return(tibble())
    
    x$slopes %>%
      rename(level = group) %>%              # <-- rename the slope's group column first
      mutate(
        biomarker = x$biomarker,
        family = assign_family(x$biomarker),
        biomarker_label = pretty_bio(x$biomarker),
        n_id = x$n_id, n_rows = x$n_rows,
        note = x$note,
        group = group_name                   # <-- now 'group' is the moderator name (A_bl / AT_bl)
      ) %>%
      relocate(group, biomarker, family, biomarker_label, level, n_id, n_rows)
  }) %>%
    left_join(tbl_main %>% dplyr::select(group, biomarker, p_3way, p_3way_fdr),
              by = c("group","biomarker")) %>%
    group_by(family, group, level) %>%
    mutate(p_within_fdr = p.adjust(p.value, method="fdr")) %>%
    ungroup()
  
  list(main=tbl_main, slopes=tbl_slopes)
}

tbl_A  <- make_tables(res_A,  "A_bl")
tbl_AT <- make_tables(res_AT, "AT_bl")

## ---- 7.3 Save tables ----
file_main_A   <- file.path(results_dir, paste0("bio_main_A_bl_", analysis_tag, ".csv"))
file_slopes_A <- file.path(results_dir, paste0("bio_slopes_A_bl_", analysis_tag, ".csv"))
file_main_AT  <- file.path(results_dir, paste0("bio_main_AT_bl_", analysis_tag, ".csv"))
file_slopes_AT<- file.path(results_dir, paste0("bio_slopes_AT_bl_", analysis_tag, ".csv"))

write_csv(tbl_A$main,   file_main_A)
write_csv(tbl_A$slopes, file_slopes_A)
write_csv(tbl_AT$main,  file_main_AT)
write_csv(tbl_AT$slopes,file_slopes_AT)

cat("\nSaved tables:\n",
    file_main_A, "\n", file_slopes_A, "\n",
    file_main_AT, "\n", file_slopes_AT, "\n", sep="")

## ---- 7.4 Create forest plots ----
cols_AT <- c(
  "A-T-" = "#2A9D8F",  # green / amyloid– tau–
  "A+T-" = "#F28E2B",  # orange / amyloid+ tau–
  "A+T+" = "#7A5195"   # purple / amyloid+ tau+
)

plot_forest <- function(slopes_df,
                            title = NULL,
                            order_level = "A+T+",
                            cols = cols_AT) {
  
  df <- slopes_df %>%
    filter(is.na(note),
           !is.na(estimate),
           !is.na(conf.low),
           !is.na(conf.high)) %>%
    mutate(
      level = as.character(level),
      
      # Etiquetas de significancia
      sig_label = ifelse(!is.na(p_3way_fdr) & p_3way_fdr < 0.05, "*", "")
    )
  
  if (nrow(df) == 0) return(invisible(NULL))
  
  sig_biomarkers <- df %>%
    group_by(biomarker_label) %>%
    summarise(any_sig = any(sig_label == "*"), .groups = "drop")
  
  df <- df %>%
    left_join(sig_biomarkers, by = "biomarker_label") %>%
    mutate(
      biomarker_label_fmt = ifelse(any_sig,
                                   paste0("**", "*", biomarker_label, "**"),
                                   biomarker_label)
    )
  
  # order biomarkers by estimate in one chosen level (e.g., A+T+)
  if (!is.null(order_level) && order_level %in% df$level) {
    ord <- df %>%
      filter(level == order_level) %>%
      distinct(biomarker_label_fmt, estimate, any_sig) %>%
      arrange(desc(any_sig), estimate) %>%
      pull(biomarker_label_fmt)
    
    df <- df %>%
      mutate(
        biomarker_label_fmt = factor(biomarker_label_fmt, levels = rev(ord))
      )
  }
  
  # keep only colors for levels that exist in df
  cols_use <- cols[names(cols) %in% unique(df$level)]
  
  dodge_w <- 0.65
  
  ggplot(df, aes(x = estimate, y = biomarker_label_fmt)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.7, alpha = 0.6, colour = "black") +
    
    geom_errorbarh(
      aes(xmin = conf.low, xmax = conf.high, colour = level),
      height = 0.18, linewidth = 0.7,
      position = position_dodge(width = dodge_w)
    ) +
    
    geom_point(
      aes(colour = level),
      position = position_dodge(width = dodge_w)
    ) +
    
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16)) +
    
    facet_wrap(~ family, scales = "free_y", ncol = 1, strip.position = "top") +
    
    scale_colour_manual(values = cols_use) +
    
    labs(
      title = title,
      x = expression(Delta~"yearly A-IADL-Q slope per +1 SD biomarker (95% CI)"),
      y = NULL
    ) +
    
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey85"),
      strip.background = element_rect(fill = "grey95", linewidth = 0.2),
      strip.text = element_text(face = "bold", size = 11),
      axis.title.x = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_markdown(size = 10),
    ) +
    scale_colour_manual(
      values = cols_use,
      name = "AT group (*3-way FDR p < 0.05)"
    )
  
}

p_AT <- plot_forest (
  tbl_AT$slopes,
  title = paste0("Biomarker effects on A-IADL-Q decline by AT status"),
  order_level = "A+T+",
  cols = cols_AT
)

print(p_AT)

ggsave(
  filename = file.path(results_dir, paste0("forest_AT_A4_landscape_", analysis_tag, ".pdf")),
  plot = p_AT_pub,
  width = 11.69, height = 8.27, units = "in"
)


ggsave(file.path(results_dir, paste0("forest_AT_pub_ (", bio_mode, ")", analysis_tag, ".png")),
       p_AT, width = 11.69, height = 8.27, dpi = 600)


# ============================================================================= #
# 8. -2.2 SELECTION OF OPTIMAL INTERVAL FOR DECLINE CLASIFICATION ----
# ============================================================================= #

## ---- 8.1 Create comparison dataset for comparison of V3-V1 vs V2-V1 ----
# Followup of participants
fu_summary <- alfa_study %>%
  filter(visit %in% c("V1", "V2", "V3")) %>%
  group_by(alfa_id) %>%
  summarise(
    # Baseline date (V1)
    v1_date = min(date[visit == "V1"], na.rm = TRUE),
    
    # Number of visits with non-missing IADL across V1–V3
    n_visits = n_distinct(visit[visit %in% c("V1", "V2", "V3") & !is.na(iadl_score)]),
    
    # Last follow-up date with IADL (V2 or V3)
    fu_last_date = max(date[visit %in% c("V2", "V3") & !is.na(iadl_score)], na.rm = TRUE),
    fu_last_date = if_else(is.infinite(fu_last_date), as.Date(NA), fu_last_date),
    
    # Follow-up time in years (V1 → last FU)
    fu_time_years = if_else(
      !is.na(fu_last_date),
      as.numeric(fu_last_date - v1_date) / 365.25,
      0
    ),
    .groups = "drop"
  ) %>%
  select(alfa_id, n_visits, fu_time_years)

# IADL scores at each visit in wide format
iadl_wide <- alfa_study %>%
  filter(visit %in% c("V1", "V2", "V3")) %>%
  select(alfa_id, visit, iadl_score) %>%
  tidyr::pivot_wider(
    names_from  = visit,
    values_from = iadl_score,
    names_prefix = "iadl_"
  )

# Baseline dataset: 1 row per participant at V1
baseline_tbl <- alfa_study %>%
  filter(visit == "V1") %>%
  distinct(alfa_id, .keep_all = TRUE) %>%
  left_join(fu_summary, by = "alfa_id") %>%
  left_join(iadl_wide,  by = "alfa_id") %>%
  mutate(
    # sex is already text in your dataset
    sex = case_when(
      sex %in% c("Men", "Man", "M", "1") ~ "Men",
      sex %in% c("Women", "Woman", "W", "2") ~ "Women",
      TRUE ~ NA_character_     # anything else becomes NA (no "Unknown")
    ),
    sex = factor(sex, levels = c("Men", "Women")),
    AT_bl = droplevels(AT_bl)
  )

## ---- 8.2 Create correlation table between V3 and V2 decliners ----
iadl_tbl <- alfa_study %>%
  select(alfa_id, visit, iadl_score) %>%
  pivot_wider(
    names_from = visit,
    values_from = iadl_score,
    names_prefix = "iadl_"
  ) %>%
  mutate(
    iadl_change_V2 = iadl_V2 - iadl_V1,
    iadl_change_V3 = iadl_V3 - iadl_V1,
    
    iadl_decline_V2 = case_when(
      !is.na(iadl_change_V2) & iadl_change_V2 <= -22 ~ "Decliner V2",
      !is.na(iadl_change_V2) ~ "Non-decliner",
      TRUE ~ NA_character_
    ),
    
    iadl_decline_V3 = case_when(
      !is.na(iadl_change_V3) & iadl_change_V3 <= -22 ~ "Decliner V3",
      !is.na(iadl_change_V3) ~ "Non-decliner",
      TRUE ~ NA_character_
    )
  )

tbl_overlap <-
  iadl_tbl %>%
  select(iadl_decline_V2, iadl_decline_V3) %>%
  tbl_cross(
    row = iadl_decline_V2,
    col = iadl_decline_V3
  ) %>%
  modify_header(label = "**V2–V1 decline**") %>%
  modify_spanning_header(all_stat_cols() ~ "**V3–V1 decline**")

tbl_overlap

## ---- 8.3 Identify transient decliners ----
transient_decliners <- iadl_tbl %>%
  filter(
    iadl_decline_V2 == "Decliner V2",
    iadl_decline_V3 == "Non-decliner"
  )

# IDs of interest
ids_transient <- transient_decliners$alfa_id

length(ids_transient)

# Add group label to baseline dataset
baseline_tbl2 <- baseline_tbl %>%
  mutate(
    transient_decliner = case_when(
      alfa_id %in% ids_transient ~ "Decliner V2-V1 only",
      TRUE ~ "Others"
    ),
    transient_decliner = factor(
      transient_decliner,
      levels = c("Others", "Decliner V2-V1 only")
    )
  )

transient_decliners %>%
  select(alfa_id, iadl_V1, iadl_V2, iadl_V3,
         iadl_change_V2, iadl_change_V3)

## ---- 8.4 Create dataset for LOESS trajectory plots ----
plot_iadl <- alfa_study %>%
  filter(visit %in% c("V1","V2","V3")) %>%
  mutate(
    transient_decliner = case_when(
      alfa_id %in% ids_transient ~ "Decliner V2-V1 only",
      TRUE ~ "Others"
    ),
    transient_decliner = factor(
      transient_decliner,
      levels = c("Others","Decliner V2-V1 only")
    )
  )

## ---- 8.5 LOESS trajectory plots only between V2-V1 and others ----
ggplot(plot_iadl,
       aes(x = visit, y = iadl_score,
           group = alfa_id,
           color = transient_decliner)) +
  
  geom_line(alpha = 0.25) +
  geom_point(size = 1.5, alpha = 0.7) +
  
  stat_summary(
    aes(group = transient_decliner),
    fun = mean,
    geom = "line",
    linewidth = 1.5
  ) +
  
  stat_summary(
    aes(group = transient_decliner),
    fun = mean,
    geom = "point",
    size = 3
  ) +
  
  scale_color_manual(values = c(
    "Others" = "grey50",
    "Decliner V2-V1 only" = "#D55E00"
  )) +
  
  labs(
    x = "Visit",
    y = "A-IADL-Q score",
    color = "Group"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(face = "bold")
  )

plot_iadl %>%
  filter(transient_decliner == "Decliner V2-V1 only") %>%
  ggplot(aes(visit, iadl_score, group = alfa_id)) +
  
  geom_line(alpha = 0.5) +
  geom_point(size = 2) +
  
  labs(
    title = "Participants declining at V2 but not at V3",
    x = "Visit",
    y = "A-IADL-Q score"
  ) +
  
  theme_bw()

plot_time <- dt0 %>%
  mutate(
    transient_decliner = case_when(
      alfa_id %in% ids_transient ~ "Decliner V2-V1 only",
      TRUE ~ "Others"
    )
  )

final_plot <- ggplot(plot_time,
                     aes(time_years, iadl_score,
                         group = alfa_id,
                         color = transient_decliner)) +
  
  geom_line(alpha = 0.25) +
  geom_point(alpha = 0.5) +
  
  stat_smooth(
    aes(group = transient_decliner),
    method = "loess",
    se = TRUE,
    linewidth = 1.2
  ) +
  
  theme_bw() +
  labs(
    x = "Years from baseline",
    y = "A-IADL-Q score"
  )

final_plot

## ---- 8.6 LOESS trajectory plots for all decliners groups ----
iadl_groups <- iadl_tbl %>%
  mutate(
    decline_group = case_when(
      iadl_decline_V2 == "Decliner V2" & iadl_decline_V3 == "Decliner V3" ~ "Decliner V2 & V3",
      iadl_decline_V2 == "Decliner V2" & iadl_decline_V3 == "Non-decliner" ~ "Decliner V2 only",
      iadl_decline_V2 == "Non-decliner" & iadl_decline_V3 == "Decliner V3" ~ "Decliner V3 only",
      iadl_decline_V2 == "Non-decliner" & iadl_decline_V3 == "Non-decliner" ~ "Non-decliner",
      TRUE ~ NA_character_
    ),
    decline_group = factor(
      decline_group,
      levels = c(
        "Non-decliner",
        "Decliner V2 only",
        "Decliner V3 only",
        "Decliner V2 & V3"
      )
    )
  )

plot_time2 <- dt0 %>%
  left_join(
    iadl_groups %>% select(alfa_id, decline_group),
    by = "alfa_id"
  )

Four_group_trajectory_plot <- ggplot(plot_time2,
                                     aes(time_years, iadl_score,
                                         group = alfa_id,
                                         color = decline_group)) +
  
  geom_line(alpha = 0.05) +
  geom_point(alpha = 0.1) +
  
  stat_smooth(
    aes(group = decline_group),
    method = "loess",
    span = 0.9,
    se = FALSE,
    linewidth = 1.3
  ) +
  
  scale_color_manual(values = c(
    "Non-decliner" = "black",
    "Decliner V2 only" = "#E69F00",
    "Decliner V3 only" = "#56B4E9",
    "Decliner V2 & V3" = "#D55E00"
  )) +
  
  labs(
    x = "Years from baseline",
    y = "A-IADL-Q score",
    color = "Group"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(face = "bold")
  )

Four_group_trajectory_plot

# ============================================================================= #
# 9. DECLINER VS NO-DECLINER COMPARISON (-2.2 POINTS FROM V3-V1) ----
# ============================================================================= #

## ---- 9.1 SUCOG overall cohort check ----
demo_bl_by_AT_SUCOG <- baseline_tbl %>%
  select(
    AT_bl,
    Age_V1,
    YearsEducation,
    sex,
    iadl_V1,
    iadl_V2,
    iadl_V3,
    n_visits,
    fu_time_years,
    sucog_cogprob
  ) %>%
  tbl_summary(
    by = AT_bl,
    missing = "no",
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous() ~ 1,
      fu_time_years ~ 2
    ),
    label = list(
      #AT_bl        ~ "AT group (baseline)",
      Age_V1       ~ "Age at baseline (years)",
      YearsEducation ~ "Education (years)",
      sex          ~ "Sex",
      iadl_V1      ~ "A-IADL-Q score at V1 (baseline)",
      iadl_V2      ~ "A-IADL-Q score at V2",
      iadl_V3      ~ "A-IADL-Q score at V3",
      n_visits     ~ "Number of visits with A-IADL-Q",
      fu_time_years ~ "Follow-up time (years)"
    )
  ) %>%
  add_overall() %>%
  add_p(test = list(
    all_continuous()  ~ "kruskal.test",
    all_categorical() ~ "chisq.test"
  )) %>%
  add_n()

demo_bl_by_AT_SUCOG

wide_SUCOG <- alfa_long %>%
  tidyr::pivot_wider(
    id_cols   = alfa_id,
    names_from = visit,
    values_from = sucog_cogprob
  ) %>%
  select(V1, V2, V3)
all_patterns <- c("0,0,0","0,0,1","0,1,0","0,1,1",
                  "1,0,0","1,0,1","1,1,0","1,1,1")
wide_SUCOG$pattern <- paste(wide_SUCOG$V1, wide_SUCOG$V2, wide_SUCOG$V3, sep = ",")
table(wide_SUCOG$pattern)

wide_filtered <- wide_SUCOG %>%
  dplyr::filter(rowSums(is.na(dplyr::across(c(V1, V2, V3)))) <= 1)

table_full <- table(factor(wide_filtered$pattern, levels = all_patterns))
table_full

# Table with all patterns (including NA patterns)
pattern_all <- as.data.frame(table(wide_filtered$pattern)) %>%
  rename(pattern = Var1, n = Freq)

# Table with only valid complete patterns
pattern_complete <- as.data.frame(table_full) %>%
  rename(pattern = Var1, n = Freq)

pattern_complete %>%
  mutate(percent = scales::percent(n / sum(n))) %>%
  gt() %>%
  tab_header(
    title = "SCD Patterns Across Visits",
    subtitle = "Complete cases only"
  ) %>%
  cols_label(
    pattern = "SCD Pattern (V1,V2,V3)",
    n = "N",
    percent = "%"
  )
pattern_all %>%
  gt() %>%
  tab_header(
    title = "All Observed SCD Patterns",
    subtitle = "Including incomplete visit data"
  ) %>%
  cols_label(
    pattern = "SCD Pattern (V1,V2,V3)",
    n = "N"
  )

## ---- 9.2 Correlation with SCD ----
wide <- alfa_long %>%
  select(alfa_id, visit, sucog_cogprob, iadl_score, ME_Total_SUCOG, LE_Total_SUCOG, FE_Total_SUCOG, Total_SUCOG) %>%
  pivot_wider(
    id_cols = alfa_id,
    names_from = visit,
    values_from = c(sucog_cogprob, iadl_score, ME_Total_SUCOG, LE_Total_SUCOG, FE_Total_SUCOG, Total_SUCOG),
    names_sep = "_"
  )

# Only using 3 visits SCD decline

wide <- wide %>%
  mutate(
    pattern = paste(sucog_cogprob_V1, sucog_cogprob_V2, sucog_cogprob_V3, sep = ","),
    
    SUCOG_group = case_when(
      pattern == "0,0,0" ~ "SCD-",
      pattern == "1,1,1" ~ "SCD+",
      pattern %in% c("0,0,1", "0,1,1") ~ "Decliner",
      TRUE ~ NA_character_
    )
  )

# Using 2 and 3 visits SCD decline
wide <- wide %>%
  mutate(
    pattern = paste(sucog_cogprob_V1, sucog_cogprob_V2, sucog_cogprob_V3, sep = ","),
    
    SUCOG_group = case_when(
      pattern %in% c("0,0,0", "0,0,NA", "0,NA,0") ~ "SCD-",
      pattern %in% c("1,1,1", "1,1,NA", "1,NA,1") ~ "SCD+",
      pattern %in% c("0,0,1", "0,1,1", "0,1,NA", "0,NA,1") ~ "Decliner",
      TRUE ~ NA_character_
    )
  )

# AIDL decline
wide <- wide %>%
  mutate(
    iadl_diff = iadl_score_V1 - iadl_score_V3,
    
    AIDL_decline = case_when(
      !is.na(iadl_diff) & iadl_diff >= 22 ~ "Decliner",
      !is.na(iadl_diff) ~ "Non decliner",
      TRUE ~ NA_character_
    )
  )

# Keep patients with non NA results (SUCOG 1/0)
analysis_df <- wide %>%
  filter(!is.na(SUCOG_group), !is.na(AIDL_decline))

# Keep patients with non NA results (SUCOG continuos)
analysis_df <- wide %>%
  filter(!is.na(AIDL_decline))

# Create correlation table
summary_table <- analysis_df %>%
  count(SUCOG_group, AIDL_decline) %>%
  group_by(SUCOG_group) %>%
  mutate(
    percent = n / sum(n),
    label = paste0(n, " (", percent(percent, accuracy = 0.1), ")")
  ) %>%
  ungroup() %>%
  select(SUCOG_group, AIDL_decline, label) %>%
  tidyr::pivot_wider(
    names_from = AIDL_decline,
    values_from = label,
    values_fill = "0 (0%)"
  )

totals <- analysis_df %>%
  count(SUCOG_group, name = "Total")

summary_table <- summary_table %>%
  left_join(totals, by = "SUCOG_group")

summary_table %>%
  gt() %>%
  tab_header(
    title = "Association Between SCD Status and IADL Functional Decline"
    ) %>%
  cols_label(
    SUCOG_group = "SCD Group",
    Total = "Total N"
  ) %>%
  tab_spanner(
    label = "IADL Functional Status n (%)*",
    columns = c(`Decliner`, `Non decliner`)
  ) %>%
  cols_align(
    align = "center",
    -SUCOG_group
  ) %>%
  tab_source_note(
    source_note = md("*AIDL decline defined as V1–V3 difference > 2.2 points.")
  )

# Create summary mean table

library(dplyr)
library(gt)

# Select SUCOG variables for all visits
sucog_vars <- c(
  "Total_SUCOG_V1", "Total_SUCOG_V2", "Total_SUCOG_V3",
  "ME_Total_SUCOG_V1", "ME_Total_SUCOG_V2", "ME_Total_SUCOG_V3",
  "LE_Total_SUCOG_V1", "LE_Total_SUCOG_V2", "LE_Total_SUCOG_V3",
  "FE_Total_SUCOG_V1", "FE_Total_SUCOG_V2", "FE_Total_SUCOG_V3"
)

# Compute mean scores per IADL group
summary_scores_visits <- analysis_df %>%
  group_by(AIDL_decline) %>%
  summarise(across(all_of(sucog_vars), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

# Optional: pivot longer for a cleaner GT table with Visit columns
summary_scores_long <- summary_scores_visits %>%
  pivot_longer(
    cols = -AIDL_decline,
    names_to = c("Domain", "Visit"),
    names_pattern = "(.*)_V(\\d)"
  ) %>%
  pivot_wider(
    names_from = Visit,
    values_from = value,
    names_prefix = "V"
  )

# Create GT table
summary_scores_long %>%
  gt() %>%
  tab_header(
    title = "Mean SUCOG Scores by IADL Functional Decline and Visit"
  ) %>%
  cols_label(
    AIDL_decline = "IADL Group",
    Domain = "SUCOG Domain",
    V1 = "Visit 1",
    V2 = "Visit 2",
    V3 = "Visit 3"
  ) %>%
  cols_align(
    align = "center",
    everything()
  ) %>%
  tab_source_note(
    source_note = md("Mean SUCOG scores per visit, stratified by IADL decline (-2.2 threshold) group.")
  )

#PLOT teh differences

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)

# Prepare data: long format for plotting
plot_data <- analysis_df %>%
  select(
    AIDL_decline,
    Total_SUCOG_V1, Total_SUCOG_V2, Total_SUCOG_V3,
    ME_Total_SUCOG_V1, ME_Total_SUCOG_V2, ME_Total_SUCOG_V3,
    LE_Total_SUCOG_V1, LE_Total_SUCOG_V2, LE_Total_SUCOG_V3,
    FE_Total_SUCOG_V1, FE_Total_SUCOG_V2, FE_Total_SUCOG_V3
  ) %>%
  pivot_longer(
    cols = -AIDL_decline,
    names_to = c("Domain", "Visit"),
    names_pattern = "(.*)_V(\\d)"
  ) %>%
  mutate(
    Visit = as.integer(Visit)
  ) %>%
  group_by(AIDL_decline, Domain, Visit) %>%
  summarise(
    mean_score = mean(value, na.rm = TRUE),
    sd_score = sd(value, na.rm = TRUE),
    n = n(),
    se = sd_score / sqrt(n)  # standard error for error bars
  ) %>%
  ungroup()

# Plot: line plot with error bars
ggplot(plot_data, aes(x = Visit, y = mean_score, color = AIDL_decline, group = AIDL_decline)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se), width = 0.2) +
  facet_wrap(~ Domain, scales = "free_y") +  # one plot per SUCOG domain
  scale_x_continuous(breaks = 1:3, labels = paste0("V", 1:3)) +
  scale_color_manual(values = c("Decliner" = "#D55E00", "Non decliner" = "#0072B2")) +
  labs(
    title = "Trajectory of SUCOG Scores by IADL Decline Group",
    subtitle = "Mean scores with standard error across three visits",
    x = "Visit",
    y = "Mean SUCOG Score",
    color = "IADL Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12)
  )

## ---- 9.3 Datset creation for comparison table ----

# 1) Baseline table and id_table as reference
id_tbl <- alfa_study %>%
  distinct(alfa_id)

baseline_tbl <- alfa_study %>%
  filter(visit == "V1") %>%
  distinct(alfa_id, .keep_all = TRUE) %>%
  select(
    alfa_id,
    Age_V1,
    sex,
    YearsEducation,
    APOE_binary,
    AT_bl
  )


# 2) Follow-up summary (n_visits & fu_time_years)
fu_summary <- alfa_study %>%
  filter(visit %in% c("V1", "V2", "V3")) %>%
  group_by(alfa_id) %>%
  summarise(
    # Baseline date (V1)
    v1_date = min(date[visit == "V1"], na.rm = TRUE),
    
    # Number of visits with non-missing IADL across V1–V3
    n_visits = n_distinct(visit[visit %in% c("V1", "V2", "V3") & !is.na(iadl_score)]),
    
    # Last follow-up date with IADL (V2 or V3)
    fu_last_date = max(date[visit %in% c("V2", "V3") & !is.na(iadl_score)], na.rm = TRUE),
    fu_last_date = if_else(is.infinite(fu_last_date), as.Date(NA), fu_last_date),
    
    # Follow-up time in years (V1 → last FU)
    fu_time_years = if_else(
      !is.na(fu_last_date),
      as.numeric(fu_last_date - v1_date) / 365.25,
      0
    ),
    .groups = "drop"
  ) %>%
  select(alfa_id, n_visits, fu_time_years)

fu_tbl <- fu_summary


# 3) IADL score change and threshold division
iadl_tbl <- alfa_study %>%
  select(alfa_id, visit, iadl_score) %>%
  pivot_wider(
    names_from = visit,
    values_from = iadl_score,
    names_prefix = "iadl_"
  ) %>%
  mutate(
    iadl_change = iadl_V3 - iadl_V1,
    iadl_decline = case_when(
      !is.na(iadl_change) & iadl_change <= -22 ~ "Decliner",
      !is.na(iadl_change) ~ "Non-decliner",
      TRUE ~ NA_character_
    ),
    iadl_decline = factor(iadl_decline,
                          levels = c("Non-decliner","Decliner"))
  ) %>%
  select(alfa_id, iadl_decline, iadl_V1, iadl_change)


# 4) Biomarkers in wide format and extract V1
biomarkers <- c(
  "BIO_Abeta1.40_CSF_Elecsys_UGOT_Roche",
  "BIO_Abeta1.40_PLASMA_Elecsys_UGOT_Roche",
  "BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche",
  "BIO_Abeta1.42_PLASMA_Elecsys_UGOT_Roche",
  "BIO_GFAP_CSF_Elecsys_UGOT_Roche",
  "BIO_GFAP_PLASMA_Elecsys_UGOT_Roche",
  "BIO_NFL_CSF_Elecsys_UGOT_Roche",
  "BIO_NFL_PLASMA_Elecsys_UGOT_Roche",
  "BIO_pTau181_CSF_Elecsys_UGOT_Roche",
  "BIO_pTau181_PLASMA_Elecsys_UGOT_Roche",
  "BIO_tTau_CSF_Elecsys_UGOT_Roche",
  "BIO_pTau217_PLASMA_Elecsys_UGOT_Roche",
  "BIO_pTau217_PLASMA_MSD_LILLY_in.house",
  "BIO_pTau217_CSF_MSD_LILLY_in.house",
  "BIO_pTau231_CSF_ELISA_UGOT_ADx.NeuroSciences"
)

biomarker_tbl <- alfa_study %>%
  filter(visit == "V1") %>%
  select(alfa_id, all_of(biomarkers))


# 5) SUCOG V1 and yearly change extraction
sucog_tbl <- alfa_study %>%
  mutate(
    across(
      c(sucog_decr2yrs, sucog_cogprob, Total_SUCOG, ME_Total_SUCOG, FE_Total_SUCOG, LE_Total_SUCOG),
      ~ as.numeric(as.character(.))
    )
  ) %>%
  select(
    alfa_id, visit, sucog_decr2yrs, sucog_cogprob,
    Total_SUCOG, ME_Total_SUCOG,
    FE_Total_SUCOG, LE_Total_SUCOG,
    FechaResultado_SUCOG
  ) %>%
  pivot_wider(
    names_from = visit,
    values_from = c(
      sucog_decr2yrs, sucog_cogprob, Total_SUCOG, ME_Total_SUCOG,
      FE_Total_SUCOG, LE_Total_SUCOG,
      FechaResultado_SUCOG
    ),
    names_glue = "{.value}_{visit}"
  ) %>%
  mutate(
    sucog_tot_slope =
      (Total_SUCOG_V3 - Total_SUCOG_V1) /
      as.numeric(FechaResultado_SUCOG_V3 - FechaResultado_SUCOG_V1) * 365.25,
    
    sucog_me_slope =
      (ME_Total_SUCOG_V3 - ME_Total_SUCOG_V1) /
      as.numeric(FechaResultado_SUCOG_V3 - FechaResultado_SUCOG_V1) * 365.25,
    
    sucog_fe_slope =
      (FE_Total_SUCOG_V3 - FE_Total_SUCOG_V1) /
      as.numeric(FechaResultado_SUCOG_V3 - FechaResultado_SUCOG_V1) * 365.25,
    
    sucog_le_slope =
      (LE_Total_SUCOG_V3 - LE_Total_SUCOG_V1) /
      as.numeric(FechaResultado_SUCOG_V3 - FechaResultado_SUCOG_V1) * 365.25
  ) %>%
  select(
    alfa_id,
    sucog_decr2yrs_V1,
    sucog_decr2yrs_V2,
    sucog_decr2yrs_V3,
    sucog_cogprob_V1,
    sucog_cogprob_V2,
    sucog_cogprob_V3,
    Total_SUCOG_V1,
    ME_Total_SUCOG_V1,
    FE_Total_SUCOG_V1,
    LE_Total_SUCOG_V1,
    sucog_tot_slope,
    sucog_me_slope,
    sucog_fe_slope,
    sucog_le_slope
  )

# 6) PACC yearly change
pacc_tbl <- alfa_study %>%
  select(alfa_id, visit, PACC, FechaResultado, PACC_attention, PACC_memory, PACC_executive, PACC_language, PACC_visuospatial) %>%
  pivot_wider(
    names_from = visit,
    values_from = c(PACC,PACC_attention, PACC_memory, PACC_executive, PACC_language, PACC_visuospatial, FechaResultado),
    names_glue = "{.value}_{visit}"
  ) %>%
  mutate(pacc_slope = (PACC_V3 - PACC_V1) / 
           as.numeric(FechaResultado_V3 - FechaResultado_V1) * 365.25,
         pacc_attention_slope = (PACC_attention_V3 - PACC_attention_V1) / 
           as.numeric(FechaResultado_V3 - FechaResultado_V1) * 365.25,
         pacc_memory_slope = (PACC_memory_V3 - PACC_memory_V1) / 
           as.numeric(FechaResultado_V3 - FechaResultado_V1) * 365.25,
         pacc_executive_slope = (PACC_executive_V3 - PACC_executive_V1) / 
           as.numeric(FechaResultado_V3 - FechaResultado_V1) * 365.25,
         pacc_language_slope = (PACC_language_V3 - PACC_language_V1) / 
           as.numeric(FechaResultado_V3 - FechaResultado_V1) * 365.25,
         pacc_visuospatial_slope = (PACC_visuospatial_V3 - PACC_visuospatial_V1) / 
           as.numeric(FechaResultado_V3 - FechaResultado_V1) * 365.25
  ) %>%
  select(alfa_id, 
         PACC_V1, 
         PACC_attention_V1,
         PACC_memory_V1,
         PACC_executive_V1, 
         PACC_language_V1, 
         PACC_visuospatial_V1, 
         pacc_slope,
         pacc_attention_slope,
         pacc_memory_slope,
         pacc_executive_slope, 
         pacc_language_slope, 
         pacc_visuospatial_slope)


# 7) Frailty Index & SCD
frailtyscd_tbl <- alfa_study %>%
  select(alfa_id, visit, SCD_def, Frailty_Index) %>%
  pivot_wider(
    names_from = visit,
    values_from = c(SCD_def, Frailty_Index),
    names_glue = "{.value}_{visit}"
  ) %>%
  select(alfa_id, Frailty_Index_V1, SCD_def_V1, SCD_def_V2, SCD_def_V3)


# 8) Diagnostic
diagnostic_tbl <- alfa_study %>%
  select(alfa_id, visit, diagnostic) %>%
  pivot_wider(
    names_from = visit,
    values_from = c(diagnostic),
    names_glue = "{.value}_{visit}"
  ) %>%
  select(alfa_id, diagnostic_V1, diagnostic_V2, diagnostic_V3)


# 8) TMTB
TMT_tbl <- alfa_study %>%
  select(alfa_id, visit, TMTB_completion, TMTB_time) %>%
  pivot_wider(
    names_from = visit,
    values_from = c(TMTB_completion, TMTB_time),
    names_glue = "{.value}_{visit}"
  ) %>%
  select(alfa_id,
         TMTB_completion_V1,
         TMTB_time_V1,
         TMTB_completion_V2,
         TMTB_time_V2,
         TMTB_completion_V3,
         TMTB_time_V3)

# 9) Build final dataset
decline_tbl <- id_tbl %>%
  left_join(baseline_tbl, by="alfa_id") %>%
  left_join(iadl_tbl, by="alfa_id") %>%
  left_join(biomarker_tbl, by="alfa_id") %>%
  left_join(sucog_tbl, by="alfa_id") %>%
  left_join(pacc_tbl, by="alfa_id") %>%
  left_join(fu_tbl, by="alfa_id") %>%
  left_join(frailtyscd_tbl, by="alfa_id") %>%
  left_join(diagnostic_tbl, by="alfa_id") %>%
  left_join(TMT_tbl, by="alfa_id") %>%
  mutate(
    sex = factor(sex, levels=c("Men","Women"), labels = c("Men", "Women")),
    APOEe4 = factor(APOE_binary, levels=c(0,1),
                    labels=c("Non-carrier","Carrier")),
    AT_bl = droplevels(AT_bl),
    diagnostic_V1 = factor(diagnostic_V1, levels=c("1","2", "3", "4"), labels = c("Cognitively Unimpaired", "SCD", "MCI", "Dementia")),
    diagnostic_V3 = factor(diagnostic_V3, levels=c("1","2", "3", "4"), labels = c("Cognitively Unimpaired", "SCD", "MCI", "Dementia"))
  )

alfa_study_decliners <- alfa_study %>%
  left_join(iadl_tbl, by="alfa_id")
out_file_rdata <- file.path("C:/Users/U272674/Desktop", "full_alfa_study_decliners.RData")
save(alfa_study_decliners, file = out_file_rdata)
cat("Saved RData to:\n  ", out_file_rdata, "\n", sep="")

## ---- 9.4 Create final comparison table

# Main table
iadl_decline_table <-
  decline_tbl %>%
  select(
    iadl_decline,
    Age_V1,
    sex,
    YearsEducation,
    APOEe4,
    AT_bl,
    Frailty_Index_V1,
    
    # biomarkers
    starts_with("BIO"),
    
    # baseline measures
    sucog_decr2yrs_V1,
    sucog_cogprob_V1,
    SCD_def_V1,
    diagnostic_V1,
    Total_SUCOG_V1,
    ME_Total_SUCOG_V1,
    FE_Total_SUCOG_V1,
    LE_Total_SUCOG_V1,
    PACC_V1,
    PACC_attention_V1,
    PACC_memory_V1,
    PACC_executive_V1,
    PACC_language_V1,
    PACC_visuospatial_V1,
    TMTB_completion_V1,
    TMTB_time_V1,
    
    # Evolution measures
    sucog_decr2yrs_V2,
    sucog_decr2yrs_V3,
    sucog_cogprob_V2,
    sucog_cogprob_V3,
    SCD_def_V2,
    SCD_def_V3,
    diagnostic_V2,
    diagnostic_V3,
    TMTB_completion_V2,
    TMTB_completion_V3,
    TMTB_time_V2,
    TMTB_time_V3,
    
    # slopes
    sucog_tot_slope,
    sucog_me_slope,
    sucog_fe_slope,
    sucog_le_slope,
    pacc_slope,
    pacc_attention_slope,
    pacc_memory_slope,
    pacc_executive_slope, 
    pacc_language_slope, 
    pacc_visuospatial_slope,
    
    # IADL + FU
    iadl_V1,
    fu_time_years
  ) %>%
  
  tbl_summary(
    by = iadl_decline,
    missing = "no",
    
    type = list(
      Total_SUCOG_V1 ~ "continuous",
      ME_Total_SUCOG_V1 ~ "continuous",
      FE_Total_SUCOG_V1 ~ "continuous",
      LE_Total_SUCOG_V1 ~ "continuous",
      sucog_decr2yrs_V1 ~ "categorical",
      sucog_decr2yrs_V2 ~ "categorical",
      sucog_decr2yrs_V3 ~ "categorical",
      sucog_cogprob_V1 ~ "categorical",
      sucog_cogprob_V2 ~ "categorical",
      sucog_cogprob_V3 ~ "categorical",
      SCD_def_V1 ~ "categorical",
      SCD_def_V2 ~ "categorical",
      SCD_def_V3 ~ "categorical",
      TMTB_completion_V2 ~ "categorical",
      TMTB_completion_V3 ~ "categorical",
      diagnostic_V1 ~ "categorical",
      diagnostic_V2 ~ "categorical",
      diagnostic_V3 ~ "categorical",
      sucog_tot_slope ~ "continuous",
      sucog_me_slope ~ "continuous",
      sucog_fe_slope ~ "continuous",
      sucog_le_slope ~ "continuous",
      pacc_slope ~ "continuous",
      pacc_attention_slope ~ "continuous",
      pacc_memory_slope ~ "continuous", 
      pacc_executive_slope ~ "continuous", 
      pacc_language_slope ~ "continuous", 
      pacc_visuospatial_slope ~ "continuous"
    ),
    
    label = list(
      AT_bl ~ "AT group (baseline)",
      Age_V1 ~ "Age at baseline (years)",
      YearsEducation ~ "Education (years)",
      sex ~ "Sex",
      Frailty_Index_V1 ~ "Frailty Index",
      
      BIO_Abeta1.40_CSF_Elecsys_UGOT_Roche ~ "CSF Amyloid β40",
      BIO_Abeta1.40_PLASMA_Elecsys_UGOT_Roche ~ "Plasma Amyloid β40",
      BIO_Abeta1.42_CSF_Elecsys_UGOT_Roche ~ "CSF Amyloid β42",
      BIO_Abeta1.42_PLASMA_Elecsys_UGOT_Roche ~ "Plasma Amyloid β42",
      BIO_GFAP_CSF_Elecsys_UGOT_Roche ~ "CSF GFAP",
      BIO_GFAP_PLASMA_Elecsys_UGOT_Roche ~ "Plasma GFAP",
      BIO_NFL_CSF_Elecsys_UGOT_Roche ~ "CSF NfL",
      BIO_NFL_PLASMA_Elecsys_UGOT_Roche ~ "Plasma NfL",
      BIO_pTau181_CSF_Elecsys_UGOT_Roche ~ "CSF pTau181",
      BIO_pTau181_PLASMA_Elecsys_UGOT_Roche ~ "Plasma pTau181",
      BIO_tTau_CSF_Elecsys_UGOT_Roche ~ "CSF tTau",
      BIO_pTau217_PLASMA_Elecsys_UGOT_Roche ~ "Plasma pTau217 (Roche)",
      BIO_pTau217_PLASMA_MSD_LILLY_in.house ~ "Plasma pTau217 (Lilly)",
      BIO_pTau217_CSF_MSD_LILLY_in.house ~ "CSF pTau217",
      BIO_pTau231_CSF_ELISA_UGOT_ADx.NeuroSciences ~ "CSF pTau231",
      
      sucog_decr2yrs_V1 ~ "SUCOG 2-yr decline V1 (baseline)",
      sucog_decr2yrs_V2 ~ "SUCOG 2-yr decline V2",
      sucog_decr2yrs_V3 ~ "SUCOG 2-yr decline V3",
      SCD_def_V1 ~ "SCD patient V1 (baseline)",
      SCD_def_V2 ~ "SCD patient V2",
      SCD_def_V3 ~ "SCD patient V3",
      diagnostic_V1 ~ "Diagnostic (baseline)",
      diagnostic_V2 ~ "Diagnostic V2",
      diagnostic_V3 ~ "Diagnostic V3",
      Total_SUCOG_V1 ~ "SUCOG Total Score V1 (baseline)",
      ME_Total_SUCOG_V1 ~ "SUCOG Memory Score V1 (baseline)",
      FE_Total_SUCOG_V1 ~ "SUCOG Executive Score V1 (baseline)",
      LE_Total_SUCOG_V1 ~ "SUCOG Language Score V1 (baseline)",
      TMTB_completion_V1 ~ "TMTB completion (baseline)",
      TMTB_completion_V2 ~ "TMTB completion V2",
      TMTB_completion_V3 ~ "TMTB completion V3",
      TMTB_time_V1 ~ "TMTB time employed (baseline)",
      TMTB_time_V2 ~ "TMTB time employed V2",
      TMTB_time_V3 ~ "TMTB time employed V3",
      
      
      sucog_tot_slope ~ "Yearly change in SUCOG (Total)",
      sucog_me_slope ~ "Yearly change in SUCOG (Memory)",
      sucog_fe_slope ~ "Yearly change in SUCOG (Executive)",
      sucog_le_slope ~ "Yearly change in SUCOG (Language)",
      pacc_slope ~ "Yearly change in PACC",
      pacc_attention_slope ~ "Yearly change in PACC (Attention)",
      pacc_memory_slope ~ "Yearly change in PACC (Memory)",
      pacc_executive_slope ~ "Yearly change in PACC (Executive)", 
      pacc_language_slope ~ "Yearly change in PACC (Language)", 
      pacc_visuospatial_slope ~ "Yearly change in PACC (Visuospatial)",
      
      iadl_V1 ~ "A-IADL-Q score at baseline",
      fu_time_years ~ "Follow-up time (years)"
    ),
    
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    )
  ) %>%
  
  add_p() %>%
  bold_p(t = 0.05) %>%

# Define subsections of the table
modify_table_body(
  ~ .x %>%
    mutate(
      group = case_when(
        variable %in% c("Age_V1","sex","YearsEducation","APOEe4","AT_bl", "Frailty_Index")
        ~ "Demographics",
        
        grepl("^BIO", variable)
        ~ "Biomarkers",
        
        variable %in% c("sucog_decr2yrs_V1","sucog_cogprob_V1","SCD_def_V1",
                        "diagnostic_V1", "Total_SUCOG_V1","ME_Total_SUCOG_V1","FE_Total_SUCOG_V1",
                        "LE_Total_SUCOG_V1", "PACC_V1", "PACC_attention_V1", 
                        "PACC_memory_V1", "PACC_executive_V1", "PACC_language_V1", 
                        "PACC_visuospatial_V1", "TMTB_time_V1", "TMTB_completion_V1")
        ~ "Baseline cognition",
        
        variable %in% c("sucog_decr2yrs_V2", "sucog_decr2yrs_V3", "sucog_cogprob_V2",
                        "sucog_cogprob_V3", "SCD_def_V2", "SCD_def_V3", "diagnostic_V2",
                        "diagnostic_V3", "TMTB_time_V2", "TMTB_time_V3", 
                        "TMTB_completion_V2", "TMTB_completion_V3")
        ~ "Evolution cognition",
        
        grepl("slope", variable)
        ~ "Cognitive decline slopes",
        
        variable %in% c("iadl_V1","fu_time_years")
        ~ "Functional measures",
        
        TRUE ~ NA_character_
      )
    )
) %>%
  
  as_gt() %>%

# Improve Style
  
# zebra rows
gt::opt_row_striping() %>%
  
  # make section labels appear
  gt::tab_row_group(
    label = "Demographics",
    rows = group == "Demographics"
  ) %>%
  gt::tab_row_group(
    label = "Biomarkers",
    rows = group == "Biomarkers"
  ) %>%
  gt::tab_row_group(
    label = "Baseline cognition",
    rows = group == "Baseline cognition"
  ) %>%
  gt::tab_row_group(
    label = "Evolution cognition",
    rows = group == "Evolution cognition"
  ) %>%
  gt::tab_row_group(
    label = "Cognitive decline slopes",
    rows = group == "Cognitive decline slopes"
  ) %>%
  gt::tab_row_group(
    label = "Functional measures",
    rows = group == "Functional measures"
  ) %>%
  
  # bold only section headers
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_row_groups()
  ) %>%
  
  # nicer title
  gt::tab_header(
    title = md("**Baseline characteristics by IADL decline status**"),
    subtitle = "Mean (SD) or N (%)"
  )

iadl_decline_table

## ---- 9.4 Save decliners comparison table
gtsave(iadl_decline_table, "IADL_decline_comparison_20260505.png")
