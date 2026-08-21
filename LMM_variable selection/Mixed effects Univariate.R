#Mixed model univariate analysis

library(lme4)
library(lmerTest)
library(broom.mixed)
library(dplyr)
library(purrr)
library(gt)
library(gtsummary)
library(ggeffects) 
library(performance) 
library(marginaleffects)



# Variables to test
vars = c('Gender_v1', 'AgeatDiagnosis_z', 'num_anti', 'DKA', 'bmi_calc_v1_z',
         'autoimmune','GAD_v1_z','IA2_v1_z','ZNT8_v1_z',
         'wh_ratio_v1_z','GAD_bin_v1',
         'IA2_bin_v1','ZNT8_bin_v1','HbA1c_at_diagnosis_v1_z',
         'famhisinsdiab',
         'T1DGRS2_z')

#loop function
lmm_results = map_df(vars, function(v){
  
  form <- as.formula(
    paste0("scale(log(UCPCR)) ~ dur_diab_months_z * ", v, " + (dur_diab_months_z || Study_ID)")
  )
  
  model <- lmer(form, data = df_UCPCR_clean)
  
  tidy(model, effects = "fixed") %>%
    mutate(variable = v)
})
lmm_results


lmm_table = lmm_results %>%
  filter(term != "(Intercept)" &
           (grepl("^dur_diab_months_z:", term) |
              !grepl("^dur_diab_months_z$", term))) %>%
  mutate(
    effect_type = if_else(
      grepl("^dur_diab_months_z:", term),
      "Interaction with time",
      "Main effect"
    ),
    category = term
  ) %>%
  select(variable, category, effect_type,
         estimate, std.error, p.value)


lmm_table %>%
  mutate(
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    p.value = signif(p.value, 3)
  ) %>%
  arrange(p.value) %>%
  gt() %>%
  tab_header(title = "Bivariate Linear Mixed Effects Models") %>%
  gt::gtsave("Bivariate_LMM_grouped.png")


