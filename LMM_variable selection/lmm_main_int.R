
#18 variables analysis
complete_lmm = df_UCPCR_clean %>% select(UCPCR, Study_ID, dur_diab_months_z,Gender_v1, AgeatDiagnosis_z, num_anti,
                                         DKA, bmi_calc_v1_z,
                                         autoimmune,GAD_v1_z,IA2_v1_z,ZNT8_v1_z,
                                         wh_ratio_v1_z,GAD_bin_v1,
                                         IA2_bin_v1,ZNT8_bin_v1,HbA1c_at_diagnosis_v1_z,
                                         Smoker_v1,famhisinsdiab,
                                         T1DGRS2_z,T2DGRS_z) %>% drop_na()
# Variables to test
vars = c('Gender_v1', 'AgeatDiagnosis_z', 'num_anti', 'DKA', 'bmi_calc_v1_z',
         'autoimmune','GAD_v1_z','IA2_v1_z','ZNT8_v1_z',
         'wh_ratio_v1_z','GAD_bin_v1',
         'IA2_bin_v1','ZNT8_bin_v1','HbA1c_at_diagnosis_v1_z',
         'Smoker_v1','famhisinsdiab',
         'T1DGRS2_z','T2DGRS_z')

#loop function
lmm_results = map_df(vars, function(v){
  
  form <- as.formula(
    paste0("scale(log(UCPCR)) ~ dur_diab_months_z * ", v, " + (dur_diab_months_z || Study_ID)")
  )
  
  model <- lmer(form, data = complete_lmm)
  
  tidy(model, effects = "fixed") %>%
    mutate(variable = v)
})
lmm_results


lmm_table <- lmm_results %>%
  filter(term %in% vars | grepl("dur_diab_months_z:", term)) %>%
  mutate(effect_type = case_when(
    term %in% vars ~ "Main effect",
    TRUE ~ "Interaction with time"
  )) %>%
  select(variable, effect_type, estimate, std.error, p.value)


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


#1 Gender
l1 = lmer(scale(log(UCPCR)) ~ dur_diab_months_z +
            (dur_diab_months_z || Study_ID), 
          data = complete_lmm)
gen_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + Gender_v1 + 
                 dur_diab_months_z * Gender_v1 + 
                 (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(l1,gen_lmm)


#2
age_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + AgeatDiagnosis_z + 
                 dur_diab_months_z * AgeatDiagnosis_z + 
                 (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(l1,age_lmm)


#3 DKA
dka_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + DKA + 
                 dur_diab_months_z * DKA + 
                 (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(l1,dka_lmm)


#4 BMI
#completes cases for LRT Comparison
bmi_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + bmi_calc_v1_z + 
                 dur_diab_months_z * bmi_calc_v1_z + 
                 (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(l1,bmi_lmm)

  

#5 ZNT8 STATUS
zn_stat_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + ZNT8_bin_v1 + 
                     dur_diab_months_z * ZNT8_bin_v1 + 
                     (dur_diab_months_z || Study_ID), data =complete_lmm)
anova(l1,zn_stat_lmm)

#6 ZNT8 TITRE
zn_titre_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + ZNT8_v1_z + 
                     dur_diab_months_z * ZNT8_v1_z + 
                     (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(l1,zn_titre_lmm)


#7 GAD STATUS
gadS_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + GAD_bin_v1 + 
                  dur_diab_months_z * GAD_bin_v1 + 
                  (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(l1,gadS_lmm)

#8 GAD Titre
gadT_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + GAD_v1_z + 
                  dur_diab_months_z * GAD_v1_z + 
                  (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(l1,gadT_lmm)



#9 IA2 Status
ia2S_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + IA2_bin_v1 + 
                  dur_diab_months_z * IA2_bin_v1 + 
                  (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(l1,ia2S_lmm)

  
#10 IA2 Titre
ia2T_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + IA2_v1_z + 
                  dur_diab_months_z * IA2_v1_z + 
                  (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(l1,ia2T_lmm)


#11 History of additional auto-immune disease
auto_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + autoimmune + 
                  dur_diab_months_z * autoimmune + 
                  (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(l1,auto_lmm)

#12 Number of antibodies
num_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + num_anti + 
                  dur_diab_months_z * num_anti + 
                  (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(l1,num_lmm)

#13 waist hip ratio
whr_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + wh_ratio_v1_z + 
                 dur_diab_months_z * wh_ratio_v1_z + 
                 (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(l1,whr_lmm)


#14 T1DGRS2
t1d_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + T1DGRS2_z + 
                dur_diab_months_z * T1DGRS2_z + 
                (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(l1,t1d_lmm)

  

#15 Famiy history of insulin treated diabetes
famhis_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + famhisinsdiab + 
                 dur_diab_months_z * famhisinsdiab + 
                 (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(l1,famhis_lmm)
  

#16 HBA1C
 hb_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + HbA1c_at_diagnosis_v1_z + 
                    dur_diab_months_z * HbA1c_at_diagnosis_v1_z + 
                    (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(l1,hb_lmm)



