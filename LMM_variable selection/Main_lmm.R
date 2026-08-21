#LMM with only main effects
library(lme4test)

#1 Gender
gen_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + Gender_v1 +
                  (dur_diab_months_z || Study_ID), 
                data = df_UCPCR_clean) %>% summary()

#2
age_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z +AgeatDiagnosis_z +
                  (dur_diab_months_z || Study_ID), 
                data = df_UCPCR_clean)  %>% summary()

#3 DKA
dka_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z +DKA +
                  (dur_diab_months_z || Study_ID), 
                data = df_UCPCR_clean) %>% summary()

#4 BMI
bmi_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + bmi_calc_v1_z +
                  (dur_diab_months_z || Study_ID), 
                data = df_UCPCR_clean) %>% summary()


#5 ZNT8 STATUS
zn_stat_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + ZNT8_bin_v1 +
                      (dur_diab_months_z || Study_ID), 
                    data = df_UCPCR_clean) %>% summary()

#6 ZNT8 TITRE
zn_titre_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + ZNT8_v1_z +
                       (dur_diab_months_z || Study_ID), 
                     data = df_UCPCR_clean) %>% summary()

#7 GAD STATUS
gadS_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + GAD_bin_v1 +
                   (dur_diab_months_z || Study_ID), 
                 data = df_UCPCR_clean) %>% summary()


#8 GAD Titre
gadT_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + GAD_v1_z +
                   (dur_diab_months_z || Study_ID), 
                 data =df_UCPCR_clean) %>% summary()

#9 IA2 Status
ia2S_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + IA2_bin_v1 +
                   (dur_diab_months_z || Study_ID), 
                 data = df_UCPCR_clean) %>% summary()

#10 IA2 Titre
ia2T_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + IA2_v1_z +
                   (dur_diab_months_z || Study_ID), 
                 data = df_UCPCR_clean) %>% summary()

#11 History of additional auto-immune disease
auto_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + autoimmune +
                   (dur_diab_months_z || Study_ID), 
                 data = df_UCPCR_clean) %>% summary()

#12 Number of antibodies
num_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + num_anti +
                  (dur_diab_months_z || Study_ID), 
                data = df_UCPCR_clean) %>% summary()

#13 waist hip ratio
whr_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + wh_ratio_v1_z +
                  (dur_diab_months_z || Study_ID), 
                data = df_UCPCR_clean) %>% summary()

#14 T1DGRS2
t1d_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + T1DGRS2_z +
                  (dur_diab_months_z || Study_ID), 
                data = df_UCPCR_clean) %>% summary()


#15 Famiy history of insulin treated diabetes
famhis_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + famhisinsdiab +
                     (dur_diab_months_z || Study_ID), 
                   data = df_UCPCR_clean) %>% summary()

#16 HBA1C
hb_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + HbA1c_at_diagnosis_v1_z +
                 (dur_diab_months_z || Study_ID), 
               data = df_UCPCR_clean) %>% summary()

