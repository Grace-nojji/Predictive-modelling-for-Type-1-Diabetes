#16 Linear mixed model comparing models with and withou interaction term using LRT
#ggpredict used to visualization
library(marginaleffects)
library(ggeffects)
complete_lmm = df_UCPCR_clean %>% select(UCPCR, Study_ID, dur_diab_months_z,Gender_v1, AgeatDiagnosis_z, num_anti,
                                         DKA, bmi_calc_v1_z,wh_ratio_v1,
                                         autoimmune,GAD_v1_z,IA2_v1_z,ZNT8_v1_z,
                                         wh_ratio_v1_z,GAD_bin_v1,
                                         IA2_bin_v1,ZNT8_bin_v1,HbA1c_at_diagnosis_v1_z,
                                         famhisinsdiab,
                                         T1DGRS2_z) %>% drop_na()



#1 Gender
gen_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + Gender_v1 +
            (dur_diab_months_z || Study_ID), 
          data = complete_lmm)
gen_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + Gender_v1 + 
                 dur_diab_months_z * Gender_v1 + 
                 (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(gen_main,gen_lmm)

 #gt_summary table
gen_tabs = gen_lmm %>% tbl_regression(
  label = list(
    dur_diab_months_z ~ "Duration of Diabetes (Scaled)",
    Gender_v1 ~ "Gender",
    `dur_diab_months_z:Gender_v1` ~ "Duration : Gender Interaction"
  )
) %>% 
  bold_labels() %>%      
  italicize_levels() %>% 
  modify_caption('**LMM analysis with Gender**')

gen_tabs %>%
  as_gt() %>%
  gt::gtsave(filename = "LMM_Gender_Table.png",vheight=200, vwidth=350) 
   

#2
age_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z +AgeatDiagnosis_z +
                  (dur_diab_months_z || Study_ID), 
                data = complete_lmm)
age_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + AgeatDiagnosis_z + 
                 dur_diab_months_z * AgeatDiagnosis_z + 
                 (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(age_main,age_lmm)


#3 DKA
dka_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z +DKA +
                  (dur_diab_months_z || Study_ID), 
                data = complete_lmm)
dka_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + DKA + 
                 dur_diab_months_z * DKA + 
                 (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(dka_main,dka_lmm)
  #gt_summary table
  dka_tabs = dka_lmm %>% tbl_regression(
  label = list(
    dur_diab_months_z ~ "Duration of Diabetes (Scaled)",
    `dur_diab_months_z:DKA` ~ "Duration : DKA"
  )
) %>% 
  bold_labels() %>%      
  italicize_levels() %>% 
  modify_caption('**LMM analysis with DKA**') %>%
  modify_footnote(everything() ~ "DKA: Diabetic Ketoacidosis")

dka_tabs %>%
  as_gt() %>%
  gt::gtsave(filename = "LMM_DKA_Table.png",vheight=200, vwidth=350) 


#4 BMI
bmi_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + bmi_calc_v1_z +
                  (dur_diab_months_z || Study_ID), 
                data = complete_lmm)
bmi_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + bmi_calc_v1_z + 
                 dur_diab_months_z * bmi_calc_v1_z + 
                 (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(bmi_main ,bmi_lmm)

#5 ZNT8 STATUS
zn_stat_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + ZNT8_bin_v1 +
                  (dur_diab_months_z || Study_ID), 
                data = complete_lmm)
zn_stat_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + ZNT8_bin_v1 + 
                     dur_diab_months_z * ZNT8_bin_v1 + 
                     (dur_diab_months_z || Study_ID), data =complete_lmm)
anova(zn_stat_main,zn_stat_lmm)

#6 ZNT8 TITRE
zn_titre_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + ZNT8_v1_z +
                      (dur_diab_months_z || Study_ID), 
                    data = complete_lmm)
zn_titre_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + ZNT8_v1_z + 
                      dur_diab_months_z * ZNT8_v1_z + 
                      (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(zn_titre_main,zn_titre_lmm)


#7 GAD STATUS
gadS_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + GAD_bin_v1 +
                       (dur_diab_months_z || Study_ID), 
                     data = complete_lmm)
gadS_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + GAD_bin_v1 + 
                  dur_diab_months_z * GAD_bin_v1 + 
                  (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(gadS_main,gadS_lmm)



#8 GAD Titre
gadT_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + GAD_v1_z +
                   (dur_diab_months_z || Study_ID), 
                 data = complete_lmm)
gadT_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + GAD_v1_z + 
                  dur_diab_months_z * GAD_v1_z + 
                  (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(gadT_main,gadT_lmm)
 
   #gt_summary table
   gad_tabs = gadT_lmm %>% tbl_regression(
   label = list(
    dur_diab_months_z ~ "Duration of Diabetes (Scaled)",
    GAD_v1_z ~ 'GAD Titre',
    `dur_diab_months_z:GAD_v1_z` ~ "Duration : GAD Titre"
    )
  ) %>% 
  bold_labels() %>%      
  italicize_levels() %>% 
  modify_caption('**LMM analysis with GAD Titre**') %>%
  modify_footnote(
    everything() ~ "GAD: Glutamic acid decarboxylase")

  gad_tabs %>%
  as_gt() %>%
  gt::gtsave(filename = "LMM_gad_Table.png",vheight=200, vwidth=350) 


#9 IA2 Status
ia2S_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + IA2_bin_v1 +
                   (dur_diab_months_z || Study_ID), 
                 data = complete_lmm)
ia2S_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + IA2_bin_v1 + 
                  dur_diab_months_z * IA2_bin_v1 + 
                  (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(ia2S_main,ia2S_lmm)

  #gt_summary table
  ia2_tabs = ia2S_lmm %>% tbl_regression(
  label = list(
    dur_diab_months_z ~ "Duration of Diabetes (Scaled)",
    IA2_bin_v1 ~ 'IA2 Positivity',
    `dur_diab_months_z:IA2_bin_v1` ~ "Duration : IA2 Positivity"
  )
  ) %>% 
  bold_labels() %>%      
  italicize_levels() %>% 
  modify_caption('**LMM analysis with IA2 Positivity**') %>%
  modify_footnote(
    everything() ~ "IA2: Islet Antigen 2")

  ia2_tabs %>%
  as_gt() %>%
  gt::gtsave(filename = "LMM_ia2_Table.png",vheight=200, vwidth=350) 

 
#10 IA2 Titre
ia2T_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + IA2_v1_z +
                   (dur_diab_months_z || Study_ID), 
                 data = complete_lmm)
ia2T_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + IA2_v1_z + 
                  dur_diab_months_z * IA2_v1_z + 
                  (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(ia2T_main,ia2T_lmm)

#gt_summary table
ia2T_tabs = ia2T_lmm %>% tbl_regression(
  label = list(
    dur_diab_months_z ~ "Duration of Diabetes (Scaled)",
    IA2_v1_z ~ 'IA2 Titre',
    `dur_diab_months_z:IA2_v1_z` ~ "Duration : IA2 Titre"
  )
) %>% 
  bold_labels() %>%      
  italicize_levels() %>% 
  modify_caption('**LMM analysis with IA2 Titre**') %>%
  modify_footnote(
    everything() ~ "IA2: Islet antigen 2")

ia2T_tabs %>%
  as_gt() %>%
  gt::gtsave(filename = "LMM_ia2T_Table.png",vheight=200, vwidth=350) 


#11 History of additional auto-immune disease
auto_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + autoimmune +
                   (dur_diab_months_z || Study_ID), 
                 data = complete_lmm)
auto_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + autoimmune + 
                  dur_diab_months_z * autoimmune + 
                  (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(auto_main,auto_lmm)


#12 Number of antibodies
num_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + num_anti +
                   (dur_diab_months_z || Study_ID), 
                 data = complete_lmm)
num_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + num_anti + 
                 dur_diab_months_z * num_anti + 
                 (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(l1,num_lmm)

num_tabs = num_lmm %>% tbl_regression(
  label = list(
    dur_diab_months_z ~ "Duration of Diabetes (Scaled)",
    num_anti ~ 'Number of antibodies',
    `dur_diab_months_z:num_anti` ~ "Duration : Number of antibodies"
  )
) %>% 
  bold_labels() %>%      
  italicize_levels() %>% 
  modify_caption('**LMM analysis with Number of antibodies**')

num_tabs %>%
  as_gt() %>%
  gt::gtsave(filename = "LMM_num_Table.png",vheight=200, vwidth=350) 


#13 waist hip ratio
whr_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + wh_ratio_v1_z +
                  (dur_diab_months_z || Study_ID), 
                data = complete_lmm)
whr_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + wh_ratio_v1_z + 
                 dur_diab_months_z * wh_ratio_v1_z + 
                 (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(whr_main,whr_lmm)

whr_tabs = whr_lmm %>% tbl_regression(
  label = list(
    dur_diab_months_z ~ "Duration of Diabetes (Scaled)",
    wh_ratio_v1_z ~ 'Waist hip ratio',
    `dur_diab_months_z:wh_ratio_v1_z` ~ "Duration : Waist hip ratio"
  )
) %>% 
  bold_labels() %>%      
  italicize_levels() %>% 
  modify_caption('**LMM analysis with Waist-hip ratio**')

whr_tabs %>%
  as_gt() %>%
  gt::gtsave(filename = "LMM_whr_Table.png",vheight=200, vwidth=300) 


#14 T1DGRS2
t1d_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + T1DGRS2_z +
                  (dur_diab_months_z || Study_ID), 
                data = complete_lmm)
t1d_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + T1DGRS2_z + 
                 dur_diab_months_z * T1DGRS2_z + 
                 (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(t1d_main,t1d_lmm)


#15 Famiy history of insulin treated diabetes
famhis_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + famhisinsdiab +
                  (dur_diab_months_z || Study_ID), 
                data = complete_lmm)
famhis_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + famhisinsdiab + 
                    dur_diab_months_z * famhisinsdiab + 
                    (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(famhis_main,famhis_lmm)




#16 HBA1C
hb_main = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + HbA1c_at_diagnosis_v1_z +
                  (dur_diab_months_z || Study_ID), 
                data = complete_lmm)
hb_lmm = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + HbA1c_at_diagnosis_v1_z + 
                dur_diab_months_z * HbA1c_at_diagnosis_v1_z + 
                (dur_diab_months_z || Study_ID), data = complete_lmm)
anova(hb_main,hb_lmm)

