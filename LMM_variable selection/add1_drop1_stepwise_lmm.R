#VARIABLE SELECTION PRT 3

library(lme4)
library(lmerTest)

#Using drop function to find significant variables
complete_lmm = df_UCPCR_clean %>% select(UCPCR, Study_ID, dur_diab_months_z,Gender_v1, AgeatDiagnosis_z, num_anti, DKA, bmi_calc_v1_z,
                                         autoimmune,GAD_v1_z,IA2_v1_z,ZNT8_v1_z,
                                         wh_ratio_v1_z,GAD_bin_v1,
                                         IA2_bin_v1,ZNT8_bin_v1,HbA1c_at_diagnosis_v1_z,
                                         famhisinsdiab,
                                         T1DGRS2_z) %>% drop_na()

#Forward selection, comparison baseline with model + main effects using add1 function
basic_mod = lmer(scale(log(UCPCR)) ~ dur_diab_months_z +
                          (dur_diab_months_z || Study_ID),
                           data = complete_lmm)

add1(
  basic_mod, 
  scope = ~ dur_diab_months_z + Gender_v1 + AgeatDiagnosis_z + num_anti + DKA + bmi_calc_v1_z +
               autoimmune + GAD_v1_z + IA2_v1_z + ZNT8_v1_z +
               wh_ratio_v1_z + GAD_bin_v1 +
               IA2_bin_v1 + ZNT8_bin_v1 + HbA1c_at_diagnosis_v1_z +
               famhisinsdiab +
               T1DGRS2_z ,
  test = 'Chisq')

#Backward selection, comparing with or without interaction term
#1) 
gen_mod = lmer(scale(log(UCPCR)) ~ dur_diab_months_z + Gender_v1 +
                            dur_diab_months_z * Gender_v1 +
                   (dur_diab_months_z || Study_ID),
                 data = df_UCPCR_clean)
drop1(gen_mod,  test='Chisq')
gen_mod1 = update(gen_mod,.~.-dur_diab_months_z:Gender_v1) 
drop1(gen_mod,  test='Chisq')

#2) 
DKA_mod = lmer(scale(log(UCPCR)) ~ dur_diab_months_z*DKA +
                 (dur_diab_months_z||Study_ID),
               data = df_UCPCR_clean)
drop1(DKA_mod,  test='Chisq')
DKA_mod1 = update(DKA_mod,.~.-dur_diab_months_z:DKA) 
drop1(DKA_mod1,  test='Chisq')

#3)
anti_mod = lmer(scale(log(UCPCR)) ~ dur_diab_months_z*num_anti +
                 (dur_diab_months_z||Study_ID),
               data = df_UCPCR_clean)
drop1(anti_mod,  test='Chisq')
anti_mod1 = update(anti_mod,.~.-dur_diab_months_z:num_anti) 
drop1(anti_mod1,  test='Chisq')

#4)
gad_mod = lmer(scale(log(UCPCR)) ~ dur_diab_months_z*GAD_bin_v1 +
                 (dur_diab_months_z||Study_ID),
               data = df_UCPCR_clean)
drop1(gad_mod,  test='Chisq')
gad_mod1 = update(gad_mod,.~.-dur_diab_months_z:GAD_bin_v1) 
drop1(gad_mod1,  test='Chisq')

#5)
znt8_mod = lmer(scale(log(UCPCR)) ~ dur_diab_months_z*ZNT8_bin_v1 +
                 (dur_diab_months_z||Study_ID),
               data = df_UCPCR_clean)
drop1(znt8_mod,  test='Chisq')
znt8_mod1 = update(znt8_mod,.~.-dur_diab_months_z:ZNT8_bin_v1) 
drop1(znt8_mod1,  test='Chisq')

#6)
IA2_mod = lmer(scale(log(UCPCR)) ~ dur_diab_months_z*IA2_bin_v1 +
                  (dur_diab_months_z||Study_ID),
                data = df_UCPCR_clean)
drop1(IA2_mod,  test='Chisq')
IA2_mod1 = update(IA2_mod,.~.-dur_diab_months_z:IA2_bin_v1) 
drop1(IA2_mod1 ,  test='Chisq')

#7)
hba1c_mod = lmer(scale(log(UCPCR)) ~ dur_diab_months_z*HbA1c_at_diagnosis_v1_z +
                 (dur_diab_months_z||Study_ID),
               data = df_UCPCR_clean)
drop1(hba1c_mod,  test='Chisq')
hba1c_mod1 = update(hba1c_mod,.~.-dur_diab_months_z:HbA1c_at_diagnosis_v1_z ) 
drop1(hba1c_mod1,  test='Chisq')

#8)
age_mod = lmer(scale(log(UCPCR)) ~ dur_diab_months_z*AgeatDiagnosis_z +
                   (dur_diab_months_z||Study_ID),
                 data = df_UCPCR_clean)
drop1(age_mod,  test='Chisq')
age_mod1 = update(age_mod,.~.-dur_diab_months_z:AgeatDiagnosis_z ) 
drop1(age_mod1,  test='Chisq')

#9)
t1dgrs_mod = lmer(scale(log(UCPCR)) ~ dur_diab_months_z*T1DGRS2_z+
                 (dur_diab_months_z||Study_ID),
               data = df_UCPCR_clean)
drop1(t1dgrs_mod ,  test='Chisq')
t1dgrs_mod1 = update(t1dgrs_mod,.~.-dur_diab_months_z:T1DGRS2_z) 
drop1(t1dgrs_mod1,  test='Chisq')

#10)
famins_mod = lmer(scale(log(UCPCR)) ~ dur_diab_months_z*famhisinsdiab+
                    (dur_diab_months_z||Study_ID),
                  data = df_UCPCR_clean)
drop1(famins_mod ,  test='Chisq')
famins_mod1 = update(famins_mod,.~.-dur_diab_months_z:famhisinsdiab) 
drop1(famins_mod1,  test='Chisq')

#11)
auto_mod = lmer(scale(log(UCPCR)) ~ dur_diab_months_z*autoimmune+
                    (dur_diab_months_z||Study_ID),
                  data = df_UCPCR_clean)
drop1(auto_mod ,  test='Chisq')
auto_mod1 = update(auto_mod,.~.-dur_diab_months_z:autoimmune) 
drop1(auto_mod1,  test='Chisq')

#12)
whr_mod = lmer(scale(log(UCPCR)) ~ dur_diab_months_z*wh_ratio_v1_z+
                    (dur_diab_months_z||Study_ID),
                  data = df_UCPCR_clean)
drop1(whr_mod ,  test='Chisq')
whr_mod1 = update(whr_mod,.~.-dur_diab_months_z:wh_ratio_v1_z) 
drop1(whr_mod1,  test='Chisq')

#14)
GadT_mod = lmer(scale(log(UCPCR)) ~ dur_diab_months_z*GAD_v1_z+
                 (dur_diab_months_z||Study_ID),
               data = df_UCPCR_clean)
drop1(GadT_mod,  test='Chisq')
GadT_mod1 = update(GadT_mod,.~.-dur_diab_months_z:GAD_v1_z_v1_z) 
drop1(GadT_mod1,  test='Chisq')

#15)
ZNT8T_mod = lmer(scale(log(UCPCR)) ~ dur_diab_months_z*ZNT8_v1_z+
                 (dur_diab_months_z||Study_ID),
               data = df_UCPCR_clean)
drop1(ZNT8T_mod ,  test='Chisq')
ZNT8T_mod1 = update(ZNT8T_mod,.~.-dur_diab_months_z:ZNT8_v1_z) 
drop1(ZNT8_v1_z1,  test='Chisq')

#16)
IA2T_mod = lmer(scale(log(UCPCR)) ~ dur_diab_months_z*IA2_v1_z+
                 (dur_diab_months_z||Study_ID),
               data = df_UCPCR_clean)
drop1(IA2T_mod  ,  test='Chisq')
IA2T_mod1 = update(IA2T_mod,.~.-dur_diab_months_z:IA2_v1_z) 
drop1(IA2T_mod1,  test='Chisq')



#Variables selected for model building: Gender,DKA,hx of additional auto-immune disease,/
 #ZNT8,IA2,HBA1C,T1DGRS, Family history of insulin treated diabetes,and wh_ratio

#Model building using stepwise elimination by stepwise

model1= (UCPCR_log ~ dur_diab_months_z * (Gender_v1 + DKA +
                      autoimmune + ZNT8_bin_v1 +
                      T1DGRS2_z + IA2_bin_v1 +famhisinsdiab + 
                      wh_ratio_v1_z + HbA1c_at_diagnosis_v1_z) +
                      (1 | Study_ID))

#forward selection
model1_for = buildmer(model1, data=df_UCPCR_clean, buildmerControl=buildmerControl(direction='order',
                                               args=list(control=lmerControl(optimizer='bobyqa'))))

(q <- formula(model1_for@model))

summary(model1_for)
#backward selection
model1_bac = buildmer(model1,data = df_UCPCR_clean, buildmerControl=list(direction='backward',
                                                 args=list(control=lmerControl(optimizer='bobyqa'))))
summary(j)


