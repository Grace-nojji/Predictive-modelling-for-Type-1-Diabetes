#Forward selection log reg using step-wise

vars_needed = c('Study_ID', 'c_peptide_v1_z', 'c_peptide_group', 'Gender_v1', 'AgeatDiagnosis_z', 
                                        'num_anti', 'DKA', 'bmi_calc_v1_z',
                                         'autoimmune', 'GAD_v1_z', 'IA2_v1_z', 'ZNT8_v1_z',
                                         'wh_ratio_v1_z', 'GAD_bin_v1',
                                         'IA2_bin_v1', 'ZNT8_bin_v1', 'HbA1c_at_diagnosis_v1_z',
                                         'Smoker_v1', 'famhisinsdiab',
                                         'T1DGRS2_z','T2DGRS_z')

#removing all missing values
df_complete <- df_define_filtered[complete.cases(df_define_filtered[, vars_needed]), ]


formular_log = c_peptide_group ~ c_peptide_v1_z + Gender_v1 + DKA + IA2_bin_v1 +
  num_anti + GAD_bin_v1 + GAD_v1_z 

#Forward selection
forward_step = stepwise(formula = formular_log,
                data = df_complete,
                type = "logit",
                include = c('c_peptide_v1_z'),
                strategy = "forward",
                metric = "AIC",
)
summary(forward_step$forward$AIC)

#Backward selection
back_step = stepwise(formula = formular_log,
                        data = df_complete,
                        type = "logit",
                        include = c('c_peptide_v1_z'),
                        strategy = "backward",
                        metric = "AIC",
)
summary(back_step$backward$AIC)

#Bidirectional selection
bi_step = stepwise(formula = formular_log,
                     data = df_complete,
                     type = "logit",
                     include = c('c_peptide_v1_z'),
                     strategy = "bidirection",
                     metric = "AIC",
)
summary(bi_step$bidirection$AIC)

#Forward using Add1 function

#Forward selection, comparison baseline with model + main effects using add1 function
basic_mod = glm(c_peptide_group ~ c_peptide_v1_z,
                family = 'binomial',
                data = df_complete) 

add1(
  basic_mod, 
  scope = ~ c_peptide_v1_z + Gender_v1 + AgeatDiagnosis_z + num_anti + DKA + bmi_calc_v1_z +
    autoimmune + GAD_v1_z + IA2_v1_z + ZNT8_v1_z +
    wh_ratio_v1_z + GAD_bin_v1 +
    IA2_bin_v1 + ZNT8_bin_v1 + HbA1c_at_diagnosis_v1_z +
    Smoker_v1 + famhisinsdiab +
    T1DGRS2_z + T2DGRS_z  ,
  test = 'Chisq')


