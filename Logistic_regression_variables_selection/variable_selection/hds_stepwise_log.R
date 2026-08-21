library(StepReg)
library(dplyr)

#Selected variables Gender, DKA, number of antibodies, GAD status, IA2 Status, GAD Titre

vars_needed <- c('c_peptide_group', 'c_peptide_v1', 'Gender_v1', 'DKA',
                 'num_anti', 'GAD_bin_v1', 'IA2_bin_v1', 
                 'GAD_v1')

#removing all missing values
df_complete = df_define_filtered[complete.cases(df_define_filtered[, vars_needed]), ]
df_complete = df_complete %>% mutate(c_peptide_v1_z = scale(c_peptide_v1),
                                     GAD_v1_z = scale(GAD_v1))
#Stepwise regression
formular_log = c_peptide_group ~ c_peptide_v1_z + Gender_v1 + DKA + IA2_bin_v1 +
  num_anti + GAD_bin_v1 + GAD_v1_z 

res4 = stepwise(formula = formular_log,
                data = df_complete,
                type = "logit",
                include = c('c_peptide_v1_z'),
                strategy = "bidirection",
                metric = "AIC",
)
res4
summary(res4$bidirection$AIC)

#Variables selected from stepwise regression, Gender, GAD titre, IA2 status, dka

#Multivariate logistic regression
df_define_filtered = df_define_filtered %>% mutate(c_peptide_v1_z = scale(c_peptide_v1),
                                     GAD_v1_z = scale(GAD_v1))
multiple_logistic_reg = glm(c_peptide_group ~ c_peptide_v1_z + Gender_v1 + IA2_bin_v1 + DKA +
                    GAD_v1_z,
                  data = df_define_filtered,
                  family = binomial) 
summary(multiple_logistic_reg)
