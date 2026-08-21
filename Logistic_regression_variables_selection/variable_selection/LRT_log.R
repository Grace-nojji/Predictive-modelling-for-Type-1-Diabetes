#LRT for logistic reg
#Comparing baseline models, to those with variables above and beyond
library(pROC)
library(reliabilitydiag)
library(gtsummary)
library(patchwork)
library(ggplot2)

#Categorizing c-peptide in <200 and >=200
df_define = df_define  %>% mutate(c_peptide_group = ifelse(c_peptide_v4 < 200,'1','0'))
df_define $c_peptide_group = as.factor(df_define $c_peptide_group)

#filtering out participants with no c_peptide measure at v4
df_define_filtered = df_define %>% filter(!is.na(c_peptide_group))

#dropping all missing variables, to allow a fair comparison
complete_cases = df_define_filtered %>% select(c_peptide_group, Study_ID,c_peptide_v1_z,
                                         Gender_v1, AgeatDiagnosis_z, num_anti,
                                         DKA, bmi_calc_v1_z,
                                         autoimmune,GAD_v1_z,IA2_v1_z,ZNT8_v1_z,
                                         wh_ratio_v1_z,GAD_bin_v1,
                                         IA2_bin_v1,ZNT8_bin_v1,HbA1c_at_diagnosis_v1_z,
                                         Smoker_v1,famhisinsdiab,
                                         T1DGRS2_z,T2DGRS_z) %>% drop_na()

#Using drop function to find significant variables
#0) 
cpep_log = glm(c_peptide_group ~ c_peptide_v1_z,
               family = binomial,
               data = complete_cases) 

#1) 
gen_log = glm(c_peptide_group ~ Gender_v1 + c_peptide_v1_z,
              family = binomial,
              data = complete_cases) 
anova(cpep_log,gen_log)

#2) 
dka_log = glm(c_peptide_group ~ DKA + c_peptide_v1_z,
              family = binomial,
              data = complete_cases) 
anova(cpep_log,dka_log )

#3)
num_log = glm(c_peptide_group ~ num_anti + c_peptide_v1_z,
              family = binomial,
              data = complete_cases) 
anova(cpep_log,num_log )

#4)
gad_log = glm(c_peptide_group ~ GAD_bin_v1 + c_peptide_v1_z,
              family = binomial,
              data = complete_cases)
anova(cpep_log,gad_log )

#5)
znt8_log = glm(c_peptide_group ~ ZNT8_bin_v1 + c_peptide_v1_z,
               family = binomial,
               data = complete_cases) 
anova(cpep_log,znt8_log )

#6)
ia2_log = glm(c_peptide_group ~ IA2_bin_v1 + c_peptide_v1_z,
              family = binomial,
              data = complete_cases) 
anova(cpep_log,ia2_log )

#7)
hba1c_log = glm(c_peptide_group ~ HbA1c_at_diagnosis_v1_z + c_peptide_v1_z,
                family = binomial,
                data = complete_cases)
anova(cpep_log,hba1c_log )

#8)
age_log = glm(c_peptide_group ~ AgeatDiagnosis_z + c_peptide_v1_z,
              family = binomial,
              data = complete_cases) 
anova(cpep_log,age_log )

#9)
t1dgrs_log = glm(c_peptide_group ~ T1DGRS2_z + c_peptide_v1_z,
                 family = binomial,
                 data = complete_cases)
anova(cpep_log,t1dgrs_log)

#10)
t2d_log = glm(c_peptide_group ~ T2DGRS_z + c_peptide_v1_z,
              family = binomial,
              data = complete_cases) 
anova(cpep_log,t2d_log)

#11)
famins_log = glm(c_peptide_group ~ famhisinsdiab + c_peptide_v1_z,
                 family = binomial,
                 data = complete_cases) 
anova(cpep_log,famins_log)

#12)
auto_log = glm(c_peptide_group ~ autoimmune + c_peptide_v1_z,
               family = binomial,
               data = complete_cases)
anova(cpep_log,auto_log)

#13)
smk_log = glm(c_peptide_group ~ Smoker_v1 + c_peptide_v1_z,
              family = binomial,
              data = complete_cases) 
anova(cpep_log,smk_log)

#14)
whr_log = glm(c_peptide_group ~ wh_ratio_v1_z + c_peptide_v1_z,
              family = binomial,
              data = complete_cases) 
anova(cpep_log,whr_log)

#15)
titre_ia2_log = glm(c_peptide_group ~ IA2_v1_z + c_peptide_v1_z,
                    family = binomial,
                    data = complete_cases)
anova(cpep_log,titre_ia2_log)

#16)
titre_gad_log = glm(c_peptide_group ~ GAD_v1_z + c_peptide_v1_z,
                    family = binomial,
                    data = complete_cases) 
anova(cpep_log,titre_gad_log)

#17)
titre_znt8_log = glm(c_peptide_group ~ ZNT8_v1_z + c_peptide_v1_z,
                     family = binomial,
                     data = complete_cases) 
anova(cpep_log,titre_znt8_log)

#18)
bmi_log = glm(c_peptide_group ~ bmi_calc_v1_z + c_peptide_v1_z,
              family = binomial,
              data = complete_cases) 
anova(cpep_log,bmi_log)

