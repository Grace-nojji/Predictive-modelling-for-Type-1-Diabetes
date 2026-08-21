#Topic : Developing prediction models for beta cell function decline in type 1 diabetes
#Installing necessary packages 

library(ggplot2)
library(tidyverse)
library(gtsummary)
library(dplyr)

#Importing dataset
#Renaming dataset 
df = SR_SSU_2026

#Extracting the duration of diabetes in years for v1 and v2
df = df %>% mutate(home_urine_sample_received_v1 = as.Date(home_urine_sample_received_v1),
                   V2DateUrineReceive_v2         = as.Date(V2DateUrineReceive_v2),
                   date_of_diagnosis_v1          = as.Date(date_of_diagnosis_v1),
                   dur_diab_sample_yr_v1 = as.numeric(home_urine_sample_received_v1 - date_of_diagnosis_v1)/365.25,
                   dur_diab_sample_yr_v2 = as.numeric(V2DateUrineReceive_v2 - date_of_diagnosis_v1)/365.25 
                   )

#Recoding the antibody titres to Zero for the negatives
df = df %>% mutate(GAD_v1 = ifelse(GAD_v1 == "negative", 10, GAD_v1))
df <- df %>% mutate(IA2_v1 = ifelse(IA2_v1 == "negative", 7, IA2_v1))
df <- df %>%
  mutate(
    ZNT8_v1 = case_when(
      ZNT8_v1 == "negative" & AgeatDiagnosis < 30 ~ 64,
      ZNT8_v1 == "negative" & AgeatDiagnosis >= 30 ~ 10,
      TRUE ~ as.numeric(ZNT8_v1)
    )
  )
#Defining Type1 DM ( 1+ anti pos and type1 or 2+ anti pos)
df_define = df %>% filter((clinical_diagnosis_v1=='Type 1' & num_anti == 1) | (num_anti >= 2))


#Correcting data formats
df_define$Gender_v1 = as.factor(df_define$Gender_v1)
df_define$DKA = as.factor( df_define$DKA)
df_define$num_anti = as.factor(df_define$num_anti)
df_define$Smoker_v1 = as.factor(df_define$Smoker_v1 )
df_define$autoimmune = as.factor( df_define$autoimmune )
df_define$famhisdiab = as.factor(df_define$famhisdiab )
df_define$famhisinsdiab = as.factor( df_define$famhisinsdiab )
df_define$GAD_bin_v1 = as.factor( df_define$GAD_bin_v1 )
df_define$IA2_bin_v1 = as.factor( df_define$IA2_bin_v1 )
df_define$ZNT8_bin_v1 = as.factor( df_define$ZNT8_bin_v1)
df_define$famhisauto = as.factor( df_define$famhisauto)
df_define$famhisnoninsdiab = as.factor( df_define$famhisnoninsdiab)
df_define$Mother_diabetes_v1 = as.factor( df_define$Mother_diabetes_v1)
df_define$Father_diabetes_v1 = as.factor(df_define$Father_diabetes_v1)
df_define$GAD_v1 = as.numeric(df_define$GAD_v1)
df_define$ZNT8_v1 = as.numeric(df_define$ZNT8_v1)
df_define$IA2_v1 = as.numeric(df_define$IA2_v1)

#scaling numeric variables
df_define = df_define %>%
  mutate(
    AgeatDiagnosis_z = as.numeric(scale(AgeatDiagnosis)),
    HbA1c_at_diagnosis_v1_z = as.numeric(scale(HbA1c_at_diagnosis_v1)),
    T1DGRS2_z = as.numeric(scale(T1DGRS2)),
    T2DGRS_z = as.numeric(scale(T2DGRS)),
    GAD_v1_z = as.numeric(scale(GAD_v1)),
    ZNT8_v1_z = as.numeric(scale(ZNT8_v1)),
    IA2_v1_z = as.numeric(scale(IA2_v1)),
    bmi_calc_v1_z = as.numeric(scale(bmi_calc_v1)),
    wh_ratio_v1_z = as.numeric(scale(wh_ratio_v1)),
    c_peptide_v1_z = as.numeric(scale(c_peptide_v1)) 
  )


#Converting data from wide format to long format for LINEAR MIXED MODEL
df_long = df_define %>%
  pivot_longer(
    cols = matches("^(UCPCR|dur_diab_sample_yr)_v\\d+$"),
    names_to = c(".value", "visit"),
    names_pattern = "(.*)_v(\\d+)"
  )

#Converting years to month
df_long$dur_diab_months = df_long$dur_diab_sample_yr * 12

#Removing patients without atleast one recording of UCPCR
df_UCPCR_clean = df_long %>%
  group_by(Study_ID) %>%
  filter(any(!is.na(UCPCR))) %>%
  ungroup()

#Grouping fast and slow regressors for LOGISTIC REGRESSION MODEL
#Categorizing c-peptide in <200 and >=200
df_define = df_define  %>% mutate(c_peptide_group = as.factor(
  ifelse(c_peptide_v4 > 200,'slow_regressor', 'fast_regressor')))

#filtering out participants with no c_peptide measure at v4 for LOGISTIC REGRESSION MODEL
df_define_filtered = df_define %>% filter(!is.na(c_peptide_group))



