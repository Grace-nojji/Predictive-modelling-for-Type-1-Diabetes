#Variable selection for logistic regression
#18 variables explored, backward elimination using drop1 function

library(pROC)
library(reliabilitydiag)
library(gtsummary)
library(patchwork)
library(ggplot2)

#Categorizing c-peptide in <200 and >=200
df_define = df_define  %>% mutate(c_peptide_group = as.factor(
          ifelse(c_peptide_v4 < 200,'1','0')))

#filtering out participants with no c_peptide measure at v4
df_define_filtered = df_define %>% filter(!is.na(c_peptide_group))

#Using drop function to find significant variables
#0) 
cpep_log = glm(c_peptide_group ~ c_peptide_v1_z,
              family = 'binomial',
              data = df_define_filtered) %>% summary()
drop1(cpep_log ,  test='Chisq')

#ROC AUC
df_define_filtered$cpep_pred =  cpep_log %>% predict(df_define_filtered, type = 'response')
roc_obj_cpep_pred = roc(c_peptide_group ~ cpep_pred, data=df_define_filtered,  ci=TRUE) #0.7586
plot(roc_obj_cpep_pred)

#Reliability
cpep_re = reliabilitydiag(
  fitted(cpep_log),
  y = as.numeric(model.frame(cpep_log)$c_peptide_group) - 1
)
plot(cpep_re)


#1) 
gen_log = glm(c_peptide_group ~ Gender_v1 + c_peptide_v1_z,
                family = binomial,
               data = df_define_filtered) 
drop1(gen_log ,  test='Chisq')

#ROC AUC
df_define_filtered$gen_pred =  gen_log %>% predict(df_define_filtered, type = 'response')
roc_obj_gen_pred = roc(c_peptide_group ~ gen_pred, data=df_define_filtered, ci=TRUE) #0.763
plot(roc_obj_gen_pred)

#Reliability
gen_re = reliabilitydiag(
  fitted(gen_log),
  y = as.numeric(model.frame(gen_log)$c_peptide_group) - 1
)

#2) 
dka_log = glm(c_peptide_group ~ DKA + c_peptide_v1_z,
              family = binomial,
              data = df_define_filtered)
drop1(dka_log ,  test='Chisq')

#ROC AUC
df_define_filtered$dka_pred =  dka_log %>% predict(df_define_filtered, type = 'response')
roc_obj_dka_pred = roc(c_peptide_group ~ dka_pred, data=df_define_filtered,  ci=TRUE) #0.760
plot(roc_obj_dka_pred)

#Reliability
dka_re = reliabilitydiag(
  fitted(dka_log),
  y = as.numeric(model.frame(dka_log)$c_peptide_group) - 1
)

#3)
num_log = glm(c_peptide_group ~ num_anti + c_peptide_v1_z,
              family = binomial,
              data = df_define_filtered) 
drop1(num_log,  test='Chisq')

#ROC AUC
df_define_filtered$num_pred =  num_log %>% predict(df_define_filtered, type = 'response')
roc_obj_num_pred = roc(c_peptide_group ~ num_pred, data=df_define_filtered,  ci=TRUE) #0.768
plot(roc_obj_num_pred)

#Reliability
num_re = reliabilitydiag(
  fitted(num_log),
  y = as.numeric(model.frame(num_log)$c_peptide_group) - 1
)

#4)
gad_log = glm(c_peptide_group ~ GAD_bin_v1 + c_peptide_v1_z,
              family = binomial,
              data = df_define_filtered) 
drop1(gad_log,  test='Chisq')

#ROC AUC
df_define_filtered$gad_pred =  gad_log %>% predict(df_define_filtered, type = 'response')
roc_obj_gad_pred= roc(c_peptide_group ~ gad_pred, data=df_define_filtered, ci=TRUE) #0.766
plot(roc_obj_gad_pred)

#Reliability
gad_re = reliabilitydiag(
  fitted(gad_log),
  y = as.numeric(model.frame(gad_log)$c_peptide_group) - 1
)

#5)
znt8_log = glm(c_peptide_group ~ ZNT8_bin_v1 + c_peptide_v1_z,
              family = binomial,
              data = df_define_filtered)
drop1(znt8_log,  test='Chisq')

#ROC AUC
df_define_filtered$znt8_pred =  znt8_log %>% predict(df_define_filtered, type = 'response')
roc_obj_znt8_pred= roc(c_peptide_group ~ znt8_pred, data=df_define_filtered, ci=TRUE) #0.760
plot(roc_obj_znt8_pred)

#Reliability
zn_re = reliabilitydiag(
  fitted(znt8_log),
  y = as.numeric(model.frame(znt8_log)$c_peptide_group) - 1
)

#6)
ia2_log = glm(c_peptide_group ~ IA2_bin_v1 + c_peptide_v1_z,
               family = binomial,
               data = df_define_filtered) 
drop1(ia2_log,  test='Chisq')

#ROC AUC
df_define_filtered$ia2_pred = ia2_log %>% predict(df_define_filtered, type = 'response')
roc_obj_ia2_pred = roc(c_peptide_group ~ ia2_pred, data=df_define_filtered, ci=TRUE) #0.763
plot(roc_obj_ia2_pred)

#Reliability
ia2_re = reliabilitydiag(
  fitted(ia2_log),
  y = as.numeric(model.frame(ia2_log)$c_peptide_group) - 1
)

#7)
hba1c_log = glm(c_peptide_group ~ HbA1c_at_diagnosis_v1_z + c_peptide_v1_z,
              family = binomial,
              data = df_define_filtered) 
drop1(hba1c_log,  test='Chisq')

#ROC AUC
df_define_filtered$hba1c_pred = hba1c_log %>% predict(df_define_filtered, type = 'response')
roc_obj_hba1c_pred = roc(c_peptide_group ~ hba1c_pred, data=df_define_filtered,  ci=TRUE ) #0.750
plot(roc_obj_hba1c_pred)

#8)
age_log = glm(c_peptide_group ~ AgeatDiagnosis_z + c_peptide_v1_z,
                family = binomial,
                data = df_define_filtered) 
drop1(age_log,  test='Chisq')
df_define_filtered$age_pred = age_log %>% predict(df_define_filtered, type = 'response')
roc_obj_age = roc(c_peptide_group ~ age_pred, data=df_define_filtered, ci=TRUE) #0.76(0.71-0.81)

#9)
t1dgrs_log = glm(c_peptide_group ~ T1DGRS2_z + c_peptide_v1_z,
              family = binomial,
              data = df_define_filtered) 
drop1(t1dgrs_log,  test='Chisq')
df_define_filtered$t1dgrs_pred = t1dgrs_log  %>% predict(df_define_filtered, type = 'response')
roc_obj_t1dgrs = roc(c_peptide_group ~ t1dgrs_pred, data=df_define_filtered, ci=TRUE) #0.76(0.71-0.81)


#10)
t2d_log = glm(c_peptide_group ~ T2DGRS_z + c_peptide_v1_z,
                 family = binomial,
                 data = df_define_filtered) 
drop1(t2d_log,  test='Chisq')
df_define_filtered$t2d_pred = t2d_log %>% predict(df_define_filtered, type = 'response')
roc_obj_t2d = roc(c_peptide_group ~ t2d_pred, data=df_define_filtered, ci=TRUE) #0.76(0.72-0.81)


#11)
famins_log = glm(c_peptide_group ~ famhisinsdiab + c_peptide_v1_z,
              family = binomial,
              data = df_define_filtered) 
drop1(famins_log,  test='Chisq')
df_define_filtered$famins_pred = famins_log %>% predict(df_define_filtered, type = 'response')
roc_obj_famins = roc(c_peptide_group ~ famins_pred, data=df_define_filtered, ci=TRUE) #0.76(0.71-0.81)

#12)
auto_log = glm(c_peptide_group ~ autoimmune + c_peptide_v1_z,
                 family = binomial,
                 data = df_define_filtered) 
drop1(auto_log,  test='Chisq')
df_define_filtered$auto_pred = auto_log %>% predict(df_define_filtered, type = 'response')
roc_obj_auto = roc(c_peptide_group ~ auto_pred, data=df_define_filtered, ci=TRUE) #0.76(0.71-0.81)

#13)
smk_log = glm(c_peptide_group ~ Smoker_v1 + c_peptide_v1_z,
               family = binomial,
               data = df_define_filtered) 
drop1(smk_log,  test='Chisq')
df_define_filtered$smk_pred = smk_log %>% predict(df_define_filtered, type = 'response')
roc_obj_smk = roc(c_peptide_group ~smk_pred, data=df_define_filtered, ci=TRUE) #0.76(0.71-0.81)

#14)
whr_log = glm(c_peptide_group ~ wh_ratio_v1_z + c_peptide_v1_z,
              family = binomial,
              data = df_define_filtered)
drop1(whr_log,  test='Chisq')
      df_define_filtered$whr_pred = whr_log %>% predict(df_define_filtered, type = 'response')
      roc_obj_whr = roc(c_peptide_group ~ whr_pred, data=df_define_filtered, ci=TRUE) #0.75(0.71-0.80)
      
#15)
titre_ia2_log = glm(c_peptide_group ~ IA2_v1_z + c_peptide_v1_z,
                    family = binomial,
                    data = df_define_filtered) 
      drop1(titre_ia2_log,  test='Chisq') 
      
#ROC AUC
df_define_filtered$ti_ia2_pred = titre_ia2_log %>% predict(df_define_filtered, type = 'response')
roc_obj_ti_ia2_pred = roc(c_peptide_group ~ ti_ia2_pred, data=df_define_filtered, ci=TRUE) #0.76
 plot(roc_obj_ti_ia2_pred)      
      
#16)
titre_gad_log = glm(c_peptide_group ~ GAD_v1_z + c_peptide_v1_z,
                          family = binomial,
                          data = df_define_filtered)
      drop1(titre_gad_log,  test='Chisq') 
      
#ROC AUC
      df_define_filtered$ti_gad_pred = titre_gad_log %>% predict(df_define_filtered, type = 'response')
      roc_obj_ti_gad_pred = roc(c_peptide_group ~ ti_gad_pred, data=df_define_filtered, ci=TRUE) #0.766
      plot(roc_obj_ti_gad_pred) 
      
      #Reliability
      gad_titre_re = reliabilitydiag(
        fitted(titre_gad_log),
        y = as.numeric(model.frame(titre_gad_log)$c_peptide_group) - 1
      )
      
#17)
titre_znt8_log = glm(c_peptide_group ~ ZNT8_v1_z + c_peptide_v1_z,
                          family = binomial,
                          data = df_define_filtered) 
      drop1(titre_znt8_log,  test='Chisq') 
      #ROC AUC
      df_define_filtered$titre_znt8_pred = titre_znt8_log %>% predict(df_define_filtered, type = 'response')
      roc_obj_titre_znt8 = roc(c_peptide_group ~ titre_znt8_pred, data=df_define_filtered, ci=TRUE) #0.75(0.71-0.80)
      plot(roc_obj_titre_znt8) 
      
#18)
      bmi_log = glm(c_peptide_group ~ bmi_calc_v1_z + c_peptide_v1_z,
                           family = binomial,
                           data = df_define_filtered) 
      drop1(bmi_log,  test='Chisq') 
      #ROC AUC
      df_define_filtered$bmi_pred = bmi_log %>% predict(df_define_filtered, type = 'response')
      roc_obj_bmi = roc(c_peptide_group ~ bmi_pred, data=df_define_filtered, ci=TRUE) #0.76(0.71-0.81)
      plot(roc_obj_bmi) 
      
      
          
#18)
first_log = glm(c_peptide_group ~  c_peptide_v1_z + 
                                   Gender_v1 + DKA + num_anti_z + GAD_bin_v1 +
                                   ZNT8_bin_v1 +  GAD_v1_z,
                           family = binomial,
                           data = df_define_filtered) %>% summary()
#ROC AUC
df_define_filtered$firstlog_pred = first_log %>% predict(df_define_filtered, type = 'response')
roc_obj_firstlog_pred  = roc(c_peptide_group ~ firstlog_pred , data=df_define_filtered, ci=TRUE) #0.79 (0.74-0.83)
plot(roc_obj_ti_gad_pred)  

first_log_table = tbl_regression(first_log,
                                   label = list(
                                     c_peptide_v1_z ~ "c_peptide at visit 1",
                                     Gender_v1               ~ "Gender",
                                     DKA                     ~ "DKA",
                                     num_anti_z              ~ "Number of antibodies",
                                     GAD_bin_v1             ~ "GAD antibody positivity",
                                     GAD_v1_z               ~ "GAD Titre",
                                     ZNT8_bin_v1           ~ 'ZNT8 antibody positivity'
                                   )
                                 ) %>%
                                   modify_caption("**Multivariate Logistic-A**") %>%
                                   add_global_p() %>%
                                   modify_footnote(
                                     everything() ~ NA  # clears default footnotes
                                   ) %>%
                                   modify_table_styling(
                                     columns = label,
                                     footnote = "DKA: Diabetic Ketoacidosis; 
                                     ZNT8:Zinc Transporter 8; 
                                     GAD: Glutamic acid decarboxylase"
                                   )
first_log_table %>% 
  as_gt() %>%
  #gt::tab_options(table.width = gt::px(100)) %>%
  gt::gtsave("first_log_table.png", vwidth = 400, vheight = 100)
                                 

#Reliability
first_re = reliabilitydiag(
  fitted(first_log ),
  y = as.numeric(model.frame(first_log )$c_peptide_group) - 1
)
plot(first_re)     

#plotting Reliabilty diagram of variables of the variables
p1 <- plot(cpep_re) + labs(title = "C-peptide at first visit")
p2 <- plot(gen_re)  + labs(title = "Gender")
p3 <- plot(dka_re)  + labs(title = "DKA")
p4 <- plot(gad_re)  + labs(title = "GAD status")
p5 <- plot(gad_titre_re)  + labs(title = "GAD Titre")
p6 <- plot(zn_re)  + labs(title = "ZNT8 Status")
p7 <- plot(first_re)  + labs(title = "All Variables")


Cali_plots = patchwork::wrap_plots(
  p1,p2,p3,p4,p5,p6,p7,p8,
  nrow = 2, ncol=4
) + patchwork::plot_annotation(tag_levels = 'a')
Cairo::CairoPNG('Cairo_cali_plots.png', width = 40, height = 15, units = "in", dpi = 300)
print(Cali_plots)
dev.off()


Cairo::CairoPNG('Cairo_cali_plots.png', width = 35, height = 20)
Cali_plots
dev.off()
ggsave('cali_plots.pdf', Cali_plots, width = 35, height = 20 )

#Selected variables gender, dka, number of antibodies, antibodies(GAD,ZNT8), gad titres

vars_needed <- c('c_peptide_group', 'c_peptide_v1_z', 'Gender_v1', 'DKA',
                 'num_anti_z', 'GAD_bin_v1', 'IA2_bin_v1', 
                 'GAD_v1_z')

#removing all missing values
df_complete <- df_define_filtered[complete.cases(df_define_filtered[, vars_needed]), ]


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
plot(res4, strategy = "bidirection", process = "overview")


#Multivariate logistic regression
multi_log_2 = glm(c_peptide_group ~ c_peptide_v1_z + Gender_v1 + num_anti_z +
                  GAD_v1_z,
                 data = df_define_filtered,
                family = binomial) %>% summary()
#regression table
log_table = tbl_regression(multi_log_2 ,
               label = list(
                 c_peptide_v1_z ~ "c_peptide at visit 1",
                 Gender_v1               ~ "Gender",
                num_anti_z              ~ "Number of antibodies",
                 GAD_v1_z             ~ "GAD antibody positivity"
               )
) %>%
  modify_caption("**Multivariate Logistic_B**") %>%
  add_global_p() %>%
  modify_footnote(
    everything() ~ NA  # clears default footnotes
  ) %>%
  modify_table_styling(
    columns = label,
    footnote = "DKA: Diabetic Ketoacidosis; 
    GAD: Glutamic acid decarboxylase"
  )
log_table  %>% 
  as_gt() %>%
  gt::gtsave("log_table.png", vwidth = 350, vheight = 400)

#Prediction for outcome probabilities
df_define_filtered$ multipred2 =  multi_log_2 %>% predict(df_define_filtered, type = 'response')#predicted probabilities..similar
#choosing thresholds
roc_obj_multi2 = roc(c_peptide_group ~ multipred2, data=df_define_filtered, ci=TRUE) #0.78 (0.73-0.83)
plot(roc_obj_multi2, ci=TRUE)

#Discrimination
#ROC 
png("roc_log.png", width = 6, height = 4, units = "in", res = 300)
auc_multi_2= roc(c_peptide_group ~ multipred2, data=df_define_filtered, plot = TRUE, print.auc = TRUE, ci=TRUE)
dev.off()

#Calibration
mf = model.frame(multi_log_2)

png("roc_calibration_log.png", width = 3, height = 3, units = "in", res = 300)
cal = reliabilitydiag(
  fitted(multi_log_2),
  y = as.numeric(mf$c_peptide_group) - 1
)
plot(cal)
dev.off()

par(mfrow = c(1, 2))

plot(roc_obj_firstlog_pred, main = " A",
     print.auc = TRUE, print.auc.ci = TRUE)
plot(roc_obj_multi2, main = "B ",
     print.auc = TRUE, print.auc.ci = TRUE)

par(mfrow = c(1, 1))






