#Logistic regression model, to assess clinical features that predict c-peptide levels <200


library(broom.mixed)
library(dplyr)
library(purrr)
library(gt)
library(pROC)
library(yardstick)
library(Cairo)


#Categorizing c-peptide in <200 and >=200
df_define = df_define  %>% mutate(c_peptide_group = ifelse(c_peptide_v4 < 200,'1','0'))
df_define $c_peptide_group = as.factor(df_define $c_peptide_group)

#filtering out participants with no c_peptide measure at v4
df_define_filtered = df_define %>% filter(!is.na(c_peptide_group))


# Variables to test
vars = c('Gender_v1', "AgeatDiagnosis_z", "num_anti_z", "DKA", "bmi_calc_v1_z", 
         "autoimmune",'GAD_v1_z','IA2_v1_z','ZNT8_v1_z',
         'wh_ratio_v1_z','GAD_bin_v1',
         'IA2_bin_v1','ZNT8_bin_v1','HbA1c_at_diagnosis_v1_z',
         'Smoker_v1','famhisinsdiab',
         'T1DGRS2_z','T2DGRS_z')

#loop function
log_results = map_df(vars, function(v){
  
  form <- as.formula(
    paste0("c_peptide_group ~ c_peptide_v1_z + ",v )
  )
  
  model <- glm(form, data = df_define_filtered,
               family = binomial)
  
  tidy(model) %>%
    mutate(variable = v)
})
  
log_results

#Tabulating results
log_table <- log_results %>%
  filter(term != "(Intercept)") %>%
  mutate(
    OR = exp(estimate)
  ) %>%
  select(variable, term, OR, p.value)

log_table %>%
  mutate(
    OR = round(OR, 3),
    p.value = signif(p.value, 3)
  ) %>%
  arrange(p.value) %>%
  gt() %>%
  tab_header(title = "Univariate Logistic Regression Screening") %>%
  gt::gtsave("univariate_Log_reg.png")


#log reg for GAD titre
GAD_v1 = glm(c_peptide_group ~ GAD_v1, data = df_define_filtered, family=binomial) %>% summary()
tbl_regression(GAD_v1, exponentiate = TRUE) %>% modify_caption('Univariate logistic reg with GAD titres')
gad_v1=ggplot( df_define_filtered, aes(x=GAD_v1, y=as.numeric(c_peptide_group)))+
  geom_point()+
  geom_smooth(method = "glm")+
  ylab('c_peptide group')+
  xlab('GAD Titres') 
  ggsave('GAD_titre_log.png', gad_v1, height=5, width=5)

  #Prediction
  df_define_filtered$GADpredlink = GAD_v1 %>% predict(df_define_filtered, type = 'link')
  df_define_filtered$GADpred = GAD_v1 %>% predict(df_define_filtered, type = 'response')#predicted probabilities..similar
  #choosing thresholds
  roc_obj_gad = roc(c_peptide_group ~ GADpred, data=df_define_filtered)
  coords(roc_obj_gad, "best", best.method = "youden", ret = "threshold")
  df_define_filtered$GADpredbin = df_define_filtered$GADpred > 0.619
  df_cal <- df_define_filtered %>%
    filter(!is.na(GADpred))
   #predicted probabilities
  p1 = ggplot(df_define_filtered, aes(x=GAD_v1)) +
    geom_line(aes(y=GADpred))
  #predicted log odd
  p2 = ggplot(df_define_filtered, aes(x=GAD_v1)) +
    geom_line(aes(y=GADpredlink))
  p = p1 + p2

  #Discrimination
  #ROC for gad titre
  auc_gadT = roc(c_peptide_group ~ GADpred, data=df_define_filtered,
                      plot = TRUE, print.auc = TRUE, ci=TRUE)
  
  dis_GAD = df_define_filtered %>% mutate(dis_predbin = ifelse(GADpred>0.619, "1", "0")) %>%
    mutate(dis_predbin = factor(dis_predbin, levels = levels(c_peptide_group)))
  
  sens(dis_GAD, dis_predbin, c_peptide_group, event_level = 'second' )
  spec(dis_GAD, dis_predbin, c_peptide_group, event_level = 'second' )
  
  #log reg for GAD_bin
  Gad_bin = glm(c_peptide_group ~ GAD_bin_v1, data = df_define_filtered, family=binomial) %>% summary()
  #regression table
  tbl_regression( Gad_bin, exponentiate = TRUE) %>% modify_caption('Univariate logistic reg with GAD status')
  
  #Prediction for outcome probabilities
  df_define_filtered$GadSpred = Gad_bin %>% predict(df_define_filtered, type = 'response')#predicted probabilities..similar
  #choosing thresholds
  roc_obj_gadS = roc(c_peptide_group ~ GadSpred, data=df_define_filtered)
  coords(roc_obj_gadS , "best", best.method = "youden", ret = "threshold")
  df_define_filtered$GadSpredbin = df_define_filtered$GadSpred > 0.541
  
  #Discrimination
  #ROC 
  auc_gadS = roc(c_peptide_group ~ GadSpred, data=df_define_filtered, plot = TRUE, print.auc = TRUE, ci=TRUE)
  
  dis_GadS = df_define_filtered %>% mutate(dis_gadSpredbin = ifelse(GadSpred>0.541, "1", "0")) %>%
    mutate(dis_gadSpredbin = factor(dis_gadSpredbin, levels = levels(c_peptide_group)))
  #sensitivity
  sens(dis_GadS, dis_gadSpredbin, c_peptide_group, event_level = 'second' ) #0.631
  #specificity
  spec(dis_GadS, dis_gadSpredbin, c_peptide_group, event_level = 'second' ) #0.548
  
  
  
  #log reg for GENDER
  Gender = glm(c_peptide_group ~ Gender_v1, data = df_define_filtered, family=binomial) %>% summary()
  #regression table
  tbl_regression(Gender, exponentiate = TRUE) %>% modify_caption('Univariate logistic reg with Gender')
  
  #Prediction for outcome probabilities
  df_define_filtered$Genpred =  Gender %>% predict(df_define_filtered, type = 'response')#predicted probabilities..similar
  #choosing thresholds
  gen_obj_gad = roc(c_peptide_group ~ Genpred, data=df_define_filtered)
  coords(gen_obj_gad, "best", best.method = "youden", ret = "threshold")
  df_define_filtered$Genpredbin = df_define_filtered$Genpred > 0.613
  
  #Discrimination
  #ROC 
  auc_gen = roc(c_peptide_group ~ Genpred, data=df_define_filtered, plot = TRUE, print.auc = TRUE, ci=TRUE)
  
  dis_Gen = df_define_filtered %>% mutate(dis_genpredbin = ifelse(Genpred>0.613, "1", "0")) %>%
    mutate(dis_genpredbin = factor(dis_genpredbin, levels = levels(c_peptide_group)))
  #sensitivity
  sens(dis_Gen, dis_genpredbin, c_peptide_group, event_level = 'second' ) ...#0.724
  #specificity
  spec(dis_Gen, dis_genpredbin, c_peptide_group, event_level = 'second' ) #0.498
  
  
  #log reg for DKA
  dka = glm(c_peptide_group ~ DKA, data = df_define_filtered, family=binomial) %>% summary()
  #regression table
  tbl_regression(dka, exponentiate = TRUE) %>% modify_caption('Univariate logistic reg with DKA')
  
  #Prediction for outcome probabilities
  df_define_filtered$dkapred = dka %>% predict(df_define_filtered, type = 'response')#predicted probabilities..similar
  #choosing thresholds
  roc_obj_dka = roc(c_peptide_group ~ dkapred, data=df_define_filtered)
  coords(roc_obj_dka, "best", best.method = "youden", ret = "threshold")
  df_define_filtered$dkapredbin = df_define_filtered$dkapred > 0.677
  
  #Discrimination
  #ROC 
  auc_dka = roc(c_peptide_group ~ dkapred, data=df_define_filtered, plot = TRUE, print.auc = TRUE, ci=TRUE)
  
  dis_dka = df_define_filtered %>% mutate(dis_dkapredbin = ifelse(dkapred>0.677, "1", "0")) %>%
    mutate(dis_dkapredbin = factor(dis_dkapredbin, levels = levels(c_peptide_group)))
  #sensitivity
  sens(dis_dka, dis_dkapredbin, c_peptide_group, event_level = 'second' ) #0.779
  #specificity
  spec(dis_dka, dis_dkapredbin, c_peptide_group, event_level = 'second' ) #0.426
  
  
  #log reg for number of antibodies
  num = glm(c_peptide_group ~ num_anti, data = df_define_filtered, family=binomial) %>% summary()
  #regression table
  tbl_regression(num, exponentiate = TRUE) %>% modify_caption('Univariate logistic reg with number of antibodies')
  
  #Prediction for outcome probabilities
  df_define_filtered$numpred = num %>% predict(df_define_filtered, type = 'response')#predicted probabilities..similar
  #choosing thresholds
  roc_obj_num = roc(c_peptide_group ~ numpred, data=df_define_filtered)
  coords(roc_obj_num, "best", best.method = "youden", ret = "threshold")
  
  #Discrimination
  #ROC 
  auc_num = roc(c_peptide_group ~ numpred, data=df_define_filtered, plot = TRUE, print.auc = TRUE, ci=TRUE)
  
  dis_num = df_define_filtered %>% mutate(dis_numpredbin = ifelse(numpred>0.567, "1", "0")) %>%
    mutate(dis_numpredbin = factor(dis_numpredbin, levels = levels(c_peptide_group)))
  #sensitivity
  sens(dis_num, dis_numpredbin, c_peptide_group, event_level = 'second' ) #0.665
  #specificity
  spec(dis_num, dis_numpredbin, c_peptide_group, event_level = 'second' ) #0.483
  
  
  #log reg for IA2 status
  ia2 = glm(c_peptide_group ~ IA2_bin_v1, data = df_define_filtered, family=binomial) %>% summary()
  #regression table
  tbl_regression(ia2, exponentiate = TRUE) %>% modify_caption('Univariate logistic reg with IA2_status')
  
  #Prediction for outcome probabilities
  df_define_filtered$ia2pred = ia2 %>% predict(df_define_filtered, type = 'response')#predicted probabilities..similar
  #choosing thresholds
  roc_obj_ia2 = roc(c_peptide_group ~ ia2pred, data=df_define_filtered)
  coords(roc_obj_ia2, "best", best.method = "youden", ret = "threshold")
  
  #Discrimination
  #ROC 
  auc_ia2 = roc(c_peptide_group ~ ia2pred, data=df_define_filtered, plot = TRUE, print.auc = TRUE, ci=TRUE)
  
  dis_ia2 = df_define_filtered %>% mutate(dis_ia2predbin = ifelse(ia2pred>0.607, "1", "0")) %>%
    mutate(dis_ia2predbin = factor(dis_ia2predbin, levels = levels(c_peptide_group)))
  #sensitivity
  sens(dis_ia2, dis_ia2predbin, c_peptide_group, event_level = 'second' ) #0.673
  #specificity
  spec(dis_ia2, dis_ia2predbin, c_peptide_group, event_level = 'second' ) #0.459
  
  
  #log reg for wh_ratio
  whr = glm(c_peptide_group ~ wh_ratio_v1, data = df_define_filtered, family=binomial) %>% summary()
  #regression table
  tbl_regression(whr, exponentiate = TRUE) %>% modify_caption('Univariate logistic reg with waist hip ratio')
  
  #Prediction for outcome probabilities
  df_define_filtered$whrpred = whr %>% predict(df_define_filtered, type = 'response')#predicted probabilities..similar
  #choosing thresholds
  roc_obj_whr = roc(c_peptide_group ~ whrpred, data=df_define_filtered)
  coords(roc_obj_whr, "best", best.method = "youden", ret = "threshold")
  
  #Discrimination
  #ROC 
  auc_whr = roc(c_peptide_group ~ whrpred, data=df_define_filtered, plot = TRUE, print.auc = TRUE, ci=TRUE)
  
  dis_whr = df_define_filtered %>% mutate(dis_whrpredbin = ifelse(whrpred>0.605, "1", "0")) %>%
    mutate(dis_whrpredbin = factor(dis_whrpredbin, levels = levels(c_peptide_group)))
  #sensitivity
  sens(dis_whr, dis_whrpredbin, c_peptide_group, event_level = 'second' ) #0.680
  #specificity
  spec(dis_whr, dis_whrpredbin, c_peptide_group, event_level = 'second' ) #0.477
  
  
  
  #Multivariate logistic regression
  multi_log = glm(c_peptide_group ~ Gender_v1 + DKA + num_anti +
                                   IA2_bin_v1 + wh_ratio_v1 + GAD_v1 + GAD_bin_v1
                                   , data = df_define_filtered,
                  family = binomial)
  #regression table
  tbl_regression( multi_log, exponentiate = TRUE) %>% modify_caption('Multi-variate logistic reg ')
  
  #Prediction for outcome probabilities
  df_define_filtered$ multipred =  multi_log%>% predict(df_define_filtered, type = 'response')#predicted probabilities..similar
  #choosing thresholds
  roc_obj_multi = roc(c_peptide_group ~ multipred, data=df_define_filtered)
  coords(roc_obj_multi, "best", best.method = "youden", ret = "threshold")
  
  #Discrimination
  #ROC 
  auc_multi= roc(c_peptide_group ~ multipred, data=df_define_filtered, plot = TRUE, print.auc = TRUE, ci=TRUE)
  
  multi_log2 = glm(c_peptide_group ~ Gender_v1 + DKA 
                  , data = df_define_filtered,
                  family = binomial)
  #regression table
  tbl_regression( multi_log2, exponentiate = TRUE) %>% modify_caption('Multi-variate logistic reg ')
  
  #Prediction for outcome probabilities
  df_define_filtered$ multi2pred =  multi_log2 %>% predict(df_define_filtered, type = 'response')#predicted probabilities..similar
  #choosing thresholds
  roc_obj_multi2 = roc(c_peptide_group ~ multi2pred, data=df_define_filtered)
  coords(roc_obj_multi2, "best", best.method = "youden", ret = "threshold")
  
  #Discrimination
  #ROC 
  auc_multi2= roc(c_peptide_group ~ multi2pred, data=df_define_filtered, plot = TRUE, print.auc = TRUE, ci=TRUE)

    
  #plotting ROC and AUC of the variables
  par(mfrow = c(3, 3))
  
  plot(auc_dka, main = " DKA",
       print.auc = TRUE, print.auc.ci = TRUE)
  
  plot(auc_gadS, main = " GAD_status",
       print.auc = TRUE, print.auc.ci = TRUE)
  
  plot(auc_gadT, main = "GAD_Titre",
       print.auc = TRUE, print.auc.ci = TRUE)
  
  plot(auc_gen, main = "Gender",
       print.auc = TRUE, print.auc.ci = TRUE)
  
  plot(auc_ia2, main = "IA2_status",
       print.auc = TRUE, print.auc.ci = TRUE)
  
  plot(auc_num, main = "Number of Antibodies",
       print.auc = TRUE, print.auc.ci = TRUE)
  
  plot(auc_whr, main = "Waist-hip ratio",
       print.auc = TRUE, print.auc.ci = TRUE)
  
  plot(auc_multi, main = 'All variables',
       print.auc = TRUE, print.auc.ci = TRUE)
  
  plot(auc_multi2, main = "Gender+DKA+Num_antibody",
       print.auc = TRUE, print.auc.ci = TRUE)

  
  par(mfrow = c(1, 1))