library('lme4')
library('lmerTest')
library('buildmer')
library('ggeffects')
library('marginaleffects')


df_UCPCR_clean$scaled_logUCPCR = scale(log(df_UCPCR_clean$UCPCR))


rirs_step = lmer(scale(log(UCPCR)) ~ 
                   dur_diab_months_z * Gender_v1 +
                   dur_diab_months_z * num_anti +
                   dur_diab_months_z * DKA +
                   (dur_diab_months_z || Study_ID), data = df_UCPCR_clean)
rirs_step  %>% summary()


full_model_lmm_table = tbl_regression(
  rirs_step ,
  label = list(
    dur_diab_months_z       ~ "Duration of diabetes",
    Gender_v1               ~ "Gender",
    DKA                     ~ "DKA",
    num_anti ~ 'Number of antibodies'
  )
) %>%
  modify_caption("**Multivariable Linear Mixed Model**") %>%
  add_global_p() %>%
  modify_footnote(
    everything() ~ NA  # clears default footnotes
  ) %>%
  modify_table_styling(
    columns = label,
    footnote = "DKA: Diabetic Ketoacidosis"
  )
full_model_lmm_table  %>% 
  as_gt() %>%
  gt::gtsave("full_lmm_tables.png", vwidth = 1000, vheight = 400)


############
#VISUALIZING
############
complete_lmm = df_UCPCR_clean %>% select(UCPCR, Study_ID, dur_diab_months_z,Gender_v1, AgeatDiagnosis_z, num_anti, DKA, bmi_calc_v1_z,
                                         autoimmune,GAD_v1_z,IA2_v1_z,ZNT8_v1_z,
                                         wh_ratio_v1_z,GAD_bin_v1,
                                         IA2_bin_v1,ZNT8_bin_v1,HbA1c_at_diagnosis_v1_z,
                                         Smoker_v1,famhisinsdiab,
                                         T1DGRS2_z,T2DGRS_z) %>% drop_na()
#DKA in multi-model
dka_model = predictions(
  rirs_step,
  newdata = datagrid(Study_ID = NA,
                     DKA = unique(df_UCPCR_clean$DKA), 
                     dur_diab_months_z = seq(min(complete_lmm$dur_diab_months_z),
                                             max(complete_lmm$dur_diab_months_z),
                                             length.out = 100)),
  re.form = NA)


model_dka = ggplot(dka_model, aes(x = dur_diab_months_z, y = estimate, 
                             color = DKA, fill = DKA,
                             ymin = conf.low, ymax = conf.high)) +
  geom_ribbon(alpha = .15, color = NA) + 
  geom_line(linewidth = 1) +
  xlab('Duration of Diabetes (scaled)') +
  ylab("log(UCPCR)") 
  ggsave('model_dka.png',model_dka, width = 6, height = 6)
  
#Gender
  gen_model = predictions(
    rirs_step,
    newdata = datagrid(Study_ID = NA,
                       Gender_v1 = unique(complete_lmm$Gender_v1), 
                       dur_diab_months_z = seq(min(complete_lmm$dur_diab_months_z),
                                               max(complete_lmm$dur_diab_months_z),
                                               length.out = 100)),
    re.form = NA)
  
  
  model_gen = ggplot(gen_model, aes(x = dur_diab_months_z, y = estimate, 
                                    color = Gender_v1, fill = Gender_v1,
                                    ymin = conf.low, ymax = conf.high)) +
    geom_ribbon(alpha = .15, color = NA) + 
    geom_line(linewidth = 1) +
    xlab('Duration of Diabetes (scaled)') +
    ylab("log(UCPCR)") +
    labs(color = "Gender", fill = "Gender")
  ggsave('model_gen.png',model_gen, width = 6, height = 6)
  
#Number of antibodies
  num_model = predictions(
    rirs_step,
    newdata = datagrid(Study_ID = NA,
                       num_anti = unique(complete_lmm$num_anti), 
                       dur_diab_months_z = seq(min(complete_lmm$dur_diab_months_z),
                                               max(complete_lmm$dur_diab_months_z),
                                               length.out = 100)),
    re.form = NA)
  
  
  model_num = ggplot(num_model, aes(x = dur_diab_months_z, y = estimate, 
                                    color = num_anti, fill = num_anti,
                                    ymin = conf.low, ymax = conf.high)) +
    geom_ribbon(alpha = .15, color = NA) + 
    geom_line(linewidth = 1) +
    xlab('Duration of Diabetes (scaled)') +
    ylab("log(UCPCR)") +
    labs(color = "Number of antibodies", fill = "Number of antibodies")
  ggsave('model_num.png',model_num, width = 6, height = 6)
  
  
#Visualizing multi-variables models
model_plot = patchwork::wrap_plots(
             model_num, 
             model_dka,
             model_gen,
             nrow = 2, ncol=2
) + patchwork::plot_annotation(tag_levels = 'a') 
ggsave('model_lmm.png',model_plot, width = 15, height = 6)






