#loading the package
library(rms)
library(pROC)

#Original model
original_model = glm(c_peptide_group ~ c_peptide_v1_z + Gender_v1 + DKA + IA2_bin_v1 +
      GAD_bin_v1 + GAD_v1_z, family = 'binomial',
      data = df_define_filtered)
pred_origin = predict(original_model,
                      df_define_filtered,
                      type = 'response')

#ROC AUC for original model
roc_obj_original =  roc(df_define_filtered$c_peptide_group, pred_origin)
auc_value_origin = auc(roc_obj_original)

# Logistic Regression Model
fit <- lrm(c_peptide_group ~ c_peptide_v1_z + Gender_v1 + DKA + IA2_bin_v1 +
             GAD_bin_v1 + GAD_v1_z, data = df_define_filtered,
             x = TRUE, y = TRUE)

# Setting a seed 
set.seed(123)
# Running with 200 bootstrap samples
val_results <- validate(fit, method = "boot", B = 200)
print(val_results)

#Calculating ROC AUC from dxy
corrected_dxy = val_results["Dxy", "index.corrected"]

# 4. Convert it to optimism-corrected ROC AUC
corrected_auc = (corrected_dxy + 1) / 2
print(corrected_auc)

#Boot strapping using resample
B = 200
auc_values_train = numeric(B)
auc_values_test = numeric(B)

for (i in 1:B){
  # Bootstrap sample
   boot_idx = sample(1:nrow(df_define_filtered),
              size = nrow(df_define_filtered),
              replace = TRUE)
   boot_data = df_define_filtered[boot_idx, ]
   
   #Fitting logistic regressionoutcome
   model = glm(c_peptide_group ~ c_peptide_v1_z + Gender_v1 + DKA + IA2_bin_v1 +
                      GAD_bin_v1 + GAD_v1_z, family = 'binomial',
                     data = boot_data)
   #Predicted probability (sampled model and sampled data)
   pred_prob_train = predict(model,
                             boot_data,
                             type = 'response')
   #Predicted probability (sampled model and original data)
   pred_prob_test =  predict(model,
                            df_define_filtered,
                            type = 'response')
   # ROC and AUC
   roc_obj_train =  roc(boot_data$c_peptide_group,pred_prob_train)
   roc_obj_test =  roc(df_define_filtered$c_peptide_group,pred_prob_test)
  
   auc_values_train[i] = auc(roc_obj_train)
   auc_values_test[i] = auc(roc_obj_test)
   
   optimism = mean(auc_values_train) - mean(auc_values_test)
    }

#Corrected Model AUC = Original Model AUC - Optimism
Corrected_AUC = auc_value_origin - optimism #0.77
#Average Optimism value = 0.012

