library(performance)
library(Metrics)

#Diagnostics
qqnorm(residuals(rirs_step))
qqline(residuals(rirs_step))

diag_lmm = plot(fitted(rirs_step), resid(rirs_step))
abline(h = 0)

# QQ plot
png("qqplot_rirs.png", width = 800, height = 800, res = 200)
qqnorm(residuals(rirs_step))
qqline(residuals(rirs_step))
dev.off()

# Residuals vs fitted
png("resid_fitted_rirs.png", width = 800, height = 800, res = 200)
plot(fitted(rirs_step), resid(rirs_step),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)
dev.off()


#Performance
# Marginal & conditional R² 
r2_vals = r2(rirs_step)
r2_vals
r2_vals$R2_marginal    #0.25
r2_vals$R2_conditional #0.74

#Using bootmer for parametric sampling
metrics <- function(model) {
  r2_out <- r2(model)
  c(R2 = r2_out$R2_conditional)
}
boo_r2 <- bootMer(rirs_step, FUN = metrics, nsim = 200, use.u = FALSE, type = "parametric")
boo_r2
confint(boo_r2)







# RMSE and MAE 
df_fitted = model.frame(rirs_step)
df_fitted$scaled_log_ucpcr = df_fitted$`scale(log(UCPCR))`

rmse(
  df_fitted$scaled_log_ucpcr,
  df_fitted$pred_UCPCR
) #0.38

mae(df_fitted$scaled_log_ucpcr,
    df_fitted$pred_UCPCR
) #0.38

