## ----Load necessary packages--------------------------------------------------
library(rplec)

## ----Load our example data----------------------------------------------------
beta_values_case <- download_beta_values_case()
beta_values_control <- download_beta_values_control()
data(ga)
data(phenotype)

## ----Normalize DNA methylation values per sample------------------------------
norm_beta_values_case <- bmiq_norm_450k(beta_values_case)
norm_beta_values_control <- bmiq_norm_450k(beta_values_control)

## ----Estimate DNA-methylation-based gestational age---------------------------
dnam_ga_case <- plec(norm_beta_values_case)
dnam_ga_control <- plec(norm_beta_values_control)
dnam_ga <- rbind(dnam_ga_case, dnam_ga_control)

## ----Perform quality control--------------------------------------------------
plec_qc <- qc(dnam_ga, ga, phenotype)

## ----figure-1, echo=FALSE, fig.height=5, fig.width=6--------------------------
plec_qc

## ----Estimate placental aging for either case or control----------------------
aging_case <- plec(norm_beta_values_case, type = "residual")
aging_control <- plec(norm_beta_values_control, type = "residual")
aging <- rbind(aging_case, aging_control)

## ----Compare case and control to identify placental aging---------------------
ipla_results <- ipla(aging, ga, phenotype)

## ----figure-2, echo=FALSE, fig.height=5, fig.width=7--------------------------
ipla_results

## ----Conduct statistical test-------------------------------------------------
ipla_stats <- ipla(aging, ga, phenotype, method = "Mann-Whitney U")

## ----figure-3, echo=FALSE, fig.height=5, fig.width=7--------------------------
ipla_stats

## ----Conduct statistical test for a specific range of GA----------------------
ipla_stats_5_20 <-
  ipla(aging, ga, phenotype, method = "Mann-Whitney U", from = 5, to = 20)

## ----figure-4, echo=FALSE, fig.height=5, fig.width=7--------------------------
ipla_stats_5_20

## ----session-info, echo=TRUE--------------------------------------------------
sessionInfo()

