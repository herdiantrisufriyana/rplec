# Use this to avoid R CMD check notes about undefined global variables
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "x", "probe_info_450k", "probe_id", "GA", "Y", "Yhat", 
    "aging_difference", "metric", "value", "ci", "lb", "ub",
    "plec_int_coef", "plec_scaler_mean", "plec_scaler_scale",
    "predictor", "prob", "est", "avg", "report", "results", "ga_min",
    "ga_max", "lower_ci", "upper_ci", "estimates", "beta_values_case",
    "beta_values_control"
  ))
}