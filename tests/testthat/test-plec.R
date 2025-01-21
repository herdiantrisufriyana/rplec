library(testthat)

test_that("plec validates inputs correctly", {
  # Create mock data for valid input
  beta <- download_beta_values_case()
  beta <- beta[, 1, drop = FALSE]
  norm_beta <- bmiq_norm_450k(beta)
  
  valid_types <- c(
    "stack", "normal", "residual", "condition", "trimester",
    "ga_est", "ga_res_comb_pr_est", "ga_res_comb_tb_est", "ga_res_comb_ta_est",
    "ga_res_conds_pred_est", 
    paste0(
      "ga_res_conds_",
      c("fgr", "pe", "pe_onset", "preterm", "anencephaly",
        "spina_bifida", "gdm", "diandric_triploid", "miscarriage",
        "lga", "subfertility", "hellp", "chorioamnionitis"
      ),
      "_est"
    ),
    paste0(
      c("fgr", "pe", "pe_onset", "preterm", "anencephaly", "spina_bifida",
        "gdm", "diandric_triploid", "miscarriage", "lga", "subfertility", 
        "hellp", "chorioamnionitis"
      ),
      "_pred"
    )
  )
  
  # Test valid input
  expect_error(plec(norm_beta, type = "stack"), NA)
  
  # Test norm_beta is not a data frame
  expect_error(
    plec(as.matrix(norm_beta), type = "stack"),
    "'norm_beta' must be a data frame."
  )
  
  # Test invalid 'type' argument
  expect_error(
    plec(norm_beta, type = "invalid_type"),
    "'type' must be one of:"
  )
  
  # Test valid 'type' argument
  for (valid_type in valid_types) {
    expect_error(plec(norm_beta, type = valid_type), NA)
  }
  
  # Test verbose is not a logical scalar
  expect_error(
    plec(norm_beta, type = "stack", verbose = "yes"),
    "'verbose' must be a logical scalar."
  )
  expect_error(
    plec(norm_beta, type = "stack", verbose = c(TRUE, FALSE)),
    "'verbose' must be a logical scalar."
  )
})