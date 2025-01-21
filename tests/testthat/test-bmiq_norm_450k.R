library(testthat)

test_that("bmiq_norm_450k validates inputs correctly", {
  # Load required probes from probe_info_450k
  probe_info_450k <- get("probe_info_450k", envir = environment())
  required_probes <- probe_info_450k$probeID
  
  # Create mock data for valid input
  beta <- download_beta_values_case()
  beta <- beta[, 1, drop = FALSE]
  
  # Test valid input
  expect_error(bmiq_norm_450k(beta), NA)
  
  # Test 'beta' is not a data frame
  expect_error(
    bmiq_norm_450k(as.matrix(beta)),
    "'beta' must be a data frame."
  )
  
   # Test invalid 'cores' argument
  expect_error(
    bmiq_norm_450k(beta, cores = 0),
    "'cores' must be a positive integer."
  )
  expect_error(
    bmiq_norm_450k(beta, cores = 1.5),
    "'cores' must be a positive integer."
  )
  expect_error(
    bmiq_norm_450k(beta, cores = -1),
    "'cores' must be a positive integer."
  )
  
  # Test invalid 'verbose' argument
  expect_error(
    bmiq_norm_450k(beta, verbose = "yes"),
    "'verbose' must be a logical scalar."
  )
  expect_error(
    bmiq_norm_450k(beta, verbose = c(TRUE, FALSE)),
    "'verbose' must be a logical scalar."
  )
  
  # Test missing required probes
  beta_missing_probes <- beta[-1, , drop = FALSE] # Remove one required probe
  expect_error(
    bmiq_norm_450k(beta_missing_probes),
    "The following required probes are missing from 'beta':"
  )
})