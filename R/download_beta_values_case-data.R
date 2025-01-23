#' Access Beta values at individual CpG sites for case group
#'
#' Downloads and loads beta values for the case group. Data contains beta 
#' values for 5 samples and 452,453 probes.
#'
#' @return A data frame with 452,453 rows and 5 columns:
#' \describe{
#'   \item{GSM1931565}{Beta values for sample GSM1931565.}
#'   \item{GSM5114811}{Beta values for sample GSM5114811.}
#'   \item{GSM2589558}{Beta values for sample GSM2589558.}
#'   \item{GSM1842848}{Beta values for sample GSM1842848.}
#'   \item{GSM1843045}{Beta values for sample GSM1843045.}
#' }
#' @source Derived from the 2024 Placental Clock DREAM Challenge.
#' @keywords dataset
#' @export
#'
#' @examples
#' 
#' beta_values_case <- download_beta_values_case()
#' head(beta_values_case)

download_beta_values_case <- function() {
  url <-
    paste0(
      "https://raw.githubusercontent.com/herdiantrisufriyana/rplec/master/data/"
      , "beta_values_case.rda"
    )
  temp <- tempfile(fileext = ".rda")
  utils::download.file(url, temp, mode = "wb")
  load(temp)
  unlink(temp)  # Clean up the temporary file
  return(beta_values_case)
}