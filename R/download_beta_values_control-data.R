#' Access Beta values at individual CpG sites for control group
#'
#' Downloads and loads beta values for the control group. Data contains beta 
#' values for 5 samples and 452,453 probes.
#'
#' @return A data frame with 452,453 rows and 5 columns:
#' \describe{
#'   \item{GSM7115144}{Beta values for sample GSM7115144.}
#'   \item{GSM1702248}{Beta values for sample GSM1702248.}
#'   \item{GSM3179749}{Beta values for sample GSM3179749.}
#'   \item{GSM4281756}{Beta values for sample GSM4281756.}
#'   \item{GSM5210472}{Beta values for sample GSM5210472.}
#' }
#' @source Derived from the 2024 Placental Clock DREAM Challenge.
#' @keywords dataset
#' @export
#'
#' @examples
#' \dontrun{
#'   beta_values_control <- download_beta_values_control()
#'   head(beta_values_control)
#' }

download_beta_values_control <- function() {
  url <-
    paste0(
      "https://raw.githubusercontent.com/herdiantrisufriyana/rplec/master/data/"
      , "beta_values_control.rda"
    )
  temp <- tempfile(fileext = ".rda")
  utils::download.file(url, temp, mode = "wb")
  load(temp)
  unlink(temp)  # Clean up the temporary file
  return(beta_values_control)
}