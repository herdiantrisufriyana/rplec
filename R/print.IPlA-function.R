#' Print method for an IPlA object
#'
#' This function provides a custom print method for an object of class `IPlA`. 
#' It is the output of \code{\link{ipla}} function. When this object is 
#' printed, it automatically displays the plot stored in the `plot` component 
#' of the object.
#'
#' @param x An object of class `IPlA`.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns `x`.
#'
#' @keywords placental-epigenetic-clock, placental-aging, gestational-age
#'
#' @export
#'
#' @examples
#'
#' # Prepare data
#' data(aging)
#' data(ga)
#' data(phenotype)
#' 
#' # Identify placental aging
#' ipla_results <- ipla(aging, ga, phenotype)
#' 
#' # Displays the plot
#' ipla_results

print.IPlA <- function(x, ...) {
  x$plot
}