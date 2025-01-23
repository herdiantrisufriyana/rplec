#' Indices for random beta values
#'
#' A list of indices for `rand.idx` in \code{\link{bmiq_norm_450k}} function. 
#' The input probes must be already filtered and ordered the same way to the 
#' that when we developed our placental epigenetic clock. Run 
#' \code{data(probe_info_450k)} and find the required probes in 
#' `prob_info_450k$probeID`.
#'
#' @format A list of 2 elements where each has a length of nfit of 10000:
#' \describe{
#'   \item{beta1.v}{An integer indicating the selected indices.}
#'   \item{beta2.v}{An integer indicating the selected indices.}
#' }
#'
#' @source Derived from `ChAMP` R package.
#' @keywords dataset
#' @name beta_v_indices
'beta_v_indices'
