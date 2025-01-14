#' Probe info for 450K
#'
#' A list of 450K probe information for \code{\link{bmiq_norm}} function. The 
#' probes are already filtered and ordered the same way to the input when we 
#' developed our placental epigenetic clock.
#'
#' @format A list of 2 elements where each has a length of 346407:
#' \describe{
#'   \item{Design}{An integer indicating design type 1 or 2.}
#'   \item{probeID}{A character for each probe identifier.}
#' }
#'
#' @source Derived from `ChAMP` R package.
#' @keywords dataset
#' @name probe_info_450k
'probe_info_450k'
