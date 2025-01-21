#' Perform quality control
#'
#' This function evaluates the precision of DNA-methylation-based (DNAm)
#' gestational age (GA) based on calibration, root mean square error (RMSE), 
#' mean absolute error (MAE), and Pearson's correlation coefficient (r). The 
#' sample identifiers (IDs) are automatically matched among the DNAm-GA, GA, 
#' and phenotype (optional). Only GA from 5 to 44 weeks' gestation are shown in 
#' the calibration plot.
#'
#' @param dnam_ga A data frame of DNA-methylation-based GA. This data frame 
#' must be the output of \code{\link{plec}} function. The rows must be named 
#' according to the sample IDs.
#' @param ga A data frame of GA. There is only one column, i.e., `GA`, as shown 
#' in `ga`. Use \code{data(ga)} to load this data frame. The rows must be named 
#' according to the sample IDs.
#' @param phenotype A data frame of phenotype (optional). There is only one 
#' column, i.e., `phenotype`, as shown in `phenotype`. Use 
#' \code{data(phenotype)} to load this data frame. The rows must be named 
#' according to the sample IDs.
#'
#' @return A ggplot object of calibration plot with RMSE, MAE, and r.
#'
#' @keywords placental-epigenetic-clock, quality-control, gestational-age
#'
#' @export
#'
#' @importFrom dplyr group_by mutate mutate_at n pull rename select slice 
#' summarize
#' @importFrom ggplot2 aes annotate coord_equal element_blank geom_abline 
#' geom_point ggplot scale_x_continuous scale_y_continuous theme theme_minimal
#' @importFrom purrr reduce
#' @importFrom stats cor.test qnorm sd
#' @importFrom tidyr gather unite
#'
#' @examples
#'
#' \dontrun{
#'   data(beta_values_case)
#'   norm_beta_values_case <- bmiq_norm_450k(beta_values_case)
#'   dnam_ga_case <- plec(norm_beta_values_case)
#'   
#'   data(ga)
#'   ga_case <- ga[phenotype$phenotype == "Case", , drop = FALSE]
#'   qc(dnam_ga_case, ga_case)
#' }

qc <- function(dnam_ga, ga, phenotype = NULL){
  
  # Core sub-function for RMSE
  rmse <- function(actual, predicted, na.rm = FALSE) {
    sqrt(mean((predicted - actual) ^ 2, na.rm = na.rm))
  }
  
  # Core sub-function for MAE
  mae <- function(actual, predicted, na.rm = FALSE) {
    mean(abs(predicted - actual), na.rm = na.rm)
  }
  
  # Implement all sub-functions
  ## Input validation
  
  ### Check that dnam_ga is a data frame
  if (!is.data.frame(dnam_ga)) {
    stop("'dnam_ga' must be a data frame.")
  }
  
  ### Check that dnam_ga has row names
  if (is.null(rownames(dnam_ga)) || any(is.na(rownames(dnam_ga)))) {
    stop("'dnam_ga' must have row names representing sample IDs.")
  }
  
  ### Check that dnam_ga has a column named 'output'
  if (!"output" %in% colnames(dnam_ga)) {
    stop("'dnam_ga' must have a column named 'output'.")
  }
  
  ### Check that ga is a data frame
  if (!is.data.frame(ga)) {
    stop("'ga' must be a data frame.")
  }
  
  ### Check that ga has row names
  if (is.null(rownames(ga)) || any(is.na(rownames(ga)))) {
    stop("'ga' must have row names representing sample IDs.")
  }
  
  ### Check that ga has a column named 'GA'
  if (!"GA" %in% colnames(ga)) {
    stop("'ga' must have a column named 'GA'.")
  }
  
  ### If phenotype is provided, validate it
  if (!is.null(phenotype)) {
    
    #### Check that phenotype is a data frame
    if (!is.data.frame(phenotype)) {
      stop("'phenotype' must be a data frame.")
    }
    
    #### Check that phenotype has row names
    if (is.null(rownames(phenotype)) || any(is.na(rownames(phenotype)))) {
      stop("'phenotype' must have row names representing sample IDs.")
    }
    
    #### Check that phenotype has a column named 'phenotype'
    if (!"phenotype" %in% colnames(phenotype)) {
      stop("'phenotype' must have a column named 'phenotype'.")
    }
  }
  
  ## Main codes
  dnam_ga_order <- match(rownames(ga), rownames(dnam_ga))
  dnam_ga <- rename(slice(dnam_ga, dnam_ga_order), Yhat = output)
  
  if(!is.null(phenotype)){
    phenotype_order <- match(rownames(ga), rownames(phenotype))
    phenotype <- slice(phenotype, phenotype_order)
    
    data <- cbind(dnam_ga, select(ga, Y = GA), phenotype)
  }else{
    data <- cbind(dnam_ga, select(ga, Y = GA))
  }
  
  seed <- 2025-01-10
  boot_eval <- list()
  
  for(b in seq(30)){
    set.seed(seed + b)
    
    i <-
      data |>
      nrow() |>
      seq() |>
      sample(size = 10 * nrow(data), replace = TRUE)
    
    boot_eval[[b]] <-
      data[i, , drop = FALSE] |>
      summarize(
        RMSE = rmse(Y, Yhat, na.rm =TRUE)
        ,MAE = mae(Y, Yhat, na.rm =TRUE)
        ,r = cor.test(Y, Yhat)$estimate
        ,.groups = "drop"
      )
  }
  
  boot_eval <-
    boot_eval |>
    reduce(rbind) |>
    mutate_at(c("RMSE", "MAE", "r"), as.numeric) |>
    gather(metric, value) |>
    mutate_at("metric", \(x) factor(x, unique(x))) |>
    group_by(metric) |>
    summarize(
      avg = mean(value)
      ,lb = mean(value) - qnorm(0.975) * sd(value) / sqrt(n())
      ,ub = mean(value) + qnorm(0.975) * sd(value) / sqrt(n())
      ,.groups = "drop"
    ) |>
    mutate_at(c("avg", "lb", "ub"), format, 3, digits = 3) |>
    unite(ci, lb, ub, sep = ", ") |>
    mutate(ci = paste0("(", ci, ")")) |>
    unite(estimates, avg, ci, sep = " ") |>
    unite(report, metric, estimates, sep = " = ") |>
    pull(report) |>
    paste0(collapse = "\n")
  
  if(!is.null(phenotype)){
    output <-
      data |>
      ggplot(aes(Y, Yhat, color = phenotype))
  }else{
    output <-
      data |>
      ggplot(aes(Y, Yhat))
  }
  
  output <-
    output +
    geom_abline(intercept = 0, slope = 1, lty = 2) +
    geom_point(na.rm = TRUE) +
    annotate(
      geom = "label"
      , x = 5
      , y = 44
      , label = boot_eval
      , size = 3, hjust = 0, vjust = 1
    ) +
    coord_equal() +
    scale_x_continuous("Gestational age (weeks)", limits = c(5, 44)) +
    scale_y_continuous(
      "DNA-methylation-based\ngestational age (weeks)"
      , limits = c(5, 44)
    ) +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  return(output)
}
