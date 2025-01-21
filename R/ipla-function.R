#' Identify placental aging
#'
#' This function identifies placental aging based on the case-control aging 
#' difference. Placental aging is defined as the residual DNA-methylation-based 
#' (DNAm) gestational ages (GA). Only GA from 5 to 44 weeks' gestation are 
#' shown in the placental aging plot.
#'
#' @param aging A data frame of residual DNA-methylation-based GA. This data 
#' frame must be the output of \code{\link{plec}} function with argument 
#' `type="residual"`. The rows must be named according to the sample IDs.
#' @param ga A data frame of GA. There is only one column, i.e., `GA`, as shown 
#' in `ga`. Use \code{data(ga)} to load this data frame. The rows must be named 
#' according to the sample IDs.
#' @param phenotype A data frame of phenotype (optional. There is only one 
#' column, i.e., `phenotype`, as shown in `phenotype`. Use 
#' \code{data(phenotype)} to load this data frame. The rows must be named 
#' according to the sample IDs.
#' @param case A character of of case name in `phenotype` (default="Case").
#' @param control A character of of case name  in `phenotype` 
#' (default="Control").
#' @param method A character of of the method of statistical test (optional), 
#' i.e., "Mann-Whitney U" or "Permutation".
#' @param from An integer from 5 to 44 indicating minimum GA (weeks) to be 
#' included in the statistical test. If it is undefined, the minimum GA in 
#' either case or control is applied.
#' @param to An integer from 5 to 44 indicating maximum GA (weeks) to be 
#' included in the statistical test. If it is undefined, the maximum GA in 
#' either case or control is applied.
#'
#' @return An ggplot object consisting the aging plot without or with 
#' statistical test results.
#'
#' @keywords placental-epigenetic-clock placental-aging gestational-age
#'
#' @export
#'
#' @importFrom dplyr filter rename select slice
#' @importFrom ggplot2 aes coord_equal geom_errorbar geom_hline geom_label 
#' geom_line geom_ribbon ggplot scale_x_continuous scale_y_continuous theme 
#' theme_minimal
#' @importFrom purrr reduce
#' @importFrom stats aggregate approx qnorm sd wilcox.test
#'
#' @examples
#'
#' \dontrun{
#'   # Prepare data
#'   data(aging)
#'   data(ga)
#'   data(phenotype)
#'   
#'   # Identify placental aging
#'   ipla(aging, ga, phenotype)
#'   
#'   ## Conduct statistical test
#'   ipla(aging, ga, phenotype, method = "Mann-Whitney U")
#'   
#'   ## Conduct statistical test for a specific range of GA
#'   ipla(aging, ga, phenotype, method = "Mann-Whitney U", from = 5, to = 20)
#' }

ipla <- 
  function(
    aging
    , ga
    , phenotype
    , case = "Case"
    , control = "Control"
    , method = NULL
    , from = NULL
    , to = NULL
  ){
    
    # Core sub-function for bootstrapping interpolation
    bootstrap_interpolation <- function(x, y, common_x, n_bootstrap) {
      # Store bootstrapped interpolations
      bootstrap_results <-
        matrix(NA, nrow = n_bootstrap, ncol = length(common_x))
      
      for (i in 1:n_bootstrap) {
        # Resample x and y with replacement
        set.seed(seed + i)
        resample_indices <- sample(seq_along(x), replace = TRUE)
        x_resampled <- x[resample_indices]
        y_resampled <- y[resample_indices]
        
        # Aggregate duplicates by averaging
        unique_data <- aggregate(y_resampled ~ x_resampled, FUN = mean)
        colnames(unique_data) <- c("x", "y")  # Ensure column names are clear
        
        # Interpolate resampled y-values at common_x
        bootstrap_results[i, ] <-
          approx(unique_data$x, unique_data$y, xout = common_x, rule = 2)$y
      }
      
      return(bootstrap_results)
    }
    
    # Implement all sub-functions
    ## Input validation
    
    ### Check that aging is a data frame
    if (!is.data.frame(aging)) {
      stop("'aging' must be a data frame.")
    }
    
    ### Check that aging has row names
    if (is.null(rownames(aging)) || any(is.na(rownames(aging)))) {
      stop("'aging' must have row names representing sample IDs.")
    }
    
    ### Check that aging has a column named 'output'
    if (!"output" %in% colnames(aging)) {
      stop("'aging' must have a column named 'output'.")
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
    
    ### Check that phenotype is a data frame
    if (!is.null(phenotype)) {
      if (!is.data.frame(phenotype)) {
        stop("'phenotype' must be a data frame if provided.")
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
    
    ### Check case and control are valid strings
    if (!is.character(case) || length(case) != 1) {
      stop("'case' must be a single character string.")
    }
    if (!is.character(control) || length(control) != 1) {
      stop("'control' must be a single character string.")
    }
    
    ### Check method is valid if provided
    if (!is.null(method)) {
      if (!method %in% c("Mann-Whitney U", "Permutation")) {
        stop(
          paste0(
            "'method' must be either 'Mann-Whitney U' or 'Permutation' if "
            , "provided."
          )
        )
      }
    }
    
    ### Check from and to are valid integers
    if (!is.null(from)) {
      if (!is.numeric(from) || length(from) != 1 || from < 5 || from > 44) {
        stop("'from' must be a single integer between 5 and 44.")
      }
    }
    if (!is.null(to)) {
      if (!is.numeric(to) || length(to) != 1 || to < 5 || to > 44) {
        stop("'to' must be a single integer between 5 and 44.")
      }
    }
    
    ### Additional logic to ensure valid range for from and to
    if (!is.null(from) && !is.null(to) && from > to) {
      stop("'from' cannot be greater than 'to'.")
    }
    
    ## Main codes
    seed <- 2025-01-10
    
    aging_order <- match(rownames(ga), rownames(aging))
    aging <- rename(slice(aging, aging_order), aging = output)
    phenotype_order <- match(rownames(ga), rownames(phenotype))
    phenotype <- slice(phenotype, phenotype_order)
    data <- cbind(aging, select(ga, ga = GA), phenotype)
    
    ## Filter data into separate x and y for Control and Case
    x1 <- filter(data, phenotype == control)$ga
    y1 <- filter(data, phenotype == control)$aging
    x2 <- filter(data, phenotype == case)$ga
    y2 <- filter(data, phenotype == case)$aging
    
    ## Step 1: Define a common x-axis range for interpolation
    common_x <- seq(min(c(x1, x2)), max(c(x1, x2)), length.out = 100)
    
    ## Step 2a: Interpolate y-values for both lines, allowing extrapolation
    y1_interp <- approx(x1, y1, xout = common_x, rule = 2)$y
    y2_interp <- approx(x2, y2, xout = common_x, rule = 2)$y
    
    ## Step 2b: Bootstrapping function to generate y1 and y2 interpolations
    
    ### Number of bootstrap iterations
    n_bootstrap <- 30
    
    ### Perform bootstrapping for y1 and y2
    y1_bootstrap <- bootstrap_interpolation(x1, y1, common_x, n_bootstrap)
    y2_bootstrap <- bootstrap_interpolation(x2, y2, common_x, n_bootstrap)
    y_diff_bootstrap <- y2_bootstrap - y1_bootstrap
    
    ## Step 3: Calculate 95% CIs for each value in common_x
    ### Calculate mean and standard error for each value in common_x
    y1_mean <- y1_interp
    y2_mean <- y2_interp
    y_diff_mean <- y2_interp - y1_interp
    
    y1_sd <- apply(y1_bootstrap, 2, sd, na.rm = TRUE)
    y2_sd <- apply(y2_bootstrap, 2, sd, na.rm = TRUE)
    y_diff_sd <- apply(y_diff_bootstrap, 2, sd, na.rm = TRUE)
    
    ### Compute 95% CI using mean +/- z * (sd / sqrt(n))
    z_value <- qnorm(0.975)
    y1_ci_lower <- y1_mean - z_value * (y1_sd / sqrt(n_bootstrap))
    y1_ci_upper <- y1_mean + z_value * (y1_sd / sqrt(n_bootstrap))
    
    y2_ci_lower <- y2_mean - z_value * (y2_sd / sqrt(n_bootstrap))
    y2_ci_upper <- y2_mean + z_value * (y2_sd / sqrt(n_bootstrap))
    
    y_diff_ci_lower <- y_diff_mean - z_value * (y_diff_sd / sqrt(n_bootstrap))
    y_diff_ci_upper <- y_diff_mean + z_value * (y_diff_sd / sqrt(n_bootstrap))
    
    y1_ci <-
      rbind(matrix(y1_ci_lower, nrow = 1), matrix(y1_ci_upper, nrow = 1))
    y2_ci <-
      rbind(matrix(y1_ci_lower, nrow = 1), matrix(y1_ci_upper, nrow = 1))
    y_diff_ci <-
      rbind(matrix(y_diff_ci_lower, nrow = 1), matrix(y_diff_ci_upper, nrow = 1))
    
    ### Combine into a data frame for plotting
    bootstrap_results <- data.frame(
      ga = rep(common_x, 2),
      mean = c(y1_mean, y2_mean),
      lower_ci = c(y1_ci[1, ], y2_ci[1, ]),
      upper_ci = c(y1_ci[2, ], y2_ci[2, ]),
      phenotype = rep(c(control, case), each = length(common_x))
    )
    
    ## Step 4: Statistical Test (Mann-Whitney U test and permutation)
    
    ### Subset interpolated values within the min and max of common_x
    if(!is.null(method)){
      
      if(is.null(from)){
        range_min <- min(common_x)
      }else{
        range_min <- from
      }
      
      if(is.null(to)){
        range_max <- max(common_x)
      }else{
        range_max <- to
      }
      
      valid_idx <- common_x >= range_min & common_x <= range_max
      y1_interp_subset <- y1_interp[valid_idx]
      y2_interp_subset <- y2_interp[valid_idx]
      
      if(method == "Mann-Whitney U"){
        ### Mann-Whitney U Test within range
        wilcox_test <- 
          suppressWarnings(
            wilcox.test(
              abs(y1_interp_subset), abs(y2_interp_subset)
              , alternative = "two.sided"
            )
          )
        
        stat_test_obj <- wilcox_test
        stat_test <- "Mann-Whitney U Test"
        p_value <- wilcox_test$p.value
      }else{
        ### Permutation Test within range
        set.seed(seed)
        n_permutations <- 1000
        distance1 <- mean(abs(y1_interp_subset))
        distance2 <- mean(abs(y2_interp_subset))
        observed_diff <- distance1 - distance2
        combined_y <- c(y1_interp_subset, y2_interp_subset)
        perm_diffs <- replicate(n_permutations, {
          shuffled <- sample(combined_y)
          perm_group1 <- shuffled[1:length(y1_interp_subset)]
          perm_group2 <-
            shuffled[(length(y1_interp_subset) + 1):length(shuffled)]
          mean(abs(perm_group1)) - mean(abs(perm_group2))
        })
        p_value_perm <- mean(abs(perm_diffs) >= abs(observed_diff))
        
        stat_test_obj <-
          list(
            distance1 = distance1, distance2 = distance2, p.value = p_value_perm
          )
        stat_test <- "Permutation Test"
        p_value <- p_value_perm
      }
      
      p_value <-
        ifelse(
          p_value <= 0.001
          , "<0.001"
          , ifelse(
            p_value > 0.05
            , ">0.05"
            , paste0("=", format(p_value, digits = 3))
          )
        )
      
      ### Create data frame for statistical test results
      stat_results_data <-
        data.frame(
          ga = range_min + (range_max - range_min) / 2
          , ga_min = range_min
          , ga_max = range_max
          , aging_difference = 10
          , results = paste0(stat_test, ":\np-value", p_value)
        )
    }
    
    ## Step 5: Plot with 95% CI
    
    ### Compute the difference
    difference_data <-
      data.frame(
        ga = common_x
        , aging_difference = y2_interp - y1_interp
        , lower_ci = y_diff_ci[1, ]
        , upper_ci = y_diff_ci[2, ]
      )
    
    ### Plot the aging difference
    aging_difference_plot <-
      difference_data |>
      ggplot(aes(x = ga, y = aging_difference)) +
      geom_hline(yintercept = 0, lty = 2)
    
    if(!is.null(method)){
      aging_difference_plot <-
        aging_difference_plot +
        geom_errorbar(
          data = stat_results_data, aes(xmin = ga_min, xmax = ga_max)
        )
    }
    
    aging_difference_plot <-
      aging_difference_plot +
      geom_line(na.rm = TRUE) +
      geom_ribbon(
        aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, color = NA
      )
    
    if(!is.null(method)){
      aging_difference_plot <-
        aging_difference_plot +
        geom_label(data = stat_results_data, aes(label = results), size = 3)
    }
    
    aging_difference_plot <-
      aging_difference_plot +
      scale_x_continuous("Gestational age (weeks)", limits = c(5, 44)) +
      scale_y_continuous(
        "Case-control aging difference (weeks)"
        , limits = c(-10, 10)
      ) +
      theme_minimal()
    
    # Conclude the output
    output <- aging_difference_plot
    
    return(output)
  }
