#' Estimate DNA-methylation-based gestational age
#'
#' This function estimate gestational age (GA) using BMIQ-normalized beta 
#' values. The estimated GA is a sum of normal and residual GAs. The latter 
#' is a sum of condition- and trimester-specific, residual GAs.
#'
#' @param norm_beta A data frame of normalized beta values where each column 
#' represents a sample and each row represent a probe. This data frame must be 
#' the output of \code{\link{bmiq_norm_450k}} function. The rows must be named 
#' according to the probe IDs. Meanwhile, the columns must be named according to 
#' the sample IDs.
#' @param type An character indicating the type of outputs which are primarily: 
#' (1) "stack" (default) for the estimated GA; (2) "normal" for the estimated 
#' normal GA; (3) "residual" for the estimated residual GA; (4) 
#' "condition" for the condition-specific, estimated residual GA; and (5) 
#' "trimester" for the trimester-specific, estimated residual GA. In addition, 
#' a user can obtain the output of a single submodel using the column name 
#' (except `predictor`) in `plec_int_coef`. Use \code{data(plec_int_coef)} to 
#' load this data frame.
#' @param verbose A logical scalar indicating whether to show a progress bar.
#'
#' @return A data frame of the estimated GA.
#'
#' @keywords placental-epigenetic-clock, predict, gestational-age
#'
#' @export
#'
#' @importFrom dplyr filter mutate rename_at
#' @importFrom pbapply pblapply
#' @importFrom purrr imap reduce
#' @importFrom stringr str_detect
#' @importFrom tibble column_to_rownames
#' @importFrom utils capture.output data
#'
#' @examples
#'
#' \dontrun{
#'   data(beta_values_case)
#'   norm_beta_values_case <- bmiq_norm_450k(beta_values_case)
#'   dnam_ga_case <- plec(norm_beta_values_case)
#' }

plec <- function(norm_beta, type = "stack", verbose = TRUE){
  
  # Wrapper sub-function for running 1 sample (message output is hidden)
  plec_df1 <- function(x, data, ...){
    suppressMessages(invisible(capture.output({
      results <-
        data[, x, drop = TRUE] |>
        `names<-`(rownames(data)) |>
        plec_req(...)
    })))
    
    results
  }
  
  # Core sub-function for the requested placental epigenetic clock
  plec_req <- function(norm_beta, type, sampleID){
    
    intercepts <-
      plec_int_coef |>
      filter(predictor == "(Intercept)")
    
    cpg_coefs <-
      plec_int_coef |>
      filter(str_detect(predictor, "^cg"))
    
    non_cpg_coefs <-
      plec_int_coef |>
      filter(predictor != "(Intercept)" & !str_detect(predictor, "^cg"))
    
    cpg_scaler_mean <-
      plec_scaler_mean |>
      filter(str_detect(predictor, "^cg"))
    
    non_cpg_scaler_mean <-
      plec_scaler_mean |>
      filter(!str_detect(predictor, "^cg"))
    
    cpg_scaler_scale <-
      plec_scaler_scale |>
      filter(str_detect(predictor, "^cg"))
    
    non_cpg_scaler_scale <-
      plec_scaler_scale |>
      filter(!str_detect(predictor, "^cg"))
    
    if(type %in% c("ga_est", "stack", "normal", "residual", "trimester")){
      ga_est <-
        plec_vec(
          norm_beta, intercepts, cpg_coefs, cpg_scaler_mean, cpg_scaler_scale
          , "ga_est"
        )
    }
    
    conditions <-
      c("fgr", "pe", "pe_onset", "preterm", "anencephaly", "spina_bifida"
        , "gdm", "diandric_triploid", "miscarriage", "lga", "subfertility"
        , "hellp", "chorioamnionitis"
      )
    
    pred_prob <- list()
    ga_res_conds_est <- list()
    
    for(condition in conditions){
      if(
        type
        %in% c(
          paste0(condition, "_pred")
          , "ga_res_conds_pred_est"
          , "stack", "residual", "condition", "trimester"
        )
      ){
        pred_prob[[condition]] <-
          plec_vec(
            norm_beta, intercepts, cpg_coefs, cpg_scaler_mean, cpg_scaler_scale
            , paste0(condition, "_pred")
          )
        }
      
      if(
        type
        %in% c(
          paste0("ga_res_conds_", condition, "_est")
          , "ga_res_conds_pred_est"
          , "stack", "residual", "condition", "trimester"
        )
        ){
          ga_res_conds_est[[condition]] <-
            plec_vec(
              norm_beta, intercepts, cpg_coefs
              , cpg_scaler_mean, cpg_scaler_scale
              , paste0("ga_res_conds_", condition, "_est")
            )
        }
    }
    
    if(
      type
      %in% c(
        "ga_res_conds_pred_est", "stack", "residual", "condition", "trimester"
      )
    ){
      pred_prob <-
        pred_prob |>
        imap(
          ~ data.frame(cond = .y, prob = .x) |>
          column_to_rownames(var = "cond")
        ) |>
        reduce(rbind)
      
      ga_res_conds_est <-
        ga_res_conds_est |>
        imap(
          ~ data.frame(cond = .y, est = .x) |>
          column_to_rownames(var = "cond")
        ) |>
        reduce(rbind)
      
      ga_res_conds_pred_est <-
        cbind(pred_prob, ga_res_conds_est) |>
        mutate(residual = prob * est)
      
      ga_res_conds_pred_est <-
        ga_res_conds_pred_est$residual |>
        `names<-`(rownames(ga_res_conds_pred_est)) |>
        plec_vec(
          intercepts, non_cpg_coefs, non_cpg_scaler_mean, non_cpg_scaler_scale
          , "ga_res_conds_pred_est"
        )
    }
    
    if(type %in% c("ga_res_comb_pr_est", "stack", "residual", "trimester")){
      ga_res_comb_pr_est <-
        plec_vec(
          norm_beta, intercepts, cpg_coefs, cpg_scaler_mean, cpg_scaler_scale
          , "ga_res_comb_pr_est"
        )
    }
    
    if(type %in% c("ga_res_comb_tb_est", "stack", "residual", "trimester")){
      ga_res_comb_tb_est <-
        plec_vec(
          norm_beta, intercepts, cpg_coefs, cpg_scaler_mean, cpg_scaler_scale
          , "ga_res_comb_tb_est"
        )
    }
    
    if(type %in% c("ga_res_comb_ta_est", "stack", "residual", "trimester")){
      ga_res_comb_ta_est <-
        plec_vec(
          norm_beta, intercepts, cpg_coefs, cpg_scaler_mean, cpg_scaler_scale
          , "ga_res_comb_ta_est"
        )
    }
    
    if(type %in% c("stack", "residual", "trimester")){
      ga_res_comb_est <-
        ga_est + ga_res_conds_pred_est
      
      ga_res_comb_est <-
        ifelse(
          ceiling(ga_res_comb_est) <= 36
          , ifelse(is.na(ga_res_comb_pr_est), 0, ga_res_comb_pr_est)
          , ifelse(
            ceiling(ga_res_comb_est) >= 37 & ceiling(ga_res_comb_est) <= 40
            , ifelse(is.na(ga_res_comb_tb_est), 0, ga_res_comb_tb_est)
            , ifelse(is.na(ga_res_comb_ta_est), 0, ga_res_comb_ta_est)
          )
        )
    }
    
    if(type %in% c("stack", "residual")){
      ga_residual <- ga_res_conds_pred_est + ga_res_comb_est
    }
    
    if(type == "normal"){
      output <- ga_est
    }else if(type == "residual"){
      output <- ga_residual
    }else if(type %in% c("condition", "ga_res_conds_pred_est")){
      output <- ga_res_conds_pred_est
    }else if(type == "trimester"){
      output <- ga_res_comb_est
    }else if(type == "ga_est"){
      output <- ga_est
    }else if(str_detect(type, "_pred$")){
      for(condition in conditions){
        if(type == paste0(condition, "_pred")) output <- pred_prob[[condition]]
      }
    }else if(str_detect(type, "^ga_res_conds_")){
      for(condition in conditions){
        if(type == paste0("ga_res_conds_", condition, "_est")){
          output <- ga_res_conds_est[[condition]]
        }
      }
    }else if(type == "ga_res_comb_pr_est"){
      output <- ga_res_comb_pr_est
    }else if(type == "ga_res_comb_tb_est"){
      output <- ga_res_comb_tb_est
    }else if(type == "ga_res_comb_ta_est"){
      output <- ga_res_comb_ta_est
    }else{
      output <- ga_est + ga_residual
    }
    
    output <-
      data.frame(output = output) |>
      rename_at("output", \(x) sampleID)
    
    return(output)
  }
  
  # Core sub-function for individual placental epigenetic clock
  plec_vec <- function(data, int, coefs, scaler_mean, scaler_scale, colname){
    
    data <- data[match(coefs$predictor, names(data))]
    scaled_data <- data - scaler_mean[, colname, drop = TRUE]
    scaled_data <- scaled_data / scaler_scale[, colname, drop = TRUE]
    output <- int[, colname, drop = TRUE]
    output <- output + sum(scaled_data * coefs[, colname, drop = TRUE])
    
    if(str_detect(colname, "_pred$")) output <- 1 / (1 + exp(-output))
    
    return(output)
  }
  
  # Implement all sub-functions
  data("plec_int_coef")
  data("plec_scaler_mean")
  data("plec_scaler_scale")
  
  if(verbose){
    looping_fn <- pblapply
  }else{
    looping_fn <- lapply
  }
  
  colnames(norm_beta) |>
    looping_fn(
      \(x)
      plec_df1(
        x = x
        , data = norm_beta
        , type = type
        , sampleID = x
      )
    ) |>
    reduce(cbind) |>
    `rownames<-`("output") |>
    t() |>
    as.data.frame()
  
}
