#' Calibrate risk 
#'
#' Calculate the re-calibrated risk within each group.
#'
#' @param data a data.frame containing the risk score, observed response and group label
#' @param covariate a data.frame containing the covariates Z
#' @param controlvar the binary indicator of being in control set
#' @param validatevar the binary indicator of being in validation set
#' @param groupvar the group label, categorical
#' @param response the observed outcome or response variable, should be binary 0/1
#' @param risk the (original) risk score under evaluation, scalar
#' @param transform apply a logit transformation to the risk score, converts probability to a linear predictor
#' @param method method to be used to calibrate the risks, options: 'logit' for logistic regression
#' @param quietly suppress output messages, default = TRUE
#'
#' @return list with data.frame with the calibrated risk results and MSPE if cv option selected
#' 
#' @import dplyr
#' @importFrom stats glm predict poly
#'
#' @export


calibrateRisk_new_validate <- function(data 
                            , covariate
                            , controlvar
                            , validatevar
                            , groupvar
                            , response
                            , risk
                            , transform = TRUE
                            , method = "logit"
                            , quietly = TRUE
){
  
  # construct dataset
  s <- data %>% dplyr::select({{groupvar}}); names(s) ='s'
  y <- data %>% dplyr::select({{response}}); names(y) = 'y'
  gX <- data %>% dplyr::select({{risk}}); names(gX) = 'gX'
  control <- data %>% dplyr::select({{controlvar}}); names(control) = 'control'
  validate <- data %>% dplyr::select({{validatevar}}); names(validate) = 'validate'
  
  if(transform == TRUE){
    # get linear predictor version of risk using logit transfomation
    lp.gX <- logit(gX)   ; names(lp.gX) = 'lp.gX'
  }else{
    lp.gX <- gX; names(lp.gX) = 'lp.gX'
  }
  
  df <- data.frame(s, y, gX, lp.gX, control, validate)
  
  # Determine group levels 
  slist <-  df %>%
    dplyr::arrange(.data$s) %>%
    dplyr::pull(.data$s) %>%
    unique()
  
  ns = length(slist)
  
  if(quietly != TRUE){
    cat("There are ",ns," groups: ",slist,"\n")
  }
  
  # loop over groups
  for(i.s in 1:ns){
    
    # Get group subset
    df.s <- df %>% dplyr::filter(.data$s == slist[i.s])
    validate <- df[df$validate == 1, ]
    covariate.validate <- covariate[df$validate == 1, ]
    covariate.s <- covariate[df$s == slist[i.s],]
    # Calibrate risk according to method
    df.s <-  calibrate_new_validate(train = df.s
                           , covariate=covariate.s
                           , validate=validate
                           , covariate.validate=covariate.validate
                           , method='logit')
    
    # stack results by s level
    if(i.s == 1){
      caldf <- df.s
    }else{
      caldf <- caldf %>%
        dplyr::bind_rows(df.s)
    }
  } # end loop over s
  
  out <- list(caldf= caldf)
  
  return(out)
}
