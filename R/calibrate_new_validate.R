#' calibrate
#'
#' Inner function for re-calibrated risk within a group.
#'
#' @param train a data.frame containing the risk score, observed response and group label
#' @param covariate a data.frame containing the covariates Z
#' @param method method to be used to calibrate the risks, options: 'logit' for logistic regression
#' 
#' @return data.frame with the calibrated risk results
#' 
#' @import dplyr
#' @importFrom stats glm predict poly
#'
#' @export

calibrate_new_validate <- function(train
                      ,covariate
                      ,validate
                      ,covariate.validate
                      ,method
){
  
  if(tolower(method) == 'logit'){
    Y=train$y
    p=dim(covariate)[2]
    X <- covariate
    gx=train$gX
    logit.gx=train$lp.gX
    current_s <- train$s[1]
    
    total.validate=cbind(validate,covariate.validate)
    #pheno1=glm(y~ as.factor(s)+logit.gX+Z1+Z2+Z3,family=binomial, data=total.validate)
    pheno1=glm(y~ as.factor(s)+Z1+Z2+Z3,family=binomial, data=total.validate)
    
    tau_hat=pheno1$coefficients
    
    sample1=which(train$control==1)
    Y_N0=Y[sample1]
    X_N0=as.matrix(X[sample1,])
    gx_N0=gx[sample1]
    logit.gx_N0=logit.gx[sample1]
    
    sample2=which(train$validate==1)
    Y_n1r1=Y[sample2]
    X_n1r1=as.matrix(X[sample2,])
    gx_n1r1=gx[sample2]
    logit.gx_n1r1=logit.gx[sample2]
    
    sample3=setdiff(1:length(train$y), union(sample1, sample2))#unknown
    Y_n1r0=Y[sample3]
    X_n1r0=as.matrix(X[sample3,])
    gx_n1r0=gx[sample3]
    logit.gx_n1r0=logit.gx[sample3]
    
    #perform phenotyping on the validated data
    #tau_hat=c(1,rep(0.5,p))
    
    
    objective_function <- function(params) {
      t0 <- params[1]
      t1 <- params[2]
      a  <- params[3]
      b  <- params[4]
      
      # if (current_s == 1) {
      #   Y_n1r1_estimate <- hf(a + b * (tau_hat[1] + logit.gx_n1r1 * tau_hat[3] + X_n1r1 %*% tau_hat[4:6]))
      #   Y_n1r0_estimate <- hf(a + b * (tau_hat[1] + logit.gx_n1r0 * tau_hat[3] + X_n1r0 %*% tau_hat[4:6]))
      # } else if (current_s == 2) {
      #   Y_n1r1_estimate <- hf(a + b * (tau_hat[1] + tau_hat[2] + logit.gx_n1r1 * tau_hat[3] + X_n1r1 %*% tau_hat[4:6]))
      #   Y_n1r0_estimate <- hf(a + b * (tau_hat[1] + tau_hat[2] + logit.gx_n1r0 * tau_hat[3] + X_n1r0 %*% tau_hat[4:6]))
      # }
      if (current_s == 1) {
        Y_n1r1_estimate <- hf(a + b * (tau_hat[1] +  X_n1r1 %*% tau_hat[3:5]))
        Y_n1r0_estimate <- hf(a + b * (tau_hat[1] +  X_n1r0 %*% tau_hat[3:5]))
      } else if (current_s == 2) {
        Y_n1r1_estimate <- hf(a + b * (tau_hat[1] + tau_hat[2] + X_n1r1 %*% tau_hat[3:5]))
        Y_n1r0_estimate <- hf(a + b * (tau_hat[1] + tau_hat[2] + X_n1r0 %*% tau_hat[3:5]))
      }
      # if (current_s == 1) {
      #   Y_n1r1_estimate <- hf(tau_hat[1] +  X_n1r1 %*% tau_hat[3:5])
      #   Y_n1r0_estimate <- hf(tau_hat[1] +  X_n1r0 %*% tau_hat[3:5])
      # } else if (current_s == 2) {
      #   Y_n1r1_estimate <- hf(tau_hat[1] + tau_hat[2] + X_n1r1 %*% tau_hat[3:5])
      #   Y_n1r0_estimate <- hf(tau_hat[1] + tau_hat[2] + X_n1r0 %*% tau_hat[3:5])
      # }
      
      
      
      # Compute residuals for the four equations
      # eq1_1 <- 2 * sum((Y_n1r1 - Y_n1r1_estimate) * logit.gx_n1r1 * bh(t0 + t1 * logit.gx_n1r1)) +
      #   sum((Y_n1r1_estimate - Y_n1r1) * logit.gx_n1r1)
      # eq1_2 <- 2 * sum((Y_n1r1 - Y_n1r1_estimate) * bh(t0 + t1 * logit.gx_n1r1)) +
      #   sum((Y_n1r1_estimate - Y_n1r1) * logit.gx_n1r1)
      eq1_1 <- sum((Y_n1r1 - Y_n1r1_estimate) * logit.gx_n1r1 * bh(t0 + t1 * logit.gx_n1r1)) -
        sum((1 - Y_n1r1 - (1-Y_n1r1_estimate)) * logit.gx_n1r1*hf(t0 + t1 * logit.gx_n1r1))
      eq1_2 <-sum((Y_n1r1 - Y_n1r1_estimate) * bh(t0 + t1 * logit.gx_n1r1)) -
        sum((1 - Y_n1r1 - (1-Y_n1r1_estimate)) * hf(t0 + t1 * logit.gx_n1r1))
      
      eq2_1 <- sum((Y_n1r1 == 1) * logit.gx_n1r1 * bh(t0 + t1 * logit.gx_n1r1) - 
                     (Y_n1r1 == 0) * logit.gx_n1r1 * hf(t0 + t1 * logit.gx_n1r1)) +
        sum(Y_n1r0_estimate * logit.gx_n1r0 * bh(t0 + t1 * logit.gx_n1r0) - 
              (1 - Y_n1r0_estimate) * logit.gx_n1r0 * hf(t0 + t1 * logit.gx_n1r0)) -
        sum(logit.gx_N0 * hf(t0 + t1 * logit.gx_N0))
      
      eq2_2 <- sum((Y_n1r1 == 1) * bh(t0 + t1 * logit.gx_n1r1) - 
                     (Y_n1r1 == 0) * hf(t0 + t1 * logit.gx_n1r1)) +
        sum(Y_n1r0_estimate * bh(t0 + t1 * logit.gx_n1r0) - 
              (1 - Y_n1r0_estimate) * hf(t0 + t1 * logit.gx_n1r0)) -
        sum(hf(t0 + t1 * logit.gx_N0))
      
      # Return the sum of squared residuals
      #sum(eq1_1^2 + eq1_2^2 + eq2_1^2 + eq2_2^2)
      c(eq1_1, eq1_2, eq2_1, eq2_2)
      #c(eq2_1, eq2_2)
    }
    
    init1=glm(Y_n1r1~logit.gx_n1r1,family=binomial)
    initial_params <- c(init1$coefficients[1], init1$coefficients[2],0,1)  # Replace with reasonable starting values
    #result <- optim(initial_params, objective_function, method = "BFGS")
    #print(initial_params)
    #-----try-catch-----#
    # max_iter <- 100   # maximum number of attempts to avoid an infinite loop
    # iter <- 1
    # result <- NULL
    # 
    # repeat {
    #   if (iter > max_iter) {
    #     stop("Max iterations reached without a valid result.")
    #   }
    #   
    #   result <- tryCatch({
    #     #optim(initial_params, objective_function, method = "BFGS")
    #     nleqslv(initial_params, objective_function)
    #   }, error = function(e) {
    #     message("Error during optimization: ", e$message)
    #     NULL
    #   })
    #   
    #   # Check if result is valid and parameters are within bounds
    #   if (!is.null(result) && !any(result$x > 100)) {
    #     break  # Valid result found, exit the loop
    #   } else {
    #     iter <- iter + 1
    #   }
    # }
    
    #print(result$par)
    
    #rs.gX <- hf(result$par[1] + result$par[2] * logit(train$gX))
    
    result <- nleqslv(initial_params, objective_function)
    
    print(result$x)
    
    rs.gX <- hf(result$x[1] + result$x[2] * logit(train$gX))
    train <- train %>%
      bind_cols(rs.gX = rs.gX)
    
  } else {stop("unknown calibration method specified")}
  
  
  return(train)
}
