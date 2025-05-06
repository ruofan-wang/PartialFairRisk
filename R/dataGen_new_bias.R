dataGen_new_bias <- function(ngps = 2 # number of groups
                        , n = c(1000,1000) # sample sizes for group 1, 2
                        , N0 = c(200,200)
                        , n1r1 = c(200,200)
                        , alpha = c(2,4) # beta distribution shape parameter
                        , beta = c(8,8) # beta distribution scale parameter
                        , sigeps = c(0.01,0.01) # noise
                        , theta_a = c(0,0) # calibration intercept
                        , theta_b = c(1,1) # calibration slope
                        , true_control = c(0.9,0.8)
                        , true_case = c(1,0.9)
                        , h = c("logit","logit")){
  
  # ------------ Data generation ------------- #
  # group label: s
  s <- rep(c(1:ngps),times=n)
  
  # true risk: pi |s ~ Beta(alpha_s,beta_s)
  alpha_s <- rep(alpha,times=n)
  beta_s <- rep(beta,times=n)
  pi_s <- rbeta(sum(n),alpha_s,beta_s)
  
  # y: y | pi ~ Bern(pi)
  y <- rbinom(sum(n),1,pi_s)
  
  # gX = h(pi_s,theta)
  a_s <- rep(theta_a,times=n)
  b_s <- rep(theta_b,times=n)
  
  sigeps_s <- rep(sigeps,times=n)
  epsilon <- rnorm(sum(n),0,sigeps_s)
  
  if(h[1] == "logit"){
    h.pi_s <- log(pi_s/(1-pi_s))
  } else if(h[1] == "cloglog"){
    h.pi_s <- log(-log(1-pi_s))
  } else{
    print("choice of h is unknown")
    stop()
  }
  
  logit.gX <- a_s + b_s * h.pi_s + epsilon
  gX <- exp(logit.gX)/(1+exp(logit.gX))
  
  #Z1=h.pi_s+rnorm(length(gX),0,0.1)
  Z1=h.pi_s+rnorm(length(gX),0,0.01)
  Z2=Z1+rnorm(length(gX),0,1)
  Z3=Z2+rnorm(length(gX),0,5)
  
  #Z1=h.pi_s+rnorm(length(gX),0,0.01)
  #Z2=5+0.5*logit.gX+rnorm(length(gX),0,1)
  #Z3=gX+rnorm(length(gX),0,0.09)
  
  print(cor(Z1,logit.gX))
  print(cor(Z1,h.pi_s))
  
  validate <- integer(sum(n))  # all 0
  control  <- integer(sum(n))  # all 0
  
  for(i in 1:ngps) {
    # Indices for group i
    idx_i <- which(s == i)
    # Subset of group i where y=0
    idx_i_y0 <- idx_i[y[idx_i] == 0]
    
    # Randomly choose N0[i] among y=0 for control=1
    idx_control <- sample(idx_i_y0, size = N0[i], replace = FALSE)
    control[idx_control] <- 1
    
    # From the remaining in group i, choose n1r1[i] for validate=1
    idx_remaining <- setdiff(idx_i, idx_control)
    idx_validate <- sample(idx_remaining, size = n1r1[i], replace = FALSE)
    validate[idx_validate] <- 1
  }
  
  y_new <- y  # initialize with the original values
  
  for(i in 1:ngps) {
    idx_group <- which(s == i)
    idx_y1 <- idx_group[y[idx_group] == 1]
    idx_y0 <- idx_group[y[idx_group] == 0]
    # For these, simulate misclassification: they are identified as cases
    # with probability true_case[i]
    y_new[idx_y1] <- rbinom(length(idx_y1), 1, true_case[i])
    y_new[idx_y0] <- rbinom(length(idx_y0), 1, 1-true_control[i])
  }
  #y_new: observed y value with misspecification
  #y: true y value
  #y_final: y value with misspecification if validate=1 then the true value
  #y_final <- ifelse(validate == 1, y, y_new)
  # Create dataframe 
  #simData <- data.frame(s, pi_s, y,logit.gX, gX, control, validate)
  #simData <- data.frame(s, pi_s, y = y_final, y_true = y, y_new=y_new, logit.gX, gX, control, validate)
  simData <- data.frame(s, pi_s, y = y, y_new=y_new, logit.gX, gX, control, validate)
  Z <- data.frame(Z1,Z2,Z3)
  return(list(simData = simData, Z = Z))
}
