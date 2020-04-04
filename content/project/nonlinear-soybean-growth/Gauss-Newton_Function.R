gauss_newton <- function(y = y, x = x, func = expression(), theta = c(), alpha = 0.05, lambda = 0.5, epsilon = 0.00001, maxiter = 20){
  require(Deriv)
  require(MASS)
  
  p <- length(theta)
  #counts number of iterations
  i = 0
  #initial SS
  S_old = 1000000000
  
  repeat{
    i = i + 1
    
    if (p >= 1){
      a  <- theta[1];
      da <- Deriv(func, c("a"))
    }
    if (p >= 2){
      b  <- theta[2];
      db <- Deriv(func, c("b"))
    }
    if (p >= 3){
      c  <- theta[3];
      dc <- Deriv(func, c("c"))
    }
    if (p >= 4){
      d  <- theta[4]; 
      dd <- Deriv(func, c("d"))
    }
    if (p >= 5){
      e  <- theta[5]; 
      de <- Deriv(func, c("e"))
    }
    if (p >= 6){
      f  <- theta[6]; 
      df <- Deriv(func, c("f"))
    }
    
    if (p == 1) {
      if(eval(da)[1] == 1){da <- rep(1,length(y))}
      V <- matrix(eval(da), ncol = 1, nrow = length(y))
    }
    if (p == 2) {
      if(eval(da)[1] == 1){da <- rep(1,length(y))}
      if(eval(db)[1] == 1){db <- rep(1,length(y))}
      V <- matrix(c(eval(da), eval(db)), ncol = 2, nrow = length(y))
    }
    if (p == 3) {
      if(eval(da)[1] == 1){da <- rep(1,length(y))}
      if(eval(db)[1] == 1){db <- rep(1,length(y))}
      if(eval(dc)[1] == 1){dc <- rep(1,length(y))}
      V <- matrix(c(eval(da), eval(db), eval(dc)), ncol = 3, nrow = length(y))
    }
    if (p == 4) {
      if(eval(da)[1] == 1){da <- rep(1,length(y))}
      if(eval(db)[1] == 1){db <- rep(1,length(y))}
      if(eval(dc)[1] == 1){dc <- rep(1,length(y))}
      if(eval(dd)[1] == 1){dd <- rep(1,length(y))}
      V <- matrix(c(eval(da), eval(db), eval(dc), eval(dd)), ncol = 4, nrow = length(y))
    }
    if (p == 5) {
      if(eval(da)[1] == 1){da <- rep(1,length(y))}
      if(eval(db)[1] == 1){db <- rep(1,length(y))}
      if(eval(dc)[1] == 1){dc <- rep(1,length(y))}
      if(eval(dd)[1] == 1){dd <- rep(1,length(y))}
      if(eval(de)[1] == 1){de <- rep(1,length(y))}
      V <- matrix(c(eval(da), eval(db), eval(dc), eval(dd), eval(de)), ncol = 5, nrow = length(y))
    }
    if (p == 6) {
      if(eval(da)[1] == 1){da <- rep(1,length(y))}
      if(eval(db)[1] == 1){db <- rep(1,length(y))}
      if(eval(dc)[1] == 1){dc <- rep(1,length(y))}
      if(eval(dd)[1] == 1){dd <- rep(1,length(y))}
      if(eval(de)[1] == 1){de <- rep(1,length(y))}
      if(eval(df)[1] == 1){df <- rep(1,length(y))}
      V <- matrix(c(eval(da), eval(db), eval(dc), eval(dd), eval(de), eval(df)), ncol = 6, nrow = length(y))
    }
    
    yhat = eval(func)
    S <- t(y-yhat)%*%(y-yhat)
    delta <- ginv(t(V)%*%V)%*%t(V)%*%(y-yhat)
    offset <- sqrt((t(y-yhat)%*%V%*%ginv(t(V)%*%V)%*%t(V)%*%(y-yhat))/(t(y-yhat)%*%(y-yhat)))
    
    #store old values
    S_old <- S
    theta_old <- theta
    #update values
    theta <- theta_old + delta
    #check for step necesity
    if (p >= 1){a  <- theta[1]}
    if (p >= 2){b  <- theta[2]}
    if (p >= 3){c  <- theta[3]}
    if (p >= 4){d  <- theta[4]}
    if (p >= 5){e  <- theta[5]}
    yhat = eval(func)
    S <- t(y-yhat)%*%(y-yhat)
    if (S > S_old){
      theta = theta_old + lambda*delta
    }
    
    if (offset < epsilon | i > maxiter) break
  }
  
  e           <- y - yhat
  dgf          <- (length(y)-length(theta))
  sigma_sq    <- as.numeric(S)/dgf
  se_theta    <- sqrt(diag(as.numeric(S)/dgf*ginv(t(V)%*%V)))
  lower_theta <- theta - qt(1-alpha/2, dgf)*se_theta 
  upper_theta <- theta + qt(1-alpha/2, dgf)*se_theta 
  theta_est   <- cbind(theta, se_theta, lower_theta, upper_theta)
  colnames(theta_est) <- c("Estimate", "Std Error", "Lower", "Upper")
  if (p == 1){rownames(theta_est) <- c("Theta 1")}
  if (p == 2){rownames(theta_est) <- c("Theta 1", "Theta 2")}
  if (p == 3){rownames(theta_est) <- c("Theta 1", "Theta 2", "Theta 3")}
  if (p == 4){rownames(theta_est) <- c("Theta 1", "Theta 2", "Theta 3","Theta 4")}
  if (p == 5){rownames(theta_est) <- c("Theta 1", "Theta 2", "Theta 3","Theta 4","Theta 5")}
  
  
  return(list(theta = round(theta,6), sigma_sq = round(sigma_sq,3), yhat = yhat, e = e, SSE = S, V = V, df = dgf, se_theta = se_theta, lower_theta = lower_theta, upper_theta = upper_theta, theta_est = theta_est))
}