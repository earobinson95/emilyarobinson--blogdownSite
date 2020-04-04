
# 3 Step GLS Function

WLS <- function(y = y, x = x, func = expression(), theta0 = c(), psi0 = 0.5, tol = 10^{-8}, maxiter = 20){
  require(Deriv)
  require(MASS)
  
  iteration <- 0
  p <- length(theta0)
  W <- diag(rep(1,length(y)))
  theta_hat_old <- rep(10000,length(theta0))
  
  psi_old       <- psi0
  dummy_data <- rep(0,length(y))
  Fpsi_func <- expression(e*(geo_mean/yhat)^psi_old)
  dpsi      <- Deriv(Fpsi_func, c("psi_old"))
  
  # initial run
  if (p >= 1){
    a  <- theta0[1];
    da <- Deriv(func, c("a"))
  }
  if (p >= 2){
    b  <- theta0[2];
    db <- Deriv(func, c("b"))
  }
  if (p >= 3){
    c  <- theta0[3];
    dc <- Deriv(func, c("c"))
  }
  if (p >= 4){
    d  <- theta0[4]; 
    dd <- Deriv(func, c("d"))
  }
  if (p >= 5){
    e  <- theta0[5]; 
    de <- Deriv(func, c("e"))
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
  
  yhat <- eval(func)
  e    <- y - yhat
  
  delta     <- solve(t(V)%*%W%*%V)%*%t(V)%*%W%*%e
  theta_new <- theta0 + delta
  
  repeat{
    
    repeat{
      # converge to theta hats
      theta_old <- theta_new
      
      if (p >= 1){
        a  <- theta_old[1];
        da <- Deriv(func, c("a"))
      }
      if (p >= 2){
        b  <- theta_old[2];
        db <- Deriv(func, c("b"))
      }
      if (p >= 3){
        c  <- theta_old[3];
        dc <- Deriv(func, c("c"))
      }
      if (p >= 4){
        d  <- theta_old[4]; 
        dd <- Deriv(func, c("d"))
      }
      if (p >= 5){
        e  <- theta_old[5]; 
        de <- Deriv(func, c("e"))
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
      
      yhat <- eval(func)
      e    <- y - yhat
      
      delta     <-  solve(t(V)%*%W%*%V)%*%t(V)%*%W%*%e
      theta_new <- theta_old + delta
      
      
      # check theta convergence
      max_relative_change_ols <- max(abs(theta_new - theta_old)/abs(theta_old))
      
      if (max_relative_change_ols < tol) break
    }
    
    # check theta hat convergence
    theta_hat_new <- theta_new
    iteration <- iteration + 1
    max_relative_change_gls <- max(abs(theta_hat_new - theta_hat_old)/abs(theta_hat_old))
    theta_hat_old <- theta_new
    
    yhat <- eval(func)
    e    <- y - yhat
    
    # estimate psi -------------------------------------------------------
    geo_mean   <- prod(yhat)^{1/length(yhat)}
    
    repeat{
      
      Vpsi      <- matrix(eval(dpsi), ncol = 1, nrow = length(y))
      
      Fi       <- eval(Fpsi_func)
      e_psi    <- dummy_data - Fi
      
      delta_psi <-  solve(t(Vpsi)%*%Vpsi)%*%t(Vpsi)%*%e_psi
      psi_new   <- psi_old + delta_psi/2
      psi_new
      psi_old <- as.numeric(psi_new)
      
      
      # check psi convergence
      max_relative_change_psi <- max(abs(psi_new - psi_old)/abs(psi_old))
      max_relative_change_psi
      psi_old <- as.numeric(psi_new)
      
      if (max_relative_change_psi < 10^{-6}) break
    }
    
    # create weights -----------------------------------------------------
    W <- diag(1/eval(func)^{2*psi_new})
    delta     <- solve(t(V)%*%W%*%V)%*%t(V)%*%W%*%e
    theta_new <- theta_hat_old + delta
    
    #if (max_relative_change_gls < tol | iteration == maxiter) break
    if (iteration == maxiter) break
  }
  
  yhat              <- eval(func)
  OLS_resid         <- y - yhat
  sigma_sq          <- (t(e)%*%W%*%e)/(length(y)-length(theta0))
  
  ystar             <- diag(sqrt(W)*y)
  Vstar             <- sqrt(W)%*%V
  wr                <- ystar - diag(sqrt(W)*yhat)
  hat               <- Vstar%*%solve(t(Vstar)%*%Vstar)%*%t(Vstar)
  studentized_resid <- wr/(sqrt(sigma_sq)*sqrt(rep(1, length(y)) - diag(hat)))
  
  return(list(C = iteration, theta = theta_new, psi = psi_new, V = V, Vstar = Vstar, W = W, OLS_resid = OLS_resid, wr = wr, studentized_resid = studentized_resid, sigma_sq = sigma_sq, ystar = ystar, yhat = yhat))
  
}

# x     <- c(2,14,18)
# y     <- c(13,273,551)
# trial <- GLS_unkownpsi(y = y, x = x, func = expression(exp(a+b*x)),  psi0 = 0.5, theta0 = c(2,0.25), tol = 10^{-8}, maxiter = 20)
# trial

