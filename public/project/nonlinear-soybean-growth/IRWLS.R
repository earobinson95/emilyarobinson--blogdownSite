
# IRWLS Function

IRWLS <- function(y = y, x = x, func = expression(), theta0 = c(), psi = 0.5, tol = 10^{-8}, maxiter = 20){
  require(Deriv)
  
  iteration <- 0
  p <- length(theta0)
  theta_hat_old <- rep(10000,length(theta0))
  
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
  W <- diag(rep(1,length(y)))
  delta     <- solve(t(V)%*%W%*%V)%*%t(V)%*%W%*%e
  theta_new <- theta0 + delta
  
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
      
      W <- diag(1/eval(func)^{2*psi})
      delta     <-  solve(t(V)%*%W%*%V)%*%t(V)%*%W%*%e
      theta_new <- theta_old + delta
      
      iteration <- iteration + 1
      # check theta convergence
      max_relative_change <- max(abs(theta_new - theta_old)/abs(theta_old))
      
      if (max_relative_change < tol | iteration == maxiter) break
  }
  if(iteration == maxiter){print("Did not converge")}
  
  yhat <- eval(func)
  e    <- y - yhat
  sigma_sq <- (t(e)%*%W%*%e)/(length(y)-length(theta0))
  
  return(list(C = iteration, theta = theta_new, sigma_sq = sigma_sq, yhat = yhat))
}

# x     <- c(2,14,18)
# y     <- c(13,273,551)
# trial <- IRWLS(y = y, x = x, func = expression(exp(a+b*x)),  psi = 0.5, theta0 = c(2,0.25), tol = 10^{-8}, maxiter = 20)
# trial
