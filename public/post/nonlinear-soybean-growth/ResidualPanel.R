resid_panel <- function(data = data, yhat = yhat, observed = y){
  require(ggplot2)
  require(gridExtra)
  e = observed - yhat
  resid1 <- ggplot(data, aes(x = yhat, y = e)) +
    geom_point() +
    geom_line(aes(y = 0)) +
    theme_minimal() +
    ggtitle("Residual vs Predicted")
  
  resid2 <- ggplot(data, aes(x = observed, y = e)) +
    geom_point() +
    geom_line(aes(y = 0)) +
    theme_minimal() +
    ggtitle("Residual vs Observed")
  
  resid3 <- ggplot(data, aes(x = yhat, y = e^2)) +
    geom_point() +
    theme_minimal() +
    ggtitle("Sq. Residual vs Observed")
  
  resid4 <- ggplot(data, aes(x = yhat, y = abs(e))) +
    geom_point() +
    theme_minimal() +
    ggtitle("Abs. Residual vs Observed")
  
  grid.arrange(resid1, resid2, resid3, resid4, ncol = 2)
}