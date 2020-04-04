
# READ IN DATA ---------------------------------------------------------------------------
soybean_data <- read.csv("C:/Users/EmilyARobinson/Dropbox/Nonlinear/Soybean Growth/Data/soybean_data.csv")
cols <- c(1,3)
soybean_data[cols]  <- lapply(soybean_data[cols], factor)

library(fastDummies)
soybean_data2 <- soybean_data[order(soybean_data$Genotype),]
soybean_data2 <- dummy_cols(soybean_data2, select_columns = "Genotype", remove_first_dummy = TRUE)

# 3 Parameter Logistic OLS ------------------------------------------------------------------
ols3p_nls <- nls(Leaf_Weight ~ a/(1+b*exp(-c*Days)), 
                 data = soybean_data, 
                 start = list(a = 20, b = 700,c = 0.125))
round(summary(ols3p_nls)$coeff[,1:2],3)


# 4 Parameter Logistic OLS ------------------------------------------------------------------
ols4p_nls <- nls(Leaf_Weight ~ a + (b-a)/(1+exp((c-Days)/d)), 
                 data = soybean_data, 
                 start = list(a = 0.2, b = 20,c = 50, d = 8))
round(summary(ols4p_nls)$coeff[,1:2],3)

# 3 Parameter vs 4 Parameter ----------------------------------------------------------------
anova(ols3p_nls, ols4p_nls)

# Indicator 3 Parameter Logistic OLS --------------------------------------------------------

# Full
ols3p_ind0 <- nls(Leaf_Weight ~ (a+ap*Genotype_P)/(1+(b+bp*Genotype_P)*exp(-(c+cp*Genotype_P)*Days)), 
                 data = soybean_data2, 
                 start = list(a = 16, ap = 4.78, b = 1035, bp = -490, c = 0.125, cp = -0.01))
round(summary(ols3p_ind0)$coeff[,1:2],3)

# Take away ap? No
ols3p_ind1 <- nls(Leaf_Weight ~ (a)/(1+(b+bp*Genotype_P)*exp(-(c+cp*Genotype_P)*Days)), 
                  data = soybean_data2, 
                  start = list(a = 20, b = 1035, bp = -490, c = 0.125, cp = -0.01))
anova(ols3p_ind1, ols3p_ind0)

# Take away bp? Yes
ols3p_ind2 <- nls(Leaf_Weight ~ (a + ap*Genotype_P)/(1+(b)*exp(-(c+cp*Genotype_P)*Days)), 
                  data = soybean_data2, 
                  start = list(a = 16, ap = 4.78, b = 700, c = 0.125, cp = -0.01))
anova(ols3p_ind2, ols3p_ind0)


# Take away cp? Yes
ols3p_ind3 <- nls(Leaf_Weight ~ (a + ap*Genotype_P)/(1+(b)*exp(-(c)*Days)), 
                  data = soybean_data2, 
                  start = list(a = 16, ap = 4.78, b = 700, c = 0.125))
anova(ols3p_ind3, ols3p_ind2)

library(ggplot2)
ggplot(soybean_data2, aes(x = Days, y = Leaf_Weight, group = Genotype, color = Genotype)) +
  geom_point() +
  geom_line(aes(y = fitted(ols3p_ind2), group = Genotype, color = Genotype)) +
  theme_minimal() +
  xlab("Days after Planting") +
  ylab("Average Leaf Weight/Plant (g)")

source("ResidualPanel.R")
resid_panel(data = soybean_data2, yhat = fitted(ols3p_ind2), observed = soybean_data$Leaf_Weight)

soybean_data2$OLS_yhat <- fitted(ols3p_ind2)
soybean_data2 <- soybean_data2[order(soybean_data2$Days),]

#Does not appear to be correlation
plot(soybean_data2$Days,soybean_data2$Leaf_Weight - soybean_data2$OLS_yhat, xlab = "Days", ylab = "Residuals")


# Nonconstant variance issue
source("3step_GLS_unknownpsi.R")

GLS_unknownpsi <- GLS_unkownpsi(y = soybean_data2$Leaf_Weight, x = soybean_data2$Days, func = expression((a + b*soybean_data2$Genotype_P)/(1+c*exp(-d*x))), theta0 = c(16, 4.78, 700, 0.125, -0.01), psi0 = 0.5, tol = 10^{-8}, maxiter = 20)







wfct <- function(expr)
{
  expr <- deparse(substitute(expr))
  
  ## create new environment
  newEnv <- new.env()
  
  ## get call
  mc <- sys.calls()[[1]]
  mcL <- as.list(mc)
  
  ## get data and write to newEnv
  DATA <- mcL[["data"]]
  DATA <- eval(DATA)
  DATA <- as.list(DATA)
  NAMES <- names(DATA)
  for (i in 1:length(DATA)) assign(NAMES[i], DATA[[i]], envir = newEnv)
  
  ## get parameter, response and predictor names
  formula <- as.formula(mcL[[2]])
  VARS <- all.vars(formula)
  RESP <- VARS[1]
  RHS <- VARS[-1]
  PRED <- match(RHS, names(DATA))
  PRED <- names(DATA)[na.omit(PRED)]
  
  ## calculate variances for response values if "error" is in expression
  ## and write to newEnv
  if (length(grep("error", expr)) > 0) {
    y <- DATA[[RESP]]
    x <- DATA[[PRED]]
    ## test for replication
    if (!any(duplicated(x))) stop("No replicates available to calculate error from!")
    ## calculate error
    error <- tapply(y, x, function(e) var(e, na.rm = TRUE))
    error <- as.numeric(sqrt(error))
    ## convert to original repititions
    error <- rep(error, as.numeric(table(x)))
    assign("error", error, envir = newEnv)
  }
  
  ## calculate fitted or residual values if "fitted"/"resid" is in expression
  ## and write to newEnv
  if (length(grep("fitted", expr)) > 0 || length(grep("resid", expr)) > 0) {
    mc2 <- mc
    mc2$weights <- NULL
    MODEL <- eval(mc2)
    fitted <- fitted(MODEL)
    resid <- residuals(MODEL)
    assign("fitted", fitted, newEnv)
    assign("resid", resid, newEnv)
  }
  
  ## return evaluation in newEnv: vector of weights
  OUT <- eval(parse(text = expr), envir = newEnv)
  return(OUT)
}

ols3p_GLS <- nls(Leaf_Weight ~ (a + ap*Genotype_P)/(1+(b)*exp(-(c)*Days)), 
                  data = soybean_data2, 
                  start = list(a = 16, ap = 4.78, b = 700, c = 0.125),
                  weights = wfct(1/fitted^2*1))
summary(ols3p_GLS)$coeff

ggplot(soybean_data2, aes(x = Days, y = Leaf_Weight, group = Genotype, color = Genotype)) +
  geom_point() +
  geom_line(aes(y = fitted(ols3p_GLS), group = Genotype, color = Genotype)) +
  theme_minimal() +
  xlab("Days after Planting") +
  ylab("Average Leaf Weight/Plant (g)")

resid_panel(data = soybean_data2, yhat = fitted(ols3p_GLS), observed = soybean_data2$Leaf_Weight)


```{r mod7, echo=F, message=FALSE, warning=FALSE}
wfct <- function(expr)
{
  expr <- deparse(substitute(expr))
  
  ## create new environment
  newEnv <- new.env()
  
  ## get call
  mc <- sys.calls()[[1]]
  mcL <- as.list(mc)
  
  ## get data and write to newEnv
  DATA <- mcL[["data"]]
  DATA <- eval(DATA)
  DATA <- as.list(DATA)
  NAMES <- names(DATA)
  for (i in 1:length(DATA)) assign(NAMES[i], DATA[[i]], envir = newEnv)
  
  ## get parameter, response and predictor names
  formula <- as.formula(mcL[[2]])
  VARS <- all.vars(formula)
  RESP <- VARS[1]
  RHS <- VARS[-1]
  PRED <- match(RHS, names(DATA))
  PRED <- names(DATA)[na.omit(PRED)]
  
  ## calculate variances for response values if "error" is in expression
  ## and write to newEnv
  if (length(grep("error", expr)) > 0) {
    y <- DATA[[RESP]]
    x <- DATA[[PRED]]
    ## test for replication
    if (!any(duplicated(x))) stop("No replicates available to calculate error from!")
    ## calculate error
    error <- tapply(y, x, function(e) var(e, na.rm = TRUE))
    error <- as.numeric(sqrt(error))
    ## convert to original repititions
    error <- rep(error, as.numeric(table(x)))
    assign("error", error, envir = newEnv)
  }
  
  ## calculate fitted or residual values if "fitted"/"resid" is in expression
  ## and write to newEnv
  if (length(grep("fitted", expr)) > 0 || length(grep("resid", expr)) > 0) {
    mc2 <- mc
    mc2$weights <- NULL
    MODEL <- eval(mc2)
    fitted <- fitted(MODEL)
    resid <- residuals(MODEL)
    assign("fitted", fitted, newEnv)
    assign("resid", resid, newEnv)
  }
  
  ## return evaluation in newEnv: vector of weights
  OUT <- eval(parse(text = expr), envir = newEnv)
  return(OUT)
}
mod7 <- nls(Leaf_Weight ~ (a + ap*Genotype_P)/(1+(b)*exp(-(c)*Days)), 
            data = soybean_data2, 
            start = list(a = 16, ap = 4.78, b = 700, c = 0.125),
            weights = wfct(1/fitted^2*0.5))
kable(round(summary(mod7)$coeff[,1:2],3), caption = "Estimated parameters for Model (7)")
```



source("IRWLS.R")
irwls.mod <- IRWLS(y = soybean_data2$Leaf_Weight, x = soybean_data2$Days, 
                   func = expression((a + b*soybean_data2$Genotype_P)/(1+c*exp(-d*x))), 
                   theta0 = c(16, 4.78, 700, 0.125), 
                   psi = 1, 
                   tol = 10^{-8}, 
                   maxiter = 20)
irwls.mod$theta
ggplot(soybean_data2, aes(x = Days, y = Leaf_Weight, group = Genotype, color = Genotype)) +
  geom_point() +
  geom_line(aes(y = irwls.mod$yhat, group = Genotype, color = Genotype)) +
  theme_minimal() +
  xlab("Days after Planting") +
  ylab("Average Leaf Weight/Plant (g)")
resid_panel(data = soybean_data2, yhat = irwls.mod$yhat, observed = soybean_data2$Leaf_Weight)

