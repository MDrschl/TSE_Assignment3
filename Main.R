# Time Series Econometrics: Assignment 2

# Maximilian Droschl (21-616-412)
# Mohamed Ayadi (18-814-996)
# Aram Jivan Chetirian (19-607-258)

rm(list = ls(all = TRUE))

# Load packages
library(polynom)
library(stats)
library(tseries)
library(xts)
library(zoo)
library(forecast)
library(urca)

### Task 1

### 1.a) Define the polynomials
Theta <- c(1, -4.75, 7.375, -5, 1.5625)
Phi <- c(1, -1, 0.5)

### 1.b) Calculate roots using polyroot
ar_roots <- polyroot(Theta)
print(ar_roots)

ma_roots <- polyroot(Phi)
print(ma_roots)

### 1.c) Find reduced MA(2) model through polynomial division
# Perform polynomial division (Theta / Phi)
Q <- polynom::polynomial(Theta)/polynom::polynomial(Phi)
print(Q)

# Extracting coefficients of reduced MA(2)
Phi_red <- as.numeric(coef(Q))

### 1.d)
# Roots of reduced MA(2)
ma_roots_red <- polyroot(Phi_red)
print(ma_roots_red)

### 1.f)
# Roots of new MA(2), which now is invertible
Phi_new <- c(1, -1.2, 0.32)
ma_roots_new <- polyroot(Phi_new)
print(ma_roots_new)

rm(list = ls(all = TRUE))

### Task 2
##############################################################################
# a) Time Series Plot
##############################################################################
# ------------------------------------------------
# 1) Read the data and create a time-series object
# ------------------------------------------------
file_path <- file.choose()
print(file_path)
raw_vals <- scan(file_path, skip = 7)
data <- read.table(
  file = file_path,
  skip = 7,
  header = FALSE,
  sep = "\t",
  strip.white = TRUE
)

# Assign column names manually
colnames(data) <- c("money", "gdp", "cpi")
head(data)

# Extract the GDP column
gdp <- as.numeric(data$gdp)
quarters <- seq(as.yearqtr("1970 Q1"), by = 0.25, length.out = length(gdp))
date_index <- as.Date(quarters)
gdp <- xts(gdp, order.by = date_index)

head(gdp)

# ------------------------------------------------
# 2) Plot the real GDP time series
# ------------------------------------------------
plot(gdp,
     main = "Real Euro GDP (1970Q1 - 2014Q4)",
     xlab = "Time (Quarterly)")
# ------------------------------------------------
# 3) Compute and plot the quarterly growth rate of real GDP
# ------------------------------------------------
growth <- diff(gdp)/lag(gdp, 1) * 100
plot(growth,
     main = "Quarterly Growth Rate of Real Euro GDP (log diff * 100)",
     xlab = "Time (Quarterly)",
     ylab = "Percent")
# ------------------------------------------------
# 4) Perform Augmented Dickey-Fuller (ADF) tests
# ------------------------------------------------
# In this section, we follow the procedure suggested during the lecture,
# using the urca package.

# ADF test for GDP
# Checking for significant lags
adf10 = ur.df(gdp,type="trend",lags=10)
summary(adf10) # --> significant lag at 1

# Rerunning regression with one lag:
adf1 = ur.df(gdp,type="trend",lags=1)
summary(adf1)


# ADF test for GDP growth rate
growth <- na.omit(growth)
# Checking for significant lags
adf10 = ur.df(growth,type="trend",lags=10)
summary(adf10) # --> significant lag at 1

# Rerunning regression with one lag:
adf1 = ur.df(growth,type="trend",lags=1)
summary(adf1)

##############################################################################
# b) Box-Jenkins Approach for Real GDP (Growth) Model Selection
##############################################################################


acf(growth, main = "ACF of GDP Growth")
pacf(growth, main = "PACF of GDP Growth")

# Fit AR(p) models for p = 0..9, record AIC and BIC
pmax <- 9
n <- length(growth)

results <- data.frame(
  p   = 0:pmax,
  AIC = NA_real_,
  BIC = NA_real_
)

for (p in 0:pmax) {
  fit <- Arima(growth, order = c(p, 0, 0), include.mean = TRUE)
  
  results$AIC[p + 1] <- fit$aic
  
  k <- p + 1
  results$BIC[p + 1] <- fit$aic + (2 * k - 2) * log(n)
}

print(results)

# 4. Identify which p minimizes AIC and BIC
best_aic_p <- results$p[which.min(results$AIC)]
best_bic_p <- results$p[which.min(results$BIC)]
cat("Best p by AIC:", best_aic_p, "\n")
cat("Best p by BIC:", best_bic_p, "\n")

# Check back with internal function from the forecast package
optorder = auto.arima(growth,d=0)
optorder # --> gives similar results

##############################################################################
# c) Estimate and Forecast the Chosen AR Models, Check Diagnostics
##############################################################################

# best model according to AIC 

# Fit the chosen AR(p) model
ar_fit <- Arima(growth, order = c(best_aic_p, 0, 0), include.mean = TRUE)

# Display a summary
summary(ar_fit)

# Forecast h=8 quarters ahead (2 years), for example
ar_forecast <- forecast(ar_fit, h = 8)
plot(ar_forecast, main = paste("AR(", best_aic_p, ") Forecast of GDP Growth", sep=""))

# -------------------------
# Residual Diagnostics
# -------------------------
# 1) Standardized residuals
std_resid <- residuals(ar_fit) / sd(residuals(ar_fit))
plot(std_resid, type = "l", main = "Standardized Residuals: AR(2)")
abline(h = 0, col = "red")

# 2) ACF/PACF of residuals
acf(std_resid, main = "ACF of Standardized Residuals: AR(2)")
pacf(std_resid, main = "PACF of Standardized Residuals: AR(2)")

# 3) Ljung-Box test for autocorrelation
Box.test(std_resid, lag = 12, type = "Ljung-Box")

# 4) Shapiro-Wilk test for normality (small samples) 
shapiro.test(std_resid)

############################################
# best model according to BIC 

# Fit the chosen AR(p) model
ar_fit_bic <- Arima(growth, order = c(best_bic_p, 0, 0), include.mean = TRUE)

# Display a summary
summary(ar_fit_bic)

# Forecast h=8 quarters ahead (2 years), for example
ar_forecast_bic <- forecast(ar_fit_bic, h = 8)
plot(ar_forecast_bic, main = paste("AR(", best_bic_p, ") Forecast of GDP Growth", sep=""))

# -------------------------
# Residual Diagnostics
# -------------------------
# 1) Standardized residuals
std_resid_bic <- residuals(ar_fit_bic) / sd(residuals(ar_fit_bic))
plot(std_resid_bic, type = "l", main = "Standardized Residuals: AR(1)")
abline(h = 0, col = "red")

# 2) ACF/PACF of residuals
acf(std_resid_bic, main = "ACF of Standardized Residuals: AR(1)")
pacf(std_resid_bic, main = "PACF of Standardized Residuals: AR(1)")

# 3) Ljung-Box test for autocorrelation
Box.test(std_resid_bic, lag = 12, type = "Ljung-Box")

# 4) Shapiro-Wilk test for normality (small samples) 
shapiro.test(std_resid_bic)

##############################################################################
# d) Use the Growth Rate of Real GDP from 1972Q2–2007Q4 and Estimate ARIMA
##############################################################################

growth_vec <- as.numeric(growth)

growth_ts <- ts(
  growth_vec,
  start = c(1970, 1), # year, quarter
  frequency = 4         # quarterly
)

growth_sub <- window(growth_ts, start = c(1972, 2), end = c(2007, 4))
head(growth_sub)

ar2_model <- Arima(growth_sub, order = c(2, 0, 0))

summary(ar2_model)

# Forecast h=8 quarters ahead (2 years)
ar_forecast_2 <- forecast(ar2_model, h = 8)

# Create forecast plot
plot(ar_forecast_2, main = "AR(2) 2-Year GDP Growth Forecast", 
     xlab = "Year", ylab = "GDP Growth Rate", ylim = c(-3, 2))
actual_forecast_period <- window(growth_ts, 
                                 start = c(2008, 1), 
                                 end = c(2009, 4))
lines(actual_forecast_period, col = "red", lwd = 2)
legend("topleft", 
       legend = c("Forecast", "Actual"), 
       col = c("blue", "red"), 
       lty = 1,
       lwd = 2,
       bty = "n")
