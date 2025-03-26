# Time Series Econometrics: Assignment 2

# Maximilian Droschl (21-616-412)
# Mohamed Ayadi (18-814-996)
# Aram Jivan Chetirian (19-607-258)

# Load packages
library(polynom)
library(stats)

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

#PART 2
library(tseries)
library(xts)
library(zoo)
library(forecast)
##############################################################################
# a) Time Series Plot
##############################################################################
# ------------------------------------------------
# 1) Read the data and create a time-series object
# ------------------------------------------------
file_path <- file.choose()
print(file_path)
raw_vals <- scan(file_path, skip = 7)
df <- read.table(
  file = file_path,
  skip = 7,       # skip lines 1–7
  header = FALSE, # don't use line 8 as header
  sep = "\t",     # data are tab-separated from line 8 onward
  strip.white = TRUE
)

# Assign column names manually
colnames(df) <- c("money", "gdp", "cpi")
head(df)

# Extract the GDP column
gdp <- df$GDP
gdp <- as.numeric(df$gdp)
quarters <- seq(as.yearqtr("1970 Q1"), by = 0.25, length.out = length(gdp))
date_index <- as.Date(quarters)
gdp_xts <- xts(gdp, order.by = date_index)

head(gdp_xts)

# ------------------------------------------------
# 2) Plot the real GDP time series
# ------------------------------------------------
plot(gdp_xts,
     main = "Real Euro GDP (1970Q1 - 2014Q4)",
     xlab = "Time (Quarterly)",
     ylab = "GDP (in dataset units)")
# ------------------------------------------------
# 3) Compute and plot the quarterly growth rate of real GDP
# ------------------------------------------------
gdp_growth <- diff(log(gdp_xts)) * 100
plot(gdp_growth,
     main = "Quarterly Growth Rate of Real Euro GDP (log diff * 100)",
     xlab = "Time (Quarterly)",
     ylab = "Percent")
# ------------------------------------------------
# 4) Perform Augmented Dickey-Fuller (ADF) tests
# ------------------------------------------------
adf_level <- adf.test(as.numeric(gdp_xts))
print(adf_level)
gdp_growth <- diff(log(gdp_xts)) * 100
gdp_growth <- na.omit(gdp_growth)

# Augmented Dickey-Fuller test 
adf_growth <- adf.test(as.numeric(gdp_growth))
print(adf_growth)

##############################################################################
# b) Box-Jenkins Approach for Real GDP (Growth) Model Selection
##############################################################################


gdp_growth <- na.omit(gdp_growth)
acf(gdp_growth, main = "ACF of GDP Growth")
pacf(gdp_growth, main = "PACF of GDP Growth")

# Fit AR(p) models for p = 0..9, record AIC and BIC
pmax <- 9
n <- length(gdp_growth)

results <- data.frame(
  p   = 0:pmax,
  AIC = NA_real_,
  BIC = NA_real_
)

for (p in 0:pmax) {
  fit <- Arima(gdp_growth, order = c(p, 0, 0), include.mean = TRUE)
  
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

##############################################################################
# c) Estimate and Forecast the Chosen AR Models, Check Diagnostics
##############################################################################

chosen_p <- best_aic_p
ar_fit <- Arima(gdp_growth, order = c(chosen_p, 0, 0), include.mean = TRUE)
summary(ar_fit)

# Forecast h=8 quarters ahead (2 years), for example
ar_forecast <- forecast(ar_fit, h = 8)
plot(ar_forecast, main = paste("AR(", chosen_p, ") Forecast of GDP Growth", sep=""))

# 1) Standardized residuals
std_resid <- residuals(ar_fit) / sd(residuals(ar_fit))
plot(std_resid, type = "l", main = "Standardized Residuals")
abline(h = 0, col = "red")

# 2) ACF/PACF of residuals
acf(std_resid, main = "ACF of Standardized Residuals")
pacf(std_resid, main = "PACF of Standardized Residuals")

# 3) Ljung-Box test for autocorrelatio
Box.test(std_resid, lag = 12, type = "Ljung-Box")

# 4) Shapiro-Wilk test for normality
shapiro.test(std_resid)

##############################################################################
# d) Use the Growth Rate of Real GDP from 1972Q2–2007Q4 and Estimate ARIMA
##############################################################################

gdp_growth_vec <- as.numeric(gdp_growth)
gdp_growth_ts <- ts(
  gdp_growth_vec,
  start = c(1970, 1),   
  frequency = 4         
)

gdp_growth_sub <- window(gdp_growth_ts, start = c(1972, 2), end = c(2007, 4))
head(gdp_growth_sub)

model_sub <- auto.arima(gdp_growth_sub)
summary(model_sub)

fc_sub <- forecast(model_sub, h = 4)
plot(fc_sub, main = "ARIMA Forecast (1972Q2–2007Q4 Subsample)")


