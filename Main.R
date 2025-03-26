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
#a)
# ------------------------------------------------
# 1) Read the data and create a time-series object
# ------------------------------------------------

# Read from a text file. Adjust path as necessary.
# Assuming the file has a header with one column named 'GDP'.
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

# Assign meaningful column names manually
colnames(df) <- c("money", "gdp", "cpi")

# Check result
head(df)

# Extract the GDP column (adjust the name 'GDP' if different).
gdp <- df$GDP
gdp <- as.numeric(df$gdp)

quarters <- seq(as.yearqtr("1970 Q1"), by = 0.25, length.out = length(gdp))
date_index <- as.Date(quarters)  # Convert to Date (this will be the first day of each quarter)

# Create the xts object with the GDP data and the time index
gdp_xts <- xts(gdp, order.by = date_index)

# Display the first few observations
head(gdp_xts)

# ------------------------------------------------
# 2) Plot the real GDP time series
# ------------------------------------------------
plot(gdp_xts,
     main = "Real Euro GDP (1970Q1 - 2014Q4)",
     xlab = "Time (Quarterly)",
     ylab = "GDP (in dataset units)")

# Brief Comment:
# The series trends upward over time, indicating nonstationarity with an increasing mean,
# which is typical of real GDP data.

# ------------------------------------------------
# 3) Compute and plot the quarterly growth rate of real GDP
# ------------------------------------------------
# Compute quarterly growth rate as 100 * log differences
gdp_growth <- diff(log(gdp_xts)) * 100

plot(gdp_growth,
     main = "Quarterly Growth Rate of Real Euro GDP (log diff * 100)",
     xlab = "Time (Quarterly)",
     ylab = "Percent")

# Brief Comment:
# The growth rate series appears more stationary, fluctuating around zero.

# ------------------------------------------------
# 4) Perform Augmented Dickey-Fuller (ADF) tests
# ------------------------------------------------


# ADF test on levels (real GDP)
# Convert the xts object to a numeric vector for the test
adf_level <- adf.test(as.numeric(gdp_xts))
print(adf_level)

# Compute quarterly growth rate as 100 * log differences
gdp_growth <- diff(log(gdp_xts)) * 100

# Remove NA values (often the first element is NA after differencing)
gdp_growth <- na.omit(gdp_growth)

# Load the tseries package and run the Augmented Dickey-Fuller test on the cleaned series
adf_growth <- adf.test(as.numeric(gdp_growth))
print(adf_growth)
