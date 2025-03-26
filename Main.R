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