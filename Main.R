# Time Series Econometrics: Assignment 2

# Maximilian Droschl (21-616-412)
# Mohamed Ayadi (18-814-996)
# Aram Jivan Chetirian (19-607-258)

# Load packages
library(polynom)

### Task 1

### a)

# Define the polynomials
Theta <- c(1, -4.75, 7.375, -5, 1.5625)
Phi <- c(1, -1, 0.5)

ar_roots <- polyroot(Theta)
print(ar_roots)


ma_roots <- polyroot(Phi)
print(ma_roots)

# Perform polynomial division (Theta / Phi)
Q <- polynom::polynomial(Theta) / 
polynom::polynomial(Phi)
print(Q)

mod_phi <- Mod(phi_roots)
ar_part_roots  <- phi_roots[mod_phi > 1]
ma_part_roots  <- phi_roots[mod_phi < 1]

## Reconstruct polynomials from the chosen roots:
poly_from_roots <- function(r) {
  p <- c(1)
  for (rt in r) {
    p <- convolve(p, c(1, -rt), type="open")
  }
  p
}

ar_poly  <- poly_from_roots(ar_part_roots)
ma_poly  <- poly_from_roots(ma_part_roots)

cat("Reduced AR polynomial:\n");  print(ar_poly)
cat("Reduced MA polynomial:\n");  print(ma_poly)

## (c) The coefficients in 'ma_poly' give the final MA(2) part, i.e. 1 + theta1 B + theta2 B^2.

## (d) Check stationarity and invertibility:
is_stationary   <- all(Mod(ar_part_roots) > 1)
is_invertible   <- all(Mod(ma_part_roots) > 1)

cat("Stationary?  ", is_stationary, "\n")
cat("Invertible?  ", is_invertible, "\n")
