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