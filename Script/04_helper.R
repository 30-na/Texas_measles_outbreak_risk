library(sf)
library(ggplot2)
library(reshape2)
library(lubridate)
library(dplyr)
library(tidyr)
library(readr)


#------------------ Load datasets -------------------
mmr_county <- readRDS("ProcessedData/map_county.rds")



# Load County boundries
tx_counties <- counties(
  state = "TX",
  year = 2024
) %>%
  mutate(
    county = toupper(NAME)
  )




# Set vaccine efficacy and basic reproduction number
efficacy <- 0.97
R0 <- 18


# --------------Compare PLOS and Lin Formula--------------------
plot_plos <- function(Vj, efficacy = 0.97, R0 = 18) {
  a <- 1 - efficacy * Vj
  X_vals <- seq(0, 1, length.out = 10000)
  f_X <- X_vals - a * (1 - exp(-a * R0 * X_vals))
  
  plot(X_vals, f_X,
       type = "l",
       xlab = "X ",
       ylab = "f(X)",
       main = paste("MMR =", Vj, ", R0 = ", R0, ",  e = 0.97"))
  abline(h = 0, col = "red")
  return(f_X)
}

plot_plos(.9)

plot_lin <- function(v, R0 = 18) {
  a <- 1 - v
  X_vals <- seq(0, 1, length.out = 10000)
  f_X <- X_vals - a * (1 - exp(-R0 * X_vals))
  
  plot(X_vals, f_X, type = "l",
       xlab = "X",
       ylab = "f(X)",
       main = paste("MMR =", v, " R0 = ", R0))
  
  abline(h = 0, col = "red", lwd = 2)
}


plot_lin(0.9)


#---------------------Calculate Probability of infection for county and District

# PLOS
find_internal_infection_PLOS <- function(Vj, efficacy_rate = efficacy, R0_value = R0, tol = 1e-4) {
  if (is.na(Vj)) return(NA_real_)
  a <- 1 - efficacy_rate * Vj
  X_vals <- seq(0, 1, length.out = 10000)
  f_X <- X_vals - a * (1 - exp(-a * R0_value * X_vals))
  close_to_zero <- abs(f_X) < tol
  return(max(X_vals[close_to_zero]))
}

# Interface (Lin. )
find_internal_infection_lin <- function(v, R0_value = R0, tol = 1e-4) {
  if (is.na(v)) return(NA_real_)
  a <- 1 - v
  X_vals <- seq(0, 1, length.out = 10000)
  f_X <- X_vals - a * (1 - exp(-R0_value * X_vals))
  close_to_zero <- abs(f_X) < tol
  return(max(X_vals[close_to_zero]))
}

