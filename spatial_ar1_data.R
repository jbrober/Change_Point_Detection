#'---------------------------------------------------------------------------------------------------------
#' Project Name: spatial_ar1_data.R
#' 
#' Short Description: One change point data generation and functions.
#'
#' Author: Yi Huang 
#'
#' Data Authored: 11/28/2024
#' Edited 11/29/2024 by James Robertson
#'---------------------------------------------------------------------------------------------------------

library(MASS)

#'---------------------------------------------------------------------------------------------------------
#' FUNCTION: Matern
#' PURPOSE:  Matern covariance function
#' @param d
#' @param rho
#' @param nu
#'---------------------------------------------------------------------------------------------------------
Matern <- function(d, rho, nu=0.5) {
  if(nu==0.5){ # nu=0.5 is exponential model
    result <- exp(-d/rho)
  }else{
    d1 <- sqrt(2*nu)*d/rho
    result <- 2^(1-nu)/gamma(nu)*d1^nu*besselK(d1, nu)
  }
  result[d==0] <- 1
  result
}


#'---------------------------------------------------------------------------------------------------------
#' FUNCTION: sim.y
#' PURPOSE:  Simulate one spatio-temporal segment
#' @param theta
#' @param S.dist
#' @param TT
#' @param T.burn
#'---------------------------------------------------------------------------------------------------------
sim.y <- function(theta, S.dist, TT, T.burn = 100) {
  phi <- theta[1]      # Temporal autoregressive coefficient
  rho <- theta[2]      # Range parameter for spatial correlation
  sigma2 <- theta[3]   # Variance of errors
  mu <- ifelse(length(theta) > 3, theta[4], 0) # Optional mean parameter
  
  # Spatial covariance matrix
  Sigma <- sigma2 * Matern(S.dist, rho = rho)
  
  # Initialize data and random errors
  n_locs <- nrow(S.dist)
  y <- matrix(0, nrow = TT + T.burn, ncol = n_locs)
  e <- mvrnorm(n = TT + T.burn, mu = rep(0, n_locs), Sigma = Sigma)
  
  # Simulate temporal autoregressive process
  y[1, ] <- e[1, ]
  for (t in 2:(TT + T.burn)) {
    y[t, ] <- phi * y[t - 1, ] + e[t, ]
  }
  
  # Discard burn-in period and add mean
  y <- y[-(1:T.burn), ] + mu
  return(y)
}


#'---------------------------------------------------------------------------------------------------------
#' FUNCTION: gen_one_change_point_data
#' PURPOSE:  Generate one-change-point spatio-temporal data
#' @param theta
#' @param S.dist
#' @param TT
#' @param T.burn
#'---------------------------------------------------------------------------------------------------------
gen_one_change_point_data <- function(S, TT, theta1, theta2, change_point_ratio = 0.5) {
  # Create a regular 2D grid of spatial locations
  min.dist <- 1
  max.lati <- max.long <- sqrt(S)
  n.lati <- max.lati / min.dist
  n.long <- max.long / min.dist
  mat.lattice <- cbind(rep(seq(1, max.lati, by = min.dist), n.long),
                       rep(seq(1, max.long, by = min.dist), each = n.lati))
  
  # Compute spatial distance matrix
  S.dist <- as.matrix(dist(mat.lattice))
  
  # Calculate segment lengths based on change point ratio
  T1 <- floor(TT * change_point_ratio)  # Length of first segment
  T2 <- TT - T1                        # Length of second segment
  
  # Generate data for each segment
  y1 <- sim.y(theta = theta1, S.dist = S.dist, TT = T1)
  y2 <- sim.y(theta = theta2, S.dist = S.dist, TT = T2)
  
  # Combine segments to create one change-point data
  y <- rbind(y1, y2)
  
  return(list(data = y, spatial_coords = mat.lattice, distance_matrix = S.dist))
}




#'---------------------------------------------------------------------------------------------------------
#' PURPOSE:  Example usage
#' @param S: 3^2
#' @param TT: 15
#' @param theta1
#' @param theta2
#' @param change_point_ratio
#' @param n=S*TT=135 # Note sample size n
#'---------------------------------------------------------------------------------------------------------
if (samp=='small') {
# simple size 60
# set.seed(1128)

S <- 2^2                  # Number of spatial locations (e.g., 6x6 grid)
TT <- 15                  # Total temporal length, 15 time points
theta1 <- c(-0.4, 0.5, 1) # Parameters for first segment
theta2 <- c(-0.2, 0.6, 1) # Parameters for second segment (post-change)
change_point_ratio <- 0.5 # Change occurs at 50% of the time

result <- gen_one_change_point_data(S, TT, theta1, theta2, change_point_ratio)

y <- result$data
spatial_coords <- result$spatial_coords
distance_matrix <- result$distance_matrix

data1 <- data.frame(
  time = rep(1:TT, each = S),
  location = rep(1:S, TT),
  value = as.vector(y)
)

# save(data1, file = "data_n60.RData")
} else {
# simple size 135
# set.seed(1128)

S <- 3^2                  # Number of spatial locations (e.g., 6x6 grid)
TT <- 15                  # Total temporal length, 15 time points
theta1 <- c(-0.4, 0.5, 1) # Parameters for first segment
theta2 <- c(-0.2, 0.6, 1) # Parameters for second segment (post-change)
change_point_ratio <- 0.5 # Change occurs at 50% of the time

result <- gen_one_change_point_data(S, TT, theta1, theta2, change_point_ratio)

y <- result$data
spatial_coords <- result$spatial_coords
distance_matrix <- result$distance_matrix

data2 <- data.frame(
  time = rep(1:TT, each = S),
  location = rep(1:S, TT),
  value = as.vector(y)
)

# save(data1, file = "data_n135.RData")
}

