##-----------------------------------------##
## Simulation on Selection Bias Estimators ##
## Written by: Santiago Gómez-Echeverry    ##
## Last update: 16/02/2023                 ##
##-----------------------------------------##

# To-Do list
# Check COS and Meng_rZS

#### - (I) Working space and packages - ####
    
rm(list = ls()) # - Remove all the objects in the environment
sys <- Sys.info()
options(digts = 3)
# We will assign each one of the different folders into an object so that we can access them easily
fold_code <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/1 - Selection Bias Simulations/2 - Code")
fold_data <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/1 - Selection Bias Simulations/3 - Data")
fold_graphs <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/1 - Selection Bias Simulations/4 - Graphs & Tables")
setwd(fold_code) # For the moment, let us work on the code folder
# getwd() # - To check that we are in the right folder
    
# Below are all the packages that we will use in the analyses. We will check if they are installed, install them if they
# are not, and finally load them.
packages <- c('ggplot2', 'ggpubr', 'forecast', 'reshape2', 'MASS', 'corpcor', 'ggridges', 'ltm', 'stringr', 'knitr', 'kableExtra', 'tidyr', 'dplyr', 'utils')
    
# Installing the packages
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)){
  install.packages(packages[!installed_packages]) 
}
    
# Loading the packages 
invisible(lapply(packages, library, character.only = TRUE))
    
# Besides setting the working space, I'll define in here some basic functions that are needed like the population std.
p.sd <- function(x){
 sqrt(sum((x-mean(x))^2)/length(x))}  
    
#### - (II) Selection Bias Estimators - ####

# This section creates a function for each one of the estimators that we will test on the simulated data sets. It's
# good to keep in mind that this functions will take Y, S, and in most cases Z as an input, and will return a selection
# estimate.

# (1) Meng_rYS
    
Meng_rYS <- function(Y,S){
  t_total <- ncol(Y)
  y <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
  S_Meng <- matrix(data = NA, nrow = t_total, ncol = 1)
  D_I <- D_U <- D_O <- numeric(t_total) 
  for (t in 1:t_total){
    N_t <- length(Y[,t])
    n_t <- sum(S[,t])
    y[,t] <- Y[,t]*S[,t]
    y[,t][y[,t]==0] <- NA 
    # Threefold decomposition
    r_ys <- cor.test(Y[,t],S[,t])$estimate
    sigma_t_Y <- p.sd(Y[,t])
    if (t==1){
      r_ys_hat <- NA
      sigma_t_Y_hat <- NA
    } else if (t>1) {
      r_ys_hat <- cor.test(Y[,t-1],S[,t-1])$estimate
      sigma_t_Y_hat <- p.sd(Y[,t-1])*mean(y[,t], na.rm = TRUE)/mean(y[,t-1], na.rm = TRUE)
    }
    D_I[t] <- r_ys_hat
    D_U[t] <- sigma_t_Y_hat
    D_O[t] <- sqrt((N_t-n_t)/n_t)
    S_Meng[t] <- D_I[t]*D_U[t]*D_O[t]
  }
  return(S_Meng)
}
    
# (2) Meng_rZS
    
Meng_rZS <- function(Y,S,Z){
  t_total <- ncol(Y)
  y <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
  S_Meng <- matrix(data = NA, nrow = t_total, ncol = 1) 
  D_I <- D_U <- D_O <- numeric(t_total) 
  for (t in 1:t_total){
    N_t <- length(Y[,t])
    n_t <- sum(S[,t])
    y[,t] <- Y[,t]*S[,t]
    y[,t][y[,t]==0] <- NA 
    # Threefold decomposition
    r_ys <- cor.test(Y[,t],S[,t])$estimate
    sigma_t_Y <- p.sd(Y[,t])
    r_ys_hat <- cor.test(Z[,t],S[,t])$estimate
    if (t==1){
      sigma_t_Y_hat <- NA
    } else if (t>1){
      sigma_t_Y_hat <- p.sd(Y[,t-1])*mean(y[,t], na.rm = TRUE)/mean(y[,t-1], na.rm = TRUE)
    }
    D_I[t] <- r_ys_hat
    D_U[t] <- sigma_t_Y_hat
    D_O[t] <- sqrt((N_t-n_t)/n_t)
    S_Meng[t] <- D_I[t]*D_U[t]*D_O[t]
  }
  return(S_Meng)
}
    
# (3) MUB_hist
    
MUB_hist <- function(Y,S,Z){
  t_total <- ncol(Y)
  y <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
  z <- matrix(data = NA, nrow = length(Z[,1]), ncol = t_total)
  Delta_Y <- Delta_Z <- rho_zy <- phi <- phi_hat <- phi_hat_norm <- numeric(t_total)
  g_hat_phi <- c_obs <- Y_bar_phi <- numeric(t_total)
  S_MUB <- matrix(data= NA, nrow = t_total, ncol = 1)
  for (t in 1:t_total){
    y[,t] <- Y[,t]*S[,t]
    y[,t][y[,t]==0] <- NA 
    z[,t] <- Z[,t]*S[,t]
    z[,t][z[,t]==0] <- NA 
    # Elements of phi
    Delta_Y[t] <- (mean(y[,t], na.rm = TRUE) - mean(Y[,t]))/sd(y[,t], na.rm = TRUE)
    Delta_Z[t] <- (mean(z[,t], na.rm = TRUE) - mean(Z[,t]))/sd(z[,t], na.rm = TRUE)
    rho_zy[t] <- cor(y[,t], z[,t], use = "complete.obs")
    phi[t] <- (Delta_Y[t] - rho_zy[t]*Delta_Z[t])/
      ((Delta_Y[t]-rho_zy[t]*Delta_Z[t])+(Delta_Z[t]+rho_zy[t]*Delta_Y[t]))
    if (t==1){
      phi_hat[t] <- NA
    } else if (t>1){
      phi_hat[t] <- phi[t-1]
    }
    # Truncation to [0,1]
    phi_hat_norm[t] <- ifelse(phi_hat[t]>0,ifelse(phi_hat[t]<1,phi_hat[t],1),0)
    # Remaining elements of Y(phi)
    g_hat_phi[t] <- (phi_hat_norm[t]+(1-phi_hat_norm[t])*rho_zy[t])/
      (phi_hat_norm[t]*rho_zy[t]+(1-phi_hat_norm[t]))
    c_obs[t] <- sd(y[,t], na.rm = TRUE)/sd(z[,t], na.rm = TRUE)
    # Y(phi) and MUB estimator 
    Y_bar_phi[t] <- mean(y[,t], na.rm = TRUE) - g_hat_phi[t]*c_obs[t]*(mean(z[,t], na.rm = TRUE) - mean(Z[,t]))
    S_MUB[t] <- mean(y[,t], na.rm = TRUE) - Y_bar_phi[t]
  }
  return(S_MUB)
}
    
# (4) MUB_robust
    
MUB_robust <- function(Y,S,Z){
  t_total <- ncol(Y)
  y <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
  z <- matrix(data = NA, nrow = length(Z[,1]), ncol = t_total)
  Delta_Y <- Delta_Z <- rho_zy <- phi <- phi_hat <- phi_hat_norm <- numeric(t_total)
  g_hat_phi <- c_obs <- Y_bar_phi <- numeric(t_total)
  S_MUB <- matrix(data= NA, nrow = t_total, ncol = 1)
  for (t in 1:t_total){
    y[,t] <- Y[,t]*S[,t]
    y[,t][y[,t]==0] <- NA 
    z[,t] <- Z[,t]*S[,t]
    z[,t][z[,t]==0] <- NA 
    # Elements of phi
    Delta_Y[t] <- (median(y[,t], na.rm = TRUE) - mean(Y[,t]))/sd(y[,t], na.rm = TRUE)
    Delta_Z[t] <- (median(z[,t], na.rm = TRUE) - mean(Z[,t]))/sd(z[,t], na.rm = TRUE)
    rho_zy[t] <- cor(y[,t], z[,t], use = "complete.obs")
    phi[t] <- (Delta_Y[t] - rho_zy[t]*Delta_Z[t])/
      ((Delta_Y[t]-rho_zy[t]*Delta_Z[t])+(Delta_Z[t]+rho_zy[t]*Delta_Y[t]))
    if (t==1){
      phi_hat[t] <- NA
    } else if (t>1){
      phi_hat[t] <- phi[t-1]
    }
    # Truncation to [0,1]
    phi_hat_norm[t] <- ifelse(phi_hat[t]>0,ifelse(phi_hat[t]<1,phi_hat[t],1),0)
    # Remaining elements of Y(phi)
    g_hat_phi[t] <- (phi_hat_norm[t]+(1-phi_hat_norm[t])*rho_zy[t])/
      (phi_hat_norm[t]*rho_zy[t]+(1-phi_hat_norm[t]))
    c_obs[t] <- sd(y[,t], na.rm = TRUE)/sd(z[,t], na.rm = TRUE)
    # Y(phi) and MUB estimator 
    Y_bar_phi[t] <- mean(y[,t], na.rm = TRUE) - g_hat_phi[t]*c_obs[t]*(mean(z[,t], na.rm = TRUE) - mean(Z[,t]))
    S_MUB[t] <- mean(y[,t], na.rm = TRUE) - Y_bar_phi[t]
  }
  return(S_MUB)
}
    
# (5) MUB_half
    
MUB_half <- function(Y,S,Z){
  t_total <- ncol(Y)
  y <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
  z <- matrix(data = NA, nrow = length(Z[,1]), ncol = t_total)
  Delta_Y <- Delta_Z <- rho_zy <- phi <- phi_hat <- phi_hat_norm <- numeric(t_total)
  g_hat_phi <- c_obs <- Y_bar_phi <- numeric(t_total)
  S_MUB <- matrix(data= NA, nrow = t_total, ncol = 1)
  for (t in 1:t_total){
    y[,t] <- Y[,t]*S[,t]
    y[,t][y[,t]==0] <- NA 
    z[,t] <- Z[,t]*S[,t]
    z[,t][z[,t]==0] <- NA 
    # Elements of phi
    Delta_Y[t] <- (mean(y[,t], na.rm = TRUE) - mean(Y[,t]))/sd(y[,t], na.rm = TRUE)
    Delta_Z[t] <- (mean(z[,t], na.rm = TRUE) - mean(Z[,t]))/sd(z[,t], na.rm = TRUE)
    rho_zy[t] <- cor(y[,t], z[,t], use = "complete.obs")
    phi[t] <- (Delta_Y[t] - rho_zy[t]*Delta_Z[t])/
      ((Delta_Y[t]-rho_zy[t]*Delta_Z[t])+(Delta_Z[t]+rho_zy[t]*Delta_Y[t]))
    phi_hat[t] <- 0.5
    phi_hat_norm[t] <- phi_hat[t]
    # Remaining elements of Y(phi)
    g_hat_phi[t] <- (phi_hat_norm[t]+(1-phi_hat_norm[t])*rho_zy[t])/
      (phi_hat_norm[t]*rho_zy[t]+(1-phi_hat_norm[t]))
    c_obs[t] <- sd(y[,t], na.rm = TRUE)/sd(z[,t], na.rm = TRUE)
    # Y(phi) and MUB estimator 
    Y_bar_phi[t] <- mean(y[,t], na.rm = TRUE) - g_hat_phi[t]*c_obs[t]*(mean(z[,t], na.rm = TRUE) - mean(Z[,t]))
    S_MUB[t] <- mean(y[,t], na.rm = TRUE) - Y_bar_phi[t]
  }
    return(S_MUB)
}
    
# (6) MUB_abs
    
MUB_abs <- function(Y,S,Z){
  t_total <- ncol(Y)
  y <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
  z <- matrix(data = NA, nrow = length(Z[,1]), ncol = t_total)
  Delta_Y <- Delta_Z <- rho_zy <- phi <- phi_hat <- phi_hat_norm <- numeric(t_total)
  g_hat_phi <- c_obs <- Y_bar_phi <- numeric(t_total)
  S_MUB <- matrix(data= NA, nrow = t_total, ncol = 1)
  for (t in 1:t_total){
    y[,t] <- Y[,t]*S[,t]
    y[,t][y[,t]==0] <- NA 
    z[,t] <- Z[,t]*S[,t]
    z[,t][z[,t]==0] <- NA 
    # Elements of phi
    Delta_Y[t] <- (mean(y[,t], na.rm = TRUE) - mean(Y[,t]))/sd(y[,t], na.rm = TRUE)
    Delta_Z[t] <- (mean(z[,t], na.rm = TRUE) - mean(Z[,t]))/sd(z[,t], na.rm = TRUE)
    rho_zy[t] <- cor(y[,t], z[,t], use = "complete.obs")
    phi[t] <- (Delta_Y[t] - rho_zy[t]*Delta_Z[t])/
      ((Delta_Y[t]-rho_zy[t]*Delta_Z[t])+(Delta_Z[t]+rho_zy[t]*Delta_Y[t]))
    if (t==1){
      phi_hat[t] <- NA
    } else if (t>1){
      phi_hat[t] <- abs(Delta_Y[t-1] - rho_zy[t-1]*Delta_Z[t-1])/
        (abs(Delta_Y[t-1]-rho_zy[t-1]*Delta_Z[t-1])+
           abs(Delta_Z[t-1]+rho_zy[t-1]*Delta_Y[t-1]))
    }
    # Truncation to [0,1]
    phi_hat_norm[t] <- ifelse(phi_hat[t]>0,ifelse(phi_hat[t]<1,phi_hat[t],1),0)
    # Remaining elements of Y(phi)
    g_hat_phi[t] <- (phi_hat[t]+(1-phi_hat_norm[t])*rho_zy[t])/
      (phi_hat_norm[t]*rho_zy[t]+(1-phi_hat_norm[t]))
    c_obs[t] <- sd(y[,t], na.rm = TRUE)/sd(z[,t], na.rm = TRUE)
    # Y(phi) and MUB estimator 
    Y_bar_phi[t] <- mean(y[,t], na.rm = TRUE) - g_hat_phi[t]*c_obs[t]*(mean(z[,t], na.rm = TRUE) - mean(Z[,t]))
    S_MUB[t] <- mean(y[,t], na.rm = TRUE) - Y_bar_phi[t]
  }
  return(S_MUB)
}
    
# (7) MUB_alt
    
MUB_alt <- function(Y,S,Z){
  t_total <- ncol(Y)
  y <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
  z <- matrix(data = NA, nrow = length(Z[,1]), ncol = t_total)
  Delta_Y <- Delta_Z <- rho_zy <- phi <- phi_hat <- phi_hat_norm <- numeric(t_total)
  g_hat_phi <- c_obs <- Y_bar_phi <- numeric(t_total)
  S_MUB <- matrix(data= NA, nrow = t_total, ncol = 1)
  for (t in 1:t_total){
    y[,t] <- Y[,t]*S[,t]
    y[,t][y[,t]==0] <- NA 
    z[,t] <- Z[,t]*S[,t]
    z[,t][z[,t]==0] <- NA 
    ratio <- sd(z[,t], na.rm = TRUE)/sd(y[,t], na.rm = TRUE)
    # Elements of phi
    Delta_Y[t] <- (mean(y[,t], na.rm = TRUE) - mean(Y[,t]))/sd(y[,t], na.rm = TRUE)
    Delta_Z[t] <- (mean(z[,t], na.rm = TRUE) - mean(Z[,t]))/sd(z[,t], na.rm = TRUE)
    rho_zy[t] <- cor(y[,t], z[,t], use = "complete.obs")
    phi[t] <- (Delta_Y[t] - rho_zy[t]*Delta_Z[t])/
      ((Delta_Y[t]-rho_zy[t]*Delta_Z[t])+(Delta_Z[t]+rho_zy[t]*Delta_Y[t]))
    if (t==1){
      phi_hat[t] <- NA
    } else if (t>1){
      phi_hat[t] <- (Delta_Y[t-1] - rho_zy[t-1]*ratio*Delta_Z[t-1])/
        ((Delta_Y[t-1]-rho_zy[t-1]*ratio*Delta_Z[t-1])+
           (Delta_Z[t-1]+rho_zy[t-1]*ratio*Delta_Y[t-1]))
    }
    # Truncation to [0,1]
    phi_hat_norm[t] <- ifelse(phi_hat[t]>0,ifelse(phi_hat[t]<1,phi_hat[t],1),0)
    # Remaining elements of Y(phi)
    g_hat_phi[t] <- (phi_hat_norm[t]+(1-phi_hat_norm[t])*rho_zy[t])/
      (phi_hat_norm[t]*rho_zy[t]+(1-phi_hat_norm[t]))
    c_obs[t] <- sd(y[,t], na.rm = TRUE)/sd(z[,t], na.rm = TRUE)
    # Y(phi) and MUB estimator 
    Y_bar_phi[t] <- mean(y[,t], na.rm = TRUE) - g_hat_phi[t]*c_obs[t]*(mean(z[,t], na.rm = TRUE) - mean(Z[,t]))
    S_MUB[t] <- mean(y[,t], na.rm = TRUE) - Y_bar_phi[t]
  }
  return(S_MUB)
}
    
# (8) MUB_alt_abs
    
MUB_alt_abs <- function(Y,S,Z){
  t_total <- ncol(Y)
  y <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
  z <- matrix(data = NA, nrow = length(Z[,1]), ncol = t_total)
  Delta_Y <- Delta_Z <- rho_zy <- phi <- phi_hat <- phi_hat_norm <- numeric(t_total)
  g_hat_phi <- c_obs <- Y_bar_phi <- numeric(t_total)
  S_MUB <- matrix(data= NA, nrow = t_total, ncol = 1)
  for (t in 1:t_total){
    y[,t] <- Y[,t]*S[,t]
    y[,t][y[,t]==0] <- NA 
    z[,t] <- Z[,t]*S[,t]
    z[,t][z[,t]==0] <- NA 
    ratio <- sd(z[,t], na.rm = TRUE)/sd(y[,t], na.rm = TRUE) 
    # Elements of phi
    Delta_Y[t] <- (mean(y[,t], na.rm = TRUE) - mean(Y[,t]))/sd(y[,t], na.rm = TRUE)
    Delta_Z[t] <- (mean(z[,t], na.rm = TRUE) - mean(Z[,t]))/sd(z[,t], na.rm = TRUE)
    rho_zy[t] <- cor(y[,t], z[,t], use = "complete.obs")
    phi[t] <- (Delta_Y[t] - rho_zy[t]*Delta_Z[t])/
      ((Delta_Y[t]-rho_zy[t]*Delta_Z[t])+(Delta_Z[t]+rho_zy[t]*Delta_Y[t]))
    if (t==1){
      phi_hat[t] <- NA
    } else if (t>1){
      phi_hat[t] <- abs(Delta_Y[t-1] - rho_zy[t-1]*ratio*Delta_Z[t-1])/
        (abs(Delta_Y[t-1]-rho_zy[t-1]*ratio*Delta_Z[t-1])+
           abs(Delta_Z[t-1]+rho_zy[t-1]*ratio*Delta_Y[t-1]))
    }
    # Truncation to [0,1]
    phi_hat_norm[t] <- ifelse(phi_hat[t]>0,ifelse(phi_hat[t]<1,phi_hat[t],1),0)
    # Remaining elements of Y(phi)
    g_hat_phi[t] <- (phi_hat_norm[t]+(1-phi_hat_norm[t])*rho_zy[t])/
      (phi_hat_norm[t]*rho_zy[t]+(1-phi_hat_norm[t]))
    c_obs[t] <- sd(y[,t], na.rm = TRUE)/sd(z[,t], na.rm = TRUE)
    # Y(phi) and MUB estimator 
    Y_bar_phi[t] <- mean(y[,t], na.rm = TRUE ) - g_hat_phi[t]*c_obs[t]*(mean(z[,t], na.rm = TRUE) - mean(Z[,t]))
    S_MUB[t] <- mean(y[,t], na.rm = TRUE) - Y_bar_phi[t]
 }
 return(S_MUB)
}
    
# (9) COS
  
COS <- function(Y,S,Z){
  t_total <- ncol(Y)
  z <- matrix(data = NA, nrow = length(Z[,1]), ncol = t_total)
  S_COS <- matrix(data = NA, nrow = t_total, ncol = 1) 
  for (t in 1:t_total){
    z[,t] <- Z[,t]*S[,t]
    z[,t][z[,t]==0] <- NA
    if (t==1){
      S_COS[t] <- NA
    } else if (t>1){
      r_ys_hat <- cor(Y[,t-1],S[,t-1])
      r_zs <- cor(Z[,t],S[,t])
      S_COS[t] <- (r_ys_hat/r_zs)*(mean(z[,t], na.rm = TRUE) - mean(Z[,t])) 
    }
  }
  return(S_COS)
}
  
# (10) Z_diff
    
Zdiff <- function(Y,S,Z){
  t_total <- ncol(Y)
  z <- matrix(data = NA, nrow = length(Z[,1]), ncol = t_total)
  S_Zdiff <- matrix(data = NA, nrow = t_total, ncol = 1) 
  for (t in 1:t_total){
    z[,t] <- Z[,t]*S[,t]
    z[,t][z[,t]==0] <- NA
    S_Zdiff[t] <- (mean(z[,t], na.rm = TRUE) - mean(Z[,t]))
  }
  return(S_Zdiff)
}
    
#### - (III) Performance functions - ####
    
# (i) Mean Absolute Differences (MAD)
    
# Note that M is the measure obtained with the different methods
    
MAD <- function(Y,S,M){
  t_total <- ncol(Y)
  RS <- RS_hat <- numeric(t_total)
  y <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
  for (t in 1:t_total){
    y[,t] <- Y[,t]*S[,t]
    y[,t][y[,t]==0] <- NA 
    RS[t] <- (mean(y[,t], na.rm = TRUE)- mean(Y[,t], na.rm = TRUE))/
      mean(Y[,t], na.rm = TRUE)
    RS_hat[t] <- M[t]/(mean(y[,t], na.rm = TRUE)-M[t])
  }
  MAD <- mean(abs(RS_hat[-1]-RS[-1]), na.rm = TRUE)
  MAD_list <- list("RS" = RS, "RS_hat"= RS_hat, "MAD" = MAD)
  return(MAD_list)
}
    
# (ii) Root Mean Square Deviation (RMSD)
    
RMSD <- function(Y,S,M){
  t_total <- ncol(Y)
  RS <- RS_hat <- numeric(t_total)
  y <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
  for (t in 1:t_total){
    y[,t] <- Y[,t]*S[,t]
    y[,t][y[,t]==0] <- NA 
    RS[t] <- (mean(y[,t], na.rm = TRUE)- mean(Y[,t], na.rm = TRUE))/
      mean(Y[,t], na.rm = TRUE)
    RS_hat[t] <- M[t]/(mean(y[,t], na.rm = TRUE)-M[t])
  }
  RMSD <- sqrt(mean((RS_hat[-1]-RS[-1])^2, na.rm = TRUE))
  RMSD_list <- list("RS" = RS, "RS_hat" = RS_hat, "RMSD" = RMSD)
  return(RMSD_list)
}
    
#### - (IV) Data generation - ####
    
# (i) Initial settings
    
set.seed(42) # We start by setting the seed so that we can reproduce the whole process every the code is run
  
# Now, we define the initial parameters that we will use to create the data

mu <- 0 # Mean of the normal variable
sigma <- 1 # Deviation of the normal variable
N <- 1000 # Number of observations
Time <- 10 # Number of periods
D <- 200 # Number of draws
const <- rnorm(n = 1, mean = 0, sd = 1)
epsilon_0 <- matrix(data = rnorm(n = N*D, mean = 0, sd = 1), nrow = N, ncol = D)

# (ii) Outcome variable assignation
  
# For all the variable's assignation we will employ three-dimensional arrays, where D1 = observation, D2 = Period, D3 = Draw.
psi <- numeric(Time)
Y <- epsilon <- array(data = NA, dim = c(N, Time, D))
    
# 1. Normal distribution
# First, let's assign the initial period (t = 1). To do this, we create a matrix in which we have the initial observations for
# all the possible draws.

Y_0 <- matrix(data = rnorm(n = N*D, mean = 0, sd = 1), nrow = N, ncol = D)  

# Second, we place the initial observation we've created in the final matrix, and assign the remaining periods
psi <- runif(n = D, min = -1, max = 1)

for (d in 1:D){
  # Place the initial values
  psi_0 <- runif(n= 1, min = -1, max = 1) # Autocorrelation coefficient
  psi[1] <- psi_0 
  Y[,1,d] <- Y_0[,d]; epsilon[,1,d] <- epsilon_0[,d]
    
  # Assigning the remaining periods
  for (p in 2:Time){
   psi[p] <- runif(n= 1, min = -1, max = 1)
   epsilon[,p,d] <- rnorm(n = N, mean = 0, sd = 1) 
   Y[,p,d] <- const + psi[p]*Y[,p-1,d] + epsilon[,p,d] # Here is the AR(1) process
  }
}
dimnames(Y)[[2]] <- paste0("Y", c(1:Time))
    
# 2. Beta distribution
    
# We will follow the same steps that we used before, starting with the variable assignment at the period t = 1. However, 
# since we now have 16 different distributions we will store the corresponding arrays in a list.
    
# Before actually assigning the values, we need to set some parameters and objects we will use
    
al_bt <- expand.grid(1:4,1:4)                               # Grid with all the possible values of alpha and beta
n_bt <- nrow(al_bt)
al_bt_n <- expand.grid(paste0("a",1:4), paste0("b",1:4))
al_bt_n <- paste(al_bt_n[,1], al_bt_n[,2], sep = "_")       # Vector with the beta parameter names
bt_0 <- array(data = matrix(data = NA, nrow = N, ncol = nrow(al_bt)), dim = c(N, nrow(al_bt), D))
  
# Now, we can assign the initial periods
    
for (b in 1:n_bt){
  bt_0[,b,] <- matrix(data = rbeta(n = N*D, shape1 = al_bt[b,1], shape2 = al_bt[b,2]), nrow = N, ncol = D) 
}

# And with the initial periods we can assign the remaining ones just like we did before
    
bt <- epsilon_bt <- vector(mode = "list", length = n_bt)
for (b in 1:n_bt){
  bt[[b]] <- epsilon_bt[[b]] <- array(data = NA, dim = c(N, Time, D))
  for (d in 1:D){
   bt[[b]][,1,d] <- bt_0[,b,d]
   for (p in 2:Time){
     epsilon_bt[[b]][,p,d] <- rnorm(n = N, mean = 0, sd = 1)
     bt[[b]][,p,d] <- const + psi[p]*bt[[b]][,p-1,d] + epsilon_bt[[b]][,p,d]
   }
  }
  dimnames(bt[[b]])[[2]] <- paste0("B", c(1:Time))
}
    
# Let's see how the variables end up being at the first period; we want to check that the distributions are indeed changing in terms of their skewness and kurtosis
    
bt_matrix <- as.data.frame(bt_0[,,1])
colnames(bt_matrix) <- al_bt_n
mbt <- melt(bt_matrix)
  
# A nice trait of the beta distribution is that we can vary the skewness and the kurtosis. Let's calculate their
# values for the cases we constructed
  
# Skewness
al_bt$skew <- (2*(al_bt[,2]-al_bt[,1])*sqrt(al_bt[,1]+al_bt[,2]+1))/
  ((al_bt[,1]*al_bt[,2]+2)*sqrt(al_bt[,1]*al_bt[,2]))
    
# Kurtosis
al_bt$kurt <- (6*((al_bt[,1]-al_bt[,2])^{2}*(al_bt[,1]+al_bt[,2]+1)-al_bt[,1]*al_bt[,2]*(al_bt[,1]+al_bt[,2]+2)))/
  (al_bt[,1]*al_bt[,2]*(al_bt[,1]+al_bt[,2]+2)*(al_bt[,1]+al_bt[,2]+3))
    
# These are the different 16 distributions for the first draw (n = 1)

mbt[,3:4] <- colsplit(mbt$variable, names = c("Alpha", "Beta"), pattern = "_")
mbt$Alpha <- gsub("[^0-9.-]", "", mbt$Alpha); mbt$Beta <- gsub("[^0-9.-]", "", mbt$Beta)
mbt <- mbt[,2:4]

ggplot(mbt, aes(value)) +
  geom_histogram(aes(y=after_stat(density)), 
                 binwidth = 0.05, color = "darkslategray4", fill = "darkslategray3") +
  geom_density(color = "darkslategray4") +
  facet_grid(Beta~Alpha, labeller = label_both) +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  ylab("Density") + xlab("Value ") + theme(text = element_text(size = 20))
ggsave(filename = "Betas.png", path = fold_graphs, width = 30, height = 30, units = "cm")
    
# (iii) Auxiliary variable assignation
    
# Since we have to create a variable that has a specific correlation with Y and we need this correlation to maintain over time, the values are not freely
# assigned. Instead, we are going to use the residuals of a regression to construct them in each period. The following function does this procedure
    
corZ <- function(Y, rho){
  set.seed(1234)
  Y.est <- X <- array(data = NA, dim = dim(Y))
  for (d in 1:D){
   for (p in 1:Time){
    X[,p,d] <- rnorm(n = length(Y[,p,d]), mean = 0, sd = 1)
    Y.est[,p,d] <- residuals(lm(X[,p,d] ~ Y[,p,d]))
    X[,p,d] <- rho*sd(Y.est[,p,d])*Y[,p,d] + Y.est[,p,d]*sd(Y[,p,d])*sqrt(1-rho^2)
   }
  }
  return(X)
}
    
cor_aux <- c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95) # These are the correlations that we will use in the analyses
n_aux <- length(cor_aux)
    
# Keeping in mind that the auxiliary variable depends on the target variable, we will need to create one for each Y we created (including all the betas)
# Correlated with the Y (normal) variable
Z <- vector(mode = "list", length = n_aux)
for (c in 1:n_aux){
  Z[[c]] <- corZ(Y,cor_aux[c])
  dimnames(Z[[c]])[[2]] <- paste0("Z", c(1:Time))
}
  
# Let's check the correlations
corYZ <- matrix(data = NA, nrow = Time, ncol = n_aux) 
    
for (d in 1:D){
  for (p in 1:Time){
    for (c in 1:n_aux){
      corYZ[p,c]<- cor(Y[,p,d],Z[[c]][,p,d], method = "pearson")
    }
  }
  print(colMeans(corYZ))
}

# Correlated with the Y (beta) variables
Z_bt <- vector(mode = "list", length = n_bt)
for (b in 1:n_bt){
 Z_bt[[b]] <- vector(mode = "list", length = n_aux)
 for (c in 1:n_aux){
   Z_bt[[b]][[c]] <- corZ(bt[[b]],cor_aux[c])
   dimnames(Z_bt[[b]][[c]])[[2]] <- paste0("Z", c(1:Time))
 }
}
names(Z_bt) <- paste0("Z_b",1:n_bt) # If we name the list within the list it will be 'easier' to access them

# Let's check the correlations with the Y (beta) variables
corBZ <- matrix(data = NA, nrow = Time, ncol = n_aux) 
for (b in 1:n_bt){
  cat("Correlations with Beta:",b)
  for (d in 1:D){
    for (p in 1:Time){
      for (c in 1:n_aux){
        corBZ[p,c]<- cor(bt[[b]][,p,d],Z_bt[[b]][[c]][,p,d], method = "pearson")
      }
    }
    print(colMeans(corBZ))
  }
}
    
# (iv) Selection variable assignation
    
corS <- function(Y, rho){
  set.seed(4321)
  Y.est <- X <- array(data = NA, dim = dim(Y))
  for (d in 1:D){
    for (p in 1:Time){
      X[,p,d] <- sample(x = c(0,1), size = length(Y[,p,d]), replace = TRUE)
      Y.est[,p,d] <- residuals(glm(X[,p,d] ~ Y[,p,d], family = "binomial"))
      X[,p,d] <- rho*sd(Y.est[,p,d])*Y[,p,d] + Y.est[,p,d]*sd(Y[,p,d])*sqrt(1-rho^2)
      X[,p,d] <- (X[,p,d] - min(X[,p,d]))/(max(X[,p,d])-min(X[,p,d]))
      X[,p,d]<- ifelse(X[,p,d]>median(X[,p,d]),1,0)
    }
  }
  return(X)
}
    
cor_sel <- c(0.513, 0.649, 0.801, 0.999)
n_sel <- length(cor_sel)
    
# Correlated with the Y (normal) variable
S <- vector(mode = "list", length = n_sel)
  
for (c in 1:n_sel){
  S[[c]] <- corS(Y,cor_sel[c])
  dimnames(S[[c]])[[2]] <- paste0("S", c(1:Time))
}
    
# Let's check the correlations
corYS <- matrix(data = NA, nrow = Time, ncol = n_sel) 
    
for (d in 1:D){
  for (p in 1:Time){
    for (c in 1:n_sel){
      corYS[p,c]<- cor.test(Y[,p,d],S[[c]][,p,d])$estimate
    }
  }
  print(colMeans(corYS))
}
    
# Correlated with the Y (beta) variables
S_bt <- vector(mode = "list", length = n_bt)
for (b in 1:n_bt){
  S_bt[[b]] <- vector(mode = "list", length = n_sel)
  for (c in 1:n_sel){
    S_bt[[b]][[c]] <- corS(bt[[b]],cor_sel[c])
    dimnames(S_bt[[b]][[c]])[[2]] <- paste0("S", c(1:Time))
  }
}
names(S_bt) <- paste0("S_b",1:n_bt) # If we name the list within the list it will be 'easier' to access them

# Let's check the correlations with the Y (beta) variables
corBS <- matrix(data = NA, nrow = Time, ncol = n_sel) 
for (b in 1:n_bt){
  cat("Correlations S with Beta:", b)
  for (d in 1:D){
    for (p in 1:Time){
      for (c in 1:n_sel){
        corBS[p,c]<- cor.test(bt[[b]][,p,d],S_bt[[b]][[c]][,p,d])$estimate
      }
    }
    print(colMeans(corBS))
  }
}

#### - (V) Estimations - ####

Vars <- c("S","Y","M","Alpha","Beta","Z", paste0(c("MAD_D","RMSD_D"), rep(1:D, each = 2)))
Performance <- data.frame(matrix(nrow = 0, ncol = length(Vars)))
colnames(Performance) <- Vars
Est <- c("Meng_rYS","Meng_rZS", "MUB_half", "MUB_hist", "MUB_robust", "MUB_alt", "MUB_abs", "MUB_alt_abs", "COS", "Zdiff")
n_Est <- length(Est)

# i. Normal distribution 

Sel <- array(data = NA, dim = c(Time, n_sel, D))            # Array that contains the true selection given different Cov(YS)
Perf <- vector(mode = "list", length = n_Est)               # List with all the estimator's performance given a specific Cov(YZ)
Sel_hat <- Est_hat <- vector(mode = "list", length = n_Est) # A list for the all the relative selection and the selection estimates

pb <- txtProgressBar(min = 0, max = n_aux, initial = 0, style = 3)    # Set the progress bar to see how the loop is advancing
for (cz in 1:n_aux){
  setTxtProgressBar(pb,cz) # Initiate the progress bar
  for (m in 1:n_Est){
    Sel_hat[[m]] <- Est_hat[[m]] <- array(data = NA, dim = c(Time, n_sel, D))     # Array with the relative selection and the selection estimates
    Perf[[m]] <- as.data.frame(matrix(data = NA, nrow = n_sel, ncol = (2*D)))     # Data frame with the estimator's performance.
    for (d in 1:D){
      for (cs in 1:n_sel){
        Y_temp <- as.matrix(Y[,,d])                                               # Temporary Y matrix 
        S_temp <- as.matrix(S[[cs]][,,d])                                         # Temporary S matrix
        Z_temp <- Z[[cz]][,,d]                                                    # Temporary Z matrix
        
        # Since the estimation for Meng_rYS is different (i.e., doesn't include Z) we will use a conditional
        if (Est[m]=="Meng_rYS"){
          Est_temp <- Est_hat[[m]][,cs,d] <- get(Est[m])(Y_temp, S_temp)          # Vector with the selection estimates
          Sel[,cs,d] <- as.vector(MAD(Y_temp, S_temp, Est_temp)[[1]])             # Vector with the true relative selection
        } else if (Est[m]!="Meng_rYS"){
          Est_temp <- Est_hat[[m]][,cs,d] <- get(Est[m])(Y_temp, S_temp, Z_temp)  # Vector with the selection estimates
        }
        Perf[[m]][cs,(d*2)-1] <- MAD(Y_temp, S_temp, Est_temp)[[3]]               # MAD constructed with Y, S and the selection estimates
        Perf[[m]][cs,(d*2)] <- RMSD(Y_temp, S_temp, Est_temp)[[3]]                # RMSD constructed with Y, S and the selection estimates
        Sel_hat[[m]][,cs,d] <- as.vector(MAD(Y_temp, S_temp, Est_temp)[[2]])      # Vector with the estimated relative selection
      }
    }
    # Let us asssign the results into a data frame 
    Perf[[m]] <- as.data.frame(cbind(as.vector(paste0("S", 1:n_sel)),
                                       t(matrix(data = rep(c("Y", Est[m], NA, NA, 
                                                             ifelse(Est[m]=="Meng_rYS", NA,  paste0("Z",cz))), n_sel), nrow = 5)), Perf[[m]]))
    colnames(Perf[[m]]) <- Vars # We need to add the proper names to the columns of the data frame
  }
  Performance <- rbind(Performance, bind_rows(Perf, .id = "column_label"))        # Bind the results from every estimator at a given Cov(YZ) with the preiviously obtained results
  assign(paste0("Sel_hat_Z",cz), Sel_hat)                                         # We want to store the selection results for each different Cov(YZ)
  close(pb) # Close the progress bar
}

# ii. Beta distribution
pb <- txtProgressBar(min = 0, max = n_bt, initial = 0, style = 3) 
for (b in 1:n_bt){
  setTxtProgressBar(pb,b)
  Est_hat <- vector(mode = "list", length = n_Est)
  Perf <- vector(mode = "list", length = n_Est)
  for (cz in 1:n_aux){
    for (m in 1:n_Est){
      Est_hat[[m]] <- array(data = NA, dim = c(Time, n_sel, D))                  # Array with the selection estimates
      Perf[[m]] <- as.data.frame(matrix(data = NA, nrow = n_sel, ncol = (2*D)))  # Data frame with the estimator's performance.
      for (d in 1:D){
        for (cs in 1:n_sel){
          B_temp <- as.matrix(bt[[b]][,,d])                                      # Temporary Y matrix 
          S_temp <- as.matrix(S_bt[[b]][[cs]][,,d])                              # Temporary S matrix 
          Z_temp <- Z_bt[[b]][[cz]][,,d]                                         # Temporary Z matrix 
          
          # Since the estimation for Meng_rYS is different (i.e., doesn't include Z) we will use a conditional
          if (Est[m]=="Meng_rYS"){
            Est_temp <- Est_hat[[m]][,cs,d] <- get(Est[m])(B_temp, S_temp)       # Vector with the selection estimates
          } else if (Est[m]!="Meng_rYS"){
            Est_temp <- Est_hat[[m]][,cs,d] <- get(Est[m])(B_temp, S_temp, Z_temp) # Vector with the selection estimates
          } 
          Perf[[m]][cs,(d*2)-1] <- MAD(B_temp, S_temp, Est_temp)[[3]]            # MAD constructed with Y, S and the selection estimates
          Perf[[m]][cs,(d*2)] <- RMSD(B_temp, S_temp, Est_temp)[[3]]             # RMSD constructed with Y, S and the selection estimates
        }
      }
      # Let us asssign the results into a data frame 
      Perf[[m]] <- as.data.frame(cbind(as.vector(paste0("S", 1:n_sel)),
                                         t(matrix(data = rep(c(paste0("beta",b), Est[m], al_bt$Var1[b], al_bt$Var2[b],
                                                               ifelse(Est[m]=="Meng_rYS", NA,  paste0("Z",cz))), n_sel), nrow = 5)), Perf[[m]]))
      colnames(Perf[[m]]) <- Vars                                                # We need to add the proper names to the columns of the data frame
    }
    Performance <- rbind(Performance, bind_rows(Perf, .id = "column_label"))     # Bind the results from every estimator at a given Cov(YZ) with the preiviously obtained results
  }
  close(pb)
}

Performance <- select(Performance, -column_label)                              # Getting rid of the additional variable created when using bind_rows()
mPerformance <- melt(Performance, id.vars = c("S","M","Y","Alpha","Beta","Z")) # This is the melted (reshaped) data set. We are using this shape since it's easier to handle.

# Before saving the results, let us arrange the variables so that it's clear what we are seeing when we handle 
# the data

# Recoding the estimators
mPerformance$M <- mPerformance$M %>%
  recode("Meng_rYS" = "a. Meng_rYS", "Meng_rZS" = "b. Meng_rZS", "MUB_half" = "c. MUB_half", 
         "MUB_hist" = "d. MUB_hist", "MUB_robust" = "e. MUB_robust", "MUB_alt" = "f. MUB_alt", 
         "MUB_abs" = "g.MUB_abs", "MUB_alt_abs" = "h. MUB_alt_abs", "COS" = "i. COS", "Zdiff" = "j. Zdiff")
# Recoding the correlations with the auxiliary variable
mPerformance$Z <- mPerformance$Z %>%
  recode("Z1" = "Cov(YZ) = 0.05", "Z2" = "Cov(YZ) = 0.2",
         "Z3" = "Cov(YZ) = 0.4", "Z4" = "Cov(YZ) = 0.6",
         "Z5" = "Cov(YZ) = 0.8", "Z6" = "Cov(YZ) = 0.95")
# Recoding the correlations with the sampling variable
mPerformance$S <- mPerformance$S %>%
  recode("S1" = "Cov(YS) = 0.2", "S2" = "Cov(YS) = 0.4",
         "S3" = "Cov(YS) = 0.6", "S4" = "Cov(YS) = 0.8")
write.table(mPerformance, file = paste0(fold_data,"/Simulated_Data_N", N, "_D", D, "_T", Time, ".txt"))

#### - (VI) Graphs - ####

mPerformance <- read.table(paste0(fold_data,"/Simulated_Data_N", N, "_D", D, "_T", Time, ".txt"))
Prf_M <- c("MAD", "RMSD")
setwd(fold_graphs) # We will store all of the graphs and tables in this folder

for (i in 1:length(Prf_M)){
  # (i) Different S, constant Z - To see how the estimators react to changes in S
  Plot_S <- mPerformance %>%
    filter((Z == "Cov(YZ) = 0.4" | is.na(Z))  & grepl("Y", Y) & grepl(Prf_M[i], variable)) %>% # Selecting only the results that we want to plot 
    ggplot(aes(x = value, y = M, fill = after_stat(x))) +                        # General aesthetic
    geom_density_ridges_gradient(rel_min_height = 0.01) +                        # The type of plot that we will use
    facet_grid(~S) + xlim(0,10) +                                                # Make several panels with specific limits
    ylab("Estimator") + xlab(Prf_M[i])  +                                        # The labels of each axis
    scale_fill_viridis_c(name = Prf_M[i], option = "D") +                        # The fill of every plot 
    theme(text = element_text(size = 20))                                        # General layout options: Used to adjust the size of the text
  ggsave(paste0("Plot_S_",Prf_M[i], ".png"), plot = Plot_S, path = fold_graphs, width = 40, height = 15, units = "cm")

  # (ii) Different Z, constant S - To see how the estimators react to changes in Z

  Plot_Z <- mPerformance %>%
    filter(S == "Cov(YS) = 0.6" & grepl("Y", Y) & grepl(Prf_M[i], variable) & Z != "Cov(YZ) = 0.05"  & Z!= "Cov(YZ) = 0.95") %>% # Selecting only the results that we want to plot
    ggplot(aes(x = value, y = M, fill = after_stat(x))) +                                                                        # General aesthetic
    geom_density_ridges_gradient(rel_min_height = 0.01) +                                                                        # The type of plot that we will use
    facet_grid(~Z) + xlim(0,5) +                                                                                                 # Make several panels with specific limits
    ylab("Estimator") + xlab(Prf_M[i])  +                                                                                        # The labels of each axis
    scale_fill_viridis_c(name = Prf_M[i], option = "D") +                                                                        # The fill of every plot 
    theme(text = element_text(size = 20))                                                                                        # General layout options: Used to adjust the size of the text
  ggsave(paste0("Plot_Z_", Prf_M[i], ".png"), plot = Plot_Z, path = fold_graphs, width = 40, height = 15, units = "cm")
  
  # (iii) Different Z extreme values, constant S - To see how the estimators react to changes in Z
  
  Plot_Z <- mPerformance %>%
    filter(S == "Cov(YS) = 0.6" & grepl("Y", Y) & grepl(Prf_M[i], variable) & Z != "Cov(YZ) = 0.2"  & Z!= "Cov(YZ) = 0.8") %>% # Selecting only the results that we want to plot
    ggplot(aes(x = value, y = M, fill = after_stat(x))) +                                                                      # General aesthetic
    geom_density_ridges_gradient(rel_min_height = 0.01) +                                                                      # The type of plot that we will use
    facet_grid(~Z) + xlim(0,5) +                                                                                               # Make several panels with specific limits
    ylab("Estimator") + xlab(Prf_M[i])  +                                                                                      # The labels of each axis
    scale_fill_viridis_c(name = Prf_M[i], option = "D") +                                                                      # The fill of every plot
    theme(text = element_text(size = 20))                                                                                      # General layout options: Used to adjust the size of the text
  ggsave(paste0("Plot_Z_", Prf_M[i], "_ext.png"), plot = Plot_Z, path = fold_graphs, width = 40, height = 15, units = "cm")
  
  # (iv) Different Y, ceteris paribus - To check how the estimators react to changes in Y distribution

  Plot_B <- mPerformance %>%
    filter(Z == "Cov(YZ) = 0.6", S == "Cov(YS) = 0.4" & !is.na(Alpha) & grepl(Prf_M[i], variable)) %>% # Selecting only the results that we want to plot
    ggplot(aes(x = value, y = M, fill = after_stat(x))) +                                              # General aesthetic
    geom_density_ridges_gradient(rel_min_height = 0.01) +                                              # The type of plot that we will use
    facet_grid(Beta~Alpha, labeller = label_both) + xlim(0,5) +                                        # Make several panels with specific limits
    ylab("Estimator") + xlab(Prf_M[i]) +                                                               # The labels of each axis
    scale_fill_viridis_c(name = Prf_M[i], option = "D") +                                              # The fill of every plot
    theme(text = element_text(size = 20))                                                              # General layout options: Used to adjust the size of the text
  ggsave(paste0("Plot_B_", Prf_M[i], ".png"), plot = Plot_B, path = fold_graphs, width = 50, height = 50, units = "cm")

  # (v) Different Y, with different S

  Plot_BS <- mPerformance %>%
    filter(Z == "Cov(YZ) = 0.6" & !is.na(Alpha) & grepl(Prf_M[i], variable) & Beta == 1) %>% # Selecting only the results that we want to plot
    ggplot(aes(x = value, y = M, fill = after_stat(x))) +                                    # General aesthetic
    geom_density_ridges_gradient(rel_min_height = 0.01) +                                    # The type of plot that we will use
    facet_grid(S~Alpha, labeller = label_both) + xlim(0,5) +                                 # Make several panels with specific limits
    ylab("Estimator") + xlab(Prf_M[i]) +                                                     # The labels of each axis
    scale_fill_viridis_c(name = Prf_M[i], option = "B") +                                    # The fill of every plot
    theme(text = element_text(size = 20))                                                    # General layout options: Used to adjust the size of the text
  ggsave(paste0("Plot_BS_", Prf_M[i], ".png"), plot = Plot_BS, path = fold_graphs, width = 50, height = 50, units = "cm")

  # (vi) Different Y, with different Z
  
  Plot_BZ <- mPerformance %>%
    filter(S == "Cov(YS) = 0.4" & !is.na(Alpha) & !is.na(Z) & grepl(Prf_M[i], variable) & Beta == 1 &
             Z != "Cov(YZ) = 0.05"  & Z!= "Cov(YZ) = 0.95") %>%                              # Selecting only the results that we want to plot
    ggplot(aes(x = value, y = M, fill = after_stat(x))) +                                    # General aesthetic
    geom_density_ridges_gradient(rel_min_height = 0.01) +                                    # The type of plot that we will use
    facet_grid(Z~Alpha, labeller = label_both) + xlim(0,5) +                                 # Make several panels with specific limits
    ylab("Estimator") + xlab(Prf_M[i]) +                                                     # The labels of each axis
    scale_fill_viridis_c(name = Prf_M[i], option = "B") +                                    # The fill of every plot
    theme(text = element_text(size = 20))                                                    # General layout options: Used to adjust the size of the text
  ggsave(paste0("Plot_BZ_", Prf_M[i], ".png"), plot = Plot_BZ, path = fold_graphs, width = 50, height = 50, units = "cm")
}

# How are the sample mean and the median of these estimates?

Tab_Av <- mPerformance %>%
  filter(grepl("Y", Y) & grepl("MAD", variable) & S == "Cov(YS) = 0.4" & !is.na(Z)) %>%
  group_by(M, Z) %>%
  summarize(MeanValue=mean(value)) %>%
  spread(Z, MeanValue) %>%
  kable(format = "latex")
save_kable(Tab_Av, file = "Table_Means_Z.pdf")

Tab_Med <- mPerformance %>%
  filter(grepl("Y", Y) & grepl("MAD", variable) & S == "Cov(YS) = 0.4" & !is.na(Z)) %>%
  group_by(M, Z) %>%
  summarize(MeanValue=median(value)) %>%
  spread(Z, MeanValue) %>%
  kable(format = "latex")
save_kable(Tab_Med, file = "Table_Means_Z.pdf")

# Graphs to check the selection

# True Selection
AvSel <- Sel[,,] %>% 
  abs() %>%            # Take the absolute value
  colMeans() %>%       # Average across period
  rowMeans() %>%       # Average across draw
  cbind(paste0("Cov(YS) = ", seq(from = 0.2, to = 0.8, by = 0.2))) %>% 
  as.data.frame() %>% 
  rename("Sel" = 1 , "S" = 2)

# Estimated Selection
TempAvSel_Z <- vector(mode = "list", length = n_Est)
AvSel_hat <- data.frame(matrix(nrow = 0, ncol = 4))
for (j in 1:n_aux){
  for (i in 1:n_Est){
    TempAvSel_Z[[i]] <- get(paste0("Sel_hat_Z",j))[[i]][,,] %>% 
      abs() %>%                        # Take the absolute value
      colMeans(na.rm = TRUE) %>%       # Average across period
      rowMeans(na.rm = TRUE) %>%       # Average across draw 
      cbind(paste0("Cov(YS) = ", seq(from = 0.2, to = 0.8, by = 0.2))) %>% 
      as.data.frame() %>% 
      rename("Sel" = 1 , "S" = 2) %>% 
      mutate(M = Est[i], Z = paste0("Cov(YZ) = ", cor_aux[j]))
  }
  AvSel_hat <- rbind(AvSel_hat, bind_rows(TempAvSel_Z, .id = "column_label"))
}

ggplot(data = AvSel_hat, aes(x = M, y = as.numeric(Sel), fill = M)) + 
  geom_col() +
  facet_grid(Z~S) +
  coord_flip() +
  geom_hline(data =AvSel, aes(yintercept=abs(as.numeric(Sel))), col = "blue", linetype = "dashed") +
  theme(strip.background = element_blank(), strip.placement = "outside", legend.position="none") +
  ylab("Selectivity") + xlab("Estimator")  + theme(text = element_text(size = 15)) +
  scale_fill_brewer(palette="Spectral")
ggsave("Plot_Sel.png", path = fold_graphs, width = 50, height = 50, units = "cm")
