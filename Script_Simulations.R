##-----------------------------------------##
## Simulation on Selection Bias Estimators ##
## Written by: Santiago Gómez-Echeverry    ##
## Last update: 17/11/2022                 ##
##-----------------------------------------##
    
#### - (I) Working space and packages - ####
    
rm(list = ls()) # - Remove all the objects in the environment
sys <- Sys.info()
memory.limit(size = 20000)
fold_code <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/1 - Selection Bias Simulations/2 - Code")
fold_data <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/1 - Selection Bias Simulations/3 - Data")
fold_graphs <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/1 - Selection Bias Simulations/4 - Graphs & Tables")
setwd(fold_code)
# getwd() # - To check that we are in the right folder
    
packages <- c('ggplot2', 'ggpubr', 'forecast', 'reshape2', 'MASS', 'corpcor', 'ggridges', 'ltm', 'stringr')
    
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
    r_ys <- cor(Y[,t],S[,t])
    sigma_t_Y <- p.sd(Y[,t])
    if (t==1){
      r_ys_hat <- NA
      sigma_t_Y_hat <- NA
    } else if (t>1) {
      r_ys_hat <- cor(Y[,t-1],S[,t-1])
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
    r_ys <- cor(Y[,t],S[,t])
    sigma_t_Y <- p.sd(Y[,t])
    if (t==1){
      r_ys_hat <- NA
      sigma_t_Y_hat <- NA
    } else if (t>1){
      r_ys_hat <- cor(Z[,t],S[,t])
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
    r_ys <- cor(Y[,t],S[,t])
    r_zs <- cor(Z[,t],S[,t])
    S_COS[t] <- (r_ys/r_zs)*(mean(z[,t], na.rm = TRUE) - mean(Z[,t]))
  }
  return(S_COS)
}
  
# (10) Zdiff
    
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
    RS[t] <- (mean(y[,1], na.rm = TRUE)- mean(Y[,t], na.rm = TRUE))/
      mean(Y[,t], na.rm = TRUE)
    RS_hat[t] <- M[t]/(mean(y[,t], na.rm = TRUE)-M[t])
  }
  MAD <- mean(abs(RS_hat-RS), na.rm = TRUE)
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
    RS[t] <- (mean(y[,1], na.rm = TRUE)- mean(Y[,t], na.rm = TRUE))/
      mean(Y[,t], na.rm = TRUE)
    RS_hat[t] <- M[t]/(mean(y[,t], na.rm = TRUE)-M[t])
  }
  RMSD <- sqrt(mean((RS_hat-RS)^2, na.rm = TRUE))
  RMSD_list <- list("RS" = RS, "RS_hat" = RS_hat, "RMSD" = RMSD)
  return(RMSD_list)
}
    
#### - (IV) Data generation - ####
    
# (i) Initial settings
    
set.seed(1234) # We start by setting the seed so that we can reproduce the whole process every the code is run
  
# Now, we define the initial parameters that we will use to create the data
    
mu <- 0 # Mean of the normal variable
sigma <- 1 # Deviation of the normal variable
N <- 1000 # Number of observations
Time <- 10 # Number of periods
D <- 500 # Number of draws
psi_0 <- runif(n= 1, min = -1, max = 1) # Autocorrelation coefficient
const <- rnorm(n = 1, mean = 0, sd = 1)
epsilon_0 <- matrix(data = NA, nrow = N, ncol = D)
for (d in 1:D){
  epsilon_0[,d] <- rnorm(n = N, mean = 0, sd = 1)
}
    
# (ii) Outcome variable assignation
  
# For all the variable's assignation we will employ three-dimensional arrays, where D1 = observation, D2 = Period, D3 = Draw.
psi <- numeric(Time)
Y <- epsilon <- array(data = NA, dim = c(N, Time, D))
    
# 1. Normal distribution
    
# First, let's assign the initial period (t = 1). To do this, we create a matrix in which we have the initial observations for
# all the possible draws.
    
Y_0 <- matrix(data = NA, nrow = N, ncol = D)
for (d in 1:D){
  Y_0[,d] <- rnorm(n = N, mean = mu, sd = sigma)
}
    
# Second, we place the initial observation we've created in the final matrix, and assign the remaining periods

psi[1] <- psi_0
for (d in 1:D){
  # Place the initial values
  Y[,1,d] <- Y_0[,d]; epsilon[,1,d] <- epsilon_0[,d]
    
  # Assigning the remaining periods
  for (p in 2:Time){
   psi[p] <- psi_0^p
   epsilon[,p,d] <- rnorm(n = N, mean = 0, sd = 1) 
   Y[,p,d] <- const + psi[p]*Y[,p-1,d] + epsilon[,p,d] # Here is the AR(1) process
  }
}
dimnames(Y)[[2]] <- paste0("Y", c(1:Time))
    
# To finalize, let's see the variable at the first period and how it evolves over time
    
## ! - Include some graphics on the normal distributed Y
    
# 2. Beta distribution
    
# We will follow the same steps that we used before, starting with the variable assignment at the period t = 1. However, 
# since we now have 16 different distributions we will store the corresponding arrays in a list.
    
# Before actually assigning the values, we need to set some parameters and objects we will use
    
alpha_beta <- expand.grid(1:4,1:4)
alpha_beta_names <- expand.grid(paste0("a",1:4), paste0("b",1:4))
alpha_beta_names <- paste(alpha_beta_names[,1], alpha_beta_names[,2], sep = "_") # Vector with the beta parameter names
betas_0 <- array(data = matrix(data = NA, nrow = N, ncol = nrow(alpha_beta)), dim = c(N, nrow(alpha_beta), D))
  
# Now, we can assign the initial periods
    
for (d in 1:D){
  for (b in 1:nrow(alpha_beta)){
   betas_0[,b,d] <- rbeta(n = N, shape1 = alpha_beta[b,1], shape2 = alpha_beta[b,2]) 
  }
}
    
# And with the initial periods we can assign the remaining ones just like we did before
    
betas <- epsilon_betas <- vector(mode = "list", length = nrow(alpha_beta))
for (b in 1:length(betas)){
  betas[[b]] <- epsilon_betas[[b]] <- array(data = NA, dim = c(N, Time, D))
  for (d in 1:D){
   betas[[b]][,1,d] <- betas_0[,b,d]
   for (p in 2:Time){
     epsilon_betas[[b]][,p,d] <- rnorm(n = N, mean = 0, sd = 1)
     betas[[b]][,p,d] <- const + psi[p]*betas[[b]][,p-1,d] + epsilon_betas[[b]][,p,d]
   }
  }
  dimnames(betas[[b]])[[2]] <- paste0("B", c(1:Time))
}
    
# Let's see how the variables end up being
    
beta_matrix <- as.data.frame(betas_0[,,1])
colnames(beta_matrix) <- alpha_beta_names
mbetas <- melt(beta_matrix)
  
# A nice trait of the beta distribution is that we can vary the skewness and the kurtosis. Let's calculate their
# values for the cases we constructed
  
# Skewness
alpha_beta$skew <- (2*(alpha_beta[,2]-alpha_beta[,1])*sqrt(alpha_beta[,1]+alpha_beta[,2]+1))/
  ((alpha_beta[,1]*alpha_beta[,2]+2)*sqrt(alpha_beta[,1]*alpha_beta[,2]))
    
# Kurtosis
alpha_beta$kurt <- (6*((alpha_beta[,1]-alpha_beta[,2])^{2}*(alpha_beta[,1]+alpha_beta[,2]+1)-alpha_beta[,1]*alpha_beta[,2]*(alpha_beta[,1]+alpha_beta[,2]+2)))/
  (alpha_beta[,1]*alpha_beta[,2]*(alpha_beta[,1]+alpha_beta[,2]+2)*(alpha_beta[,1]+alpha_beta[,2]+3))
    
# These are the different 16 distributions for the first draw (n = 1)

mbetas[,3:4] <- colsplit(mbetas$variable, names = c("Alpha", "Beta"), pattern = "_")
mbetas$Alpha <- gsub("[^0-9.-]", "", mbetas$Alpha); mbetas$Beta <- gsub("[^0-9.-]", "", mbetas$Beta)
mbetas <- mbetas[,2:4]

    
ggplot(mbetas, aes(value)) +
  geom_histogram(aes(y=..density..), 
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
    
cor_aux <- c(0.2, 0.4, 0.6, 0.8) # These are the correlations that we will use in the analyses
    
# Keeping in mind that the auxiliary variable depends on the target variable, we will need to create one for each Y we created (including all the betas)
# Correlated with the Y (normal) variable
Z <- vector(mode = "list", length = length(cor_aux))
for (c in 1:length(cor_aux)){
  Z[[c]] <- corZ(Y,cor_aux[c])
  dimnames(Z[[c]])[[2]] <- paste0("Z", c(1:Time))
}
  
# Let's check the correlations
corYZ <- matrix(data = NA, nrow = Time, ncol = length(cor_aux)) 
    
for (d in 1:D){
  for (p in 1:Time){
    for (c in 1:length(cor_aux)){
      corYZ[p,c]<- cor(Y[,p,d],Z[[c]][,p,d], method = "pearson")
    }
  }
  print(mean(corYZ[,1])); print(mean(corYZ[,2])); print(mean(corYZ[,3]));
  print(mean(corYZ[,4]))
}
    
# Correlated with the Y (beta) variables
Z_betas <- vector(mode = "list", length = length(betas))
for (b in 1:length(betas)){
 Z_betas[[b]] <- vector(mode = "list", length = length(cor_aux))
 for (c in 1:length(cor_aux)){
   Z_betas[[b]][[c]] <- corZ(betas[[b]],cor_aux[c])
   dimnames(Z_betas[[b]][[c]])[[2]] <- paste0("Z", c(1:Time))
 }
}
names(Z_betas) <- paste0("Z_b",1:16) # If we name the list within the list it will be 'easier' to access them
    
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
    
# Correlated with the Y (normal) variable
S <- vector(mode = "list", length = length(cor_sel))
  
for (c in 1:length(cor_sel)){
  S[[c]] <- corS(Y,cor_sel[c])
  dimnames(S[[c]])[[2]] <- paste0("S", c(1:Time))
}
    
# Let's check the correlations
corYS <- matrix(data = NA, nrow = Time, ncol = length(cor_sel)) 
    
for (d in 1:D){
  for (p in 1:Time){
    for (c in 1:length(cor_sel)){
      corYS[p,c]<- cor.test(Y[,p,d],S[[c]][,p,d])$estimate
    }
  }
  print(mean(corYS[,1])); print(mean(corYS[,2])); print(mean(corYS[,3]));
  print(mean(corYS[,4]))
}
    
# Correlated with the Y (beta) variables
S_betas <- vector(mode = "list", length = length(betas))
for (b in 1:length(betas)){
  S_betas[[b]] <- vector(mode = "list", length = length(cor_sel))
  for (c in 1:length(cor_aux)){
    S_betas[[b]][[c]] <- corS(betas[[b]],cor_aux[c])
    dimnames(S_betas[[b]][[c]])[[2]] <- paste0("S", c(1:Time))
  }
}
names(S_betas) <- paste0("S_b",1:16) # If we name the list within the list it will be 'easier' to access them
    
# (v) Normalization of the variables
    
# 1. Outcome variables
  
# 1.1. Normally distributed Y
nY <- array(data = NA, dim = c(N, Time, D))
    
for (d in 1:D){
  for (p in 1:Time){
    nY[,p,d] <- (Y[,p,d] - min(Y[,p,d]))/(max(Y[,p,d] - min(Y[,p,d])))
  }
}
    
# 1.2. Beta distributed Y
nbetas <- vector(mode = "list", length = length(betas))
    
for (b in 1:length(betas)){
  nbetas[[b]] <- array(data = NA, dim = c(N, Time, D))
  for (d in 1:D){
    for (p in 1:Time){
      nbetas[[b]][,p,d]<- (betas[[b]][,p,d] - min(betas[[b]][,p,d]))/(max(betas[[b]][,p,d]) - min(betas[[b]][,p,d]))  
    }
  }
}
    
# 2. Auxiliary variables
  
# 2.1. Correlated with normally distributed Y
nZ <- vector(mode = "list", length = length(cor_aux))
    
for (c in 1:length(cor_aux)){
 nZ[[c]] <- array(data = NA, dim = c(N, Time, D))
 for (d in 1:D){
   for (p in 1:Time){
     nZ[[c]][,p,d] <- (Z[[c]][,p,d] - min(Z[[c]][,p,d]))/(max(Z[[c]][,p,d]) - min(Z[[c]][,p,d]))
   }
 } 
}
    
# 2.2. Correlated with beta distributed Y
    
nZ_betas <- vector(mode = "list", length = length(Z_betas))
    
for (b in 1:length(Z_betas)){
  nZ_b <- vector(mode = "list", length = length(cor_aux))
  for (c in 1:length(cor_aux)){
    nZ_b[[c]] <- array(data = NA, dim = c(N, Time, D))
    for (d in 1:D){
      for (p in 1:Time){
        nZ_b[[c]][,p,d] <- (Z_betas[[b]][[c]][,p,d] - min(Z_betas[[b]][[c]][,p,d]))/
          (max(Z_betas[[b]][[c]][,p,d]) - min(Z_betas[[b]][[c]][,p,d]))
      }
    } 
  }
  nZ_betas[[b]] <- nZ_b
  rm(nZ_b)
}
    
#### - (V) Estimations - ####
    
Threshold <- 0
    
# (i) Meng_rYS estimation
    
# 1. Normal distribution
Selh_Meng_rYS_hat <- Sel <- Meng_rYS_hat <- array(data = NA, dim = c(Time, length(cor_sel), D))
Perf <- as.data.frame(matrix(data = NA, nrow = length(cor_sel), ncol = (2*D)))
  
for (d in 1:D){
  for (c in 1:length(cor_sel)){
   Meng_rYS_hat[,c,d] <- Meng_rYS(as.matrix(nY[,,d]), as.matrix(S[[c]][,,d]))
   Perf[c,(d*2)-1] <- MAD(as.matrix(nY[,,d]), as.matrix(S[[c]][,,d]), Meng_rYS_hat[,c,d])[[3]]
   Perf[c,(d*2)] <- RMSD(as.matrix(nY[,,d]), as.matrix(S[[c]][,,d]), Meng_rYS_hat[,c,d])[[3]]
   Sel[,c,d] <- as.vector(MAD(as.matrix(nY[,,d]), as.matrix(S[[c]][,,d]), Meng_rYS_hat[,c,d])[[1]])
   Selh_Meng_rYS_hat[,c,d] <- as.vector(MAD(as.matrix(nY[,,d]), as.matrix(S[[c]][,,d]), Meng_rYS_hat[,c,d])[[2]])
  }
}
    
Perf <- as.data.frame(cbind(as.vector(paste0("S",1:4)), as.vector(rep("a. Meng_rYS",4)),
                            as.vector(rep("Y",4)),as.vector(rep(NA,4)),Perf))
    
names <- expand.grid(c("MAD_D","RMSD_D"),1:D)
colnames(Perf) <-  c(c("S","M","Y","Z"), as.character(paste0(names$Var1,names$Var2)))
    
# 2. Beta distribution
Performance <- Perf # Final data frame in which we will append the rest of the results
for (b in 1:length(betas)){
  Meng_rYS_hat <- array(data = NA, dim = c(Time, length(cor_sel), D))
  Perf <- as.data.frame(matrix(data = NA, nrow = length(cor_sel), ncol = (2*D)))
  for (d in 1:D){
   for (c in 1:length(cor_sel)){
     Meng_rYS_hat[,c,d] <- Meng_rYS(as.matrix(nbetas[[b]][,,d]), as.matrix(S_betas[[b]][[c]][,,d]))
     Perf[c,(d*2)-1] <- MAD(as.matrix(nbetas[[b]][,,d]), as.matrix(S_betas[[b]][[c]][,,d]), Meng_rYS_hat[,c,d])[[3]]
     Perf[c,(d*2)] <- RMSD(as.matrix(nbetas[[b]][,,d]), as.matrix(S_betas[[b]][[c]][,,d]), Meng_rYS_hat[,c,d])[[3]]
   }
  }
    
  Perf <- as.data.frame(cbind(as.vector(paste0("S",1:4)), as.vector(rep("a. Meng_rYS",4)),
                              as.vector(rep(paste0("beta",b),4)),as.vector(rep(NA,4)),Perf))
      
  names <- expand.grid(c("MAD_D","RMSD_D"),1:D)
  colnames(Perf) <-  c(c("S","M","Y","Z"), as.character(paste0(names$Var1,names$Var2)))
  Performance <- rbind(Performance, Perf) 
}
    
# (ii) Remaining measures
    
Measures <- c("Meng_rZS", "MUB_half", "MUB_hist", "MUB_robust", "MUB_alt", "MUB_abs", "MUB_alt_abs", "COS", "Zdiff")
    
# 1. Normal distribution
Selh_Meas_hat <- Meas_hat <- vector(mode = "list", length = length(Measures))
Perf_M <- vector(mode = "list", length = length(Measures))
    
for (cz in 1:length(cor_aux)){
  for (m in 1:length(Measures)){
    Selh_Meas_hat[[m]] <- Meas_hat[[m]] <- array(data = NA, dim = c(Time, length(cor_sel), D))
    Perf_M[[m]] <- as.data.frame(matrix(data = NA, nrow = length(cor_sel), ncol = (2*D)))
    for (d in 1:D){
     for (cs in 1:length(cor_sel)){
      Meas_hat[[m]][,cs,d] <- get(Measures[m])(as.matrix(nY[,,d]), as.matrix(S[[cs]][,,d]), nZ[[cz]][,,d])
      Perf_M[[m]][cs,(d*2)-1] <- MAD(as.matrix(nY[,,d]), as.matrix(S[[cs]][,,d]), Meas_hat[[m]][,cs,d])[[3]]
      Perf_M[[m]][cs,(d*2)] <- RMSD(as.matrix(nY[,,d]), as.matrix(S[[cs]][,,d]), Meas_hat[[m]][,cs,d])[[3]]
      Selh_Meas_hat[[m]][,cs,d] <- as.vector(MAD(as.matrix(nY[,,d]), as.matrix(S[[cs]][,,d]), Meas_hat[[m]][,cs,d])[[2]])
     }
    }
    
    Perf_M[[m]] <- as.data.frame(cbind(as.vector(paste0("S", 1:4)), as.vector(rep(Measures[m],4)), as.vector(rep("Y", 4)),
                                       as.vector(rep(paste0("Z",cz),4)), as.matrix(Perf_M[[m]])))
    names <- expand.grid(c("MAD_D", "RMSD_D"), 1:D)
    colnames(Perf_M[[m]]) <- c(c("S", "M", "Y", "Z"), as.character(paste0(names$Var1, names$Var2)))
  }
  Performance <- rbind(Performance, Perf_M[[1]], Perf_M[[2]], Perf_M[[3]], Perf_M[[4]], 
                       Perf_M[[5]], Perf_M[[6]], Perf_M[[7]], Perf_M[[8]], Perf_M[[9]])
}
    
# 2. Beta distribution
for (b in 1:length(betas)){
  Meas_hat <- vector(mode = "list", length = length(Measures))
  Perf_M <- vector(mode = "list", length = length(Measures))
  for (cz in 1:length(cor_aux)){
    for (m in 1:length(Measures)){
      Meas_hat[[m]] <- array(data = NA, dim = c(Time, length(cor_sel), D))
      Perf_M[[m]] <- as.data.frame(matrix(data = NA, nrow = length(cor_sel), ncol = (2*D)))
      for (d in 1:D){
        for (cs in 1:length(cor_sel)){
          Meas_hat[[m]][,cs,d] <- get(Measures[m])(as.matrix(nbetas[[b]][,,d]), as.matrix(S_betas[[b]][[cs]][,,d]), nZ_betas[[b]][[cz]][,,d])
          Perf_M[[m]][cs,(d*2)-1] <- MAD(as.matrix(nbetas[[b]][,,d]), as.matrix(S_betas[[b]][[cs]][,,d]), Meas_hat[[m]][,cs,d])[[3]]
          Perf_M[[m]][cs,(d*2)] <- RMSD(as.matrix(nbetas[[b]][,,d]), as.matrix(S_betas[[b]][[cs]][,,d]), Meas_hat[[m]][,cs,d])[[3]]
        }
      }
      
      Perf_M[[m]] <- as.data.frame(cbind(as.vector(paste0("S", 1:4)), as.vector(rep(Measures[m],4)), as.vector(rep(paste0("beta",b), 4)),
                                         as.vector(rep(paste0("Z",cz),4)), as.matrix(Perf_M[[m]])))
      names <- expand.grid(c("MAD_D", "RMSD_D"), 1:D)
      colnames(Perf_M[[m]]) <- c(c("S", "M", "Y", "Z"), as.character(paste0(names$Var1, names$Var2)))
    }
    Performance <- rbind(Performance, Perf_M[[1]], Perf_M[[2]], Perf_M[[3]], Perf_M[[4]], 
                         Perf_M[[5]], Perf_M[[6]], Perf_M[[7]], Perf_M[[8]], Perf_M[[9]])
  }
}
    
mPerformance <- melt(Performance, id.vars = c("S","M","Y","Z")) # This is the melted (reshaped) data set. We are using this shape since it's easier to handle.
mPerformance$M[mPerformance$M=="Meng_rZS"] <- "b. Meng_rZS"
mPerformance$M[mPerformance$M=="MUB_half"] <- "c. MUB_half"
mPerformance$M[mPerformance$M=="MUB_hist"] <- "d. MUB_hist"
mPerformance$M[mPerformance$M=="MUB_robust"] <- "e. MUB_robust"
mPerformance$M[mPerformance$M=="MUB_alt"] <- "f. MUB_alt"
mPerformance$M[mPerformance$M=="MUB_abs"] <- "g. MUB_abs"
mPerformance$M[mPerformance$M=="MUB_alt_abs"] <- "h. MUB_alt_abs"
mPerformance$M[mPerformance$M=="COS"] <- "i. COS"
mPerformance$M[mPerformance$M=="Zdiff"] <- "j. Zdiff"
  
write.table(mPerformance, file = paste0(fold_data,"/Simulated_Data.txt"))

#### - (VI) Graphs - ####

mPerformance <- read.table(paste0(fold_data,"/Simulated_Data.txt"))
Test <- matrix(data = unlist(strsplit(mPerformance$M, split = ". ")), nrow = nrow(mPerformance), ncol = 2, byrow = TRUE)
mPerformance$L <- Test[,1]

# (i) Different S, constant Z - To see how the estimators react to changes in S
  
data <- mPerformance[which(mPerformance$Z=="Z1" | is.na(mPerformance$Z)),]
data <- data[grep("Y",data$Y),]
data$value <- as.numeric(data$value)
data$S[data$S=="S1"] <- "Cov(YS) = 0.2"
data$S[data$S=="S2"] <- "Cov(YS) = 0.4"
data$S[data$S=="S3"] <- "Cov(YS) = 0.6"
data$S[data$S=="S4"] <- "Cov(YS) = 0.8"


Plots_S_MAD <- ggplot(data = data[grep("MAD",data$variable),], aes(x = `value`, y = `M`, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) + xlim(0,0.4) +  ylab("Estimator") +
  scale_fill_viridis_c(name = "MAD", begin = 0, end = 1, option = "D") +
  theme(text = element_text(size = 20)) + xlab("MAD") + facet_grid(~S)
ggsave("Plot_S_MAD.png", path = fold_graphs, width = 40, height = 15, units = "cm")

Plots_S_RMSD <- ggplot(data = data[grep("RMSD",data$variable),], aes(x = `value`, y = `M`, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) + xlim(0,0.4) +  ylab("Estimator") +
  scale_fill_viridis_c(name = "RMSD", begin = 0, end = 1, option = "D") +
  theme(text = element_text(size = 20)) + xlab("RMSD") + facet_grid(~S)
ggsave("Plot_S_RMSD.png", path = fold_graphs, width = 40, height = 15, units = "cm")


# (ii) Different Z, constant S - To see how the estimators react to changes in Z
  
data <- mPerformance[which(mPerformance$S=="S3"),]
data <- data[grep("Y",data$Y),]
data$value <- as.numeric(data$value)
data$Z[data$Z=="Z1"] <- "Cov(YZ) = 0.2"
data$Z[data$Z=="Z2"] <- "Cov(YZ) = 0.4"
data$Z[data$Z=="Z3"] <- "Cov(YZ) = 0.6"
data$Z[data$Z=="Z4"] <- "Cov(YZ) = 0.8"
data <- data[!is.na(data$Z),]

Plots_Z_MAD <- ggplot(data = data[grep("MAD",data$variable),], aes(x = `value`, y = `M`, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) + xlim(0,0.4) +  ylab("Estimator") +
  scale_fill_viridis_c(name = "MAD", begin = 0, end = 1, option = "D") +
  theme(text = element_text(size = 20)) + xlab("MAD") + facet_grid(~Z)
Plots_Z_MAD
ggsave("Plot_Z_MAD.png", path = fold_graphs, width = 40, height = 15, units = "cm")

Plots_Z_RMSD <- ggplot(data = data[grep("RMSD",data$variable),], aes(x = `value`, y = `M`, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) + xlim(0,0.4) +  ylab("Estimator") +
  scale_fill_viridis_c(name = "RMSD", begin = 0, end = 1, option = "D") +
  theme(text = element_text(size = 20)) + xlab("RMSD") + facet_grid(~Z)
ggsave("Plot_Z_RMSD.png", path = fold_graphs, width = 40, height = 15, units = "cm")

  
# (iii) Different Y, ceteris paribus - To check how the estimators react to changes in Y distribution
  
for (b in 1:length(betas)){
  data <- mPerformance[which(mPerformance$Z=="Z1" | is.na(mPerformance$Z)),]
  data <- data[grep("S3",data$S),]
  data$Alpha <-  data$Beta <- NA
  for (i in 1:4){
    data$Alpha[grep(paste0("beta",i),data$Y)] <- data$Alpha[grep(paste0("beta",i+4),data$Y)] <- data$Alpha[grep(paste0("beta",i+8),data$Y)] <- data$Alpha[grep(paste0("beta",i+12),data$Y)] <- i
    data$Beta[grep(paste0("beta",i),data$Y)] <- 1
    data$Beta[grep(paste0("beta",i+4),data$Y)] <- 2
    data$Beta[grep(paste0("beta",i+8),data$Y)] <- 3
    data$Beta[grep(paste0("beta",i+12),data$Y)] <- 4
  }
  data <- data[!is.na(data$Alpha),]
  data <- data[!is.na(data$Beta),]
  data$value <- as.numeric(data$value)
}  

Plot_B_MAD <- ggplot(data = data[grep("MAD",data$variable),], aes(x = `value`, y = `M`, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01)  + xlim(0,0.75) + ylab("Estimator") +
  scale_fill_viridis_c(name = "MAD", option = "C") + theme(text = element_text(size = 25)) +
  xlab("MAD") + facet_grid(Beta~Alpha, labeller = label_both)
ggsave("Plot_B_MAD.png", path = fold_graphs, width = 50, height = 50, units = "cm")

Plot_B_RMSD <- ggplot(data = data[grep("RMSD",data$variable),], aes(x = `value`, y = `L`, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01)  + xlim(0,0.75) + ylab("Estimator") +
  scale_fill_viridis_c(name = "RMSD", option = "C") + theme(text = element_text(size = 25)) +
  xlab("RMSD") + facet_grid(Beta~Alpha, labeller = label_both)
ggsave("Plot_B_RMSD.png", path = fold_graphs, width = 50, height = 50, units = "cm")

# (iv) Different Y, with different S

for (b in 1:length(betas)){
  data <- mPerformance[which(mPerformance$Z=="Z1" | is.na(mPerformance$Z)),]
  data$Alpha <-  data$Beta <- NA
  for (i in 1:4){
    data$Alpha[grep(paste0("beta",i),data$Y)] <- data$Alpha[grep(paste0("beta",i+4),data$Y)] <- data$Alpha[grep(paste0("beta",i+8),data$Y)] <- data$Alpha[grep(paste0("beta",i+12),data$Y)] <- i
    data$Beta[grep(paste0("beta",i),data$Y)] <- 1
    data$Beta[grep(paste0("beta",i+4),data$Y)] <- 2
    data$Beta[grep(paste0("beta",i+8),data$Y)] <- 3
    data$Beta[grep(paste0("beta",i+12),data$Y)] <- 4
  }
  data <- data[!is.na(data$Alpha),]
  data <- data[!is.na(data$Beta),]
  data$value <- as.numeric(data$value)
}
data <- data[data$Beta==4,]
data$S[data$S=="S1"] <- "Cov(YS) = 0.2"
data$S[data$S=="S2"] <- "Cov(YS) = 0.4"
data$S[data$S=="S3"] <- "Cov(YS) = 0.6"
data$S[data$S=="S4"] <- "Cov(YS) = 0.8"
data$Betas <- paste(paste0("Alpha=",data$Alpha),paste0("Beta=",data$Beta), sep = ", ") 

Plot_BS_MAD <- ggplot(data = data[grep("MAD",data$variable),], aes(x = `value`, y = `M`, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01)  + xlim(0,0.75) + ylab("Estimator") +
  scale_fill_viridis_c(name = "MAD", option = "B") + theme(text = element_text(size = 25)) +
  xlab("MAD") + facet_grid(S~Betas)
ggsave("Plot_BS_MAD.png", path = fold_graphs, width = 50, height = 50, units = "cm")

Plot_BS_RMSD <- ggplot(data = data[grep("RMSD",data$variable),], aes(x = `value`, y = `L`, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01)  + xlim(0,0.75) + ylab("Estimator") +
  scale_fill_viridis_c(name = "RMSD", option = "B") + theme(text = element_text(size = 25)) +
  xlab("RMSD") + facet_grid(S~Betas)
ggsave("Plot_BS_RMSD.png", path = fold_graphs, width = 50, height = 50, units = "cm")

# (v) Different Y, with different Z
  
for (b in 1:length(betas)){
  data <- mPerformance[which(mPerformance$S=="S3" | is.na(mPerformance$Z)),]
  data$Alpha <-  data$Beta <- NA
  for (i in 1:4){
    data$Alpha[grep(paste0("beta",i),data$Y)] <- data$Alpha[grep(paste0("beta",i+4),data$Y)] <- data$Alpha[grep(paste0("beta",i+8),data$Y)] <- data$Alpha[grep(paste0("beta",i+12),data$Y)] <- i
    data$Beta[grep(paste0("beta",i),data$Y)] <- 1
    data$Beta[grep(paste0("beta",i+4),data$Y)] <- 2
    data$Beta[grep(paste0("beta",i+8),data$Y)] <- 3
    data$Beta[grep(paste0("beta",i+12),data$Y)] <- 4
  }
  data <- data[!is.na(data$Alpha),]
  data <- data[!is.na(data$Beta),]
  data$value <- as.numeric(data$value)
}
data <- data[data$Beta==4,]
data$Z[data$Z=="Z1"] <- "Cov(YZ) = 0.2"
data$Z[data$Z=="Z2"] <- "Cov(YZ) = 0.4"
data$Z[data$Z=="Z3"] <- "Cov(YZ) = 0.6"
data$Z[data$Z=="Z4"] <- "Cov(YZ) = 0.8"
data <- data[!is.na(data$Z),]
data$Betas <- paste(paste0("Alpha=",data$Alpha),paste0("Beta=",data$Beta), sep = ", ") 

Plot_BZ_MAD <- ggplot(data = data[grep("MAD",data$variable),], aes(x = `value`, y = `M`, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01)  + xlim(0,0.75) + ylab("Estimator") +
  scale_fill_viridis_c(name = "MAD", option = "B") + theme(text = element_text(size = 25)) +
  xlab("MAD") + facet_grid(Z~Betas)
ggsave("Plot_BZ_MAD.png", path = fold_graphs, width = 50, height = 50, units = "cm")

Plot_BZ_RMSD <- ggplot(data = data[grep("RMSD",data$variable),], aes(x = `value`, y = `L`, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01)  + xlim(0,0.75) + ylab("Estimator") +
  scale_fill_viridis_c(name = "RMSD", option = "B") + theme(text = element_text(size = 25)) +
  xlab("RMSD") + facet_grid(Z~Betas)
ggsave("Plot_BZ_RMSD.png", path = fold_graphs, width = 50, height = 50, units = "cm")

# Graphs to check the selection

plot_sel <- dat_sel <- vector(mode = "list", ncol(Sel))
cols <- c("forestgreen","sienna1", "sienna4", "lightsteelblue", "lightsteelblue1", "lightsteelblue2", "lightsteelblue3", "lightsteelblue4", "lightskyblue3", "slateblue3", "yellow4" )

for (i in 1:4){
  dat_sel[[i]] <- cbind(c(1:10), abs(Sel[,i,2]), Selh_Meng_rYS_hat[,i,2], Selh_Meas_hat[[1]][,i,2], Selh_Meas_hat[[2]][,i,2], Selh_Meas_hat[[3]][,i,2], Selh_Meas_hat[[4]][,i,2],
                        Selh_Meas_hat[[5]][,i,2], Selh_Meas_hat[[6]][,i,2], Selh_Meas_hat[[7]][,i,2], Selh_Meas_hat[[8]][,i,2], Selh_Meas_hat[[9]][,i,2])  
  colnames(dat_sel[[i]]) <- c("t","True","Meng_rYS", "Meng_rZS", "MUB_half", "MUB_hist", "MUB_robust", "MUB_alt", "MUB_abs", "MUB_alt_abs", "COS", "Zdiff" )

  dat_sel[[i]] <- melt(dat_sel[[i]], id = "t")
  dat_sel[[i]] <- dat_sel[[i]][11:120,]
  plot_sel[[i]] <- ggplot(data = dat_sel[[i]], aes(x=Var1, y = value, colour = Var2)) + geom_line() + ggtitle(paste0("S",i)) + ylim(0,0.30) +
    scale_color_manual(values = cols)
}
ggarrange(plot_sel[[1]], plot_sel[[2]], plot_sel[[3]], plot_sel[[4]],
          common.legend = TRUE, legend = "right", ncol = 2, nrow = 2)


