##-----------------------------------------##
## Simulation on Selection Bias Estimators ##
## Written by: Santiago Gómez-Echeverry    ##
## Last update: 30/08/2023                 ##
##-----------------------------------------##

#### - (I) Working space and packages - ####
    
rm(list = ls()) # - Remove all the objects in the environment
sys <- Sys.info()
options(digits = 3)
options(scipe = 10)
# We will assign each one of the different folders into an object so that we can access them easily
fold_code <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/1 - Selection Bias Simulations/2 - Code")
fold_data <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/1 - Selection Bias Simulations/3 - Data")
fold_graphs <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/1 - Selection Bias Simulations/4 - Graphs & Tables")
setwd(fold_code) # For the moment, let us work on the code folder
# getwd() # - To check that we are in the right folder
    
# Below are all the packages that we will use in the analyses. We will check if they are installed, install them if they
# are not, and finally load them. Note: dutchmasters is installed from github
packages <- c('ggplot2', 'ggpubr', 'forecast', 'reshape2', 'MASS', 'corpcor', 'ggridges', 'ltm', 'stringr', 'knitr', 'kableExtra', 'tidyr', 'dplyr', 'utils',
              'robustbase', 'forcats', 'cols4all', 'matrixStats')
    
# Installing the packages
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)){install.packages(packages[!installed_packages])}

# Loading the packages 
invisible(lapply(packages, library, character.only = TRUE))
    
# Besides setting the working space, I'll define in here some simple functions that are required across the code.
# 1. Population Standard Deviation
p.sd <- function(x){sqrt(sum((x - mean(x))^2)/length(x))}  

# 2. Median across columns of an array
colMed <- function(array){
  mat_array <- matrix(NA, nrow = dim(array)[2], ncol = dim(array)[3])
  for (i in 1:dim(array)[2]){
    mat_array[i,] <- colMedians(as.matrix(array[,i,]), na.rm = TRUE) 
  } 
  return(mat_array)
}

nice_palette <- c4a(palette = "brewer.paired", n = 11)

#### - (II) Selection Bias Estimators - ####

# This section creates a function for each one of the estimators that we will test on the simulated data sets. It's
# good to keep in mind that this functions will take Y, S, and in most cases Z as an input, and will return a selection
# estimate.

# (1) Meng_h
    
Meng_h <- function(Y,S){
  t_total <- ncol(Y)
  y <- yl <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
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
      yl[,t] <- NA
      r_ys_hat <- NA
      sigma_t_Y_hat <- NA
    } else if (t>1) {
      yl[,t] <- Y[,t-1]*S[,t]
      yl[,t][yl[,t]==0] <- NA 
      r_ys_hat <- cor.test(Y[,t-1],S[,t-1])$estimate
      sigma_t_Y_hat <- p.sd(Y[,t-1])*sd(y[,t], na.rm = TRUE)/sd(yl[,t], na.rm = TRUE)
    }
    D_I[t] <- r_ys_hat
    D_U[t] <- sigma_t_Y_hat
    D_O[t] <- sqrt((N_t-n_t)/n_t)
    S_Meng[t] <- D_I[t]*D_U[t]*D_O[t]
  }
  res <- list(S_Meng, D_I, D_U, r_ys)
  return(res)
}
    
# (2) Meng_Z1
    
Meng_Z1 <- function(Y,S,Z){
  t_total <- ncol(Y)
  y <- yl <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
  S_Meng <-  matrix(data = NA, nrow = t_total, ncol = 1) 
  D_I <- D_U <- D_O <- numeric(t_total) 
  for (t in 1:t_total){
    N_t <- length(Y[,t])
    n_t <- sum(S[,t])
    y[,t] <- Y[,t]*S[,t]
    y[,t][y[,t]==0] <- NA
    # Threefold decomposition
    r_ys <- cor.test(Y[,t],S[,t])$estimate
    r_ys_hat <- cor.test(Z[,t],S[,t])$estimate
    sigma_t_Y <- p.sd(Y[,t])
    if (t==1){
      yl[,t] <- NA
      sigma_t_Y_hat <- NA
    } else if (t>1){
      sigma_t_Y_hat <- p.sd(Z[,t])*p.sd(Y[,t-1])/p.sd(Z[,t-1])
    }
    D_I[t] <- r_ys_hat
    D_U[t] <- sigma_t_Y_hat
    D_O[t] <- sqrt((N_t-n_t)/n_t)
    S_Meng[t] <- D_I[t]*D_U[t]*D_O[t]
  }
  res <- list(S_Meng, D_I, D_U, r_ys)
  return(res)
}

# (3) Meng_Z2

Meng_Z2 <- function(Y,S,Z){
  t_total <- ncol(Y)
  y <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
  z <- matrix(data = NA, nrow = length(Z[,1]), ncol = t_total)
  S_Meng <- matrix(data = NA, nrow = t_total, ncol = 1) 
  D_I <- D_U <- D_O <- numeric(t_total) 
  for (t in 1:t_total){
    N_t <- length(Y[,t])
    n_t <- sum(S[,t])
    y[,t] <- Y[,t]*S[,t]
    y[,t][y[,t]==0] <- NA
    z[,t] <- Z[,t]*S[,t]
    z[,t][z[,t]==0] <- NA
    # Threefold decomposition
    r_ys <- cor.test(Y[,t],S[,t])$estimate
    sigma_t_Y <- p.sd(Y[,t])
    r_ys_hat <- cor.test(Z[,t],S[,t])$estimate
    if (t==1){
      sigma_t_Y_hat <- NA
    } else if (t>1){
      sigma_t_Y_hat <- (sd(y[,t], na.rm = TRUE)/sd(z[,t], na.rm = TRUE))*p.sd(Z[,t])
    }
    D_I[t] <- r_ys_hat
    D_U[t] <- sigma_t_Y_hat
    D_O[t] <- sqrt((N_t-n_t)/n_t)
    S_Meng[t] <- D_I[t]*D_U[t]*D_O[t]
  }
  res <- list(S_Meng, D_I, D_U, r_ys)
  return(res)
}

# (4) SMUB_f

MUB_f <- function(Y,S,Z){
  t_total <- ncol(Y)
  y <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
  z <- matrix(data = NA, nrow = length(Z[,1]), ncol = t_total)
  Delta_Y <- Delta_Z <- rho_zy <- phi <- phi_hat <- numeric(t_total)
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
      ((Delta_Y[t]-rho_zy[t]*Delta_Z[t])+(Delta_Z[t]-rho_zy[t]*Delta_Y[t]))
    phi_hat[t] <- 0.5
    # Remaining elements of Y(phi)
    g_hat_phi[t] <- (phi_hat[t]+(1-phi_hat[t])*rho_zy[t])/(phi_hat[t]*rho_zy[t]+(1-phi_hat[t]))
    c_obs[t] <- sd(y[,t], na.rm = TRUE)/sd(z[,t], na.rm = TRUE)
    # Y(phi) and MUB estimator 
    Y_bar_phi[t] <- mean(y[,t], na.rm = TRUE) - g_hat_phi[t]*c_obs[t]*(mean(z[,t], na.rm = TRUE) - mean(Z[,t]))
    S_MUB[t] <- mean(y[,t], na.rm = TRUE) - Y_bar_phi[t]
  }
  res <- list(S_MUB, phi_hat, phi)
  return(res)
}
    
# (5) SMUB_ht
    
MUB_ht <- function(Y,S,Z){
  t_total <- ncol(Y)
  y <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
  z <- matrix(data = NA, nrow = length(Z[,1]), ncol = t_total)
  Delta_Y <- Delta_Z <- rho_zy <- phi <- phi_hat <- phi_hat_trunc <- numeric(t_total)
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
      ((Delta_Y[t]-rho_zy[t]*Delta_Z[t])+(Delta_Z[t]-rho_zy[t]*Delta_Y[t]))
    if (t==1){
      phi_hat[t] <- NA
    } else if (t>1){
      phi_hat[t] <- phi[t-1]
    }
    # Truncation to [0,1]
    phi_hat[t] <- ifelse(phi_hat[t]>0,ifelse(phi_hat[t]<1,phi_hat[t],1),0)
    # Remaining elements of Y(phi)
    g_hat_phi[t] <- (phi_hat[t]+(1-phi_hat[t])*rho_zy[t])/(phi_hat[t]*rho_zy[t]+(1-phi_hat[t]))
    c_obs[t] <- sd(y[,t], na.rm = TRUE)/sd(z[,t], na.rm = TRUE)
    # Y(phi) and MUB estimator 
    Y_bar_phi[t] <- mean(y[,t], na.rm = TRUE) - g_hat_phi[t]*c_obs[t]*(mean(z[,t], na.rm = TRUE) - mean(Z[,t]))
    S_MUB[t] <- mean(y[,t], na.rm = TRUE) - Y_bar_phi[t]
  }
  res <- list(S_MUB, phi_hat, phi)
  return(res)
}

# (6) SMUB_ha

MUB_ha <- function(Y,S,Z){
  t_total <- ncol(Y)
  y <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
  z <- matrix(data = NA, nrow = length(Z[,1]), ncol = t_total)
  Delta_Y <- Delta_Z <- rho_zy <- phi <- phi_hat <- numeric(t_total)
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
      ((Delta_Y[t]-rho_zy[t]*Delta_Z[t])+(Delta_Z[t]-rho_zy[t]*Delta_Y[t]))
    if (t==1){
      phi_hat[t] <- NA
    } else if (t>1){
      phi_hat[t] <- abs(Delta_Y[t-1] - rho_zy[t-1]*Delta_Z[t-1])/
        (abs(Delta_Y[t-1]-rho_zy[t-1]*Delta_Z[t-1])+
           abs(Delta_Z[t-1]-rho_zy[t-1]*Delta_Y[t-1]))
    }
    # Remaining elements of Y(phi)
    g_hat_phi[t] <- (phi_hat[t]+(1-phi_hat[t])*rho_zy[t])/(phi_hat[t]*rho_zy[t]+(1-phi_hat[t]))
    c_obs[t] <- sd(y[,t], na.rm = TRUE)/sd(z[,t], na.rm = TRUE)
    # Y(phi) and SMUB estimator 
    Y_bar_phi[t] <- mean(y[,t], na.rm = TRUE) - g_hat_phi[t]*c_obs[t]*(mean(z[,t], na.rm = TRUE) - mean(Z[,t]))
    S_MUB[t] <- mean(y[,t], na.rm = TRUE) - Y_bar_phi[t]
  }
  res <- list(S_MUB, phi_hat, phi)
  return(res)
}

# (7) SMUB_wt

MUB_wt <- function(Y,S,Z){
  t_total <- ncol(Y)
  y <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
  z <- matrix(data = NA, nrow = length(Z[,1]), ncol = t_total)
  Delta_Y <- Delta_Z <- rho_zy <- phi <- phi_hat <- phi_hat_trunc <- numeric(t_total)
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
    rho_zy[t] <- cor(Y[,t], Z[,t], use = "complete.obs")
    phi[t] <- (Delta_Y[t] - rho_zy[t]*Delta_Z[t])/
      ((Delta_Y[t]-rho_zy[t]*Delta_Z[t])+(Delta_Z[t]-rho_zy[t]*Delta_Y[t]))
    ratio <- sd(Z[,t], na.rm = TRUE)/sd(y[,t], na.rm = TRUE)
    if (t==1){
      phi_hat[t] <- NA
    } else if (t>1){
      phi_hat[t] <- (Delta_Y[t-1] - rho_zy[t-1]*ratio*Delta_Z[t-1])/
        ((Delta_Y[t-1]-rho_zy[t-1]*ratio*Delta_Z[t-1])+
           (Delta_Z[t-1]-rho_zy[t-1]*(1/ratio)*Delta_Y[t-1]))
    }
    # Truncation to [0,1]
    phi_hat[t] <- ifelse(phi_hat[t]>0,ifelse(phi_hat[t]<1,phi_hat[t],1),0)
    # Remaining elements of Y(phi)
    g_hat_phi[t] <- (phi_hat[t]+(1-phi_hat[t])*ratio*rho_zy[t])/(phi_hat[t]*(1/ratio)*rho_zy[t]+(1-phi_hat[t]))
    c_obs[t] <- sd(y[,t], na.rm = TRUE)/sd(z[,t], na.rm = TRUE)
    # Y(phi) and MUB estimator 
    Y_bar_phi[t] <- mean(y[,t], na.rm = TRUE) - g_hat_phi[t]*c_obs[t]*(mean(z[,t], na.rm = TRUE) - mean(Z[,t]))
    S_MUB[t] <- mean(y[,t], na.rm = TRUE) - Y_bar_phi[t]
  }
  res <- list(S_MUB, phi_hat, phi)
  return(res)
}

# (8) SMUB_wa

MUB_wa <- function(Y,S,Z){
  t_total <- ncol(Y)
  y <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
  z <- matrix(data = NA, nrow = length(Z[,1]), ncol = t_total)
  Delta_Y <- Delta_Z <- rho_zy <- phi <- phi_hat <- numeric(t_total)
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
      ((Delta_Y[t]-rho_zy[t]*Delta_Z[t])+(Delta_Z[t]-rho_zy[t]*Delta_Y[t]))
    ratio <- sd(Z[,t], na.rm = TRUE)/sd(y[,t], na.rm = TRUE)
    if (t==1){
      phi_hat[t] <- NA
    } else if (t>1){
      phi_hat[t] <- abs(Delta_Y[t-1] - rho_zy[t-1]*ratio*Delta_Z[t-1])/
        (abs(Delta_Y[t-1]-rho_zy[t-1]*ratio*Delta_Z[t-1]) + abs(Delta_Z[t-1]-rho_zy[t-1]*(1/ratio)*Delta_Y[t-1]))
    }
    # Remaining elements of Y(phi)
    g_hat_phi[t] <- (phi_hat[t]+(1-phi_hat[t])*ratio*rho_zy[t])/(phi_hat[t]*(1/ratio)*rho_zy[t]+(1-phi_hat[t]))
    c_obs[t] <- sd(y[,t], na.rm = TRUE)/sd(z[,t], na.rm = TRUE)
    # Y(phi) and MUB estimator 
    Y_bar_phi[t] <- mean(y[,t], na.rm = TRUE ) - g_hat_phi[t]*c_obs[t]*(mean(z[,t], na.rm = TRUE) - mean(Z[,t]))
    S_MUB[t] <- mean(y[,t], na.rm = TRUE) - Y_bar_phi[t]
  }
  res <- list(S_MUB, phi_hat, phi)
  return(res)
}

# (9) SMUB_mt
    
MUB_mt <- function(Y,S,Z){
  t_total <- ncol(Y)
  y <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
  z <- matrix(data = NA, nrow = length(Z[,1]), ncol = t_total)
  Delta_Y <- Delta_Z <- rho_zy <- phi <- phi_hat <- phi_hat_trunc <- numeric(t_total)
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
      ((Delta_Y[t]-rho_zy[t]*Delta_Z[t])+(Delta_Z[t]-rho_zy[t]*Delta_Y[t]))
    if (t==1){
      phi_hat[t] <- NA
    } else if (t>1){
      phi_hat[t] <- phi[t-1]
    }
    # Truncation to [0,1]
    phi_hat[t] <- ifelse(phi_hat[t]>0,ifelse(phi_hat[t]<1,phi_hat[t],1),0)
    # Remaining elements of Y(phi)
    g_hat_phi[t] <- (phi_hat[t]+(1-phi_hat[t])*rho_zy[t])/(phi_hat[t]*rho_zy[t]+(1-phi_hat[t]))
    c_obs[t] <- sd(y[,t], na.rm = TRUE)/sd(z[,t], na.rm = TRUE)
    # Y(phi) and MUB estimator 
    Y_bar_phi[t] <- mean(y[,t], na.rm = TRUE) - g_hat_phi[t]*c_obs[t]*(mean(z[,t], na.rm = TRUE) - mean(Z[,t]))
    S_MUB[t] <- mean(y[,t], na.rm = TRUE) - Y_bar_phi[t]
  }
  res <- list(S_MUB, phi_hat, phi)
  return(res)
}
    
# (10) SMUB_ma

MUB_ma <- function(Y,S,Z){
  t_total <- ncol(Y)
  y <- matrix(data = NA, nrow = length(Y[,1]), ncol = t_total)
  z <- matrix(data = NA, nrow = length(Z[,1]), ncol = t_total)
  Delta_Y <- Delta_Z <- rho_zy <- phi <- phi_hat  <- numeric(t_total)
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
      ((Delta_Y[t]-rho_zy[t]*Delta_Z[t])+(Delta_Z[t]-rho_zy[t]*Delta_Y[t]))
    if (t==1){
      phi_hat[t] <- NA
    } else if (t>1){
      phi_hat[t] <- abs(Delta_Y[t-1] - rho_zy[t-1]*Delta_Z[t-1])/
        (abs(Delta_Y[t-1]-rho_zy[t-1]*Delta_Z[t-1]) + abs(Delta_Z[t-1]-rho_zy[t-1]*Delta_Y[t-1]))
    }
    # Remaining elements of Y(phi)
    g_hat_phi[t] <- (phi_hat[t]+(1-phi_hat[t])*rho_zy[t])/(phi_hat[t]*rho_zy[t]+(1-phi_hat[t]))
    c_obs[t] <- sd(y[,t], na.rm = TRUE)/sd(z[,t], na.rm = TRUE)
    # Y(phi) and MUB estimator 
    Y_bar_phi[t] <- mean(y[,t], na.rm = TRUE) - g_hat_phi[t]*c_obs[t]*(mean(z[,t], na.rm = TRUE) - mean(Z[,t]))
    S_MUB[t] <- mean(y[,t], na.rm = TRUE) - Y_bar_phi[t]
  }
  res <- list(S_MUB, phi_hat, phi)
  return(res)
}

# (11) COS
  
COS <- function(Y,S,Z){
  t_total <- ncol(Y)
  y <- z <- matrix(data = NA, nrow = length(Z[,1]), ncol = t_total)
  rs_hat <- rs <- my <- mY <-mz <- mZ  <- temp_sigma_ys <-  temp_sigma_zs <- temp_sigma_ys_hat <- temp_sigma_zs_hat <- S_COS <- S_COS_True <- matrix(data = NA, nrow = t_total, ncol = 1) 
  for (t in 1:t_total){
    z[,t] <- Z[,t]*S[,t]
    z[,t][z[,t]==0] <- NA
    y[,t] <- Y[,t]*S[,t]
    y[,t][y[,t]==0] <- NA
    if (t==1){
      S_COS[t] <- NA
    } else if (t>1){
      sigma_ys_hat<- cov(Y[,t-1],S[,t-1])
      sigma_zs <- cov(Z[,t],S[,t])
      S_COS[t] <- (sigma_ys_hat/sigma_zs)*(mean(z[,t], na.rm = TRUE) - mean(Z[,t]))
    }
  }
  res <- list(S_COS, sigma_ys_hat)
  return(res)
}
  
# (12) Z_diff
    
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
    #RS_hat[t] <- ifelse((mean(y[,t], na.rm = TRUE)-M[t])<0.005, NA, RS_hat[t])
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
    #RS_hat[t] <- ifelse((mean(y[,t], na.rm = TRUE)-M[t])<0.005, NA, RS_hat[t])
  }
  RMSD <- sqrt(mean((RS_hat[-1]-RS[-1])^2, na.rm = TRUE))
  RMSD_list <- list("RS" = RS, "RS_hat" = RS_hat, "RMSD" = RMSD)
  return(RMSD_list)
}
    
#### - (IV) Data generation - ####
    
# (i) Initial settings
    
set.seed(42) # We start by setting the seed so that we can reproduce the whole process every the code is run
  
# Now, we define the initial parameters that we will use to create the data

mu <- 10 # Mean of the normal variable
sigma <- 1 # Deviation of the normal variable
N <- 1000 # Number of observations
Time <- 10 # Number of periods
D <- 300 # Number of draws
const <- rnorm(n = 1, mean = 0, sd = 1)
epsilon_0 <- matrix(data = rnorm(n = N*D, mean = 0, sd = 1), nrow = N, ncol = D)

# (ii) Outcome variable assignation
  
# For all the variable's assignation we will employ three-dimensional arrays, where D1 = observation, D2 = Period, D3 = Draw.
psi <- numeric(Time)
Y <- epsilon <- array(data = NA, dim = c(N, Time, D))
    
# 1. Normal distribution
# First, let's assign the initial period (t = 1). To do this, we create a matrix in which we have the initial observations for
# all the possible draws.

Y_0 <- matrix(data = rnorm(n = N*D, mean = mu, sd = sigma), nrow = N, ncol = D)  

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
Perf_Mat <- data.frame(matrix(nrow = 0, ncol = length(Vars)))
colnames(Perf_Mat) <- Vars
Est <- c("Meng_h","Meng_Z1","Meng_Z2", "MUB_f", "MUB_ht","MUB_ha", "MUB_wt", "MUB_wa", "MUB_mt", "MUB_ma", "COS", "Zdiff")
n_Est <- length(Est)

# i. Normal distribution 
Perf <- vector(mode = "list", length = n_Est)                                    # List with all the estimator's performance (MAD, RMSD)
RS <-  array(data = NA, dim = c(Time, n_sel, D))                                 # Array that contains the true relative selectivity
RS_hat <- Est_hat <- vector(mode = "list", length = n_Est)                       # Lists for the relative selection estimates and selectivity

pb <- txtProgressBar(min = 0, max = n_aux, initial = 0, style = 3)               # Set the progress bar to see how the loop is advancing
for (cz in 1:n_aux){
  setTxtProgressBar(pb,cz)                                                       # Initiate the progress bar
  for (m in 1:n_Est){
    RS_hat[[m]] <- Est_hat[[m]] <- array(data = NA, dim = c(Time, n_sel, D))     # Array with the relative selection and the selection estimates per estimator
    Perf[[m]] <- as.data.frame(matrix(data = NA, nrow = n_sel, ncol = (2*D)))    # Data frame with the estimators' performance
    for (d in 1:D){
      for (cs in 1:n_sel){
        Y_temp <- as.matrix(Y[,,d])                                               # Temporary Y matrix 
        S_temp <- as.matrix(S[[cs]][,,d])                                         # Temporary S matrix
        Z_temp <- Z[[cz]][,,d]                                                    # Temporary Z matrix
        
        # Since the estimation for Meng_rYS is different (i.e., doesn't include Z) we will use a conditional
        if (m==1){
          Est_temp <- Est_hat[[m]][,cs,d] <- get(Est[m])(Y_temp, S_temp)[[1]]     # Vector with the selection estimates
          RS[,cs,d] <- as.vector(MAD(Y_temp, S_temp, Est_temp)[[1]])              # Vector with the true relative selection (only needed once!)
        } else if (m!=1 & m!=12){
          Est_temp <- Est_hat[[m]][,cs,d] <- get(Est[m])(Y_temp, S_temp, Z_temp)[[1]]  # Vector with the selection estimates
        } else if (m == 12){
          Est_temp <- Est_hat[[m]][,cs,d] <- get(Est[m])(Y_temp, S_temp, Z_temp)  # Vector with the selection estimates
        }
        Perf[[m]][cs,(d*2)-1] <- MAD(Y_temp, S_temp, Est_temp)[[3]]               # MAD constructed with Y, S and the selection estimates
        Perf[[m]][cs,(d*2)] <- RMSD(Y_temp, S_temp, Est_temp)[[3]]                # RMSD constructed with Y, S and the selection estimates
        RS_hat[[m]][,cs,d] <- as.vector(MAD(Y_temp, S_temp, Est_temp)[[2]])       # Vector with the estimated relative selection
      }
    }
    # Let us assign the results into a data frame 
    Perf[[m]] <- as.data.frame(cbind(as.vector(paste0("S", 1:n_sel)),
                                       t(matrix(data = rep(c("Y", Est[m], NA, NA, paste0("Z",cz)), n_sel), nrow = 5)), Perf[[m]]))
    colnames(Perf[[m]]) <- Vars # We need to add the proper names to the columns of the data frame
  }
  Perf_Mat <- rbind(Perf_Mat, bind_rows(Perf, .id = "column_label"))            # Bind the results from every estimator at a given Cov(YZ) with the preiviously obtained results
  assign(paste0("Est_hat_Z", cz), Est_hat)                                      # Store the selection per different cz
  assign(paste0("RS_hat_Z",cz), RS_hat)                                         # Store the relative selection per different cz
  close(pb) # Close the progress bar
}

# ii. Beta distribution
RS_by_B <- RS_hat_by_B <- Est_hat_by_B <- vector(mode = "list", length =n_bt)

pb <- txtProgressBar(min = 0, max = n_bt, initial = 0, style = 3) 
for (b in 1:n_bt){
  setTxtProgressBar(pb,b)
  Perf <- vector(mode = "list", length = n_Est)
  RS_by_B[[b]] <-  array(data = NA, dim = c(Time, n_sel, D))
  RS_hat <- Est_hat <- vector(mode = "list", length = n_Est)
  for (cz in 1:n_aux){
    for (m in 1:n_Est){
      RS_hat[[m]] <- array(data = NA, dim = c(Time, n_sel, D))                       # Array with the selection estimates
      Est_hat[[m]] <- array(data = NA, dim = c(Time, n_sel, D))
      Perf[[m]] <- as.data.frame(matrix(data = NA, nrow = n_sel, ncol = (2*D)))      # Data frame with the estimator's performance.
      for (d in 1:D){
        for (cs in 1:n_sel){
          B_temp <- as.matrix(bt[[b]][,,d])                                      # Temporary Y matrix 
          S_temp <- as.matrix(S_bt[[b]][[cs]][,,d])                              # Temporary S matrix 
          Z_temp <- Z_bt[[b]][[cz]][,,d]                                         # Temporary Z matrix 
          
          # Since the estimation for Meng_rYS is different (i.e., doesn't include Z) we will use a conditional
          if (m==1){
            Est_temp <- Est_hat[[m]][,cs,d] <- get(Est[m])(B_temp, S_temp)[[1]]              # Vector with the selection estimates
            RS_by_B[[b]][,cs,d] <- as.vector(MAD(Y_temp, S_temp, Est_temp)[[1]])             # Vector with the true relative selection
          } else if (m!=1 & m!=12){
            Est_temp <- Est_hat[[m]][,cs,d] <- get(Est[m])(B_temp, S_temp, Z_temp)[[1]]      # Vector with the selection estimates
          } else if (m == 12){
            Est_temp <- Est_hat[[m]][,cs,d] <- get(Est[m])(B_temp, S_temp, Z_temp)           # Vector with the selection estimates
          }
          Perf[[m]][cs,(d*2)-1] <- MAD(B_temp, S_temp, Est_temp)[[3]]            # MAD constructed with Y, S and the selection estimates
          Perf[[m]][cs,(d*2)] <- RMSD(B_temp, S_temp, Est_temp)[[3]]             # RMSD constructed with Y, S and the selection estimates
          RS_hat[[m]][,cs,d] <- as.vector(MAD(B_temp, S_temp, Est_temp)[[2]])      # Vector with the estimated relative selection
        }
      }
      # Let us assign the results into a data frame 
      Perf[[m]] <- as.data.frame(cbind(as.vector(paste0("S", 1:n_sel)),
                                         t(matrix(data = rep(c(paste0("beta",b), Est[m], al_bt$Var1[b], al_bt$Var2[b], paste0("Z",cz)), n_sel), nrow = 5)), Perf[[m]]))
      colnames(Perf[[m]]) <- Vars                                                # We need to add the proper names to the columns of the data frame
    }
    Perf_Mat <- rbind(Perf_Mat, bind_rows(Perf, .id = "column_label"))     # Bind the results from every estimator at a given Cov(YZ) with the previously obtained results
    assign(paste0("Est_B_hat_Z", cz), Est_hat)
    assign(paste0("RS_B_hat_Z",cz), RS_hat)    
  }
  Est_hat_by_B[[b]] <- list(Est_B_hat_Z1, Est_B_hat_Z2, Est_B_hat_Z3, Est_B_hat_Z4, Est_B_hat_Z5, Est_B_hat_Z6)
  RS_hat_by_B[[b]] <- list(RS_B_hat_Z1, RS_B_hat_Z2, RS_B_hat_Z3, RS_B_hat_Z4, RS_B_hat_Z5, RS_B_hat_Z6)
  close(pb)
  rm(Est_B_hat_Z1, Est_B_hat_Z2, Est_B_hat_Z3, Est_B_hat_Z4, Est_B_hat_Z5, Est_B_hat_Z6,
     RS_B_hat_Z1, RS_B_hat_Z2, RS_B_hat_Z3, RS_B_hat_Z4, RS_B_hat_Z5, RS_B_hat_Z6)
}

Perf_Mat <- select(Perf_Mat, -column_label)                              # Getting rid of the additional variable created when using bind_rows()
mPerf <- melt(Perf_Mat, id.vars = c("S","M","Y","Alpha","Beta","Z")) # This is the melted (reshaped) data set. We are using this shape since it's easier to handle.

# Before saving the results, let us arrange the variables so that it's clear what we are seeing when we handle the data

# Recoding the estimators
mPerf$M <- factor(mPerf$M, levels = Est)
# Recoding the correlations with the auxiliary variable
mPerf$Z <- mPerf$Z %>%
  recode("Z1" = "Cor(YZ) = 0.05", "Z2" = "Cor(YZ) = 0.2",
         "Z3" = "Cor(YZ) = 0.4", "Z4" = "Cor(YZ) = 0.6",
         "Z5" = "Cor(YZ) = 0.8", "Z6" = "Cor(YZ) = 0.95")
# Recoding the correlations with the sampling variable
mPerf$S <- mPerf$S %>%
  recode("S1" = "Cor(YS) = 0.2", "S2" = "Cor(YS) = 0.4",
         "S3" = "Cor(YS) = 0.6", "S4" = "Cor(YS) = 0.8")

write.table(mPerf, file = paste0(fold_data,"/Simulated_Data_N", N, "_D", D, "_T", Time, ".txt"))

# Trimmed data
# Since we will look at the relative selectivity, we might have a problem when the sample mean is close to the selectivity. This could happen due to mere chance, and it can
# affect any estimator. Thus, we are going to drop the draws for which we observe outliers in the aforementioned difference. We set an arbitrary threshold of 100 for this.
# It is relevant to note that this arrangement shouldn't take away more than a few observations, often around 1% of the draws.

# Normal Y variable
mPerfTr <- melt(Perf_Mat, id.vars = c("S","M","Y","Alpha","Beta","Z")) # Trimmed data
RS_tr <- RS

pb <- txtProgressBar(min = 0, max = n_bt, initial = 0, style = 3) 
for (cz in 1:n_aux){
  setTxtProgressBar(pb,cz)
  temp_RS_hat <-  get(paste0("RS_hat_Z",cz))
  for (d in 1:D){
    for (m in 1:n_Est){
      for (cs in 1:n_sel){
        Y_temp <- as.matrix(Y[,,d])
        S_temp <- as.matrix(S[[cs]][,,d])
        Z_temp <- Z[[cz]][,,d]
        if (m == 1){
          Est_temp <- get(Est[m])(Y_temp, S_temp)[[1]]  
        } else if (m!=1 & m != 12) {
          Est_temp <- get(Est[m])(Y_temp, S_temp, Z_temp)[[1]]  
        } else if (m==12) {
          Est_temp <- get(Est[m])(Y_temp, S_temp, Z_temp)
        }
        y_temp <- Y_temp*S_temp
        y_temp[y_temp==0] <- NA
        check <- abs(Est_temp - colMeans(y_temp, na.rm = TRUE))<0.003
        check_real <- colMeans(Y_temp, na.rm = TRUE)<0.003
        if (any(check_real == TRUE, na.rm = TRUE) == TRUE){
          RS_tr[,cs,d] <- as.vector(rep(NA, Time)) 
        } else{
        }
        if (any(check == TRUE, na.rm = TRUE) == TRUE){
          # Replace for NA in the Performance Matrix
          mPerfTr <- mPerfTr %>% 
            mutate(value = ifelse(Z == paste0("Z", cz) & Y == "Y" & S == paste0("S", cs) & M == Est[m] & grepl(as.character(d), variable), NA, value))
          # Replace for NA in the RS estimates
          temp_RS_hat[[m]][,cs,d] <- NA
        } else {
        }
      }
    }
  }
  assign(paste0("RS_hat_Z",cz,"tr"), temp_RS_hat)
  rm(Y_temp, S_temp, Z_temp, y_temp, check, temp_RS_hat)
  close(pb) # Close the progress bar
}

# Beta Y variables
RS_hat_by_Btr <- vector(mode = "list", length =n_bt)
pb <- txtProgressBar(min = 0, max = n_bt, initial = 0, style = 3)
for (b in 1:n_bt){
  setTxtProgressBar(pb,b)
  for (cz in 1:n_aux){
    temp_RS_hat <-  RS_hat_by_B[[b]][[cz]]
    for (d in 1:D){
      for (m in 1:n_Est){
        for (cs in 1:n_sel){
          Y_temp <- as.matrix(bt[[b]][,,d])
          S_temp <- as.matrix(S_bt[[b]][[cs]][,,d])
          Z_temp <- Z_bt[[b]][[cz]][,,d]
          if (m == 1){
            Est_temp <- get(Est[m])(Y_temp, S_temp)[[1]]  
          } else if (m!=1 & m!=12) {
            Est_temp <- get(Est[m])(Y_temp, S_temp, Z_temp)[[1]]  
          } else if (m==12){
            Est_temp <- get(Est[m])(Y_temp, S_temp, Z_temp)
          }
          y_temp <- Y_temp*S_temp
          y_temp[y_temp==0] <- NA
          check <- abs(Est_temp - colMeans(y_temp, na.rm = TRUE))<0.003
          if (any(check == TRUE, na.rm = TRUE) == TRUE){
            # Replace for NA in the Performance Matriz
            mPerfTr <- mPerfTr %>% 
              mutate(value = ifelse(Z == paste0("Z", cz) & Y == paste0("beta",b) & S == paste0("S", cs) & M == Est[m] & grepl(as.character(d), variable), NA, value))
            # Replace for NA in the RS estimates
            temp_RS_hat[[m]][,cs,d] <- NA
          } else {
          }
        }
      }
    }
    assign(paste0("RS_B_hat_Z",cz,"tr"), temp_RS_hat)
  }
  RS_hat_by_Btr[[b]] <- list(RS_B_hat_Z1tr, RS_B_hat_Z2tr, RS_B_hat_Z3tr, RS_B_hat_Z4tr, RS_B_hat_Z5tr, RS_B_hat_Z6tr)
  rm(RS_B_hat_Z1tr, RS_B_hat_Z2tr, RS_B_hat_Z3tr, RS_B_hat_Z4tr, RS_B_hat_Z5tr, RS_B_hat_Z6tr)
  rm(Y_temp, S_temp, Z_temp, y_temp, check, temp_RS_hat)
  close(pb) # Close the progress bar
}

# Recoding the correlations with the auxiliary variable
mPerfTr$Z <- mPerfTr$Z %>%
  recode("Z1" = "Cor(YZ) = 0.05", "Z2" = "Cor(YZ) = 0.2",
         "Z3" = "Cor(YZ) = 0.4", "Z4" = "Cor(YZ) = 0.6",
         "Z5" = "Cor(YZ) = 0.8", "Z6" = "Cor(YZ) = 0.95")
# Recoding the correlations with the sampling variable
mPerfTr$S <- mPerfTr$S %>%
  recode("S1" = "Cor(YS) = 0.2", "S2" = "Cor(YS) = 0.4",
         "S3" = "Cor(YS) = 0.6", "S4" = "Cor(YS) = 0.8")

write.table(mPerfTr, file = paste0(fold_data,"/Simulated_Data_N", N, "_D", D, "_T", Time,"_Trimmed", ".txt"))

#### - (VI) Graphs - ####

Prf_M <- c("MAD", "RMSD")
setwd(fold_graphs) # We will store all of the graphs and tables in this folder

# Since Meng_Z2 and SMUB_f are equivalent there is no point in plotting them both
Est2 <- c("Meng_h", "Meng_Z1", "Meng_Z2/MUB_f", "MUB_ht", "MUB_ha", "MUB_wt", "MUB_wa", "MUB_mt", "MUB_ma", "COS", "Zdiff")
Perf_Plot <- mPerfTr %>% 
  filter(M != "MUB_f") %>% 
  mutate(M = str_replace(M, "Meng_Z2", "Meng_Z2/MUB_f")) %>% 
  mutate(M = factor(M, levels = Est2))

# Recoding the estimators

for (i in 1:length(Prf_M)){
  # (i) Cross S x Z
  Plot_SZ <- Perf_Plot %>%
    filter(Z!= "Cor(YZ) = 0.05" & Z!= "Cor(YZ) = 0.95"  &  grepl("Y", Y) & grepl(Prf_M[i], variable)) %>%       # Selecting only the results that we want to plot 
    ggplot(aes(x = value, y = fct_rev(M), fill = M)) +                                                          # General aesthetic
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + scale_fill_manual(values = nice_palette) +                                         # The type of plot that we will use
    facet_grid(Z~S)  + xlim(0,1.5) +                                                                            # Make several panels with specific limits
    ylab("Estimator") + xlab(Prf_M[i])  +                                                                       # The labels of each axis
    theme(text = element_text(size = 30),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.key.size = unit(1.5, 'cm'))                                                                                  # General layout options: Used to adjust the size of the text
  ggsave(paste0("Plot_SZ_",Prf_M[i], ".png"), plot = Plot_SZ, path = fold_graphs, width = 50, height = 50, units = "cm")

  # (i) Cross S x Z - extreme values
  
  Plot_SZ <- Perf_Plot %>%
    filter(Z!= "Cor(YZ) = 0.2" & Z!= "Cor(YZ) = 0.8"  &  grepl("Y", Y) & grepl(Prf_M[i], variable)) %>%    # Selecting only the results that we want to plot 
    ggplot(aes(x = value, y = fct_rev(M), fill = M)) +                                                                   # General aesthetic
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + scale_fill_manual(values = nice_palette) +                                         # The type of plot that we will use
    facet_grid(Z~S)  + xlim(0,1.5) +                                                                                     # Make several panels with specific limits
    ylab("Estimator") + xlab(Prf_M[i])  +                                                                                # The labels of each axis
    theme(text = element_text(size = 30),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.key.size = unit(1.5, 'cm'))                                                                                  # General layout options: Used to adjust the size of the text
  ggsave(paste0("Plot_SZ_",Prf_M[i], "_ext.png"), plot = Plot_SZ, path = fold_graphs, width = 50, height = 50, units = "cm")
  
  # (iii) Different Y, ceteris paribus - To check how the estimators react to changes in Y distribution

  Plot_B <- Perf_Plot %>%
    filter((Z == "Cor(YZ) = 0.8" | is.na(Z)) & S == "Cor(YS) = 0.4" & !is.na(Alpha) & grepl(Prf_M[i], variable)) %>% # Selecting only the results that we want to plot
    ggplot(aes(x = value, y = fct_rev(M), fill = M)) +                                                                        # General aesthetic
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + scale_fill_manual(values = nice_palette) +                                                  # The type of plot that we will use
    facet_grid(Beta~Alpha, labeller = label_both) + xlim(0,1) +                                                     # Make several panels with specific limits
    ylab("Estimator") + xlab(Prf_M[i]) +                                                                             # The labels of each axis
    theme(text = element_text(size = 30),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.key.size = unit(1.5, 'cm'))                                                                              # General layout options: Used to adjust the size of the text
  ggsave(paste0("Plot_B_", Prf_M[i], ".png"), plot = Plot_B, path = fold_graphs, width = 50, height = 50, units = "cm")

  # (iv) Different Y, with different S

  Plot_BS <- Perf_Plot %>%
    filter((Z == "Cor(YZ) = 0.8" | is.na(Z)) & !is.na(Alpha) & grepl(Prf_M[i], variable) & Beta == 1) %>% # Selecting only the results that we want to plot
    ggplot(aes(x = value, y = fct_rev(M), fill = M)) +                                                             # General aesthetic
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + scale_fill_manual(values = nice_palette) +                                        # The type of plot that we will use
    facet_grid(S~Alpha, labeller = label_both) + xlim(0,5) +                                            # Make several panels with specific limits
    ylab("Estimator") + xlab(Prf_M[i]) +                                                                  # The labels of each axis
    theme(text = element_text(size = 30),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.key.size = unit(1.5, 'cm'))                                                                   # General layout options: Used to adjust the size of the text
  ggsave(paste0("Plot_BS_", Prf_M[i], ".png"), plot = Plot_BS, path = fold_graphs, width = 50, height = 50, units = "cm")

  # (v) Different Y, with different Z
  
  Plot_BZ <- Perf_Plot %>%
    filter(S == "Cor(YS) = 0.4" & !is.na(Alpha) & !is.na(Z) & grepl(Prf_M[i], variable) & Beta == 1 &
             Z != "Cor(YZ) = 0.05"  & Z!= "Cor(YZ) = 0.95") %>%                              # Selecting only the results that we want to plot
    ggplot(aes(x = value, y = fct_rev(M), fill = M)) +                                                # General aesthetic
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + scale_fill_manual(values  = nice_palette) +                           # The type of plot that we will use
    facet_grid(Z~Alpha, labeller = label_both) + xlim(0,4) +                                # Make several panels with specific limits
    ylab("Estimator") + xlab(Prf_M[i]) +                                                     # The labels of each axis
    theme(text = element_text(size = 30),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.key.size = unit(1.5, 'cm'))                                                      # General layout options: Used to adjust the size of the text
  ggsave(paste0("Plot_BZ_", Prf_M[i], ".png"), plot = Plot_BZ, path = fold_graphs, width = 50, height = 50, units = "cm")

  # Since we don't see a lot of variation given the kurtosis and the skewness, let's check the plots in two particular conditions: beta1 and beta13
  
  # (vi) Cross S x Z - beta 1
  Plot_B1_SZ <- Perf_Plot %>%
    filter(Z!= "Cor(YZ) = 0.05" & Z!= "Cor(YZ) = 0.95"  &  grepl("beta1", Y) & grepl(Prf_M[i], variable)) %>%       # Selecting only the results that we want to plot 
    ggplot(aes(x = value, y = fct_rev(M), fill = M)) +                                                          # General aesthetic
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + scale_fill_manual(values = nice_palette) +                                         # The type of plot that we will use
    facet_grid(Z~S)  + xlim(0,5) +                                                                            # Make several panels with specific limits
    ylab("Estimator") + xlab(Prf_M[i])  +                                                                       # The labels of each axis
    theme(text = element_text(size = 30),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.key.size = unit(2, 'cm'))                                                                                  # General layout options: Used to adjust the size of the text
  Plot_B1_SZ
  ggsave(paste0("Plot_B1_SZ_",Prf_M[i], ".png"), plot = Plot_B1_SZ, path = fold_graphs, width = 50, height = 50, units = "cm")
  
  # (vii) Cross S x Z - beta 13
  Plot_B13_SZ <- Perf_Plot %>%
    filter(Z!= "Cor(YZ) = 0.05" & Z!= "Cor(YZ) = 0.95"  &  grepl("beta16", Y) & grepl(Prf_M[i], variable)) %>%       # Selecting only the results that we want to plot 
    ggplot(aes(x = value, y = fct_rev(M), fill = M)) +                                                          # General aesthetic
    geom_boxplot(outlier.shape = NA, alpha = 0.5) + scale_fill_manual(values = nice_palette) +                                         # The type of plot that we will use
    facet_grid(Z~S)  + xlim(0,5) +                                                                            # Make several panels with specific limits
    ylab("Estimator") + xlab(Prf_M[i])  +                                                                       # The labels of each axis
    theme(text = element_text(size = 30),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.key.size = unit(1.5, 'cm'))                                                                                  # General layout options: Used to adjust the size of the text
  Plot_B13_SZ
  ggsave(paste0("Plot_B13_SZ_",Prf_M[i], ".png"), plot = Plot_B13_SZ, path = fold_graphs, width = 50, height = 50, units = "cm")
}

# How are the sample mean and the median of these estimates?

Tab_Av <- mPerf %>%
  filter(grepl("Y", Y) & grepl("MAD", variable) & S == "Cor(YS) = 0.4" & !is.na(Z)) %>%
  group_by(M, Z) %>%
  summarize(MeanValue=mean(value)) %>%
  spread(Z, MeanValue) %>%
  kable(format = "latex")
save_kable(Tab_Av, file = "Table_Means_Z.pdf")

Tab_Med <- mPerf %>%
  filter(grepl("Y", Y) & grepl("MAD", variable) & S == "Cor(YS) = 0.4" & !is.na(Z)) %>%
  group_by(M, Z) %>%
  summarize(MedianValue=median(value)) %>%
  spread(Z, MedianValue) %>%
  kable(format = "latex")
save_kable(Tab_Med, file = "Table_Medians_Z.pdf")

Tab_AvTr <- mPerfTr %>%
  filter(grepl("Y", Y) & grepl("MAD", variable) & S == "Cor(YS) = 0.4" & !is.na(Z)) %>%
  group_by(M, Z) %>%
  summarize(MeanValue=mean(value)) %>%
  spread(Z, MeanValue) %>%
  kable(format = "latex")
save_kable(Tab_AvTr, file = "Table_Means_Z_tr.pdf")

# Graphs to check the selection

# Mean
# True Selection
AvRS <- RS[,,] %>% 
  abs() %>%            # Take the absolute value
  colMeans(na.rm = T) %>%       # Average across period
  rowMeans(na.rm = T) %>%       # Average across draw
  cbind(paste0("Cor(YS) = ", seq(from = 0.2, to = 0.8, by = 0.2))) %>% 
  as.data.frame() %>% 
  rename("RS" = 1 , "S" = 2)

# Estimated Selection
TempAvRS_Z <- vector(mode = "list", length = n_Est)
AvRS_hat <- data.frame(matrix(nrow = 0, ncol = 4))
for (j in 1:n_aux){
  for (i in 1:n_Est){
    TempAvRS_Z[[i]] <- get(paste0("RS_hat_Z",j))[[i]][,,] %>% 
      abs() %>%                        # Take the absolute value
      colMeans(na.rm = TRUE) %>%       # Average across period
      rowMeans(na.rm = TRUE) %>%       # Average across draw 
      cbind(paste0("Cor(YS) = ", seq(from = 0.2, to = 0.8, by = 0.2))) %>% 
      as.data.frame() %>% 
      rename("RS" = 1 , "S" = 2) %>% 
      mutate(M = Est[i], Z = paste0("Cor(YZ) = ", cor_aux[j]))
  }
  AvRS_hat <- rbind(AvRS_hat, bind_rows(TempAvRS_Z, .id = "column_label"))
}


AvRS_hat <- AvRS_hat %>% 
  filter(M != "MUB_f") %>% 
  mutate(M = str_replace(M, "Meng_Z2", "Meng_Z2/MUB_f")) %>% 
  mutate(M = factor(M, levels = Est2))

AvRS_hat %>% 
  filter(Z != "Cor(YZ) = 0.05"  & Z!= "Cor(YZ) = 0.95") %>%
  ggplot(aes(x = fct_rev(M), y = as.numeric(RS), fill = M)) + 
  geom_col(alpha = 0.5) +
  facet_grid(Z~S) +
  coord_flip() +
  geom_hline(data =AvRS, aes(yintercept=abs(as.numeric(RS))), col = "blue", linetype = "dashed") +
  theme(text = element_text(size = 30), strip.background = element_blank(), strip.placement = "outside", legend.key.size = unit(1.5, 'cm'),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  ylab("Relative Selectivity") + xlab("Estimator")  + theme(text = element_text(size = 15)) +
  scale_fill_manual(values=nice_palette) 
ggsave("Plot_AvSel.png", path = fold_graphs, width = 50, height = 50, units = "cm")

# Estimated Selection
# True Selection
AvRS_tr <- RS_tr[,,] %>% 
  abs() %>%               # Take the absolute value
  colMeans(na.rm = T) %>% # Average across draw
  rowMeans(na.rm = T) %>% # Average across period
  cbind(paste0("Cor(YS) = ", seq(from = 0.2, to = 0.8, by = 0.2))) %>% 
  as.data.frame() %>% 
  rename("RS" = 1 , "S" = 2)

TempAvRS_Z_tr <- vector(mode = "list", length = n_Est)
AvRS_hat_tr <- data.frame(matrix(nrow = 0, ncol = 4))
for (j in 1:n_aux){
  for (i in 1:n_Est){
    TempAvRS_Z_tr[[i]] <- get(paste0("RS_hat_Z",j,"tr"))[[i]][,,] %>% 
      abs() %>%                        # Take the absolute value
      colMeans(na.rm = TRUE) %>%       # Average across period
      rowMeans(na.rm = TRUE) %>%       # Average across draw 
      cbind(paste0("Cor(YS) = ", seq(from = 0.2, to = 0.8, by = 0.2))) %>% 
      as.data.frame() %>% 
      rename("RS" = 1 , "S" = 2) %>% 
      mutate(M = Est[i], Z = paste0("Cor(YZ) = ", cor_aux[j]))
  }
  AvRS_hat_tr <- rbind(AvRS_hat_tr, bind_rows(TempAvRS_Z_tr, .id = "column_label"))
}
AvRS_hat_tr <- AvRS_hat_tr %>% 
  filter(M != "MUB_f") %>% 
  mutate(M = str_replace(M, "Meng_Z2", "Meng_Z2/MUB_f")) %>% 
  mutate(M = factor(M, levels = Est2))

TempSeRS_Z_tr <- vector(mode = "list", length = n_Est)
SeRS_hat_tr <- data.frame(matrix(nrow = 0, ncol = 4))

# Check the number of observations

for (j in 1:n_aux){
 for (i in 1:n_Est){
   TempSeRS_Z_tr[[i]] <- get(paste0("RS_hat_Z",j,"tr"))[[i]][,,] %>% 
      abs() %>%                        # Take the absolute value
      colMeans(na.rm = TRUE) %>%       # Average across period
      rowSds(na.rm = TRUE) %>%       # Average across draw 
      cbind(paste0("Cor(YS) = ", seq(from = 0.2, to = 0.8, by = 0.2))) %>% 
      as.data.frame() %>% 
      rename("RS" = 1 , "S" = 2) %>% 
      mutate(M = Est[i], Z = paste0("Cor(YZ) = ", cor_aux[j]))
  }
  SeRS_hat_tr <- rbind(SeRS_hat_tr, bind_rows(TempSeRS_Z_tr, .id = "column_label"))
}
SeRS_hat_tr <- SeRS_hat_tr %>% 
  filter(M != "MUB_f") %>% 
  mutate(M = str_replace(M, "Meng_Z2", "Meng_Z2/MUB_f")) %>% 
  mutate(M = factor(M, levels = Est2)) %>% 
  mutate(RS = as.numeric(RS)/D)

AvRS_hat_tr <- merge(x = AvRS_hat_tr, y = SeRS_hat_tr, by = c("S", "M", "Z"))

AvRS_hat_tr %>% 
  filter(Z != "Cor(YZ) = 0.05"  & Z!= "Cor(YZ) = 0.95") %>%
  ggplot(aes(x = fct_rev(M), y = as.numeric(RS.x), fill = M)) + 
  geom_col(alpha = 0.5) +
  geom_errorbar(aes(ymin =as.numeric(RS.x)-RS.y, ymax = as.numeric(RS.x)+RS.y, width = 0.5)) +
  facet_grid(Z~S) +
  coord_flip() +
  geom_hline(data = AvRS_tr, aes(yintercept=abs(as.numeric(RS))), col = "blue", linetype = "dashed") +
  theme(text = element_text(size = 30), strip.background = element_blank(), strip.placement = "outside", legend.key.size = unit(1.5, 'cm'),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  ylab("Relative Selectivity") + xlab("Estimator")  + theme() +
  scale_fill_manual(values = nice_palette)
ggsave("Plot_AvSel_Tr.png", path = fold_graphs, width = 50, height = 50, units = "cm")

# True Selection
MedRS <- RS[,,] %>% 
  abs() %>%              # Take the absolute value
  colMed() %>%           # Average across period
  rowMedians() %>%       # Average across draw
  cbind(paste0("Cor(YS) = ", seq(from = 0.2, to = 0.8, by = 0.2))) %>% 
  as.data.frame() %>% 
  rename("RS" = 1 , "S" = 2)

# Estimated Selection
TempMedRS_Z <- vector(mode = "list", length = n_Est)
MedRS_hat <- data.frame(matrix(nrow = 0, ncol = 4))
for (j in 1:n_aux){
  for (i in 1:n_Est){
    TempMedRS_Z[[i]] <- get(paste0("RS_hat_Z",j))[[i]][,,] %>% 
      abs() %>%                        # Take the absolute value
      colMed() %>%       # Average across period
      rowMeans() %>%       # Average across draw 
      cbind(paste0("Cor(YS) = ", seq(from = 0.2, to = 0.8, by = 0.2))) %>% 
      as.data.frame() %>% 
      rename("RS" = 1 , "S" = 2) %>% 
      mutate(M = Est[i], Z = paste0("Cor(YZ) = ", cor_aux[j]))
  }
  MedRS_hat <- rbind(MedRS_hat, bind_rows(TempMedRS_Z, .id = "column_label"))
}
MedRS_hat <- MedRS_hat %>% 
  filter(M != "MUB_f") %>% 
  mutate(M = str_replace(M, "Meng_Z2", "Meng_Z2/MUB_f")) %>% 
  mutate(M = factor(M, levels = Est2))

MedRS_hat %>% 
  filter(Z != "Cor(YZ) = 0.05"  & Z!= "Cor(YZ) = 0.95") %>%
  ggplot(aes(x = fct_rev(M), y = as.numeric(RS), fill = M)) + 
  geom_col(alpha = 0.5) +
  facet_grid(Z~S) +
  coord_flip() +
  geom_hline(data = MedRS, aes(yintercept=abs(as.numeric(RS))), col = "blue", linetype = "dashed") +
  theme(text = element_text(size = 30), strip.background = element_blank(), strip.placement = "outside", legend.key.size = unit(1.5, 'cm'),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  ylab("Relative Selectivity") + xlab("Estimator")  +
  scale_fill_manual(values = nice_palette)
ggsave("Plot_MedSel.png", path = fold_graphs, width = 50, height = 50, units = "cm")

# (RS_hat - RS)  vs. MAD/RMSD

AvS <- merge(AvRS, AvRS_hat_tr, by = "S") %>% 
  mutate(dif_Sel = abs(as.numeric(RS.y) - as.numeric(RS.x))) 

for (i in 1:length(Prf_M)){
 Plot_SP <- Perf_Plot %>% 
  filter(Y=="Y" & is.na(Alpha) & is.na(Alpha) & grepl(Prf_M[i], variable) & Z != "Cor(YZ) = 0.05"  & Z!= "Cor(YZ) = 0.95") %>% 
  group_by(S, Z, M) %>% 
  summarize(mean_perf = mean(value, na.rm = TRUE)) %>% 
  merge(AvS, by = c("S", "Z", "M")) %>% 
  ggplot(aes(x = dif_Sel, y = mean_perf, color = M)) + 
  geom_point(size = 3, alpha = 0.5) + ylab(Prf_M[i]) + xlab("abs(mean(RS_hat) - mean(RS_true))") +
  facet_grid(Z~S) + scale_color_manual(values = nice_palette)
ggsave(paste0("Plot_RS_", Prf_M[i], ".png"), plot = Plot_SP, path = fold_graphs, width = 30, height = 30, units = "cm")
}

# MAD vs. RMSD

MAD_d <- Perf_Plot %>% 
  filter(Y=="Y" & is.na(Alpha) & is.na(Alpha) & grepl("MAD", variable) & Z != "Cor(YZ) = 0.05"  & Z!= "Cor(YZ) = 0.95") %>% 
  mutate(draw = substr(variable, 5, 8)) 

RMSD_d <- Perf_Plot %>% 
  filter(Y=="Y" & is.na(Alpha) & is.na(Alpha) & grepl("RMSD", variable) & Z != "Cor(YZ) = 0.05"  & Z!= "Cor(YZ) = 0.95") %>% 
  mutate(draw = substr(variable, 6, 9))

Perf_d <- merge(MAD_d[,c("S", "Z", "M", "draw", "value")], RMSD_d[,c("S", "Z", "M", "draw", "value")], by = c("S", "Z", "M", "draw"), suffixes = c("_MAD", "_RMSD")) %>% 
  group_by(S, Z, M) %>% 
  summarize(MAD = mean(value_MAD, na.rm = TRUE), RMSD = mean(value_RMSD, na.rm = TRUE)) %>% 
  ggplot(aes(x = MAD, y = RMSD, color = M)) + 
  geom_point(size = 3, alpha = 0.5) + ylab("RMSD") + xlab("MAD") + scale_color_manual(values = nice_palette) +
  facet_grid(Z~S) + 
  theme(text = element_text(size = 20),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.key.size = unit(2, 'cm'))    
ggsave(paste0("MAD_RMSD.png"), plot = Perf_d, path = fold_graphs, width = 40, height = 40, units = "cm")

# Plot phi

Est_SMUB <- c("MUB_ht","MUB_ha", "MUB_wt", "MUB_wa", "MUB_mt", "MUB_ma")
n_Est_SMUB <- length(Est_SMUB)
Cor_YS_names <- paste0("Cor(YS) = ", c(0.2, 0.4, 0.6, 0.8))
Cor_YZ_names <- paste0("Cor(YZ) = ", cor_aux)

Phi_hat <- vector(mode = "list", length = D)
Phi <- vector(mode = "list", length = D)

for (d in 1:D){
  # Estimate
  phi_hat_mat <- mPerf %>% 
    filter(grepl("MUB", M) &  !grepl("MUB_f", M) &  grepl("Y", Y)) %>%
    select(S, M, Z) %>% 
    distinct()
  phi_hat_mat <- cbind(phi_hat_mat, matrix(data = NA, nrow = nrow(phi_hat_mat), ncol = Time))
  colnames(phi_hat_mat) <- c("S", "M", "Z", paste0("phi",1:10))
  # Real
  phi_mat <- phi_hat_mat %>% 
    select(S,Z) %>% 
    distinct()
  phi_mat <- cbind(phi_mat, matrix(data = NA, nrow = nrow(phi_mat), ncol = Time))
  colnames(phi_mat) <- c("S", "Z", paste0("phi",1:10))
  
  for (cs in 1:n_sel){
    for (cz in 1:n_aux){
      for (m in 1:n_Est_SMUB){
        Y_temp <- as.matrix(Y[,,d])                       # Y matrix                                               
        S_temp <- as.matrix(S[[cs]][,,d])                 # S matrix                                   
        Z_temp <- Z[[cz]][,,d]                            # Z matrix
        
        # Estimate
        phi_hat_mat[(phi_hat_mat$S==Cor_YS_names[cs] & 
                       phi_hat_mat$Z==Cor_YZ_names[cz] & 
                       phi_hat_mat$M==Est_SMUB[m]),4:13] <- get(Est_SMUB[m])(Y_temp, S_temp, Z_temp)[[2]] 
      }
      phi_mat[(phi_mat$S==Cor_YS_names[cs] & 
                     phi_mat$Z==Cor_YZ_names[cz]),3:12] <- get(Est_SMUB[1])(Y_temp, S_temp, Z_temp)[[3]] 
    }
  }
  Phi_hat[[d]] <- phi_hat_mat
  Phi[[d]] <- phi_mat
  rm(phi_hat_mat, phi_mat)
}

Phi_av <- bind_rows(Phi, .id ="D") %>% 
  group_by(S, Z) %>% 
  summarise(across(phi1:phi10, mean)) %>% 
  mutate(M = "Obs.")

Phi_av <- bind_rows(Phi_hat, .id ="D") %>% 
  group_by(S, M, Z) %>% 
  summarise(across(phi1:phi10, mean)) %>% 
  rbind(Phi_av)

Plot_phi <- Phi_av %>% 
  pivot_longer(cols = starts_with("phi"), names_to = "period", values_to = "phi") %>% 
  mutate(period = str_replace(period, "phi", "")) %>% 
  ggplot(aes(x = as.numeric(period), y = phi, group = M)) +
  geom_line(aes(col=M)) + xlab("Period") + xlim(2,10) + ylim (0,1.5) +
  facet_grid(Z~S) +
  theme(text = element_text(size = 10),
        legend.key.size = unit(1.5, 'cm'))     
Plot_phi
ggsave("Phi_hat_by_t.png", plot = Plot_phi, path = fold_graphs, width = 40, height = 40, units = "cm")

# Plot r_YS (?), sigma_Y and sigma_Z over time

#### - (VII*) Checks - ####
# This last section of the code is to be able to make some rapid checks given certain parameters. It's important to note that the 
# object Est has the name of the different estimators.

View(mPerformance[order(-mPerformance$value),])
View(mPerformanceTr[order(-mPerformanceTr$value),])
View(mPerformanceTr[is.na(mPerformanceTr$value),])

d <- 25  # Draw
m <- 10   # Estimators (1 = Meng_rYS, 2 = Meng_rZS1, 3 = Meng_rZS2, 4 = SMUB_half, 5 = SMUB_hist, 6 = SMUB_robust, 7 = SMUB_alt, 8 = SMUB_abs, 9 = SMUB_alt_abs, 10 = S_COS, 11 = S_zdiff)
cs <- 4  # Cor(YS) (1 = 0.2, 2 = 0.4, 3 = 0.6, 4 = 0.8)
cz <- 3  # Cor(YZ) (1 = 0.05, 2 = 0.2, 3 = 0.4, 4 = 0.6, 5 = 0.8, 6 = 0.95)
Y_temp <- as.matrix(Y[,,d])                       # Y matrix                                               
S_temp <- as.matrix(S[[cs]][,,d])                 # S matrix                                   
Z_temp <- Z[[cz]][,,d]                            # Z matrix
Est_temp <- get(Est[m])(Y_temp, S_temp, Z_temp)   # Estimates
RMSD(Y_temp, S_temp, Est_temp)                     # MAD and relative selectivity (True & Estimated)

# RS_hat is pretty big, why?
time <- 2

y_temp <- Y_temp*S_temp
y_temp[y_temp==0] <-NA
y_mean_temp <- mean(y_temp[,time], na.rm = TRUE) 
y_mean_temp

Est_temp[time]; y_mean_temp
Est_temp[time]/(y_mean_temp-Est_temp[time])
