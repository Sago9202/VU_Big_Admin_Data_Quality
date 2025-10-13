# --------------------------------------------- #
# Modeling Total Error using Multi-Source Data  #
# Written by: Santiago Gómez-Echeverry          #
# Script 1 - Simulation                         #
# Last modified: 2/10/2024                      #
# --------------------------------------------- #

#### - (I) Working space and packages - ####

rm(list = ls())
options(digits = 3)
options(scipe = 10)
packs  <- c('copula', 'poLCA', 'lavaan', 'progress', 'ggplot2')
ipacks <- packs %in% rownames(installed.packages())
if(any(ipacks == F)){install.packages(packs[!ipacks])}
lapply(packs, library, character.only = T)

#### - (II) Additional Functions - ####

# Population standard deviation
p.sd <- function(x){sqrt(sum(x-mean(x))^2)/length(x)}

#### - (III) Simulation - ####

SMUBM <-  function(X, Y, W, R, Z, data, phi){
  # Defining the variables
  data <-  dat
  X <- data[, X]
  Y_p1 <- data[,Y[[1]]]
  Y_p2 <- data[,Y[[2]]]
  Y_np <- data[,Y[[3]]]
  W_p1 <- data[,W[[1]]]
  W_p2 <- data[,W[[2]]]
  W_np <- data[,W[[3]]]
  R <- data[, R]
  Z <- data[, Z]
  Z_np <- Z[R==1]

  # True values
  r_err <- mean(X[R==1], na.rm = T) - mean(X)
  m_err <- mean(Y_np, na.rm = T) - mean(X[R==1], na.rm = T)
  t_err <- m_err + r_err

  # Measurement part
  require(lavaan)
  require(lavaanPlot)
  require(dplyr)
  require(semTools)
  
  mod <- 'Yl =~ Y_p1 + Y_p2 + Y_np 
          Wl =~ W_p1 + W_p2 + W_np
          Yl ~~ Wl'
          
  fit <- sem(model = mod, data = dat, se = "bootstrap", 
             bootstrap = 100, meanstructure = T, effect.coding = T)
  
  sum_m <- summary(fit)
  sum_m$pe
  #summary(fit)
  #modindices(mes)
  
  mean_X_hat <- sum_m$pe[sum_m$pe$lhs=="Yl" & sum_m$pe$op=="~ 1",5]
  lavPredict(fit)[,1]
  mean_Y_hat <- lavPredictY(fit, ynames = "Y_np", xnames = c("Y_p1", "Y_p2","W_p1","W_p2", "W_np")) %>% 
    mean()
  sigma2_Y_hat <- lavPredictY(fit, ynames = "Y_np", xnames = c("Y_p1", "Y_p2","W_p1","W_p2", "W_np")) %>% 
    var() %>% 
    as.numeric()
  
  sigma2_X_hat <- sum_m$pe[sum_m$pe$lhs=="Yl" & sum_m$pe$op=="Yl",5]
  #sigma2_Y_hat <- sum_m$pe[9,5]
  #lambda_hat <- sum_m$pe[3,5]
  lambda_hat <- reliability(fit)[1,1]
  
  X_hat <- mean_X_hat + lambda_hat*(sigma2_X_hat/sigma2_Y_hat)*(Y_np - mean_Y_hat)

  # X_hat<- lavPredict(mes)[,1]
  # sigma2_X <- var(X_hat, na.rm = T)
  # lambda2 <- lavInspect(mes, what = "start")$lambda^2
  # sigma2_e <- lavResiduals(mes)$cov[2,2]
  # #sigma2_Y_np <- var(Y_np, na.rm = T)
  # gamma <- (lambda2[2]*sigma2_X)/(lambda2[2]*sigma2_X + sigma2_e)
  # gamma <- reliability(fit)[2,]
  m_err_hat <- mean(Y_np, na.rm = T) - mean(X_hat, na.rm = T)

  # Selection part
  r_yz <- cor(na.omit(Y_np), na.omit(Z_np))
  g_hat <- (phi + (1-phi)*r_yz)/(phi*r_yz + (1-phi))
  g_hat_m <- (phi + lambda_hat*(1-phi)*r_yz)/(phi*r_yz + (1-phi))
  c_obs <- sd(Y_np, na.rm = T)/sd(Z)
  r_err_hat <- g_hat*c_obs*(mean(Z[R==1], na.rm = T) - mean(Z))
  rm_err_hat <- g_hat_m*c_obs*(mean(Z[R==1], na.rm = T) - mean(Z))
  t_err_hat <- m_err_hat + rm_err_hat
  #m_err; r_err; t_err
  #m_err_hat; r_err_hat; rm_err_hat; t_err_hat
  
  main_res <- c(m_err, r_err, t_err, m_err_hat, r_err_hat, rm_err_hat, t_err_hat)
  add_res <- list(fit, X_hat, as.numeric(lambda_hat))
  res <- list(main_res, add_res)
  return(res)
}

# Parameters
D <- 5        # Number of draws
N <- 1000       # Population size
f_p <- 0.25    # Sample fraction, probability sample
n_p <- N*f_p   # Probability sample
f_np <- 0.8    # Sample fraction, non-probability sample
n_np <- N*f_np # Non-probability sample

# Variables involved in the analyses
# X: Latent target variable
# Y: Target variable
# R: Sampling indicator
# W: Auxiliary variable - Measurement
# Z: Auxiliary variable - Selection

# Let is start by defining the correlation between the variables. We will assume that they come from a multivariate normal distribution
rho_XW <- c(0.3, 0.5, 0.8)
rho_XZ <- c(0.3, 0.5, 0.8)
rho_ZW <- 0
phi <- c(0.3, 0.5, 0.8)
param <- expand.grid(XW = rho_XW, XZ = rho_XZ,  d = 1:D, phi = phi)
param$id <- seq.int(nrow(param))
param <- param[,c(5,1,2,3,4)]

set.seed(42) # The Life, the Universe, and Everything
SimDat <- vector(mode = "list", length = nrow(param)) 
ResSim <- data.frame(matrix(nrow = nrow(param), ncol = 7))

pb <- progress_bar$new(total = nrow(param))

for(i in 1:nrow(param)){
  pb$tick()
  ## Data generation
  Cop <- normalCopula(param = c(param[i,2], param[i,3], 0.3), dim = 3, dispstr = "un")
  dat <- rCopula(N, Cop)
  # Now that we have created the correlated normal variables with the copulas, we can easily change the distribution as we want 
  dat[,1] <- qnorm(dat[,1], mean = 0, sd = 1)
  dat[,2] <- qnorm(dat[,2], mean = 0, sd = 1)
  dat[,3] <- qnorm(dat[,3], mean = 0, sd = 1)
  
  # Let's check the correlations
  #cor(dat)
  dat <- data.frame(unlist(dat))
  dat$id <- seq.int(nrow(dat))
  dat <- dat[,c(4,1,2,3)]
  colnames(dat) <- c('id', 'X', 'W', 'Z')
  
  #b_X <- param[i,4]
  #b_W <- param[i,5]
  
  # (i) Probability sample
  dat$pr_p <- runif(N)
  dat <- dat[order(dat$pr_p),] 
  dat$Y_p1 <- dat$W_p1 <- dat$Y_p2 <- dat$W_p2 <- dat$Z_p <- NA
  dat$Y_p1[1:n_p] <- dat$X[1:n_p] + rnorm(n_p, mean = 0, sd = 0.1) + 0.3*dat$W[1:n_p]
  dat$Y_p2[1:n_p] <- dat$X[1:n_p] + rnorm(n_p, mean = 0, sd = 0.1) + 0.3*dat$W[1:n_p]
  dat$W_p1[1:n_p] <- dat$W[1:n_p] + rnorm(n_p, mean = 0, sd = 0.1)
  dat$W_p2[1:n_p] <- dat$W[1:n_p] + rnorm(n_p, mean = 0, sd = 0.1)
  dat$Z_p[1:n_p] <- dat$Z[1:n_p]
  dat$R_p <- ifelse(!is.na(dat$Y_p1), 1, 0)
  
  # (ii) Non-probability sample
  phi <- param[i,4]
  dat$Zs <- dat$Z*(sd(dat$X)/sd(dat$Z))
  dat$pr_np <- pnorm(phi*dat$X + (1-phi)*dat$Zs) # Based on Little et al. (2020) and taking g()=invnormal
  dat$R_np <- rbinom(n = N, size = 1, prob = dat$pr_np)
  dat$Y_np <- dat$W_np <- dat$Z_np <- NA
  dat$Y_np <- ifelse(dat$R_np==1,dat$X + rnorm(1, mean = 0, sd = 0.2)+ 0.2*dat$W,NA) 
  dat$W_np <- ifelse(dat$R_np==1,dat$W + rnorm(1, mean = 0, sd = 0.2),NA)
  dat$Z_np <- ifelse(dat$R_np==1,dat$Z,NA)
  dat$R_gs <- rbinom(N,1,0.5)
  
  SimDat[[i]] <- dat
  
  ## Results estimation
  mod <- SMUBM(X = "X", Y = c("Y_p1","Y_p2","Y_np"), W = c("W_p1","W_p2","W_np"), R = "R_np", Z = "Z", data = dat, phi = param[i,5])
  ResSim[i,] <- mod[[1]]
}

colnames(ResSim) <- c("m_err", "r_err", "t_err", "m_err_hat", "r_err_hat", "rm_err_hat", "t_err_hat")
ResSim <- cbind(param, ResSim)

ggplot(ResSim, aes(x = r_err, y = r_err_hat, color = as.factor(phi))) + 
  geom_point() + facet_grid(XW ~ XZ)

