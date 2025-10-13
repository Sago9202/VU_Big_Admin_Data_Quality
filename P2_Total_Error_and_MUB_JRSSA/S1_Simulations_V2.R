## ------------------------------------------- ##
## Modeling Total Error using Multisource Data ##
## Script 1 - Simulation                       ##
## Written by: Santiago Gómez-Echeverry        ##
## Last modified: 21/03/2025                   ##
## ------------------------------------------- ##

#### - (I) Working space and packages - ####
rm(list = ls())
options(digits = 3)
options(scipe = 10)
packs <- c('copula', 'poLCA', 'lavaan', 'progress', 'ggplot2', 'lavaanPlot', 'dplyr', 'semTools', 'ggh4x', 'tidyr')
ipacks <- packs %in% rownames(installed.packages())
if(any(ipacks == F)){install.packages(packs[!ipacks])}
lapply(packs, library, character.only = T)
sys <- Sys.info()
fold_graphs <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/2 - Multisource Total Error/4 - Graphs & Tables")

compute_phi <- function(rho_XS, rho_XZ, rho_ZS) {
  num <- rho_XS - rho_XZ * rho_ZS
  denom <- (rho_XS - rho_XZ * rho_ZS) + (rho_ZS - rho_XZ * rho_XS)
  return(num / denom)
}

#### - (II) Data generating process - ####

set.seed(42)
N <- 5e3                                       # Population size
D <- 50                                        # Number of draws
f_p <- 0.30                                    # Sample fraction, probability sample
f_a <- 0.05                                    # Sample fraction, Audit sample
n_p <- N*f_p                                   # Probability sample size
n_a <- N*f_a                                   # Audit sample size
rho_xs <- rho_xz <- c(0.3, 0.5, 0.7)  # Possible correlations
phi <- c(0.4, 0.6, 0.8)                        # Possible phi's
rnd_e <- 0.3                                   # Random error
sys_e <- c(0.4, 0.6, 0.8)
b <- 0.1
param <- expand.grid("rho_xz" = rho_xz, "rho_xs" = rho_xs, "phi" = phi, D = 1:D, "sys_e"=sys_e)
# We need to define the correlation between zs based on the other parameters that are already set.
param$rho_zs <- (param$rho_xs*(1-param$phi+param$phi*param$rho_xz))/(param$phi - param$phi*param$rho_xz + param$rho_xz)
nsim <- nrow(param)

Dat <- vector(mode = "list", length = nrow(param))
pb <- progress_bar$new(total = nrow(param))
for (i in 1:nsim){
  pb$tick()
  Cop <- normalCopula(param = c(param[i,"rho_xz"], param[i,"rho_xs"], param[i,"rho_zs"]), dim = 3, dispstr = "un")
  dat <- rCopula(N, Cop)
  dat[,1] <- qnorm(dat[,1], mean = 0, sd = 1)
  dat[,2] <- qnorm(dat[,2], mean = 0, sd = 1)
  dat[,3] <- qbinom(dat[,3], size = 1, prob = 0.8)
  dat <- data.frame(unlist(dat))
  dat$id <- seq.int(nrow(dat))
  dat <- dat[,c(4,1,2,3)]
  colnames(dat) <- c('id', 'X', 'Z', 'S_np')
  dat$W <- rnorm(N, mean = 0, sd = 1)
  dat$pr_p1 <- runif(N)
  dat <- dat[order(dat$pr_p1),]
  dat$S_p1 <- 0
  dat$S_p1[1:n_p] <- 1
  dat$pr_p2 <- runif(N)
  dat <- dat[order(dat$pr_p2),]
  dat$S_p2 <- 0
  dat$S_p2[1:n_p] <- 1
  dat$pr_a <- NA
  dat$pr_a[dat$S_p1==1] <- runif(n_p)
  dat <- dat[order(dat$pr_a),]
  dat$S_a <- 0
  dat$S_a[1:n_a] <- 1
  dat <- dat[order(dat$id),]
  dat$Y_np <- ifelse(dat$S_np==1, b + param[i, "sys_e"]*dat$X + rnorm(N, mean = 0, sd = rnd_e), NA)
  dat$Y_p1 <- ifelse(dat$S_p1==1, b + 0.8*dat$X + rnorm(N, mean = 0, sd = rnd_e), NA)
  dat$Y_p2 <- ifelse(dat$S_p2==1, b + 0.6*dat$X + rnorm(N, mean = 0, sd = rnd_e), NA)
  dat$Y_a <- ifelse(dat$S_a==1, dat$X, NA)
  dat$W_np <- ifelse(dat$S_np==1, dat$W + rnorm(N, mean = 0, sd = rnd_e), NA)
  dat$W_p1 <- ifelse(dat$S_p1==1, dat$W + rnorm(N, mean = 0, sd = rnd_e), NA)
  dat$W_p2 <- ifelse(dat$S_p2==1, dat$W + rnorm(N, mean = 0, sd = rnd_e), NA)
  dat$W_a <- ifelse(dat$S_a==1, dat$W, NA)
  Dat[[i]] <- dat
}

#### - (III) SMUBM function - ####

# True values

smubm <- function(X, Y, W, S, Z, data, m_method = c("effects_coding", "multigroup_1", "multigroup_2"), M, r_method = c("monte_carlo", "rho_xs","calibration")){
  m_method <- match.arg(m_method)
  r_method <- match.arg(r_method)
  X <- data[, X]
  Y_np <- data[,Y[[3]]]
  S <- data[, S]
  Z <- data[, Z]
  # True values
  m_err <- abs(mean(Y_np, na.rm = T) - mean(X[S==1], na.rm = T))
  r_err <- abs(mean(X[S==1], na.rm = T) - mean(X))
  t_err <- m_err + r_err
  sX <- sd(X)
  sXn <- sd(X[S==1])
  
  if (m_method == "effects_coding"){
    # (1) Measurement Model - Effects coding
    mod <-  'Yl =~ lmbd1*Y_p1 + lmbd2*Y_p2 + lmbd3*Y_np
             Wl =~ dlt1*W_p1 + dlt2*W_p2 + dlt3*W_np
             Yl ~~ Wl
             
             # Effects coding constraint
            lmbd1 + lmbd2 + lmbd3 == 3
            '
    # Marker variable approach
    # mod <-  'Yl =~ 1*Y_a + lmbd2*Y_p1 + lmbd3*Y_np
    #          Wl =~ 1*W_a + dlt2*W_p1 + dlt3*W_np
    #          Yl ~~ Wl
    #         
    #          Y_a ~~ 0*Y_a
    #          W_a ~~ 0*W_a
    #         '
    # fit <- sem(model = mod, data = data, meanstructure = T)
    fit <- sem(model = mod, data = data, effect.coding = T, meanstructure = T)
    conv <- lavInspect(fit, "converged")
    hey <- (any(lavInspect(fit, "theta")<0) | any(diag(lavInspect(fit, "cov.lv")) < 0 ))
    sum_fit <- summary(fit) 
    sum_fit$pe
    
    if (any(diag(lavInspect(fit, "est")$theta)<0)){
      message("Issue at ", i)
    }
    # We need to get an estimate  of the Yl
    mx_hat <- mean(lavPredict(fit)[,"Yl"])
    s2x_hat <- var(lavPredict(fit)[,"Yl"])
    my_hat <- lavPredictY(fit, ynames = "Y_np", xnames = c("Y_p1", "Y_p2","W_p1","W_p2", "W_np")) %>% 
      mean()
    s2y_hat <- lavPredictY(fit, ynames = "Y_np", xnames = c("Y_p1", "Y_p2","W_p1","W_p2", "W_np")) %>% 
      var() %>% 
      as.numeric()
    rel <- reliability(fit)
    #rel <- semTools::compRelSEM(fit, tau.eq=F, obs.var=T)
    gamma_hat <- rel[rownames(rel)=="omega", colnames(rel)=="Yl"]
    
  } else if (m_method == "multigroup_1"){
    imp <- function(data){
      sub_data <- data %>% 
        dplyr:::select(Y_p1, Y_np, W_p1, W_np)
      data$group <- !(rowSums(is.na(sub_data))==0 & is.na(data$Y_a))
      # Imputation process - Following Scholtus, et al. (2015)
      # Step 1: Generate initial random values from N(0,1) for missing data
      data$Y_a[data$group==F] <- rnorm(sum(data$group==F))
      data$W_a[data$group==F] <- rnorm(sum(data$group==F))
      
      # Step 2: Regress each missing variable on all other variables (group 2 only) and keep residuals 
      
      reg_y <- lm(Y_a ~ Y_p1 + Y_np + W_p1 + W_np, data = data[data$group==F,])
      data$Y_a[data$group == F] <- residuals(reg_y)
      
      reg_w <- lm(W_a ~ Y_p1 + Y_np + W_p1 + W_np, data = data, subset = (group==F))
      data$W_a[data$group == F] <- residuals(reg_w)
      
      reg_yw <- lm(Y_a ~ W_a, data = data, subset = (group == F))
      data$W_a[data$group == F] <- residuals(reg_yw)

      # Step 3: Standardize imputed variables (ensure mean 0, variance 1 in group 2)
      for(var in c("Y_a", "W_a")){
        mean_val <- mean(data[[var]][data$group == F])
        sd_val  <- sd(data[[var]][data$group == F])
        data[[var]][data$group == F] <- (data[[var]][data$group == F] - mean_val)/sd_val
      }
      
      mean(data$Y_a[data$group==F]); mean(data$W_a[data$group==F])
      var(data$Y_a[data$group==F]); var(data$Y_a[data$group==F])
      data[data$group==F,] %>% 
        select(Y_a, W_a, Y_p1, Y_np, W_p1, W_np) %>% 
        cov(use = "complete.obs") %>% 
        round(3)
      
      return(data)
    }
    
    
    # SEM model
    
    
    mod <- '
  # ----- Structural Part -----
  # Allow the two latent variables to covary
  Yl ~~ Wl
  
  # ----- Measurement Model for Latent Construct Yl -----
  # In both groups, Yl is measured by three indicators.
  # For the gold standard indicator Y_a:
  #   In Group 1 (audit): loading is fixed to 1 (marker) and its measurement error is minimal.
  #   In Group 2 (non-audit): Y_a is treated as pure error (i.e., no true signal) by fixing its intercept to 0 
  #   and its residual variance to 1.
  Yl =~ c(1, NA)*Y_a + c(lmbd_p, lmbd_p)*Y_p1 + c(lmbd_np, lmbd_np)*Y_np
  
  # ----- Measurement Model for Latent Construct Wl -----
  # Similarly for Wl:
  Wl =~ c(1, NA)*W_a + c(dlt_p, dlt_p)*W_p1 + c(dlt_np, dlt_np)*W_np
  
  # ----- Intercepts for the Indicators -----
  # For the gold standard indicators, fix intercepts to 0 in both groups:
  Y_a ~ c(0,0)*1
  W_a ~ c(0,0)*1
  # For the other indicators, assume intercepts are equal across groups (free otherwise):
  Y_p1 ~ c(tau1, tau1)*1
  Y_np ~ c(tau2, tau2)*1
  W_p1 ~ c(tau3, tau3)*1
  W_np ~ c(tau4, tau4)*1
  
  # ----- Residual Variances -----
  # For the gold standard indicators:
  # In Group 1, allow free estimation (or set to a very small value if needed);
  # In Group 2, fix the residual variances to 1 to force a pure error model.
  Y_a ~~ c(0,1)*Y_a
  W_a ~~ c(0,1)*W_a
  
      # For the other indicators, constrain them to be equal across groups:
  Y_p1 ~~ c(psi_Yp1, psi_Yp1)*Y_p1
  Y_np ~~ c(psi_Ynp, psi_Ynp)*Y_np
  W_p1 ~~ c(psi_Wp1, psi_Wp1)*W_p1
  W_np ~~ c(psi_Wnp, psi_Wnp)*W_np
'   
    # Fit the multigroup model
    # imp_fit <- function(data){
    # 
    #   if (any(diag(lavInspect(fit, "cov.lv")$`FALSE`) < 0 ) | any(diag(lavInspect(fit, "cov.lv")$`TRUE`) < 0 )){
    #     imp_data <- imp(data)
    #     fit <- imp_fit(data)
    #   }
    #   return(fit)
    # }
    # imp_data <- imp(data)
    fit <- sem(mod, data = imp_data, group = "group", group.equal = c("loadings", "intercepts"), estimator = "MLR")
    #fit <- imp_fit(data)
    summary(fit)
    theta <- lavInspect(fit, "theta")$`FALSE`
    hey <- any(diag(theta)<0)
    conv <- lavInspect(fit, "converged")
    Y_l <- lavPredict(fit)[[2]][,"Yl"]
    mx_hat <- mean(Y_l, na.rm = T) 
    s2x_hat <- var(Y_l, na.rm = T)
    y_hat <- lavPredictY(fit, ynames = "Y_np", xnames = c("Y_a", "Y_p1", "W_a", "W_p1", "W_np"))
    y_hat <- y_hat[y_hat$group==T,1]
    my_hat <- mean(y_hat, na.rm = T)
    s2y_hat <- var(y_hat,na.rm = T)
    rel <- reliability(fit)
    gamma_hat <- rel$`TRUE`["omega","Yl"]
  } else if (m_method == "multigroup_2"){
    data$group <- (!is.na(data$Y_a))
    imp <- function(data){
      sub_data <- data %>% 
        dplyr:::select(Y_p1, Y_np, W_p1, W_np)
      data$group <- !(rowSums(is.na(sub_data))==0 & is.na(data$Y_a))
      # Imputation process - Following Scholtus, et al. (2015)
      # Step 1: Generate initial random values from N(0,1) for missing data
      data$Y_a[data$group==F] <- rnorm(sum(data$group==F))
      data$W_a[data$group==F] <- rnorm(sum(data$group==F))
      
      # Step 2: Regress each missing variable on all other variables (group 2 only) and keep residuals 
      
      reg_y <- lm(Y_a ~ Y_p1 + Y_np + W_p1 + W_np, data = data[data$group==F,])
      data$Y_a[data$group == F] <- residuals(reg_y)
      
      reg_w <- lm(W_a ~ Y_p1 + Y_np + W_p1 + W_np, data = data, subset = (group==F))
      data$W_a[data$group == F] <- residuals(reg_w)
      
      reg_yw <- lm(Y_a ~ W_a, data = data, subset = (group == F))
      data$W_a[data$group == F] <- residuals(reg_yw)
      
      # Step 3: Standardize imputed variables (ensure mean 0, variance 1 in group 2)
      for(var in c("Y_a", "W_a")){
        mean_val <- mean(data[[var]][data$group == F])
        sd_val  <- sd(data[[var]][data$group == F])
        data[[var]][data$group == F] <- (data[[var]][data$group == F] - mean_val)/sd_val
      }
      
      mean(data$Y_a[data$group==F]); mean(data$W_a[data$group==F])
      var(data$Y_a[data$group==F]); var(data$Y_a[data$group==F])
      data[data$group==F,] %>% 
        select(Y_a, W_a, Y_p1, Y_np, W_p1, W_np) %>% 
        cov(use = "complete.obs") %>% 
        round(3)
      
      return(data)
    }
    
      imp_data <- imp(data)
      # Step 1 - Regression
      lm_Y_p1 <- lm(Y_p1 ~ Y_a, data = subset(data, group == T))
      lm_Y_np <- lm(Y_np ~ Y_a, data = subset(data, group == T))
      cy_p <- coef(lm_Y_p1)[2]; cy_np <- coef(lm_Y_np)[2]
      lm_W_p1 <- lm(W_p1 ~ W_a, data = subset(data, group == T))
      lm_W_np <- lm(W_np ~ W_a, data = subset(data, group == T))
      cw_p <- coef(lm_W_p1)[2]; cw_np <- coef(lm_W_np)[2]
      
      
      # Group 2 - SEM
      
      mod <- paste0(
        'Yl =~ 0*Y_a + ', round(cy_p, 3), '*Y_p1 + ', round(cy_np, 3), '*Y_np\n',
        'Wl =~ 0*W_a + ', round(cw_p, 3), '*W_p1 + ', round(cw_np, 3), '*W_np\n',
        'Y_a ~ 0*1\n',
        'Y_a ~~ 1*Y_a\n',
        'W_a ~ 0*1\n',
        'W_a ~~ 1*W_a\n',
        'Yl ~~ Wl'
      )
      
      # 
      # mod <- paste0('
      #   # Structural model
      #   Yl ~~ Wl
      #   # Measurement model
      #   Yl =~ NA*Y_a + ', round(cy_p, 3), '*Y_p1 + ', round(cy_np, 3), '*Y_np\n',
      #   'Wl =~ NA*W_a + ', round(cw_p, 3), '*W_p1 + ', round(cw_np, 3), '*W_np\n',
      #   'Y_a ~ NA * Y_a\n',
      #   'W_a ~ NA * W_a')
      # Fit the multigroup model
      
    fit <- sem(mod, data = imp_data, estimator = "MLR")
    conv <- lavInspect(fit, "converged")
    hey <- (any(lavInspect(fit, "theta")<0) | any(diag(lavInspect(fit, "cov.lv")) < 0 ))
    Y_l <- c(data$Y_a,lavPredict(fit)[, "Yl"])
    mx_hat <- mean(Y_l, na.rm = T)
    s2x_hat <- var(Y_l, na.rm = T)
    Y_np_hat <- c(as.numeric(predict(lm_Y_np)), as.numeric(lavPredictY(fit, ynames = "Y_np", xnames = c("Y_a", "Y_p1", "W_a", "W_p1", "W_np"))))  
    my_hat <- mean(Y_np_hat) 
    s2y_hat <- var(Y_np_hat)
    rel <- reliability(fit)
    gamma_hat <- rel["omega","Yl"]
  }

  x_hat <- mx_hat + gamma_hat*(s2x_hat/s2y_hat)*(Y_np - my_hat)
  m_err_hat <- abs(mean(Y_np, na.rm = T) - mean(x_hat[S==1], na.rm = T))
  
  # (2) Representation model
  r_yz <- cor(Y_np, Z, use = "complete.obs")
  r_zs <- cor(Z, S)
  if (r_method == "monte_carlo"){
    phi_guess <- numeric(0)
    while (length(phi_guess)<M){
      r_xs_samples <- runif(M, min = -1, max = 1)
      phi_samples <- sapply(r_xs_samples, function(rho_XS) {compute_phi(rho_XS, r_yz, r_zs)})
      phi_guess <- c(phi_guess, phi_samples[phi_samples >= 0 & phi_samples <= 1])
      phi_guess <- unique(phi_guess)
    }
    phi_guess <- sample(phi_guess, M, replace = FALSE)
    g_hat <- (phi_guess + (1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess))
    g_hat_m <- (phi_guess*sqrt(gamma_hat) + gamma_hat*(1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess)*sqrt(gamma_hat))
    c_obs <- sd(Y_np, na.rm = T)/sd(Z)
    c_obs_m <- (sd(Y_np, na.rm = T)*sqrt(gamma_hat))/sd(Z)
    r_err_hat_samples <- g_hat*c_obs*(mean(Z[S == 1], na.rm = T) - mean(Z))
    rm_err_hat_samples <- g_hat_m*c_obs_m*(mean(Z[S == 1], na.rm = T) - mean(Z))
    r_err_hat <- mean(r_err_hat_samples)
    rm_err_hat <- mean(rm_err_hat_samples)
  } else if (r_method == "rho_xs"){
    r_xs <- cor(Y_p1, S, use = "complete.obs")
    phi_guess <- compute_phi(r_xs, r_yz, r_zs)
    g_hat <- (phi_guess + (1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess))
    g_hat_m <- (phi_guess*sqrt(gamma_hat) + gamma_hat*(1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess)*sqrt(gamma_hat))
    c_obs <- sd(Y_np, na.rm = T)/sd(Z)
    c_obs_m <- (sd(Y_np, na.rm = T)*sqrt(gamma_hat))/sd(Z)
    r_err_hat <- g_hat*c_obs*(mean(Z[S == 1], na.rm = T) - mean(Z))
    rm_err_hat <- g_hat_m*c_obs_m*(mean(Z[S == 1], na.rm = T) - mean(Z))
  } else if (r_method == "calibration"){
    reg_coef <- coef(glm(S_np ~ Y_p1 + Z, data = data, family = "binomial"))
    phi_guess <- reg_coef[2]/(reg_coef[2] + reg_coef[3])
    g_hat <- (phi_guess + (1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess))
    g_hat_m <- (phi_guess*sqrt(gamma_hat) + gamma_hat*(1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess)*sqrt(gamma_hat))
    c_obs <- sd(Y_np, na.rm = T)/sd(Z)
    c_obs_m <- (sd(Y_np, na.rm = T)*sqrt(gamma_hat))/sd(Z)
    r_err_hat <- g_hat*c_obs*(mean(Z[S == 1], na.rm = T) - mean(Z))
    rm_err_hat <- g_hat_m*c_obs_m*(mean(Z[S == 1], na.rm = T) - mean(Z))
  }
  t_err_hat <- m_err_hat + rm_err_hat 
  
  
  main_res <- c(m_err, r_err, t_err, m_err_hat, r_err_hat, rm_err_hat, t_err_hat, gamma_hat, conv, hey, sX, sXn)
  add_res <- list(fit, x_hat)
  res <- list(main_res, add_res)
  return(res)
}

#### - (IV) Estimation - ####

Res <- as.data.frame(matrix(,nrow =nsim, ncol = 12))
pb <- progress_bar$new(total = nrow(param))
for (i in 1:nsim){
  pb$tick()
  temp_dat <- Dat[[i]]
  smubm_res <- smubm(X = "X", Y = c("Y_p1","Y_p2","Y_np"), W = c("W_p1","W_p2","W_np"), S = "S_np", Z = "Z", data = temp_dat, m_method = "effects_coding", M = 2000, r_method = "monte_carlo")
  Res[i,] <- smubm_res[[1]]
  rm(temp_dat)
}
colnames(Res) <- c("m_err", "r_err", "t_err", "m_err_hat", "r_err_hat", "rm_err_hat", "t_err_hat", "gamma_hat","conv","hey", "sX", "sXn")
Res_ec <- cbind(param, Res)
Res_ec$stp1_m <- rep("EC", nrow(Res_ec))

Res <- as.data.frame(matrix(,nrow =nsim, ncol = 12))
pb <- progress_bar$new(total = nrow(param))
for (i in 1:nsim){
  pb$tick()
  temp_dat <- Dat[[i]]
  smubm_res <- smubm(X = "X", Y = c("Y_p1","Y_p2","Y_np"), W = c("W_p1","W_p2","W_np"), S = "S_np", Z = "Z", data = temp_dat, m_method = "multigroup_1", M = 2000, r_method = "rho_ys")
  Res[i,] <- smubm_res[[1]]
  rm(temp_dat)
}
colnames(Res) <- c("m_err", "r_err", "t_err", "m_err_hat", "r_err_hat", "rm_err_hat", "t_err_hat", "gamma_hat","conv","hey", "sX", "sXn")
Res_mg1 <- cbind(param, Res)
Res_mg1$stp1_m <- rep("MG1", nrow(Res_mg1))

Res <- as.data.frame(matrix(,nrow =nsim, ncol = 12))
pb <- progress_bar$new(total = nrow(param))
for (i in 1:nsim){
  pb$tick()
  temp_dat <- Dat[[i]]
  smubm_res <- smubm(X = "X", Y = c("Y_p1","Y_p2","Y_np"), W = c("W_p1","W_p2","W_np"), S = "S_np", Z = "Z", data = temp_dat, m_method = "multigroup_2", M = 2000, r_method = "calibration")
  Res[i,] <- smubm_res[[1]]
  rm(temp_dat)
}
colnames(Res) <- c("m_err", "r_err", "t_err", "m_err_hat", "r_err_hat", "rm_err_hat", "t_err_hat", "gamma_hat","conv","hey", "sX", "sXn")
Res_mg2 <- cbind(param, Res)
Res_mg2$stp1_m <- rep("MG2", nrow(Res_mg2))

# Exploring the results 
FRes <- rbind(Res_ec, Res_mg1, Res_mg2)
FRes <- FRes %>%
  mutate(
    TE = abs(t_err_hat - t_err)/sX,
    ME = abs(m_err_hat - m_err)/sXn,
    RE = abs(r_err_hat - r_err)/sX,
    mRE = abs(rm_err_hat - r_err)/sX
  )

tMes <- aggregate(ME ~  rho_xs + sys_e + stp1_m, data = FRes, FUN = mean, na.rm = T)
sdMes <- aggregate(ME ~ rho_xs + sys_e + stp1_m, data = FRes, FUN = sd, na.rm = T)
tMes$sdME <- sdMes$ME


ggplot(tMes, aes(x = factor(rho_xs), y = ME, group = stp1_m, color = factor(stp1_m))) +
  geom_point(size = 3) +  # Mean values as points
  geom_errorbar(aes(ymin = ME - sdME, ymax = ME + sdME), width = 0.2) +  # Error bars
  facet_grid(. ~ sys_e, labeller = label_bquote(cols = epsilon[M]: .(sys_e))) +
  labs(
    x = expression(rho[xs]), 
    y = "SME",
    color = "Meas. Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )
# phi doesn't really matter in this plot!


ggplot(mMes, aes(x = stp1_m, y = rho_xs, fill = ME)) +
  geom_tile() +
  facet_grid(phi ~ sys_e) +
  scale_fill_viridis_c(option = "magma", direction = -1) +  
  labs(x = expression(rho[xz]), y = expression(rho[xs]), fill = "Measurement Error") +
  theme_minimal(base_size = 14)



# Calculate Differences in Errors
Res <- Res_ec %>%
  mutate(
    #dif_phi = abs(phi_g - phi),
    TE = abs(t_err_hat - t_err)/sX,
    ME = abs(m_err_hat - m_err)/sXn,
    RE = abs(r_err_hat - r_err)/sX,
    mRE = abs(rm_err_hat - r_err)/sX
  )

# Filter Data for phi_g = 0.5 and Reshape
Res_long <- Res %>%
  #filter(phi_g == 0.5) %>%
  pivot_longer(cols = c(mRE, ME, TE), 
               names_to = "Error_Type", 
               values_to = "Error_Value")

# Calculate IQR-Based Limits for Outlier Control
iqr_limits <- Res_long %>%
  group_by(Error_Type) %>%
  summarize(
    Q1 = quantile(Error_Value, 0.25, na.rm = TRUE),
    Q3 = quantile(Error_Value, 0.75, na.rm = TRUE)
  ) %>%
  mutate(IQR = Q3 - Q1, ymin = Q1 - 1.5 * IQR, ymax = Q3 + 1.5 * IQR)

# Generate Boxplot with Outlier Control
ggplot(Res_long, aes(x = interaction(rho_xz, rho_xs, sep = "|"), 
                     y = Error_Value, fill = Error_Type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  ggh4x::facet_grid2(phi ~ sys_e, 
                     labeller = label_bquote(
                       rows = phi:.(phi), 
                       cols = epsilon:.(sys_e)
                     )) +
  scale_fill_brewer(palette = "Set2", name = "Error Type") +
  labs(
    x = expression(paste(rho[XZ], " | ", rho[XS])), 
    y = "Error Value") +
  theme_minimal(base_size = 14) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.text= element_text(size=20), strip.text = element_text(
          size = 20, color = "dark green"),panel.grid.minor = element_blank(),
        legend.position = "top"
  ) +
  coord_flip()
ggsave("Plot_Sim.png", path = fold_graphs, width = 40, height = 60, units = "cm")

