## ------------------------------------------- ##
## Modeling Total Error using Multisource Data ##
## Script 1 - Simulation                       ##
## Written by: Santiago Gómez-Echeverry        ##
## Last modified: 13/10/2025                   ##
## ------------------------------------------- ##

#### - (I) Working space and packages - ####
rm(list = ls())
options(digits = 3)
options(scipe = 10)
packs <- c('copula', 'poLCA', 'lavaan', 'progress', 'ggplot2', 'lavaanPlot', 'dplyr', 'semTools', 'ggh4x', 'tidyr', 'Cairo', 'flextable', 'xtable', 'cols4all', 'stringr')
ipacks <- packs %in% rownames(installed.packages())
if(any(ipacks == F)){install.packages(packs[!ipacks])}
lapply(packs, library, character.only = T)
sys <- Sys.info()
fold_data <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/2 - Continuous Multisource TE/3 - Data")
fold_graphs <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/2 - Continuous Multisource TE/4 - Graphs & Tables")

compute_phi <- function(rho_XS, rho_XZ, rho_ZS) {
  num <- rho_XS - rho_XZ * rho_ZS
  denom <- (rho_XS - rho_XZ * rho_ZS) + (rho_ZS - rho_XZ * rho_XS)
  return(num / denom)
}

nice_palette <- c4a("brewer.set3", n = 11)

#### - (II) Data generating process - ####

set.seed(42)
N <- 5e3                                       # Population size
D <- 300                                       # Number of draws
f_p <- 0.30                                    # Sample fraction, probability sample
f_a <- 0.02                                    # Sample fraction, Audit sample
n_p <- N*f_p                                   # Probability sample size
n_a <- N*f_a                                   # Audit sample size
rho_xs <- rho_xz <- c(0.3, 0.5, 0.7)           # Possible correlations
phi <- c(0.4, 0.6, 0.8)                        # Possible phi's
rnd_e <- 0.5                                   # Random error
sys_e <- c(0.1, 0.3, 0.5)
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
  #slope_bias <- sample(c(-param[i, "sys_e"], param[i, "sys_e"]),1, replace = T, prob = c(0.8, 0.2))
  dat$Y_np <- ifelse(dat$S_np==1,(1-param[i, "sys_e"])*dat$X + rnorm(N, mean = 0, sd = rnd_e), NA)
  dat$Y_p1 <- ifelse(dat$S_p1==1, 0.8*dat$X + rnorm(N, mean = 0, sd = rnd_e), NA)
  dat$Y_p2 <- ifelse(dat$S_p2==1, 0.6*dat$X + rnorm(N, mean = 0, sd = rnd_e), NA)
  dat$Y_a <- ifelse(dat$S_a==1, dat$X, NA)
  dat$W_np <- ifelse(dat$S_np==1, 0.9*dat$W + rnorm(N, mean = 0, sd = rnd_e), NA)
  dat$W_p1 <- ifelse(dat$S_p1==1, 0.8*dat$W + rnorm(N, mean = 0, sd = rnd_e), NA)
  dat$W_p2 <- ifelse(dat$S_p2==1, 0.6*dat$W + rnorm(N, mean = 0, sd = rnd_e), NA)
  dat$W_a <- ifelse(dat$S_a==1, dat$W, NA)
  Dat[[i]] <- dat
}

#### - (III) SMUBM function - ####

# True values

smubm <- function(X, Y, W, S, Z, data, m_method = c("effects_coding", "latent_std", "marker", "multigroup", "two_step"), M, r_method = c("monte_carlo", "naive","adjusted")){
  m_method <- match.arg(m_method)
  r_method <- match.arg(r_method)
  
  # ## For tests
  # data <- Dat[[1]]
  # X <- data$X
  # Y_np <- data$Y_np
  # Y_p1 <- data$Y_p1
  # S <- data$S_np
  # Z <- data$Z
  # ##

  X <- data[, X]
  Y_np <- data[,Y[[3]]]
  Y_p1 <- data[,Y[[1]]]
  S <- data[, S]
  Z <- data[, Z]
  # True values
  m_err <- mean(Y_np, na.rm = T) - mean(X[S==1], na.rm = T)
  r_err <- mean(X[S==1], na.rm = T) - mean(X)
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
    fit <- sem(model = mod, data = data, effect.coding = T, meanstructure = T)
    summary(fit, standardized = T)
    conv <- lavInspect(fit, "converged")
    hey <- (any(lavInspect(fit, "theta")<0) | any(diag(lavInspect(fit, "cov.lv")) < 0 ))
    x_hat_std <- (lavPredict(fit)[,"Yl"]- mean(lavPredict(fit)[,"Yl"]))/sd(lavPredict(fit)[,"Yl"])
    mx_hat <- mean(x_hat_std)
    s2x_hat <- var(x_hat_std)
    my_hat <- lavPredictY(fit, ynames = "Y_np", xnames = c("Y_p1", "Y_p2","W_p1","W_p2", "W_np")) %>% 
      mean()
    s2y_hat <- lavPredictY(fit, ynames = "Y_np", xnames = c("Y_p1", "Y_p2","W_p1","W_p2", "W_np")) %>% 
      var() %>% 
      as.numeric()
    rel <- reliability(fit)
    gamma_hat <- rel[rownames(rel)=="omega", colnames(rel)=="Yl"]
    parms <- as.data.frame(parameterEstimates(fit, standardized = T))
    lambda <- parms[parms$op=="=~" & parms$rhs=="Y_np", colnames(parms)=="std.all"]
    Fit <- fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic", "bic"))
  } else if (m_method == "latent_std"){
    # (1) Measurement Model - Standardization of the latent variable
    mod <-  'Yl =~ NA*Y_p1 + lmbd2*Y_p2 + lmbd3*Y_np
             Wl =~ NA*W_p1 + dlt2*W_p2 + dlt3*W_np
             Yl ~~ Wl
             
             # Standardizing the Latent Variables
             Yl ~ 0*1
             Wl ~ 0*1
             Yl ~~ 1*Yl
             Wl ~~ 1*Wl
            '
    fit <- sem(model = mod, data = data, meanstructure = T)
    summary(fit, standardized =T)
    conv <- lavInspect(fit, "converged")
    hey <- (any(lavInspect(fit, "theta")<0) | any(diag(lavInspect(fit, "cov.lv")) < 0 ))
    x_hat_std <- (lavPredict(fit)[,"Yl"]- mean(lavPredict(fit)[,"Yl"]))/sd(lavPredict(fit)[,"Yl"])
    mx_hat <- mean(x_hat_std)
    s2x_hat <- var(x_hat_std)
    my_hat <- lavPredictY(fit, ynames = "Y_np", xnames = c("Y_p1", "Y_p2","W_p1","W_p2", "W_np")) %>% 
      mean()
    s2y_hat <- lavPredictY(fit, ynames = "Y_np", xnames = c("Y_p1", "Y_p2","W_p1","W_p2", "W_np")) %>% 
      var() %>% 
      as.numeric()
    rel <- reliability(fit)
    gamma_hat <- rel[rownames(rel)=="omega", colnames(rel)=="Yl"]
    parms <- as.data.frame(parameterEstimates(fit, standardized = T))
    lambda <- parms[parms$op=="=~" & parms$rhs=="Y_np", colnames(parms)=="std.all"]
    Fit <- fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic", "bic"))
  } else if (m_method == "marker"){
    # (1) Measurement Model - Marker variable approach
    mod <-  'Yl =~ 1*Y_a + lmbd2*Y_p1 + lmbd3*Y_np
             Wl =~ 1*W_a + dlt2*W_p1 +  dlt3*W_np
             Yl ~~ Wl

             Y_a ~~ 0*Y_a
             W_a ~~ 0*W_a
            '
    fit <- sem(model = mod, data = data, meanstructure = T)
    fitMeasures(fit)
    conv <- lavInspect(fit, "converged")
    hey <- (any(lavInspect(fit, "theta")<0) | any(diag(lavInspect(fit, "cov.lv")) < 0 ))
    x_hat_std <- (lavPredict(fit)[,"Yl"]- mean(lavPredict(fit)[,"Yl"]))/sd(lavPredict(fit)[,"Yl"])
    mx_hat <- mean(x_hat_std)
    s2x_hat <- var(x_hat_std)
    my_hat <- lavPredictY(fit, ynames = "Y_np", xnames = c("Y_a", "Y_p1","W_a","W_p1", "W_np")) %>% 
      mean()
    s2y_hat <- lavPredictY(fit, ynames = "Y_np", xnames = c("Y_a", "Y_p1","W_a","W_p1", "W_np")) %>% 
      var() %>% 
      as.numeric()
    rel <- reliability(fit)
    gamma_hat <- rel[rownames(rel)=="omega", colnames(rel)=="Yl"]
    parms <- as.data.frame(parameterEstimates(fit, standardized = T))
    lambda <- parms[parms$op=="=~" & parms$rhs=="Y_np", colnames(parms)=="std.all"]
    Fit <- fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic", "bic"))
  } else if (m_method == "multigroup"){
    imp <- function(data){
      sub_data <- data %>%
        dplyr::select(Y_p1, Y_np, W_p1, W_np)
      
      data$group <- !is.na(data$Y_a)
      data$Y_a[!data$group] <- rnorm(sum(!data$group))
      data$W_a[!data$group] <- rnorm(sum(!data$group))
    
    # Step 2: Regress each missing variable on all other variables (group 2 only) and keep residuals
    
    # For Y_a regression: require complete cases on Y_p1, Y_np, W_p1, and W_np
    idx_y <- which(!data$group & complete.cases(data[, c("Y_p1", "Y_np", "W_p1", "W_np")]))
    if(length(idx_y) > 0){
      reg_y <- lm(Y_a ~ Y_p1 + Y_np + W_p1 + W_np, data = data, subset = idx_y)
      data$Y_a[idx_y] <- residuals(reg_y)
    }
    
    # For W_a regression: require complete cases on W_p1 and W_np
    idx_w <- which(!data$group & complete.cases(data[, c("W_p1", "W_np")]))
    if(length(idx_w) > 0){
      reg_w <- lm(W_a ~ W_p1 + W_np + Y_p1 + Y_np, data = data, subset = idx_w)
      data$W_a[idx_w] <- residuals(reg_w)
    }
    
    # Step 3: Standardize imputed variables in group==FALSE (to have mean 0, variance 1)
    for(var in c("Y_a", "W_a")){
      idx <- which(!data$group)
      mean_val <- mean(data[[var]][idx], na.rm = TRUE)
      sd_val  <- sd(data[[var]][idx], na.rm = TRUE)
      data[[var]][idx] <- (data[[var]][idx] - mean_val) / sd_val
    }
    return(data)
  }
  
  # SEM model specification remains unchanged
  mod <- '
        Yl =~ c(0, 1)*Y_a + c(l_p, l_p)*Y_p1 + c(l_np, l_np)*Y_np
        Wl =~ c(0, 1)*W_a + c(d_p, d_p)*W_p1 + c(d_np, d_np)*W_np

        Y_a ~ c(0, 0)*1
        W_a ~ c(0, 0)*1

        Y_a ~~ c(1, 0)*Y_a
        W_a ~~ c(1, 0)*W_a

        Yl ~~ c(r, r)*Wl
        '
  
  # Process data with the imputation function
  imp_data <- imp(data)
  
  # Fit the multi-group SEM (using group indicator as defined in imp())
  fit <- sem(mod, data = imp_data, group = "group")
  summary(fit, standardized = TRUE)
  
  theta <- lavInspect(fit, "theta")$`FALSE`
  hey <- any(diag(theta)<0)
  conv <- lavInspect(fit, "converged")
  Y_l <- lavPredict(fit)[[2]][,"Yl"]
  mx_hat <- mean(Y_l, na.rm = T)
  s2x_hat <- var(Y_l, na.rm = T)
  y_hat <- lavPredictY(fit, ynames = "Y_np", xnames = c("Y_a", "Y_p1", "W_a", "W_p1", "W_np"))
  y_hat <- y_hat[,1]
  my_hat <- mean(y_hat, na.rm = T)
  s2y_hat <- var(y_hat,na.rm = T)
  rel <- reliability(fit)
  gamma_hat <- rel$`TRUE`["omega","Yl"]
  parms <- as.data.frame(parameterEstimates(fit, standardized = T))
  lambda <- parms[parms$op=="=~" & parms$rhs=="Y_np", colnames(parms)=="std.all"][1]
  Fit <- fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic", "bic"))
} else if (m_method == "two_step"){
    imp <- function(data){
      sub_data <- data %>%
        dplyr::select(Y_p1, Y_np, W_p1, W_np)
      
      data$group <- !is.na(data$Y_a)
      data$Y_a[!data$group] <- rnorm(sum(!data$group))
      data$W_a[!data$group] <- rnorm(sum(!data$group))
      
      # Step 2: Regress each missing variable on all other variables (group 2 only) and keep residuals
      
      # For Y_a regression: require complete cases on Y_p1, Y_np, W_p1, and W_np
      idx_y <- which(!data$group & complete.cases(data[, c("Y_p1", "Y_np", "W_p1", "W_np")]))
      if(length(idx_y) > 0){
        reg_y <- lm(Y_a ~ Y_p1 + Y_np + W_p1 + W_np, data = data, subset = idx_y)
        data$Y_a[idx_y] <- residuals(reg_y)
      }
      
      # For W_a regression: require complete cases on W_p1 and W_np
      idx_w <- which(!data$group & complete.cases(data[, c("W_p1", "W_np")]))
      if(length(idx_w) > 0){
        reg_w <- lm(W_a ~ W_p1 + W_np + Y_p1 + Y_np, data = data, subset = idx_w)
        data$W_a[idx_w] <- residuals(reg_w)
      }

      # Step 3: Standardize imputed variables in group==FALSE (to have mean 0, variance 1)
      for(var in c("Y_a", "W_a")){
        idx <- which(!data$group)
        mean_val <- mean(data[[var]][idx], na.rm = TRUE)
        sd_val  <- sd(data[[var]][idx], na.rm = TRUE)
        data[[var]][idx] <- (data[[var]][idx] - mean_val) / sd_val
      }
      return(data)
    }
    
    imp_data <- imp(data)  
    # Group 1 - Regression
    lm_y_p1 <- lm(Y_p1 ~ Y_a, data = subset(imp_data, group == T))
    l_p_est <- coef(lm_y_p1)[2]
    lm_y_np <- lm(Y_np ~ Y_a, data = subset(imp_data, group == T))
    l_np_est <- coef(lm_y_np)[2]
    lm_w_p1 <- lm(W_p1 ~ W_a, data = subset(imp_data, group == T))
    d_p_est <- coef(lm_w_p1)[2]
    lm_w_np <- lm(W_np ~ W_a, data = subset(imp_data, group == T))
    d_np_est <- coef(lm_w_np)[2]
    
    # Group 2 - SEM
    mod <- paste0(
      'Yl =~ na*Y_a + ', round(l_p_est, 3), '*Y_p1 + ', round(l_np_est, 3), '*Y_np\n',
      'Wl =~ 0*W_a + ', round(d_p_est, 3), '*W_p1 + ', round(d_np_est, 3), '*W_np\n',
      'Yl ~ Wl\n',
      'Y_a ~ 0*1\n',
      'W_a ~ 0*1\n',
      'Y_a ~~ NA*Y_a\n',
      'W_a ~~ NA*W_a\n'
    )
    
    fit <- sem(mod, data = imp_data[imp_data$group==F,])
    conv <- lavInspect(fit, "converged")
    hey <- (any(lavInspect(fit, "theta")<0) | any(diag(lavInspect(fit, "cov.lv")) < 0 ))
    Y_l <- c(data$Y_a,lavPredict(fit)[, "Yl"])
    mx_hat <- mean(Y_l, na.rm = T)
    s2x_hat <- var(Y_l, na.rm = T)
    Y_np_hat <- c(as.numeric(predict(lm_y_np)), 
                  as.numeric(lavPredictY(fit, ynames = "Y_np", xnames = c("Y_a", "Y_p1", "W_a", "W_p1", "W_np"))))
    my_hat <- mean(Y_np_hat)
    s2y_hat <- var(Y_np_hat)
    rel <- reliability(fit)
    gamma_hat <- rel["omega","Yl"]
    parms <- as.data.frame(parameterEstimates(fit, standardized = T))
    lambda <- parms[parms$op=="=~" & parms$rhs=="Y_np", colnames(parms)=="std.all"][1]
    Fit <- fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic", "bic"))
  }
  x_hat <- mx_hat + gamma_hat*(s2x_hat/s2y_hat)*(Y_np - my_hat)
  m_err_hat <- mean(Y_np, na.rm = T) - mean(x_hat[S==1], na.rm = T)
  # (2) Representation model
  r_yz <- cor(Y_np, Z, use = "complete.obs")
  r_zs <- cor(Z, S)
  if (r_method == "monte_carlo"){
    grid_n    <- 10000
    rho_grid  <- seq(-1, 1, length.out = grid_n)
    psi_grid  <- sapply(rho_grid, function(rho_xs) {
      compute_phi(rho_xs, r_yz, r_zs)
    })
    psi_valid <- psi_grid[is.finite(psi_grid) & psi_grid >= 0 & psi_grid <= 1]
    psi_min <- min(psi_valid)
    psi_max <- max(psi_valid)
    phi_guess <- runif(M, min = psi_min, max = psi_max)
    # MUB
    #phi_guess <- runif(M, min =0, max = 1)
    g_hat <- (phi_guess + (1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess))
    c_obs <- sd(Y_np, na.rm = T)/sd(Z)
    r_err_hat_samples <- g_hat*c_obs*(mean(Z[S == 1], na.rm = T) - mean(Z))
    r_err_hat <- mean(r_err_hat_samples)
    # MUBM  
    g_hat_m <- (phi_guess*sqrt(gamma_hat) + gamma_hat*(1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess)*sqrt(gamma_hat))
    c_obs_m <- (sd(Y_np, na.rm = T)*sqrt(gamma_hat))/sd(Z)
    rm_err_hat_samples <- g_hat_m*c_obs_m*(mean(Z[S == 1], na.rm = T) - mean(Z))
    rm_err_hat <- mean(rm_err_hat_samples)
    # epsilon_R
    rx_err_hat <- rm_err_hat/lambda
    r_err; r_err_hat; rm_err_hat; rx_err_hat
    
  } else if (r_method == "naive"){
    r_xs <- cor(Y_p1, S, use = "complete.obs")
    phi_guess <- compute_phi(r_xs, r_yz, r_zs)
    # MUB
    g_hat <- (phi_guess + (1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess))
    c_obs <- sd(Y_np, na.rm = T)/sd(Z)
    r_err_hat <- g_hat*c_obs*(mean(Z[S == 1], na.rm = T) - mean(Z))
    # MUBM
    g_hat_m <- (phi_guess*sqrt(gamma_hat) + gamma_hat*(1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess)*sqrt(gamma_hat))
    c_obs_m <- (sd(Y_np, na.rm = T)*sqrt(gamma_hat))/sd(Z)
    rm_err_hat <- g_hat_m*c_obs_m*(mean(Z[S == 1], na.rm = T) - mean(Z))
    # epsilon_R
    rx_err_hat <- rm_err_hat/lambda
  } else if (r_method == "adjusted"){
    reg_coef <- coef(glm(S_np ~ Y_p1 + Z, data = data, family = "binomial"))
    phi_guess <- reg_coef[2]/(reg_coef[2] + reg_coef[3])
    # MUB
    g_hat <- (phi_guess + (1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess))
    c_obs <- sd(Y_np, na.rm = T)/sd(Z)
    r_err_hat <- g_hat*c_obs*(mean(Z[S == 1], na.rm = T) - mean(Z))
    # MUBM
    g_hat_m <- (phi_guess*sqrt(gamma_hat) + gamma_hat*(1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess)*sqrt(gamma_hat))
    c_obs_m <- (sd(Y_np, na.rm = T)*sqrt(gamma_hat))/sd(Z)
    rm_err_hat <- g_hat_m*c_obs_m*(mean(Z[S == 1], na.rm = T) - mean(Z))
    # epsilon_R
    rx_err_hat <- rm_err_hat/lambda
  }
  t_err_hat <- m_err_hat + rx_err_hat 
  main_res <- c(m_err, r_err, t_err, m_err_hat, r_err_hat, rm_err_hat, rx_err_hat, t_err_hat, gamma_hat, conv, hey, sX, sXn, Fit)
  add_res <- list(fit, x_hat)
  res <- list(main_res, add_res)
  return(res)
}

#### - (IV) Estimation - ####

specs <- expand.grid(c("effects_coding", "marker","multigroup", "two_step"),
                     c("naive", "monte_carlo", "adjusted"))
nspecs <- nrow(specs)

full_res <- list()

for (m in 1:nspecs){
  m_meth <- as.character(specs[m, 1])
  r_meth <- as.character(specs[m, 2])
  Res <- as.data.frame(matrix(, nrow = nsim, ncol = 24))
  pb <- progress_bar$new(format = paste0("Running Spec. ", m, " of ", nspecs, "[:bar] :percent eta: :eta"), total = nsim, clear = F, width = 120)
  for (i in 1:nsim){
    pb$tick()
    temp_dat <- Dat[[i]]
    smub_res <- smubm(X = "X", Y = c("Y_p1", "Y_p2", "Y_np"), W = c("W_p1", "W_p2", "W_np"), S = "S_np", Z = "Z", data = temp_dat, m_method = m_meth, M = 2000, r_method = r_meth)
    Res[i,] <- smub_res[[1]]
    rm(temp_dat)
  }
  colnames(Res) <- c("m_err", "r_err", "t_err", "m_err_hat", "r_err_hat", "rm_err_hat","rx_err_hat", "t_err_hat", "gamma_hat", "conv", "hey", "sX", "sXn", 
                     "chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic", "bic")
  Res_m <- cbind(param, Res)
  Res_m$stp1 <- rep(m_meth, nrow(Res_m))
  Res_m$stp2 <- rep(r_meth, nrow(Res_m))
  full_res[[m]] <- Res_m
}

FRes <- do.call(rbind, full_res)
write.table(FRes, file = paste0(fold_data, "/Results_Simulation.txt"))

#### - (V) Graphs and Tables - ####

setwd(fold_data)
FRes <- read.table("Results_Simulation.txt")

FRes$stp1[FRes$stp1=="effects_coding"] <- "Effects coding"
FRes$stp1[FRes$stp1=="marker"] <- "Marker"
FRes$stp1[FRes$stp1=="multigroup"] <- "Multigroup"
FRes$stp1[FRes$stp1=="two_step"] <- "Constrained"
FRes$stp2[FRes$stp2=="monte_carlo"] <- "Monte Carlo"
FRes$stp2[FRes$stp2=="naive"] <- "Naïve corr."
FRes$stp2[FRes$stp2=="adjusted"] <- "Adjusted corr."

metrics <- FRes %>%
  mutate(m_bias = m_err_hat - m_err,
    r_bias = rm_err_hat - r_err,
    t_bias = t_err_hat - t_err,
    
    m_abs_error = abs(m_err_hat - m_err),
    r_abs_error = abs(rm_err_hat - r_err),
    t_abs_error = abs(t_err_hat - t_err),
    
    m_sq_error = (m_err_hat - m_err)^2,
    r_sq_error = (rm_err_hat - r_err)^2,
    t_sq_error = (t_err_hat - t_err)^2,
    
    m_rae = m_abs_error / abs(m_err),
    r_rae = r_abs_error / abs(r_err),
    t_rae = t_abs_error / abs(t_err)) 

summary_tab1 <- metrics %>%
  filter(conv == 1, hey == 0, rho_xs == 0.5, phi == 0.4) %>%
  group_by(stp1, stp2, sys_e) %>%
  summarise(
    RMSE_M = sqrt(mean(m_sq_error, na.rm = T)),
    RMSE_R = sqrt(mean(r_sq_error, na.rm = T)),
    MAE_M = mean(m_abs_error, na.rm = T),
    MAE_R = mean(r_abs_error, na.rm = T),
    CV_M = sd(m_abs_error, na.rm = T) / mean(m_abs_error, na.rm = T),
    CV_R = sd(r_abs_error, na.rm = T) / mean(r_abs_error, na.rm = T)
  )
summary_tab1
tab1_sim <- xtable(summary_tab1) 
setwd(fold_graphs)
print(tab1_sim, file = "Table_1.tex", include.rownames = F)  

summary_tab2 <- metrics %>%
  group_by(stp1, stp2) %>%
  filter(conv == 1 & hey == 0 & rho_xs == 0.5 & phi == 0.4) %>% 
  summarise(
    RMSE_M = (sqrt(mean(m_sq_error, na.rm = TRUE))),
    RMSE_R = (sqrt(mean(r_sq_error, na.rm = TRUE))),
    MAE_M = mean(m_abs_error, na.rm = T),
    MAE_R = mean(r_abs_error, na.rm = T),
    CV_M = sd(m_abs_error, na.rm = T)/mean(m_abs_error, na.rm = T),
    CV_R = sd(r_abs_error, na.rm = T)/mean(r_abs_error, na.rm = T))
tab2_sim <- xtable(summary_tab2) 
print(tab2_sim, file = "Table_2.tex", include.rownames = F)  

## - Tile plots

metrics$stp1 <- factor(metrics$stp1, levels = c("Effects coding", "Marker", "Multigroup", "Constrained"))
pt <- metrics %>%
  group_by(rho_xs, rho_xz, sys_e, stp1) %>%
  filter(phi == 0.4 & hey==0 & stp2=="Naïve corr.") %>% 
  summarise(RMSE = (sqrt(mean(t_sq_error, na.rm = TRUE)))) %>% 
  ggplot(aes(x = rho_xs, y = rho_xz, fill = RMSE)) +
  geom_tile() +
  geom_text(aes(label = round(RMSE, 2)), size = 12, color = "white") +
  facet_grid(sys_e ~ stp1, labeller = label_bquote(rows = epsilon[s]:.(sys_e)))+ 
  scale_fill_viridis_c(limits = c(0, 0.8), breaks = c(0, 0.35, 0.7), option = "magma") +
  scale_x_continuous(breaks = c(0.3, 0.5, 0.7), labels = c("0.3", "0.5", "0.7"))+
  scale_y_continuous(breaks = c(0.3, 0.5, 0.7), labels = c("0.3", "0.5", "0.7"))+
  labs(x = expression(rho[XS]), y = expression(rho[XZ])) +
  theme_minimal() + theme(text = element_text(size = 40), legend.key.height = unit(1, "cm")) 
pt
ggsave("TE_Tiles_3.png", plot = pt, path = fold_graphs, width = 60, height = 50, units = "cm")

pt <- metrics %>%
  group_by(rho_xs, rho_xz, sys_e, stp1) %>%
  filter(phi == 0.6 & hey==0 & stp2=="Naïve corr.") %>% 
  summarise(RMSE = (sqrt(mean(t_sq_error, na.rm = TRUE)))) %>%
  ggplot(aes(x = rho_xs, y = rho_xz, fill = RMSE)) +
  geom_tile() +
  geom_text(aes(label = round(RMSE, 2)), size = 12, color = "white") +
  facet_grid(sys_e ~ stp1, labeller = label_bquote(rows = epsilon[s]:.(sys_e)))+ 
  scale_fill_viridis_c(limits = c(0, 0.8), breaks = c(0, 0.35, 0.7), option = "magma") +
  scale_x_continuous(breaks = c(0.3, 0.5, 0.7), labels = c("0.3", "0.5", "0.7"))+
  scale_y_continuous(breaks = c(0.3, 0.5, 0.7), labels = c("0.3", "0.5", "0.7"))+
  labs(x = expression(rho[XS]), y = expression(rho[XZ])) +
  theme_minimal() + theme(text = element_text(size = 40), legend.key.height = unit(1, "cm")) 
ggsave("TE_Tiles_5.png", plot = pt, path = fold_graphs, width = 60, height = 50, units = "cm")

pt <- metrics %>%
  group_by(rho_xs, rho_xz, sys_e, stp1) %>%
  filter(phi == 0.8 & hey==0 & stp2=="Naïve corr.") %>% 
  summarise(RMSE = (sqrt(mean(t_sq_error, na.rm = TRUE)))) %>%
  ggplot(aes(x = rho_xs, y = rho_xz, fill = RMSE)) +
  geom_tile() +
  geom_text(aes(label = round(RMSE, 2)), size = 12, color = "white") +
  facet_grid(sys_e ~ stp1, labeller = label_bquote(rows = epsilon[s]:.(sys_e)))+ 
  scale_fill_viridis_c(limits = c(0, 0.8), breaks = c(0, 0.35, 0.7), option = "magma") +
  scale_x_continuous(breaks = c(0.3, 0.5, 0.7), labels = c("0.3", "0.5", "0.7"))+
  scale_y_continuous(breaks = c(0.3, 0.5, 0.7), labels = c("0.3", "0.5", "0.7"))+
  labs(x = expression(rho[XS]), y = expression(rho[XZ])) +
  theme_minimal() + theme(text = element_text(size = 40), legend.key.height = unit(1, "cm")) 
ggsave("TE_Tiles_7.png", plot = pt, path = fold_graphs, width = 60, height = 50, units = "cm")

# Figure 6

metrics2 <- metrics %>%
  mutate(
    abs_error_naive = abs((r_err - r_err_hat)/r_err),
    abs_error_model = abs((r_err - rx_err_hat)/r_err),
    sq_error_naive = ((r_err - r_err_hat)/r_err)^2,
    sq_error_model = ((r_err - rx_err_hat)/r_err)^2,
    improvement = abs_error_naive - abs_error_model)

summary_improvement <- metrics2 %>%
  filter(stp1 == "Marker" & phi == 0.6 & hey == 0) %>% 
  group_by(stp2, sys_e, rho_xz, rho_xs) %>%
  summarise(
    mean_abs_naive     = mean(abs_error_naive, na.rm = TRUE),
    sd_abs_naive       = sd(abs_error_naive, na.rm = TRUE),
    mean_abs_model     = mean(abs_error_model, na.rm = TRUE),
    sd_abs_model       = sd(abs_error_model, na.rm = TRUE),
    mean_improvement   = mean(improvement, na.rm = TRUE),
    sd_improvement     = sd(improvement, na.rm = TRUE),
    .groups = 'drop')


imp <- ggplot(summary_improvement, aes(x = rho_xz, y = mean_improvement, color = factor(rho_xs))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_improvement - sd_improvement,
                    ymax = mean_improvement + sd_improvement),
                width = 0.02, alpha = 0.5) +
  facet_grid(sys_e ~ stp2, labeller = label_bquote(rows = epsilon[s]:.(sys_e), cols = S2:.(stp2))) +
  labs(
    x = expression(rho[xz]), y = "|MUB RE error - MUBM RE error|",
    color = expression(rho[xs])
  ) +
  theme_minimal()
imp
ggsave("Impr_RE.png", plot = imp, path = fold_graphs, units = "cm")

summary_improvement <- summary_improvement %>%
  mutate(se_improvement = sd_improvement / sqrt(n()))

imp2 <- ggplot(summary_improvement, 
                        aes(x = factor(rho_xz), y = mean_improvement, color = factor(rho_xs))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = mean_improvement - 1.96 * se_improvement,
                    ymax = mean_improvement + 1.96 * se_improvement),
                width = 0.2,
                position = position_dodge(width = 0.5)) +
  scale_color_brewer(palette = "Set1", name = expression(rho[xs])) +
  facet_grid(sys_e ~ stp2, labeller = label_bquote(rows = epsilon[s]:.(sys_e), 
                                                   cols = S2:.(stp2))) +
  labs(
    y = "Mean(MUB Rel. Err. - MUBM Rel. Err.)", 
    x = expression(rho[xz]),
    color = expression(rho[xs])
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.border = element_rect(color = "gray70", fill = NA, size = 0.5),
    strip.background = element_rect(fill = "gray90", color = NA),
    panel.spacing = unit(0.5, "lines")
  )
imp2
ggsave("Impr_RE2.png", plot = imp2, path = fold_graphs, width = 20, height = 15,units = "cm")


## Fit - Figure C1

table(FRes$hey)/nrow(FRes)
class(FRes$phi)

# We don't really care about the second step for these plots
FRes$stp1 <- factor(FRes$stp1, levels = c("Effects coding", "Marker", "Multigroup", "Constrained"))

p_rmsea <- FRes %>% 
  filter(rho_xs == 0.5 & phi == 0.6 & conv == 1 & hey == 0) %>% 
  ggplot(aes(x = as.factor(sys_e), y = rmsea, fill = as.factor(sys_e))) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "navy")  +
  facet_wrap(~stp1, scales = "free") +
  scale_fill_brewer(palette="Spectral")

p_cfi <- FRes %>% 
  filter(rho_xs == 0.5 & phi == 0.6 & conv == 1 & hey == 0) %>% 
  ggplot(aes(x = as.factor(sys_e), y = cfi, fill = as.factor(stp1))) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0.90, linetype = "dashed", color = "navy") +
  labs(y = "CFI", "Systematic Error") +
  scale_fill_brewer(palette = "Spectral") +
  scale_y_continuous(limits = quantile(FRes$cfi,c(0.1, 0.9)))

p_tli <- FRes %>% 
  filter(rho_xs == 0.5 & phi == 0.6 & conv == 1 & hey == 0) %>% 
  ggplot(aes(x = as.factor(sys_e), y = tli, fill = as.factor(stp1))) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0.90, linetype = "dashed", color = "navy") +
  labs(y = "CFI", "Systematic Error") +
  scale_fill_brewer(palette = "Spectral") +
  scale_y_continuous(limits = quantile(FRes$cfi,c(0.1, 0.9)))

p_srmr <- FRes %>% 
  filter(rho_xs == 0.5 & phi == 0.6 & conv == 1 & hey == 0) %>% 
  ggplot(aes(x = as.factor(sys_e), y = srmr, fill = as.factor(stp1))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "navy") +
  labs(y = "SRMR", "Systematic Error") +
  scale_fill_brewer(palette = "Spectral") +
  scale_y_continuous(limits = quantile(FRes$srmr, c(0.1, 0.9)))

p_chi <- FRes %>% 
  filter(rho_xs == 0.5 & phi == 0.6 & conv == 1 & hey == 0) %>% 
  ggplot(aes(x = as.factor(sys_e), y = pvalue, fill = as.factor(stp1))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "navy") +
  labs(y = "Chi2, p-value", "Systematic Error") +
  scale_fill_brewer(palette = "Spectral") 
# Convert to long format
fit_mes <- c("χ2, p-value", "RMSEA","SRMR","CFI", "TLI")
FRes_long <- FRes %>% 
  filter(rho_xs == 0.5, phi == 0.6, conv == 1, hey == 0) %>% 
  pivot_longer(cols = c(rmsea, cfi, tli, srmr, pvalue),
    names_to = "metric",
    values_to = "value") %>% 
  mutate(metric = str_to_upper(metric)) %>% 
  mutate(metric=recode(metric, "PVALUE"="χ2, p-value")) %>% 
  mutate(metric = factor(metric, levels = fit_mes))

ref_lines_df <- data.frame(metric = factor(fit_mes, levels = fit_mes),
                           yintercept = c(0.05, 0.05, 0.08, 0.95, 0.95))

# Plot
p_fit <- ggplot(FRes_long, aes(x = as.factor(sys_e), y = value, color = as.factor(sys_e))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  ggh4x::facet_grid2(stp1 ~ metric, scales = "free", independent = "all") +
  scale_color_brewer(palette = "Set1", name = expression(epsilon[s])) +
  geom_hline(data = ref_lines_df,
             aes(yintercept = yintercept),
             linetype = "dashed",
             color = "navy",
             inherit.aes = FALSE) + labs(x = "", y = "Value") +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5)) +
  theme(legend.position = "right", legend.title = element_text(size = 30),
        legend.key.size = unit(1.5, 'cm'), text = element_text(size = 20) )
ggsave("Fit.png", plot = p_fit, path = fold_graphs, width = 40, height = 26.5, units = "cm")
