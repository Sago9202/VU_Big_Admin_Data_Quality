## ------------------------------------------- ##
## Modeling Total Error using Multisource Data ##
## Script 1 - Simulation                       ##
## Written by: Santiago Gómez-Echeverry        ##
## Last modified: 07/03/2025                   ##
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
N <- 1e3                                       # Population size
D <- 50                                        # Number of draws
f_p <- 0.25                                    # Sample fraction, probability sample
f_a <- 0.05                                    # Sample fraction, Audit sample
n_p <- N*f_p                                   # Probability sample size
n_a <- N*f_a                                   # Audit sample size
rho_xs <- rho_xz <- c(0.3, 0.5, 0.7)  # Possible correlations
phi <- c(0.4, 0.6, 0.8)                        # Possible phi's
rnd_e <- 0.1                                   # Random error
sys_e <- c(0.4, 0.6, 0.8)
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
  dat[,3] <- qbinom(dat[,3], size = 1, prob = 0.6)
  dat <- data.frame(unlist(dat))
  dat$id <- seq.int(nrow(dat))
  dat <- dat[,c(4,1,2,3)]
  colnames(dat) <- c('id', 'X', 'Z', 'S_np')
  dat$W <- rnorm(N, mean = 0, sd = 1)
  dat$Y <- (1 - param[i, "sys_e"]*0.3) * dat$X +   param[i, "sys_e"]*0.3 * dat$W
  dat$Y1 <- (1 - param[i, "sys_e"]*0.2) * dat$X +   param[i, "sys_e"]*0.2 * dat$W
  dat$Y2 <- (1 - param[i, "sys_e"]*0.4) * dat$X +   param[i, "sys_e"]*0.4 * dat$W
  dat$pr_p1 <- runif(N)
  dat <- dat[order(dat$pr_p1),]
  dat$S_p1 <- 0
  dat$S_p1[1:n_p] <- 1
  dat$pr_p2 <- runif(N)
  dat <- dat[order(dat$pr_p2),]
  dat$S_p2 <- 0
  dat$S_p2[1:n_p] <- 1
  dat$pr_a <- runif(N)
  dat$pr_a[dat$S_p1==1 & dat$S_np==1] <- 0
  dat <- dat[order(dat$pr_a),]
  dat$S_a <- 0
  dat$S_a[1:n_a] <- 1
  dat <- dat[order(dat$id),]
  dat$Y_np <- ifelse(dat$S_np==1, dat$Y + rnorm(N, mean = 0, sd = rnd_e), NA)
  dat$Y_p1 <- ifelse(dat$S_p1==1, dat$Y1 + rnorm(N, mean = 0, sd = rnd_e), NA)
  dat$Y_p2 <- ifelse(dat$S_p2==1, dat$Y2 + rnorm(N, mean = 0, sd = rnd_e), NA)
  dat$Y_a <- ifelse(dat$S_a==1, dat$X, NA)
  dat$W_np <- ifelse(dat$S_np==1, dat$W + rnorm(N, mean = 0, sd = rnd_e), NA)
  dat$W_p1 <- ifelse(dat$S_p1==1, 0.8*dat$W + rnorm(N, mean = 0, sd = rnd_e), NA)
  dat$W_p2 <- ifelse(dat$S_p2==1, 0.6*dat$W + rnorm(N, mean = 0, sd = rnd_e), NA)
  dat$W_a <- ifelse(dat$S_a==1, dat$W, NA)
  Dat[[i]] <- dat
}

#### - (III) SMUBM function - ####

# True values

smubm <- function(X, Y, W, S, Z, data, N_mc, m_method = c("effects_coding", "multigroup_1", "multigroup_2")){
  m_method <- match.arg(m_method)
  X <- data[, X]
  Y_np <- data[,Y[[3]]]
  S <- data[, S]
  Z <- data[, Z]
  # True values
  m_err <- abs(mean(Y_np, na.rm = T) - mean(X[S==1], na.rm = T))
  r_err <- abs(mean(X[S==1], na.rm = T) - mean(X))
  t_err <- m_err + r_err
  
  if (m_method == "effects_coding"){
    # (1) Measurement Model - Effects coding
    mod <-  'Yl =~ lmbd1*Y_p1 + lmbd2*Y_p2 + lmbd3*Y_np
             Wl =~ dlt1*W_p1 + dlt2*W_p2 + dlt3*W_np
             Yl ~~ Wl
             
             # Effects coding constraint
            lmbd1 + lmbd2 + lmbd3 == 3
            '
    fit <- sem(model = mod, data = data, effect.coding = T, meanstructure = T)
    sum_fit <- summary(fit) 
    sum_fit$pe
  
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
        select(Y_p1, Y_np, W_p1, W_np)
      
      data$group <- (rowSums(is.na(sub_data))==0 & is.na(data$Y_a))
      # Imputation process - Following Scholtus, et al. (2015)
      # Step 1: Generate initial random values from N(0,1) for missing data
      data$Y_a[data$group==1] <- rnorm(sum(data$group==1))
      data$W_a[data$group==1] <- rnorm(sum(data$group==1))
      
      # Step 2: Regress each missing variable on all other variables (group 2 only) and keep residuals 
      
      #sub_data <- na.omit(data[data$S_a==0,])
      reg_ya <- lm(Y_a ~ Y_p1 + Y_np + W_p1 + W_np, data = data, subset = (group==1))
      data$Y_a[data$group == 1] <- residuals(reg_ya)
      
      reg_wa <- lm(W_a ~ Y_p1 + Y_np + W_p1 + W_np, data = data, subset = (group==1))
      data$W_a[data$group == 1] <- residuals(reg_wa)
      
      # Step 3: Standardize imputed variables (ensure mean 0, variance 1 in group 2)
      for(var in c("Y_a", "W_a")){
        mean_val <- mean(data[[var]][data$group == 1])
        sd_val  <- sd(data[[var]][data$group == 1])
        data[[var]][data$group == 1] <- (data[[var]][data$group==1] - mean_val)/sd_val
      }
      return(data)
    }
    
    # SEM model
    
    mod <- '
    # Structural model
      Yl ~ Wl  # Allowing the structural relation to be estimated

    # Measurement model: Fix Y_a as the reference variable in Group 2
      Yl =~ c(1,a) * Y_a + c(lmbd_p, lmbd_p) * Y_p1 + c(lmbd_np, lmbd_np) * Y_np
      Wl =~ c(1,a) * W_a + c(dlt_p, dlt_p) * W_p1 + c(dlt_np, dlt_np) * W_np

    # Fix intercepts in both groups
      Y_a ~ c(0, 0) * 1
      W_a ~ c(0, 0) * 1

    # Fix residual variances in Group 2 for stability
      Y_a ~~ c(1, NA) * Y_a
      W_a ~~ c(1, NA) * W_a
    '
    
    # Fit the multigroup model
    imp_fit  <- function(data){
      imp_data <- imp(data)
      fit <- sem(mod, data = imp_data, group = "group")
      var <- lavInspect(fit, "cov.lv")
      neg <- ifelse(sum(diag(var$`TRUE`)<0, diag(var$`FALSE`)<0)>0, T, F)
      if (neg == T){
        return(imp_fit(data))
      }
      return(fit)
    }
    fit <- imp_fit(data)
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
    sub_data <- data %>% 
      select(Y_p1, Y_np, W_p1, W_np)
    data$group <- (rowSums(is.na(sub_data))==0 & is.na(data$Y_a))
    # Imputation process - Following Scholtus, et al. (2015)
    # Step 1: Generate initial random values from N(0,1) for missing data
    data$Y_a[data$group==1] <- rnorm(sum(data$group==1))
    data$W_a[data$group==1] <- rnorm(sum(data$group==1))
    
    # Step 2: Regress each missing variable on all other variables (group 2 only) and keep residuals 
    
    #sub_data <- na.omit(data[data$S_a==0,])
    reg_ya <- lm(Y_a ~ Y_p1 + Y_np + W_p1 + W_np, data = data, subset = (group==1))
    data$Y_a[data$group == 1] <- residuals(reg_ya)
    
    reg_wa <- lm(W_a ~ Y_p1 + Y_np + W_p1 + W_np, data = data, subset = (group==1))
    data$W_a[data$group == 1] <- residuals(reg_wa)
    
    # Step 3: Standardize imputed variables (ensure mean 0, variance 1 in group 2)
    for(var in c("Y_a", "W_a")){
      mean_val <- mean(data[[var]][data$group == 1])
      sd_val  <- sd(data[[var]][data$group == 1])
      data[[var]][data$group == 1] <- (data[[var]][data$group==1] - mean_val)/sd_val
    }
    
    # Group 1 - Regression
    lm_Y_p1 <- lm(Y_p1 ~ Y_a, data = subset(data, group == 0))
    lm_Y_np <- lm(Y_np ~ Y_a, data = subset(data, group == 0))
    cy_p <- coef(lm_Y_p1)[2]
    cy_np <- coef(lm_Y_np)[2]
    lm_W_p1 <- lm(W_p1 ~ W_a, data = subset(data, group == 0))
    lm_W_np <- lm(W_np ~ W_a, data = subset(data, group == 0))
    cw_p <- coef(lm_W_p1)[2]
    cw_np <- coef(lm_W_np)[2]
    
    
    # Group 2 - SEM
    mod <- paste0('
        # Structural model
        Yl ~ Wl
        # Measurement model
        Yl =~ a*Y_a +', round(cy_p, 3), '*Y_p1 + ', round(cy_np, 3), '*Y_np\n',
        'Wl =~ a*W_a +', round(cw_p, 3), '*W_p1 + ', round(cw_np, 3), '*W_np\n',
        'Y_a ~ 0*1\n',
        'Y_a ~~ 1*Y_a\n',
        'W_a ~ 0*1\n',
        'W_a ~~ 1*W_a')
    # Fit the multigroup model
    data_g2 <- subset(data, group == T) 
    fit <- sem(mod, data = data_g2)
    summary(fit)
    
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
  r_yz <- cor(Y_np,Z, use = "complete.obs")
  r_zs <- cor(Z, S)
  r_xs_samples <- runif(N_mc, min = -1, max = 1)
  phi_samples <- sapply(r_xs_samples, function(rho_XS) {compute_phi(rho_XS, r_yz, r_zs)})
  # Filter the rho_XS values to keep only those where phi is within [0, 1]
  phi_guess <- phi_samples[phi_samples >= 0 & phi_samples <= 1]
  g_hat <- (phi_guess + (1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess))
  g_hat_m <- (phi_guess*sqrt(gamma_hat) + gamma_hat*(1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess)*sqrt(gamma_hat))
  c_obs <- sd(Y_np, na.rm = T)/sd(Z)
  c_obs_m <- (sd(Y_np, na.rm = T)*sqrt(gamma_hat))/sd(Z)
  r_err_hat_samples <- g_hat*c_obs*(mean(Z[S == 1], na.rm = T) - mean(Z))
  rm_err_hat_samples <- g_hat_m*c_obs_m*(mean(Z[S == 1], na.rm = T) - mean(Z))
  r_err_hat <- mean(r_err_hat_samples)
  rm_err_hat <- mean(rm_err_hat_samples)
  t_err_hat <- m_err_hat + rm_err_hat
  #m_err_hat; r_err_hat; rm_err_hat; t_err_hat
  
  main_res <- c(m_err, r_err, t_err, m_err_hat, r_err_hat, rm_err_hat, t_err_hat, gamma_hat)
  add_res <- list(fit, x_hat)
  res <- list(main_res, add_res)
  return(res)
}

#### - (IV) Estimation - ####

Res <- as.data.frame(matrix(,nrow =nsim, ncol = 8))
pb <- progress_bar$new(total = nrow(param))
for (i in 1:nsim){
  pb$tick()
  temp_dat <- Dat[[i]]
  smubm_res <- smubm(X = "X", Y = c("Y_p1","Y_p2","Y_np"), W = c("W_p1","W_p2","W_np"), S = "S_np", Z = "Z", data = temp_dat, N_mc = 2000, m_method = "effects_coding")
  Res[i,] <- smubm_res[[1]]
  rm(temp_dat)
}
colnames(Res) <- c("m_err", "r_err", "t_err", "m_err_hat", "r_err_hat", "rm_err_hat", "t_err_hat", "gamma_hat")
Res_ec <- cbind(param, Res)

Res <- as.data.frame(matrix(,nrow =nsim, ncol = 8))
pb <- progress_bar$new(total = nrow(param))
for (i in 1:nsim){
  pb$tick()
  temp_dat <- Dat[[i]]
  smubm_res <- smubm(X = "X", Y = c("Y_p1","Y_p2","Y_np"), W = c("W_p1","W_p2","W_np"), S = "S_np", Z = "Z", data = temp_dat, N_mc = 2000, m_method = "multigroup_1")
  Res[i,] <- smubm_res[[1]]
  rm(temp_dat)
}
colnames(Res) <- c("m_err", "r_err", "t_err", "m_err_hat", "r_err_hat", "rm_err_hat", "t_err_hat", "gamma_hat")
Res_mg1 <- cbind(param, Res)

Res <- as.data.frame(matrix(,nrow =nsim, ncol = 8))
pb <- progress_bar$new(total = nrow(param))
for (i in 1:nsim){
  pb$tick()
  temp_dat <- Dat[[i]]
  smubm_res <- smubm(X = "X", Y = c("Y_p1","Y_p2","Y_np"), W = c("W_p1","W_p2","W_np"), S = "S_np", Z = "Z", data = temp_dat, N_mc = 2000, m_method = "multigroup_2")
  Res[i,] <- smubm_res[[1]]
  rm(temp_dat)
}
colnames(Res) <- c("m_err", "r_err", "t_err", "m_err_hat", "r_err_hat", "rm_err_hat", "t_err_hat", "gamma_hat")
Res_mg2 <- cbind(param, Res)

# Exploring the results 



# Calculate Differences in Errors
Res <- Res_ec %>%
  mutate(
    #dif_phi = abs(phi_g - phi),
    TE = abs(t_err_hat - t_err),
    ME = abs(m_err_hat - m_err),
    RE = abs(r_err_hat - r_err),
    mRE = abs(rm_err_hat - r_err)
  )

# Filter Data for phi_g = 0.5 and Reshape
Res_long <- Res %>%
  #filter(phi_g == 0.5) %>%
  pivot_longer(cols = c(RE, ME, TE), 
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
    x = expression(paste(rho[yz], " | ", rho[ys])), 
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

# 
# Res %>% 
#   filter(phi_g == 0.5) %>% 
# ggplot(aes(x = t_err, y = t_err_hat, color = as.factor(phi))) + 
#   geom_line(stat = "summary", fun = "mean", size = 1) + # Lines showing the mean error for each phi_g
#   geom_abline(intercept = 0, slope = 1, size = 0.5, linetype = "dashed") +
#   ggh4x::facet_grid2(rho_ys ~ sys_e, labeller = label_bquote(rows = phi[T]:.(rho_ys), col = epsilon:.(sys_e)))  
# 
# 
# ggplot(Res, aes(x = r_err, y = r_err_hat, color = as.factor(sys_e))) + 
#   geom_point() +
#   facet_wrap(~rho_ys)  
# 
# Res %>% 
#   filter(rho_yz==0.3 & rho_ys == 0.3 & phi == 0.4 & D==1 & phi_g == 0.3) %>% 
#   View()
# 
# 
#     
# Res$dif_m_err <- abs(Res$m_err - Res$m_err_hat)
# plot(Res$m_err_hat, Res$m_err_hat)
# 
# agg <- aggregate(list("t_err"=Res$t_err, "t_err_hat"=Res$_err),
#                  by = list("rho_ys"=Res$rho_ys, "rho_yz"=Res$rho_yz, "phi"=Res$phi, "phi_g"=Res$phi_g, "sys_e"=Res$sys_e), mean)
# 
# mad <- aggregate(list(dif_t_err, dif_m_err) ~ rho_ys + rho_yz + phi + phi_g + rho_zs, data = Res, FUN = mean, na.rm = T)
# 
# Res$dif_phi <- abs(Res$phi_g - Res$phi)
# Res$dif_t_err <- abs(Res$t_err_hat - Res$t_err)
# Res$dif_m_err <- abs(Res$m_err_hat - Res$m_err)
# Res$dif_r_err <- abs(Res$r_err_hat - Res$r_err)
# Res$dif_rm_err <- abs(Res$rm_err_hat - Res$r_err)
# 
# ggplot(Res, aes(x = dif_rm_err, y = as.factor(phi_g), fill = as.factor(phi_g))) +
#   geom_boxplot(outlier.shape = NA) +
#   ggh4x::facet_grid2(rho_yz ~ rho_ys, labeller = label_bquote(rows = rho[yz]:.(rho_yz), col = rho[ys]:.(rho_ys)))  
# 
# 
# ggplot(Res, aes(x = dif_rm_err, y = as.factor(dif_phi), fill = as.factor(dif_phi))) + 
#   geom_boxplot(outlier.shape = NA) +
#   ggh4x::facet_grid2(rho_ys ~ sys_e, labeller = label_bquote(rows = rho[ys]:.(rho_ys), col = epsilon:.(sys_e)))  
# 
# 
# 
# # What happens if we include an audit sample? - With vs. Without
#   
# 
# Res %>% 
#   filter(rho_yz==0.5 & rho_ys == 0.7 & phi == 0.6 & phi_g == 0.5) %>% 
#   ggplot(aes(x = dif_m_err, y = factor(sys_e), fill = factor(sys_e))) + 
#   geom_violin() 

