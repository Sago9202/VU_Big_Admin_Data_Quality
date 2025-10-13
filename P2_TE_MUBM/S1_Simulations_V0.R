## ------------------------------------------- ##
## Modeling Total Error using Multisource Data ##
## Script 1 - Simulation                       ##
## Written by: Santiago Gómez-Echeverry        ##
## Last modified: 06/12/2024                   ##
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

#### - (II) Data generating process - ####

set.seed(42)
N <- 2e3                                       # Population size
D <- 200                                       # Number of draws
f_p <- 0.25                                    # Sample fraction, probability sample
f_a <- 0.1                                     # Sample fraction, Audit sample
n_p <- N*f_p                                   # Probability sample size
n_a <- N*f_a                                   # Audit sample size
phi_g <- rho_ys <- rho_yz <- c(0.3, 0.5, 0.7)  # Possible correlations
phi <- c(0.4, 0.6, 0.8)                        # Possible phi's
rnd_e <- 0.1                                   # Random error
sys_e <- c(0.2, 0.6, 0.8)
param <- expand.grid("rho_xz" = rho_yz, "rho_xs" = rho_ys, "phi" = phi, D = 1:D, "phi_g" = phi_g, "sys_e"=sys_e)
# We need to define the correlation between zs based on the other parameters that are already set.
param$rho_zs <- (param$rho_xs*(1-param$phi+param$phi*param$rho_xz))/(param$phi - param$phi*param$rho_xz + param$rho_xz)
nsim <- nrow(param)

Res <- as.data.frame(matrix(,nrow =nsim, ncol = 7))

Dat <- vector(mode = "list", length = nrow(param))
pb <- progress_bar$new(total = nrow(param))
for (i in 1:nsim){
  pb$tick()
  Cop <- normalCopula(param = c(param[i,1], param[i,2], param[i,7]), dim = 3, dispstr = "un")
  dat <- rCopula(N, Cop)
  dat[,1] <- qnorm(dat[,1], mean = 0, sd = 1)
  dat[,2] <- qnorm(dat[,2], mean = 0, sd = 1)
  dat[,3] <- qbinom(dat[,3], size = 1, prob = 0.6)
  dat <- data.frame(unlist(dat))
  dat$id <- seq.int(nrow(dat))
  dat <- dat[,c(4,1,2,3)]
  colnames(dat) <- c('id', 'X', 'Z', 'S_np')
  dat$W <- rnorm(N, mean = 0, sd = 1)
  dat$Y <- (1 - param[i, "sys_e"]) * dat$X +   param[i, "sys_e"] * (dat$W + rnorm(N, mean = 0, sd = rnd_e))
  dat$Y1 <- (1 - param[i, "sys_e"]*0.8) * dat$X +   param[i, "sys_e"]*0.8 * dat$W
  dat$Y2 <- (1 - param[i, "sys_e"]*0.6) * dat$X +   param[i, "sys_e"]*0.6 * dat$W
  dat$pr_p1 <- runif(N)
  dat <- dat[order(dat$pr_p1),]
  dat$S_p1 <- 0
  dat$S_p1[1:n_p] <- 1
  dat$pr_p2 <- runif(N)
  dat <- dat[order(dat$pr_p2),]
  dat$S_p2 <- 0
  dat$S_p2[1:n_p] <- 1
  dat$pr_a <- runif(N)
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

smubm <- function(X, Y, W, S, Z, data, phi_guess){
  X <- data[, X]
  Y_p1 <- data[,Y[[1]]]
  Y_p2 <- data[,Y[[2]]]
  Y_np <- data[,Y[[3]]]
  W_p1 <- data[,W[[1]]]
  W_p2 <- data[,W[[2]]]
  W_np <- data[,W[[3]]]
  S <- data[, S]
  Z <- data[, Z]
  Z_np <- Z[S==1]
  # True values
  m_err <- abs(mean(Y_np, na.rm = T) - mean(X[S==1], na.rm = T))
  r_err <- abs(mean(X[S==1], na.rm = T) - mean(X))
  t_err <- m_err + r_err
  
  # (1) Measurement Model
  mod <-  'Yl =~ Y_p1 + Y_p2 + Y_np
           Wl =~ W_p1 + W_p2 + W_np
           Yl ~~ Wl'
  fit <- sem(model = mod, data = dat, effect.coding = T, meanstructure = T)
  sum_fit <- summary(fit) 
  sum_fit$pe
  
  # We need to get an estimate  of the Yl
  mx_hat <- mean(lavPredict(fit)[[2]][,"Yl"])
  s2x_hat <- var(lavPredict(fit)[[1]][,"Yl"])
  my_hat <- lavPredictY(fit, ynames = "Y_np", xnames = c("Y_p1", "Y_p2","W_p1","W_p2", "W_np")) %>% 
    mean()
  s2y_hat <- lavPredictY(fit, ynames = "Y_np", xnames = c("Y_p1", "Y_p2","W_p1","W_p2", "W_np")) %>% 
    var() %>% 
    as.numeric()
  rel <- reliability(fit)
  #rel <- semTools::compRelSEM(fit, tau.eq=F, obs.var=T)
  lambda_hat <- rel[rownames(rel)=="omega", colnames(rel)=="Yl"]
  x_hat <- mx_hat + lambda_hat*(s2x_hat/s2y_hat)*(Y_np - my_hat)
  m_err_hat <- abs(mean(Y_np, na.rm = T) - mean(x_hat[S==1], na.rm = T))
  
  # (2) Representation model
  r_yz <- cor(Y_np,Z, use = "complete.obs")
  g_hat <- (phi_guess + (1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess))
  g_hat_m <- (phi_guess + lambda_hat*(1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess))
  c_obs <- sd(Y_np, na.rm = T)/sd(Z)
  r_err_hat <- g_hat*c_obs*(mean(Z[S == 1], na.rm = T) - mean(Z))
  rm_err_hat <- g_hat_m*c_obs*(mean(Z[S == 1], na.rm = T) - mean(Z))
  t_err_hat <- m_err_hat + rm_err_hat
  
  main_res <- c(m_err, r_err, t_err, m_err_hat, r_err_hat, rm_err_hat, t_err_hat)
  add_res <- list(fit, x_hat, lambda_hat)
  res <- list(main_res, add_res)
  return(res)
}

#### - (IV) Estimation - ####

pb <- progress_bar$new(total = nrow(param))
for (i in 1:nsim){
  pb$tick()
  temp_dat <- Dat[[i]]
  g <- param[i,"phi_g"]
  smubm_res <- smubm(X = "X", Y = c("Y_p1","Y_p2","Y_np"), W = c("W_p1","W_p2","W_np"), S = "S_np", Z = "Z", data = temp_dat, phi_guess = g)
  Res[i,] <- smubm_res[[1]]
  rm(temp_dat, g)
}
colnames(Res) <- c("m_err", "r_err", "t_err", "m_err_hat", "r_err_hat", "rm_err_hat", "t_err_hat")
Res <- cbind(param, Res)

# Calculate Differences in Errors
Res <- Res %>%
  mutate(
    dif_phi = abs(phi_g - phi),
    TE = abs(t_err_hat - t_err),
    ME = abs(m_err_hat - m_err),
    RE = abs(r_err_hat - r_err),
    mRE = abs(rm_err_hat - r_err)
  )

# Filter Data for phi_g = 0.5 and Reshape
Res_long <- Res %>%
  filter(phi_g == 0.5) %>%
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
ggplot(Res_long, aes(x = interaction(rho_yz, rho_ys, sep = "|"), 
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

