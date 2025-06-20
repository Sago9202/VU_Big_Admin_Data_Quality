##-----------------------------------------##
## Simulation on Selection Bias Estimators ##
## Written by: Santiago Gómez-Echeverry    ##
## Last update: 19/05/2025                 ##
##-----------------------------------------##

#### - (I) Working space and packages - ####

rm(list = ls()) # - Remove all the objects in the environment
sys <- Sys.info()
options(digits = 6)
options(scipe = 10)
# We will assign each one of the different folders into an object so that we can access them easily
fold_code <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/1 - Selection Bias Simulations/2 - Code")
fold_data <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/1 - Selection Bias Simulations/3 - Data")
fold_graphs <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/1 - Selection Bias Simulations/4 - Graphs & Tables")
setwd(fold_code) # For the moment, let us work on the code folder
# getwd() # - To check that we are in the right folder

# Below are all the packages that we will use in the analyses. We will check if they are installed, install them if they
# are not, and finally load them. Note: dutchmasters is installed from github
packages <- c('ggplot2', 'ggpubr', 'forecast', 'reshape2', 'MASS', 'corpcor', 'ggridges', 'ltm', 'stringr', 'knitr', 'kableExtra', 'tidyr', 'dplyr',
              'robustbase', 'forcats', 'cols4all', 'scales','matrixStats', 'progress', 'ggh4x', 'VGAM')

# Installing the packages
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)){install.packages(packages[!installed_packages])}

# Loading the packages 
invisible(lapply(packages, library, character.only = TRUE))

# Besides setting the working space, I'll define in here some simple functions that are required across the code.
p.cov <- function(x,y){mean((x-mean(x))*(y-mean(y)))}
p.sd <- function(x){sqrt(mean((x-mean(x))^2))}
p.var <- function(x){p.sd(x)^2}
cond.var <- function(y,z){(1 - cor(y,z)^{2})*p.var(y)}
cond.cov <- function(y,x,z){p.cov(x,y) - (p.cov(y,z)*p.cov(x,z))/p.var(z)}
cond.cor <- function(y,x,z){cond.cov(y,x,z)/sqrt(cond.var(y,z)*cond.var(x,z))}

colMed <- function(array){
  mat_array <- matrix(NA, nrow = dim(array)[2], ncol = dim(array)[3])
  for (i in 1:dim(array)[2]){
    mat_array[i,] <- colMedians(as.matrix(array[,i,]), na.rm = TRUE) 
  } 
  return(mat_array)
}

nice_palette <- c4a("carto.pastel", n = 11)
nice_palette <- c(nice_palette[1:6], nice_palette[11])
#nice_palette <- c4a(palette = "brewer.paired", n = 11)

compute_phi <- function(rho_XS, rho_XZ, rho_ZS) {
  num <- rho_XS - rho_XZ * rho_ZS
  denom <- (rho_XS - rho_XZ * rho_ZS) + (rho_ZS - rho_XZ * rho_XS)
  return(num / denom)
}

safe_cor <- function(x, y){
  if (sum(!is.na(x) & !is.na(y)) >= 2) {
    return(cor(x, y, use = "pairwise.complete.obs"))
  } else {
    return(NA)
  }
}

#### - (II) Selection Bias Estimators - ####

# 1. DDP_L

DDP_L <- function(Y, S, Z){
  P <- ncol(Y)
  N <- apply(Y, 2, length)
  n <- apply(S, 2, sum)  
  ys <- ifelse(S == 1, Y, NA)
  r_ys <- mapply(cor, as.data.frame(Y), as.data.frame(S))
  sigma_y <- apply(Y, 2, p.sd)
  Yl <- (cbind(rep(NA, times = nrow(Y)), Y[,-P]))
  Sl <- (cbind(rep(NA, times = nrow(S)), S[,-P]))
  ysl <- ifelse(S == 1, Yl, NA)
  r_ys_hat <- mapply(cor, as.data.frame(Yl), as.data.frame(Sl))
  sigma_y_hat <- apply(Yl, 2, p.sd)*(colSds(ys, na.rm = T)/colSds(ysl, na.rm = T))
  # Threefold decomposition
  D_I <- r_ys_hat
  D_U <- sigma_y_hat
  D_O <- sqrt((N-n)/n)
  S_Meng <- D_I*D_U*D_O
  res <- list(S_Meng = S_Meng, D_I = D_I, D_U = D_U, D_O = D_O, r_ys = r_ys)
  return(res)
}

# 2. MUB_C.5

MUB_C.5 <- function(Y, S, Z, M = 200){
  P <- ncol(Y)
  ys <- ifelse(S==1, Y, NA)
  zs <- ifelse(S==1, Z, NA)
  #True values
  dy <- (colMeans(ys, na.rm = T) - colMeans(Y))/colSds(Y)
  dz <- (colMeans(zs, na.rm = T) - colMeans(Z))/colSds(Z)
  r_yz <- mapply(safe_cor, as.data.frame(ys), as.data.frame(zs)) # We need the safecor because of the NA's
  r_zs <- mapply(cor, as.data.frame(Z), as.data.frame(S))
  phi <- (dy - r_yz*dz)/((dy - r_yz*dz) + (dz - r_yz*dy))
  phi_hat <- 0.5
  g_hat <- (phi_hat + (1-phi_hat)*r_yz)/(phi_hat*r_yz + (1-phi_hat))
  c_obs <- colSds(ys, na.rm = T)/colSds(Z)
  S_MUB <- g_hat*c_obs*(colMeans(zs, na.rm = T) - colMeans(Z))
  res <- list(S_MUB = S_MUB, phi_hat = phi_hat, phi = phi)
  return(res)
}

# 3. MUB_C

MUB_C <- function(Y, S, Z, M = 200){
  P <- ncol(Y)
  ys <- ifelse(S==1, Y, NA)
  zs <- ifelse(S==1, Z, NA)
  #True values
  dy <- (colMeans(ys, na.rm = T) - colMeans(Y))/colSds(Y)
  dz <- (colMeans(zs, na.rm = T) - colMeans(Z))/colSds(Z)
  r_yz <- mapply(safe_cor, as.data.frame(ys), as.data.frame(zs)) # We need the safecor because of the NA's
  r_zs <- mapply(cor, as.data.frame(Z), as.data.frame(S))
  phi <- (dy - r_yz*dz)/((dy - r_yz*dz) + (dz - r_yz*dy))
  phi_hat <- matrix(runif(M*P), ncol = P, nrow = M)
  g_hat <- (phi_hat + (1-phi_hat)*r_yz)/(phi_hat*r_yz + (1-phi_hat))
  c_obs <- colSds(ys, na.rm = T)/colSds(Z)
  MUB <- g_hat*c_obs*(colMeans(zs, na.rm = T) - colMeans(Z))
  S_MUB <- colMeans(MUB)
  r_ys_hat <- apply(g_hat*r_zs,2,range)
  res <- list(S_MUB = S_MUB, phi_hat = phi_hat, phi = phi, r_ys_hat = r_ys_hat)
  return(res)
}

# 3. MUB_M

MUB_M <- function(Y, S, Z){
  P <- ncol(Y)
  ys <- ifelse(S == 1, Y, NA)
  zs <- ifelse(S == 1, Z, NA)
  Yl <- (cbind(rep(NA, times = nrow(Y)), Y[,-P]))
  Zl <- (cbind(rep(NA, times = nrow(Z)), Z[,-P]))
  Sl <- (cbind(rep(NA, times = nrow(S)), S[,-P]))
  ysl <- ifelse(S == 1, Yl, NA)
  dy <- (colMeans(ys, na.rm = T) - colMeans(Y))/colSds(Y)
  dz <- (colMeans(zs, na.rm = T) - colMeans(Z))/colSds(Z)
  r_yz_hat <- mapply(cor, as.data.frame(Yl), as.data.frame(Zl))
  phi <- (dy - r_yz_hat*dz)/((dy - r_yz_hat*dz) + (dz - r_yz_hat*dy))
  phi_hat <- c(NA, phi[-P])
  g_hat <- (phi_hat + (1-phi_hat)*r_yz_hat)/(phi_hat*r_yz_hat + (1-phi_hat))
  sigma_y_hat <- apply(Yl, 2, p.sd)*(colSds(ys, na.rm = T)/colSds(ysl, na.rm = T))
  c_obs <- sigma_y_hat/colSds(Z)
  S_MUB <- g_hat*c_obs*(colMeans(zs, na.rm = T) - colMeans(Z))
  res <- list(S_MUB = S_MUB, phi_hat = phi_hat, phi = phi)
  return(res)
} 

# 4. Cov_L

Cov_L <- function(Y, S, Z){
  P <- ncol(Y)
  Yl <- (cbind(rep(NA, times = nrow(Y)), Y[,-P]))
  Sl <- (cbind(rep(NA, times = nrow(S)), S[,-P]))
  ysl <- ifelse(Sl == 1, Yl, NA)
  S_Cov <- (colMeans(ysl, na.rm = T) - colMeans(Yl))
  res <- S_Cov
  return(res)
}

# 5. Cov_C

Cov_C <- function(Y,S,Z){
  P <- ncol(Y)
  zs <- ifelse(S == 1, Z, NA)
  S_Cov <- (colMeans(zs, na.rm = T) - colMeans(Z))
  res <- S_Cov
  return(res)
}

# 6. Cov_M

Cov_M <- function(Y,S,Z){
  P <- ncol(Y)
  ys <- ifelse(S == 1, Y, NA)
  zs <- ifelse(S == 1, Z, NA)
  Yl <- (cbind(rep(NA, times = nrow(Y)), Y[,-P]))
  Sl <- (cbind(rep(NA, times = nrow(S)), S[,-P]))
  Zl <- (cbind(rep(NA, times = nrow(Z)), Z[,-P]))
  sigma_ys <- mapply(cov, as.data.frame(Y), as.data.frame(S))
  sigma_ys_hat <- mapply(cov, as.data.frame(Yl), as.data.frame(Sl))
  sigma_zs <- mapply(cov, as.data.frame(Zl), as.data.frame(Sl))
  S_Cov <- (sigma_ys_hat/sigma_zs)*(colMeans(zs, na.rm = T) - colMeans(Z))
  res <- S_Cov
  return(res)
}

#### - (III) Performance functions - ####

MAD <- function(Y, S, M){
  P <- ncol(Y)
  ys <- ifelse(S == 1, Y, NA)
  RS <- 2*abs(colMeans(ys, na.rm = T) - colMeans(Y))/(abs(colMeans(ys, na.rm = T))+abs(colMeans(Y)))
  RS_hat <- 2*abs(M)/(abs(colMeans(ys, na.rm = T))+abs(colMeans(ys, na.rm = T) - M))
  MADt <- mean(abs(RS_hat - RS), na.rm = T)
  res <- list("RS" = RS, "RS_hat" = RS_hat, "MAD" = MADt)
  return(res)
}

RMSD <- function(Y, S, M){
  P <- ncol(Y)
  ys <- ifelse(S == 1, Y, NA)
  RS <- 2*abs(colMeans(ys, na.rm = T) - colMeans(Y))/(abs(colMeans(ys, na.rm = T))+abs(colMeans(Y)))
  RS_hat <- 2*abs(M)/(abs(colMeans(ys, na.rm = T))+abs(colMeans(ys, na.rm = T) - M))
  RMSDt <- sqrt(mean((RS_hat - RS)^2, na.rm = T))
  res <- list("RS" = RS, "RS_hat" = RS_hat, "RMSD" = RMSDt)
  return(res)
}

#### - (IV) Data generation - ####

# (i) Initial settings

set.seed(42)

mu <- 0
sigma <- 1
N <- 1000
Time <- 12
D <- 300
const <- rnorm(n = 1, mean = 0, sd = 1)
epsilon_0 <- matrix(data = rnorm(n = N*D, mean = 0, sd = 1), nrow = N, ncol = D)

# (ii) Outcome variable

# 1. Normal distribution
Y <- array(data = NA, dim = c(N, Time, D))
e <- array(data = rnorm(n = N * Time * D, mean = 0, sd = sigma), dim = c(N, Time, D))
Y[, 1, ] <- rnorm(n = N * D, mean = 0, sd = sigma)  # Random initialization
psi <- runif(n = D, min = 0, max = 1)
for (p in 2:Time){
  Y[,p,] <- const + psi*Y[,p-1,] + e[,p,]
}
dimnames(Y)[[2]] <- paste0("Y", c(1:Time))
mean(Y[,1,1])
mean(Y[,6,1])

## Check the cor(Y^{t}, Y^{t-1})

Cors <- numeric(dim(Y)[3])
for (d in 1:length(Cors)){
  cor_matrix <- cor(Y[,,d],Y[,,d])
  cor_t <- diag(cor_matrix[-1, -ncol(cor_matrix)])
  Cors[d] <- mean(cor_t, na.rm = T)
}

# Print the result
cat("The average cor(Y_t, Y_t-1) is:", mean(Cors), "\n") 

# 2. Beta distribution
ab <- expand.grid(1:4, 1:4)
colnames(ab) <- c("alpha", "beta")
n_betas <- nrow(ab)

# Store results
bt <- vector(mode = "list", length = n_betas)
bt0 <- array(data = NA, dim = c(N, n_betas, D))

# Generate the Beta AR(1) Process
for (b in 1:n_betas) {
  bt[[b]] <- array(data = NA, dim = c(N, Time, D))
  
  # Initialize first period from Beta(alpha, beta)
  bt[[b]][,1,] <- bt0[,b,] <- matrix(rbeta(N * D, shape1 = ab[b,1], shape2 = ab[b,2]), nrow = N, ncol = D)
  
  # Convert first period to logit scale
  Z_t <- qlogis(bt[[b]][,1,])  # Logit transformation
  
  for (t in 2:Time) {
    # AR(1) process in logit space
    epsilon_t <- rnorm(N, mean = 0, sd = 0.1)  # Small noise in latent space
    Z_t <- psi * Z_t + epsilon_t  # AR(1) update in logit scale
    
    # Transform back to Beta space
    bt[[b]][,t,] <- plogis(Z_t)  # Inverse logit transformation
  }
}

# Compute empirical correlation cor(X_t, X_{t-1})
Cors <- numeric(D)
for (d in 1:D) {
  X_draw <- bt[[1]][,,d]  # Extract one series
  cor_values <- apply(X_draw, 1, function(b) cor(b[-1], b[-length(b)]))
  Cors[d] <- mean(cor_values, na.rm = TRUE)
}

# Print average empirical correlation
cat("The average empirical cor(X_t, X_t-1) is:", mean(Cors), "\n")

# Let's plot the variables at the first period and draw and get the skewness and kurtosis

ab$skew <- round((2 * (ab$beta - ab$alpha) * sqrt(ab$alpha + ab$beta + 1)) /
                   ((ab$alpha * ab$beta + 2) * sqrt(ab$alpha * ab$beta)), 2)

ab$kurt <- round((6 * ((ab$alpha - ab$beta)^2 * (ab$alpha + ab$beta + 1) -
                         ab$alpha * ab$beta * (ab$alpha + ab$beta + 2))) /
                   (ab$alpha * ab$beta * (ab$alpha + ab$beta + 2) * (ab$alpha + ab$beta + 3)), 2)

ab$label <- paste0("α = ", ab$alpha, 
                   ", β = ", ab$beta, 
                   "\nSkew = ", ab$skew, 
                   ", Kurt = ", ab$kurt)

# Prepare data for plotting
bt0_d1 <- as.data.frame(bt0[,,1])
colnames(bt0_d1) <- paste0("a", ab$alpha, "_b", ab$beta)

m_bt0_d1 <- melt(bt0_d1, variable.name = "combo")
m_bt0_d1 <- separate(m_bt0_d1, combo, into = c("a", "b"), sep = "_")
m_bt0_d1$a <- as.integer(gsub("a", "", m_bt0_d1$a))
m_bt0_d1$b <- as.integer(gsub("b", "", m_bt0_d1$b))
m_bt0_d1 <- left_join(m_bt0_d1, ab, by = c("a" = "alpha", "b" = "beta"))

# Plot 1: All 16 combinations
p1 <- ggplot(m_bt0_d1, aes(value)) + 
  geom_histogram(fill = "#66c2a5", alpha = 0.8, color = "black", size = 0.3, bins = 30) + 
  facet_wrap(~ label, scales = "free_y") +
  ylab("Density") + 
  xlab("Value") + 
  theme_minimal(base_size = 20) +
  theme(
    strip.text = element_text(size = 14),
    strip.background = element_blank()
  )

ggsave("Betas_1.png", p1, path = fold_graphs, width = 30, height = 30, units = "cm")

# Plot 2: Only alpha = 1, beta = 1:4
subset_ab <- filter(ab, alpha == 1)
subset_bt0 <- bt0[, ab$alpha == 1, 1]
bt0_d1_subset <- as.data.frame(subset_bt0)
colnames(bt0_d1_subset) <- paste0("alpha1_b", subset_ab$beta)

m_bt0_d1_subset <- melt(bt0_d1_subset, variable.name = "combo")
m_bt0_d1_subset$b <- as.integer(gsub("alpha1_b", "", m_bt0_d1_subset$combo))
m_bt0_d1_subset <- left_join(m_bt0_d1_subset, subset_ab, by = c("b" = "beta"))

p2 <- ggplot(m_bt0_d1_subset, aes(value, fill = label)) + 
  geom_histogram(alpha = 0.6, color = "black", size = 0.3, bins = 30) + 
  facet_wrap(~ label, scales = "free_y") +
  scale_fill_manual(values = nice_palette, guide = "none") + 
  ylab("Density") + 
  xlab("Value") + 
  theme_minimal(base_size = 20) +
  theme(
    strip.text = element_text(size = 14),
    strip.background = element_blank()
  )

ggsave("Betas_2.png", p2, path = fold_graphs, width = 30, height = 30, units = "cm")

# (iii) Auxiliary variable

corZ <- function(Y, rho){
  set.seed(1234)
  d1 <- dim(Y)[[1]]; d2 <- dim(Y)[[2]];d3 <- dim(Y)[[3]]
  X <- array(data = NA, dim = dim(Y))
  for (d in 1:d3){
   for (p in 1:d2){
    Y_temp <- Y[,p,d]
    X_temp <- rnorm(n = d1, mean = 0, sd = 1)
    Y.est <- residuals(lm(X_temp ~ Y_temp))
    X[,p,d] <- rho*sd(Y.est)*Y_temp + Y.est*sd(Y_temp)*sqrt(1-rho^2)
   }
  }
  return(X)
}

cor_aux <- c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95)
n_aux <- length(cor_aux)

# 1. Normal distribution
Z <- vector(mode = "list", length = n_aux)

for (cz in 1:n_aux){
  Z[[cz]] <- corZ(Y, cor_aux[cz])
  dimnames(Z[[cz]])[[2]] <- paste0("Z", c(1:Time))
}

# Let's check the correlations
corYZ <- matrix(data = NA, nrow = Time, ncol = n_aux) 

for (d in 1:D){
  for (p in 1:Time){
    for (cz in 1:n_aux){
      corYZ[p,cz]<- cor(Y[,p,d],Z[[cz]][,p,d], method = "pearson")
    }
  }
  print(colMeans(corYZ))
}

# 2. Beta distribution

Z_bt <- vector(mode = "list", length = n_betas)
for (b in 1:n_betas){
  Z_bt[[b]] <- vector(mode = "list", length = n_aux)
  for (cz in 1:n_aux){
    Z_bt[[b]][[cz]] <- corZ(bt[[b]],cor_aux[cz])
    dimnames(Z_bt[[b]][[cz]])[[2]] <- paste0("Z", c(1:Time))
  }
}
names(Z_bt) <- paste0("Z_b",1:n_betas) 

# Let's check the correlations with the Y (beta) variables
corBZ <- matrix(data = NA, nrow = Time, ncol = n_aux) 
for (b in 1:n_betas){
  cat("Correlations with Beta:",b)
  for (d in 1:D){
    for (p in 1:Time){
      for (cz in 1:n_aux){
        corBZ[p,cz]<- cor(bt[[b]][,p,d],Z_bt[[b]][[cz]][,p,d], method = "pearson")
      }
    }
    print(colMeans(corBZ))
  }
}

# (iv) Selection variable

corS <- function(Y, rho, sample_rate = 0.5, var_strength = 0.05, noise_sd = 0.5){
  set.seed(4321)
  X <- array(data = NA, dim = dim(Y))
  for (d in 1:D){
    for (p in 1:Time){
      Y_temp <- Y[,p,d]
      Y.est <- residuals(lm(rnorm(length(Y_temp)) ~ Y_temp))
      X_temp <- rho*sd(Y.est)*Y_temp + Y.est*sd(Y_temp)*sqrt(1-rho^2)
      X_temp <- (X_temp - min(X_temp))/(max(X_temp)-min(X_temp))
      
      # Add small variability to the threshold around sample_rate
      rate_jitter <- runif(1, 
                           min = max(sample_rate - var_strength, 0.01), 
                           max = min(sample_rate + var_strength, 0.99))
      threshold <- quantile(X_temp, probs = 1 - rate_jitter)
      X[,p,d] <- ifelse(X_temp > threshold, 1, 0)
    }
  }
  return(X)
}


cor_selt <- c(0.2, 0.4, 0.6, 0.8) 
cor_sel <- c(0.513, 0.649, 0.801, 0.999)
n_sel <- length(cor_sel)

# Correlated with the Y (normal) variable
S <- vector(mode = "list", length = n_sel)

for (cs in 1:n_sel){
  S[[cs]] <- corS(Y,cor_sel[cs])
  dimnames(S[[cs]])[[2]] <- paste0("S", c(1:Time))
}

# Let's check the correlations
corYS <- matrix(data = NA, nrow = Time, ncol = n_sel) 

for (d in 1:D){
  for (p in 1:Time){
    for (cs in 1:n_sel){
      corYS[p,cs]<- cor.test(Y[,p,d],S[[cs]][,p,d])$estimate
    }
  }
  print(colMeans(corYS))
}

# Correlated with the Y (beta) variables
S_bt <- vector(mode = "list", length = n_betas)
for (b in 1:n_betas){
  S_bt[[b]] <- vector(mode = "list", length = n_sel)
  for (cs in 1:n_sel){
    S_bt[[b]][[cs]] <- corS(bt[[b]],cor_sel[cs])
    dimnames(S_bt[[b]][[cs]])[[2]] <- paste0("S", c(1:Time))
  }
}
names(S_bt) <- paste0("S_b",1:n_betas) # If we name the list within the list it will be 'easier' to access them

# Let's check the correlations with the Y (beta) variables
corBS <- matrix(data = NA, nrow = Time, ncol = n_sel) 
for (b in 1:n_betas){
  cat("Correlations S with Beta:", b)
  for (d in 1:D){
    for (p in 1:Time){
      for (cs in 1:n_sel){
        corBS[p,cs]<- cor.test(bt[[b]][,p,d],S_bt[[b]][[cs]][,p,d])$estimate
      }
    }
    print(colMeans(corBS))
  }
}


#### - (V) Estimations - ####

Est <- c("DDP_L", "MUB_C.5", "MUB_C","MUB_M", "Cov_L", "Cov_C", "Cov_M")
mres <- c("MAD","RMSD")
perf_y <- expand.grid("S" = 1:4, "Z" = 1:6, "M"= Est, "Y"= "Y", "beta" = NA, "a" = NA, "b" = NA, "res" = mres, "d" = 1:D)
perf_b <- expand.grid("S" = 1:4, "Z" = 1:6, "M"= Est, "Y"= "Beta","beta" = 1:16, "res" = c("MAD", "RMSD"), "d" = 1:D)
perf_b <- perf_b[order(perf_b$beta),]
perf_b$a <- rep(rep(1:4, each = nrow(perf_b)/n_betas), times = 4)
perf_b$b <- rep(rep(1:4, each = nrow(perf_b)/n_betas), each = 4)
perf <- perf_tr <- rsh <- rsh_tr <-  rbind(perf_y, perf_b)
rm(perf_y, perf_b)
n_est <- nrow(perf)

v <- rsv <-  v_tr <- rsv_tr <- vector(mode = "list", length = n_est)

# Estimations
pb <- progress_bar$new(total = n_est)
for (i in 1:n_est){
  pb$tick()
  cs <- perf[i,"S"]; cz <- perf[i, "Z"]; m <- as.character(perf[i, "M"]); d <- perf[i, "d"]
  res <- as.character(perf[i, "res"]); y <- as.character(perf[i, "Y"]); b <- perf[i, "beta"]
  if (grepl("Beta", y)){
    Yt <- bt[[b]][,,d]
    St <- S_bt[[b]][[cs]][,,d]
    Zt <- Z_bt[[b]][[cz]][,,d]
  } else if (!grepl("Beta", y)){
    Yt <- Y[,,d]
    St <- S[[cs]][,,d]
    Zt <- Z[[cz]][,,d]
  }
  yt <- Yt*St; yt[yt==0] <- NA
  if (grepl("Cov", m)){
    Et <- get(m)(Yt, St, Zt)
  } else if (!grepl("Cov", m)){
    Et <- get(m)(Yt, St, Zt)[[1]]
  }

  v[[i]] <- get(res)(Yt, St, Et)[[3]]
  rsv[[i]] <- get(res)(Yt, St, Et)[[2]]
  ck <- abs(Et - colMeans(yt, na.rm = TRUE))<0.3
  if (any(ck == T, na.rm = T)){
    v_tr[[i]] <- NA; rsv_tr[[i]] <- rep(NA, Time)
  } else{
    v_tr[[i]] <- v[[i]]; rsv_tr[[i]] <- rsv[[i]]
  }
}
perf <- cbind(perf, "value"= as.vector(unlist(v)))
rsh <- cbind(rsh, matrix(unlist(rsv), ncol = Time, byrow = T))
names(rsh)[10:21] <- paste0("RS", 1:Time)
perf_tr <- cbind(perf_tr, "value" = as.vector(unlist(v_tr)))
rsh_tr <- cbind(rsh_tr, matrix(unlist(rsv_tr), ncol = Time, byrow = T))
names(rsh_tr)[10:21] <- paste0("RS", 1:Time)
rm(v, rsv, v_tr, rsv_tr)

Est2 <- c("DDP(L)", "MUB(C.5)", "MUB(C)","MUB(M)", "Cov(L)", "Cov(C)", "Cov(M)")
perf$M <- gsub("_(.*)$", "(\\1)", perf$M)
perf <- perf %>% 
  mutate(S = fct_recode(factor(S), "0.2"="1", "0.4"="2", "0.6"="3", "0.8"="4"), M = factor(M, levels = Est2),
         Z = fct_recode(factor(Z), "0.05"="1", "0.2"="2", "0.4"="3", "0.6"="4", "0.8"="5", "0.95"="6"))
write.table(perf, file = paste0(fold_data, "/Simulated_Data_N", N, "_D", D, "_T", Time, ".txt"))

rsh$M <- gsub("_(.*)$", "(\\1)", rsh$M)
rsh <- rsh %>% 
  mutate(S = fct_recode(factor(S), "0.2"="1", "0.4"="2", "0.6"="3", "0.8"="4"), M = factor(M, levels = Est2),
         Z = fct_recode(factor(Z), "0.05"="1", "0.2"="2", "0.4"="3", "0.6"="4", "0.8"="5", "0.95"="6"))

perf_tr$M <- gsub("_(.*)$", "(\\1)", perf_tr$M)
perf_tr <- perf_tr %>% 
  mutate(S = fct_recode(factor(S), "0.2"="1", "0.4"="2", "0.6"="3", "0.8"="4"), M = factor(M, levels = Est2),
         Z = fct_recode(factor(Z), "0.05"="1", "0.2"="2", "0.4"="3", "0.6"="4", "0.8"="5", "0.95"="6"))
write.table(perf_tr, file = paste0(fold_data, "/Simulated_Data_N", N, "_D", D, "_T", Time, "_tr.txt"))

rsh_tr$M <- gsub("_(.*)$", "(\\1)", rsh_tr$M)
rsh_tr <- rsh_tr %>% 
  mutate(S = fct_recode(factor(S), "0.2"="1", "0.4"="2", "0.6"="3", "0.8"="4"), M = factor(M, levels = Est2),
         Z = fct_recode(factor(Z), "0.05"="1", "0.2"="2", "0.4"="3", "0.6"="4", "0.8"="5", "0.95"="6"))

# Actual RS
rs_y <- expand.grid("S" = 1:4, "Z" = 1:6,  "Y"= "Y", "beta" = NA, "a" = NA, "b" = NA, "d" = 1:D)
rs_b <- expand.grid("S" = 1:4, "Z" = 1:6, "Y"= "Beta","beta" = 1:16, "d" = 1:D)
rs_b <- rs_b[order(rs_b$beta),]
rs_b$a <- rep(rep(1:4, each = nrow(rs_b)/n_betas), times = 4)
rs_b$b <- rep(rep(1:4, each = nrow(rs_b)/n_betas), each = 4)
rs <- rs_tr <-  rbind(rs_y, rs_b)
rm(rs_y, rs_b)
n_rs <- nrow(rs)
rsv <- rsv_tr <- vector(mode = "list", length = n_rs)

pb <- progress_bar$new(total = n_rs)
for (i in 1:n_rs){
  pb$tick()
  cs <- rs[i,"S"]; cz <- rs[i, "Z"]; d <- rs[i, "d"]; y <- as.character(rs[i, "Y"]); b <- rs[i, "beta"]
  if (grepl("Beta", y)==T){
    Yt <- bt[[b]][,,d]
    St <- S_bt[[b]][[cs]][,,d]
    Zt <- Z_bt[[b]][[cz]][,,d]
  } else if (grepl("Beta", y)==F){
    Yt <- Y[,,d]
    St <- S[[cs]][,,d]
    Zt <- Z[[cz]][,,d]
  }
  yt <- Yt*St; yt[yt==0] <- NA
  Et <- DDP_L(Yt, St, Zt)[[1]]
  rsv[[i]] <- MAD(Yt, St, Et)[[1]]
  ck <- abs(Et - colMeans(yt, na.rm = TRUE))<0.3
  if (any(ck == T, na.rm = T)){
    rsv_tr[[i]] <- rep(NA, Time)
  } else {
    rsv_tr[[i]] <- rsv[[i]]
  }
}

rs <- cbind(rs, matrix(unlist(rsv), ncol = Time, byrow = T))
rs <- rs[rs$Z==1, -2]
names(rs)[7:18] <- paste0("RS", 1:Time)
rs_tr <- cbind(rs_tr, matrix(unlist(rsv_tr), ncol = Time, byrow = T))
rs_tr <- rs_tr[rs_tr$Z==1, -2]
names(rs_tr)[7:18] <- paste0("RS", 1:Time)

rs <- rs %>% 
  mutate(S = fct_recode(factor(S), "0.2"="1", "0.4"="2", "0.6"="3", "0.8"="4"))

rs_tr <- rs_tr %>% 
  mutate(S = fct_recode(factor(S), "0.2"="1", "0.4"="2", "0.6"="3", "0.8"="4"))

#### - (VI) Graphs - ####
setwd(fold_graphs)
perf_plot <- perf %>% 
  mutate(M = factor(M, levels = Est2))
perf_plot$S <- as.numeric(as.character(perf_plot$S));perf_plot$Z <- as.numeric(as.character(perf_plot$Z))
colnames(perf_plot) <- c('S', 'Z', 'Estimator', 'Y', 'beta', 'a', 'b', 'res', 'd', 'value')

for (i in 1:length(mres)){
  # (i) Cross S x Z
  p_sz <- perf_plot %>% 
    filter(Z!=0.05 & Z!=0.95 & Y=="Y" & grepl(mres[i], res)) %>% 
    ggplot(aes(x = value, y = fct_rev(Estimator), color = Estimator)) + 
    stat_pointinterval(position = position_dodge(width = 0.9), .width = c(0.5, 1), point_size = 3, show.legend = T) +
    scale_color_manual(values = nice_palette, name = "Estimator") +
    ggh4x::facet_grid2(Z~S, scales = "free", independent = "all", labeller = label_bquote(rows = rho[YZ]:.(Z), col = rho[YS]:.(S))) + ylab("") + xlab(mres[i]) +
    theme(text = element_text(size = 30), axis.text.y = element_blank(), axis.text.x = element_text(size = 15),
          axis.ticks.y = element_blank(), legend.key.size = unit(1.5, "cm"), panel.spacing = unit(2, "lines")) +
    guides(size = "none", override.aes = list(size = 4, shape = 2))
  ggsave(paste0("Plot_SZ_", mres[i], ".png"), plot = p_sz, path = fold_graphs, width = 50, height = 50, units = "cm")
 
  # (i) Cross S x Z - extreme values
  p_sz_e <- perf_plot %>% 
    filter(Z!=0.2 & Z!=0.8 & Y=="Y" & grepl(mres[i], res)) %>% 
    ggplot(aes(x = value, y = fct_rev(Estimator), color = Estimator)) + 
    stat_pointinterval(position = position_dodge(width = 0.9), .width = c(0.5, 1), point_size = 3, show.legend = T) + 
    scale_color_manual(values = nice_palette, name = "Estimator") +
    ggh4x::facet_grid2(Z~S, scales = "free", independent = "all", labeller = label_bquote(rows = rho[YZ]:.(Z), col = rho[YS]:.(S))) + ylab("") + xlab(mres[i]) +
    theme(text = element_text(size = 30), axis.text.y = element_blank(), axis.text.x = element_text(size = 15),
          axis.ticks.y = element_blank(), legend.key.size = unit(1.5, "cm"), panel.spacing = unit(2, "lines")) +
    guides(size = "none", override.aes = list(size = 4, shape = 2))
  ggsave(paste0("Plot_SZ_", mres[i], "_ext.png"), plot = p_sz_e, path = fold_graphs, width = 50, height = 50, units = "cm")

  # (iii) Different Y, ceteris paribus - To check how the estimators react to changes in Y distribution
  p_b <- perf_plot %>% 
    filter(Z==0.8 & S==0.6 & Y=="Beta" & grepl(mres[i], res)) %>% 
    ggplot(aes(x = value, y = fct_rev(Estimator), color = Estimator)) + 
    stat_pointinterval(position = position_dodge(width = 0.9), .width = c(0.5, 1), point_size = 3, show.legend = T) + 
    scale_color_manual(values = nice_palette, name = "Estimator") +
    ggh4x::facet_grid2(b~a, scales = "free", independent = "all", labeller = label_bquote(rows = beta:.(b), col = alpha:.(a))) + ylab("") + xlab(mres[i]) +
    theme(text = element_text(size = 30), axis.text.y = element_blank(), axis.text.x = element_text(size = 15),
        axis.ticks.y = element_blank(), legend.key.size = unit(1.5, "cm"), panel.spacing = unit(2, "lines")) +
    guides(size = "none", override.aes = list(size = 4, shape = 2))
  ggsave(paste0("Plot_B_", mres[i], ".png"), plot = p_b, path = fold_graphs, width = 50, height = 50, units = "cm")

  # (iv) Different Y, with different S
  p_bs <- perf_plot %>% 
    filter(Z==0.8 & !is.na(a) & grepl(mres[i], res) & b==1) %>% 
    ggplot(aes(x = value, y = fct_rev(Estimator), color = Estimator)) + 
    stat_pointinterval(position = position_dodge(width = 0.9), .width = c(0.5, 1), point_size = 3, show.legend = T) + 
    scale_color_manual(values = nice_palette, name = "Estimator") +    
    ggh4x::facet_grid2(S~a, scales = "free", independent = "all", labeller = label_bquote(rows = rho[YS]:.(S), col = alpha:.(a))) + ylab("") + xlab(mres[i]) +
    theme(text = element_text(size = 30), axis.text.y = element_blank(), axis.text.x = element_text(size = 15),
          axis.ticks.y = element_blank(), legend.key.size = unit(1.5, "cm"), panel.spacing = unit(2, "lines")) +
    guides(size = "none", override.aes = list(size = 4, shape = 2))
  ggsave(paste0("Plot_BS_", mres[i], ".png"), plot = p_bs, path = fold_graphs, width = 50, height = 50, units = "cm")

  # (v) Different Y, with different Z
  p_bz <- perf_plot %>% 
    filter(S==0.4 & Z!=0.05 & Z!=0.95 & !is.na(a) & grepl(mres[i], res) & b==1) %>% 
    ggplot(aes(x = value, y = fct_rev(Estimator), color = Estimator)) + 
    stat_pointinterval(position = position_dodge(width = 0.9), .width = c(0.5, 1), point_size = 3, show.legend = T) + 
    scale_color_manual(values = nice_palette, name = "Estimator") +    
    ggh4x::facet_grid2(Z~a, scales = "free", independent = "all", labeller = label_bquote(rows = rho[YZ]:.(Z), col = alpha:.(a))) + ylab("") + xlab(mres[i])  +
    theme(text = element_text(size = 30), axis.text.y = element_blank(), axis.text.x = element_text(size = 15),
          axis.ticks.y = element_blank(), legend.key.size = unit(1.5, "cm"), panel.spacing = unit(2, "lines")) +
    guides(size = "none", override.aes = list(size = 4, shape = 2))
  ggsave(paste0("Plot_BZ_", mres[i], ".png"), plot = p_bz, path = fold_graphs, width = 50, height = 50, units = "cm")
  
  # (vi) Cross S x Z - beta 1
  p_b1_sz <- perf_plot %>% 
    filter(Z!=0.05 & Z!=0.95 & !is.na(a) & grepl(mres[i], res) & beta==1) %>% 
    ggplot(aes(x = value, y = fct_rev(Estimator), color = Estimator)) + 
    stat_pointinterval(position = position_dodge(width = 0.9), .width = c(0.5, 1), point_size = 3, show.legend = T) + 
    scale_color_manual(values = nice_palette, name = "Estimator") +
    ggh4x::facet_grid2(Z~S, scales = "free", independent = "all", labeller = label_bquote(rows = rho[YZ]:.(Z), col = rho[YS]:.(S))) + ylab("") + xlab(mres[i])  +
    theme(text = element_text(size = 30), axis.text.y = element_blank(), axis.text.x = element_text(size = 15),
          axis.ticks.y = element_blank(), legend.key.size = unit(1.5, "cm"), panel.spacing = unit(2, "lines")) +
    guides(size = "none", override.aes = list(size = 4, shape = 2))
  ggsave(paste0("Plot_B1_SZ_",mres[i], ".png"), plot = p_b1_sz, path = fold_graphs, width = 50, height = 50, units = "cm")
  
  # (vii) Cross S x Z - beta 13
  p_b13_sz <- perf_plot %>% 
    filter(Z!=0.05 & Z!=0.95 & !is.na(a) & grepl(mres[i], res) & beta==13) %>% 
    ggplot(aes(x = value, y = fct_rev(Estimator), color = Estimator)) + 
    stat_pointinterval(position = position_dodge(width = 0.9), .width = c(0.5, 1), point_size = 3, show.legend = T) + 
    scale_color_manual(values = nice_palette, name = "Estimator") +
    ggh4x::facet_grid2(Z~S, scales = "free", independent = "all", labeller = label_bquote(rows = rho[YZ]:.(Z), col = rho[YS]:.(S))) + ylab("") + xlab(mres[i])  +
    theme(text = element_text(size = 30), axis.text.y = element_blank(), axis.text.x = element_text(size = 15),
          axis.ticks.y = element_blank(), legend.key.size = unit(1.5, "cm"), panel.spacing = unit(2, "lines")) +
    guides(size = "none", override.aes = list(size = 4, shape = 2))
  ggsave(paste0("Plot_B13_SZ_",mres[i], ".png"), plot = p_b13_sz, path = fold_graphs, width = 50, height = 50, units = "cm")
  
  p_b1_13 <- perf_plot %>%
    filter(Z != 0.05, Z != 0.95, !is.na(a), grepl(mres[i], res), beta %in% c(1, 13)) %>%
    mutate(alphabeta = paste0("α=", a, "|β=", b), alphabeta = factor(alphabeta, levels = unique(alphabeta))) %>%
    ggplot(aes(x = alphabeta, y = value, color = Estimator)) +
    stat_pointinterval(position = position_dodge(width = -1), .width = c(0.5, 0.8), point_size = 3, show.legend = TRUE) +
    geom_vline(xintercept = 1.5, color = "gray", linewidth = 1.2) +
    scale_color_manual(values = nice_palette, name = "Estimator") +
    ggh4x::facet_grid2(Z ~ S, scales = "free", independent = "all", labeller = label_bquote(rows = rho[YZ] : .(Z), cols = rho[YS] : .(S))) +
    labs(x = NULL, y = mres[i]) +
    theme(text = element_text(size = 30), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 30, angle = 90, hjust = 0.5),
          axis.ticks.y = element_blank(),legend.key.size = unit(1.5, "cm")) +
    guides(size = "none", override.aes = list(size = 4, shape = 2)) +
    coord_flip()
  ggsave(paste0("Plot_B1_13_SZ_",mres[i], ".png"), plot = p_b1_13, path = fold_graphs, width = 50, height = 50, units = "cm")
  
  p_b1_16 <- perf_plot %>%
    filter(Z != 0.05, Z != 0.95, !is.na(a), grepl(mres[i], res), beta %in% c(1, 16)) %>%
    mutate(alphabeta = paste0("α=", a, "|β=", b), alphabeta = factor(alphabeta, levels = unique(alphabeta))) %>%
    ggplot(aes(x = alphabeta, y = value, color = Estimator)) +
    stat_pointinterval(position = position_dodge(width = -1), .width = c(0.5, 1), point_size = 3, show.legend = TRUE) +
    geom_vline(xintercept = 1.5, color = "gray", linewidth = 1.2) +
    scale_color_manual(values = nice_palette, name = "Estimator") +
    ggh4x::facet_grid2(Z ~ S, scales = "free", independent = "all", labeller = label_bquote(rows = rho[YZ] : .(Z), cols = rho[YS] : .(S))) +
    labs(x = NULL, y = mres[i]) +
    theme(text = element_text(size = 30), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 30, angle = 90, hjust = 0.5),
          axis.ticks.y = element_blank(),legend.key.size = unit(1.5, "cm")) +
    guides(size = "none", override.aes = list(size = 4, shape = 2)) +
    coord_flip()
  ggsave(paste0("Plot_B1_16_SZ_",mres[i], ".png"), plot = p_b1_16, path = fold_graphs, width = 50, height = 50, units = "cm")
  
  p_b13_16 <- perf_plot %>%
    filter(Z != 0.05, Z != 0.95, !is.na(a), grepl(mres[i], res), beta %in% c(13, 16)) %>%
    mutate(alphabeta = paste0("α=", a, "|β=", b), alphabeta = factor(alphabeta, levels = unique(alphabeta))) %>%
    ggplot(aes(x = alphabeta, y = value, color = Estimator)) +
    stat_pointinterval(position = position_dodge(width = -1), .width = c(0.5, 1), point_size = 3, show.legend = TRUE) +
    geom_vline(xintercept = 1.5, color = "gray", linewidth = 1.2) +
    scale_color_manual(values = nice_palette, name = "Estimator") +
    ggh4x::facet_grid2(Z ~ S, scales = "free", independent = "all", labeller = label_bquote(rows = rho[YZ] : .(Z), cols = rho[YS] : .(S))) +
    labs(x = NULL, y = mres[i]) +
    theme(text = element_text(size = 30), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 30, angle = 90, hjust = 0.5),
          axis.ticks.y = element_blank(),legend.key.size = unit(1.5, "cm")) +
    guides(size = "none", override.aes = list(size = 4, shape = 2)) +
    coord_flip()
  ggsave(paste0("Plot_B13_16_SZ_",mres[i], ".png"), plot = p_b13_16, path = fold_graphs, width = 50, height = 50, units = "cm")
}

Tab_Av <- perf_plot %>%
  filter(grepl("Y", Y) & grepl("MAD", res) & S == 0.4 & !is.na(Z)) %>%
  group_by(Estimator, Z) %>%
  summarize(MeanValue=mean(value)) %>%
  spread(Z, MeanValue) %>%
  kable(format = "latex")

# Graphs to check the selection
rsh$S <- as.numeric(as.character(rsh$S));rsh$Z <- as.numeric(as.character(rsh$Z))
rs$S <- as.numeric(as.character(rs$S))

# True selection
ars  <- rs %>% 
  mutate(ars = rowMeans(rs[,grepl("RS", colnames(rs))], na.rm = T))  %>% 
  filter(Y == "Y") %>% 
  group_by(S) %>% 
  summarize(rs = mean(ars, na.rm = T))

# Estimated selection
arsh <- rsh %>% 
  mutate(ars = rowMeans(rsh[, grepl("RS", colnames(rsh))], na.rm = T)) %>% 
  filter(Y == "Y") %>% 
  group_by(S, M, Z) %>% 
  summarize(rs_h = mean(ars, na.rm = T), se_rs_h = sd(ars, na.rm = T))

arsh %>% 
  filter(Z != 0.05  & Z!= 0.95) %>%
  mutate(Estimator=M) %>% 
  ggplot(aes(x = fct_rev(Estimator), y = rs_h, fill = Estimator)) + 
  geom_col(alpha = 0.5) +
  geom_errorbar(aes(ymin = rs_h-se_rs_h,  ymax = rs_h+se_rs_h, width = 0.5)) +
  facet_grid(Z~S, labeller = label_bquote(rows = rho[YZ]:.(Z), col = rho[YS]:.(S))) +
  coord_flip() +
  geom_hline(data =ars, aes(yintercept=rs), col = "blue", linetype = "dashed") +
  theme(text = element_text(size = 30), strip.background = element_blank(), strip.placement = "outside", legend.key.size = unit(1.5, 'cm'),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), , panel.spacing = unit(2, "lines")) +
  ylab("Relative Selectivity") + xlab("Estimator") +
  scale_fill_manual(values=nice_palette) 
ggsave("Plot_AvSel.png", path = fold_graphs, width = 50, height = 50, units = "cm")

# (RS_hat - RS)  vs. MAD/RMSD
rs_f <- merge(ars, arsh, by = "S") %>% 
  mutate(d_rs = abs(rs - rs_h), Estimator=M) 

for (i in 1:length(mres)){
  p_sp <- perf_plot %>% 
    filter(Y=="Y" & is.na(a) & is.na(b) & grepl(mres[i], res) & Z!= 0.05  & Z!= 0.95) %>% 
    group_by(S, Z, Estimator) %>% 
    summarize(a_perf = mean(value, na.rm = TRUE)) %>% 
    merge(rs_f, by = c("S", "Z", "Estimator")) %>% 
    ggplot(aes(x = d_rs, y = a_perf, color = Estimator)) + 
    geom_point(size = 3, alpha = 0.5) + ylab(mres[i]) + xlab(expression(abs(bar(RS)[est] - bar(RS)[real]))) +
    facet_grid(Z~S, labeller = label_bquote(rows = rho[YZ]:.(Z), col = rho[YS]:.(S))) + scale_color_manual(values = nice_palette) +
    theme(text = element_text(size = 30), strip.background = element_blank(), strip.placement = "outside", legend.key.size = unit(1.5, 'cm'),
          panel.spacing = unit(2, "lines")) 
   ggsave(paste0("Plot_RS_", mres[i], ".png"), plot = p_sp, path = fold_graphs, width = 30, height = 30, units = "cm")
}

# MAD vs. RMSD
MAD_d <- perf_plot %>% 
  filter(Y=="Y" & is.na(a) & is.na(b) & grepl("MAD", res) & Z!= 0.05  & Z!= 0.95) 

RMSD_d <- perf_plot %>% 
  filter(Y=="Y" & is.na(a) & is.na(b) & grepl("RMSD", res) & Z!= 0.05  & Z!= 0.95) 

p_res <- merge(MAD_d[,c("S", "Z", "Estimator", "d", "value")], RMSD_d[,c("S", "Z", "Estimator", "d", "value")], by = c("S", "Z", "Estimator", "d"), suffixes = c("_MAD", "_RMSD")) %>% 
  group_by(S, Z, Estimator) %>% 
  summarize(MAD = mean(value_MAD, na.rm = TRUE), RMSD = mean(value_RMSD, na.rm = TRUE)) %>% 
  ggplot(aes(x = MAD, y = RMSD, color = Estimator)) + 
  geom_point(size = 3, alpha = 0.5) + ylab("RMSD") + xlab("MAD") + scale_color_manual(values = nice_palette) +
  facet_grid(Z~S, labeller = label_bquote(rows = rho[YZ]:.(Z), col = rho[YS]:.(S))) + 
  theme(text = element_text(size = 30),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.key.size = unit(2, 'cm'), panel.spacing = unit(2, "lines"))    
ggsave(paste0("MAD_RMSD.png"), plot = p_res, path = fold_graphs, width = 50, height = 50, units = "cm")

#### - (VII) Checks - ####

# Individual check

d <- 90
cs <- 1
cz <- 6
b <- 14
Ye <- as.matrix(Y[,,d])                       # Y matrix                                               
Se <- as.matrix(S[[cs]][,,d])                 # S matrix                                   
Ze <- as.matrix(Z[[cz]][,,d])                 # Z matrix
ye <- Ye*Se; ye[ye==0] <- NA

Me <- MUB_m(Ye, Se, Ze)[[1]]
Me
RMSD(Ye, Se, Me)

abs(Me)/(1+abs(colMeans(ye, na.rm = T)-Me))
abs(Mr)/(1+abs(colMeans(ye, na.rm = T)-Mr))

Mr <- (colMeans(ye, na.rm = T)- colMeans(Ye))
1/(colMeans(ye, na.rm = T)/Mr - 1)
Mr/(1+colMeans)

# Confirm phi

check_phi <- vector(mode = "list", length = D)
check_cors <- vector(mode = "list", length = D)

for (d in 1:D){
  check_phi[[d]] <- ifelse((Phi[[d]][,3:12]>=0 & Phi[[d]][,3:12]<=1), T, F) 
  check_phi[[d]] <- cbind(Phi[[d]][,1:2], check_phi[[d]])
  
  mat_check <- check_phi[[d]]
  mat_check[,3:12] <- NA
  for (cs in 1:n_sel){
    for (cz in 1:n_aux){
      Ye <- as.matrix(Y[,,d])                       # Y matrix                                               
      Se <- as.matrix(S[[cs]][,,d])                 # S matrix                                   
      Ze <- as.matrix(Z[[cz]][,,d])                 # Z matrix
      
      for(t in 1:Time){
        r_zs <- cor(Ze[,t], Se[,t])
        r_yz <- cor(Ye[,t], Ze[,t])
        r_ys <- cor(Ye[,t], Se[,t])
        mat_check[cs+((cz-1)*4),2+t] <- (r_zs>=r_yz*r_ys & r_ys>=r_zs*r_yz)
      }
    }
  }
  check_cors[[d]] <- mat_check
}
check_phi[[20]]==check_cors[[20]]
Phi[[5]]




# Relative absolute deviation



colnames(rs) <-  c("S", "Y", "beta", "a", "b", "d",paste0("RS",1:12,"_t"))
join_rs <- left_join(rs, rsh, by = c("S", "Y", "beta", "a", "b", "d"))

rs_cols <- grep("RS", colnames(join_rs))
RAD <- (join_rs[,rs_cols[14:24]] - join_rs[,rs_cols[2:12]])/join_rs[,rs_cols[2:12]] 
RAD <- cbind(join_rs[,grep("RS",colnames(join_rs), invert = T)], RAD)
RAD %>% 
  group_by(S, Y, beta, a, b, Z, M) %>% 
  summarize(mean_RS2 = mean(RS10, na.rm = T))









