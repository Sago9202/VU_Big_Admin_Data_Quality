## --------------------------------------------------------------- ##
## Departure from Ignorable Sample Selection with Ordinal Outcomes ##
## Script 1 - Simulation                                           ##
## Written by: Santiago Gómez-Echeverry                            ##
## Last modified: 15/10/2026                                       ##
## --------------------------------------------------------------- ##

#### - (I) Working space and packages - ####

rm(list = ls())
packs <- c('MASS', 'dplyr', 'tidyr', 'purrr', 'tibble', 'ggplot2', 'MetBrewer', 'progress', 'polycor', 'paletteer')
ipacks <- packs  %in% rownames(installed.packages()) 
if(any(ipacks == F)){install.packages(packs[!ipacks])}
invisible(lapply(packs, library, character.only = T))
#library(nationalparkcolors) get using devtools
fold_data <- "C:/Users/sagom/Dropbox/PhD VU/Papers/5 - Ordinal MUB/3 - Data"
fold_graphs <- "C:/Users/sagom/Dropbox/PhD VU/Papers/5 - Ordinal MUB/4 - Graphs and Tables"
setwd(fold_data)

#### - (II) Data Generation Process (DGP) - ####

# Before the data generation we need to create some additional functions

# 1. Get the proportions (probs) in each category
get_probs <- function(ncat, bal, gamma = 2){
  if (bal == 1){
    probs <- rep(1/ncat, ncat)
  } else {
    k <- 1:ncat
    wg <- (k*(ncat+1-k))^gamma
    probs <- wg/sum(wg)
  }
  return(probs)
}

# 2. Obtain the corresponding thresholds
# ! Note that we can apply this because U ~ N(0,1)
pr_2_thr <- function(probs){
  K <- length(probs)
  cum_probs <- cumsum(probs[1:(K-1)])
  qnorm(cum_probs)   
}

## Factors affecting the performance
## Number of categories - Done
## Distribution of categories - Done
## Selection - Done
## Proxy strength - Done
## Proxy shift - Done

set.seed(42)
N <- 1e4                    # Population size
D <-  200                   # Number of replications
pi <- 0.05          
n <- N*pi                   # Observed sample
#rho_us <- 0.4
rho_ux <- c(0.2, 0.5, 0.8)  # Proxy's strength
phi <- c(0.2, 0.5, 0.8)     # Non-ignorable selection parameter
ncat <- c(3, 5, 10)
bal_vals <- c(0,1)
shift_x <- c(-0.5, 0, 0.5)
param <- expand.grid("rho_ux" =  rho_ux, "phi" =  phi, "ncat" = ncat, "bal" = bal_vals, D = 1:D, "shift_x" = shift_x, "summary" = c("arithmetic", "closed_geometric"))
#param$rho_xs <- (rho_us*(1-param$phi+param$phi*param$rho_ux))/(param$phi - param$phi*param$rho_ux + param$rho_ux)
nsim <- nrow(param)

cuts <- vector(mode = "list")
cuts$probs <- mapply(get_probs, param$ncat, param$bal, SIMPLIFY = F)
cuts$thresholds <- lapply(cuts$probs, pr_2_thr)

## Plot of the cases 

# Function to create smooth density data for a scenario
create_density_data <- function(ncat, bal, thresholds) {
  x <- seq(-4.5, 4.5, length.out = 3000)
  dens <- dnorm(x)
  
  # Assign each x to a category
  breaks <- c(-Inf, thresholds, Inf)
  y_cat <- as.numeric(cut(x, breaks = breaks, labels = 1:ncat, right = T))
  
  data.frame(
    x = x,
    dens = dens,
    y_cat = factor(y_cat),
    ncat = factor(paste0("K = ", ncat)),
    bal = factor(ifelse(bal == 1, "Balanced", "Sparse tails"), 
                 levels = c("Balanced", "Sparse tails"))
  )
}

unique_scenarios <- unique(param[, c("ncat", "bal")])
unique_scenarios <- unique_scenarios[order(unique_scenarios$ncat, unique_scenarios$bal), ]

# Generate density data for all scenarios
dens_data <- list()
for (i in 1:nrow(unique_scenarios)) {
  ncat_i <- unique_scenarios$ncat[i]
  bal_i <- unique_scenarios$bal[i]
  idx <- which(param$ncat == ncat_i & param$bal == bal_i)[1]
  thr <- cuts$thresholds[[idx]]
  
  dens_data[[i]] <- create_density_data(ncat_i, bal_i, thr)
}
dens_df <- do.call(rbind, dens_data)

# Add vertical lines for thresholds
thr_df <- data.frame()
for (i in 1:nrow(unique_scenarios)) {
  ncat_i <- unique_scenarios$ncat[i]
  bal_i <- unique_scenarios$bal[i]
  idx <- which(param$ncat == ncat_i & param$bal == bal_i)[1]
  thr <- cuts$thresholds[[idx]]
  
  thr_df <- rbind(thr_df, data.frame(
    threshold = thr,
    ncat = factor(paste0("K = ", ncat_i)),
    bal = factor(ifelse(bal_i == 1, "Balanced", "Sparse tails"), 
                 levels = c("Balanced", "Sparse tails"))
  ))
}

label_positions <- data.frame()
leader_lines <- data.frame()

for (i in 1:nrow(unique_scenarios)) {
  ncat_i <- unique_scenarios$ncat[i]
  bal_i  <- unique_scenarios$bal[i]
  
  idx   <- which(param$ncat == ncat_i & param$bal == bal_i)[1]
  probs <- cuts$probs[[idx]]
  cum_probs <- c(0, cumsum(probs))
  
  fl <- data.frame()
  for (k in 1:ncat_i) {
    p_mid <- (cum_probs[k] + cum_probs[k+1]) / 2
    x_mid <- qnorm(p_mid)
    pct <- probs[k] * 100
    label <- if (abs(pct - round(pct)) < 1e-6) sprintf("%.0f%%", pct) else sprintf("%.1f%%", pct)
    base_y <- dnorm(x_mid) + 0.02
    
    x_lab <- x_mid
    y_lab <- base_y
    
    # ONLY for K = 10, categories 5 and 6: push up + outward, different heights
    if (ncat_i == 10 && k == 5) {
      x_lab <- x_mid - 0.15
      y_lab <- base_y + 0.03
    }
    if (ncat_i == 10 && k == 6) {
      x_lab <- x_mid + 0.15
      y_lab <- base_y + 0.07
    }
    
    fl <- rbind(fl, data.frame(k = k, label = label, x_true = x_mid, y_true = base_y,
                               x = x_lab, y = y_lab, stringsAsFactors = FALSE))
  }
  
  ncat_f <- factor(paste0("K = ", ncat_i), levels = paste0("K = ", sort(unique(param$ncat))))
  bal_f  <- factor(ifelse(bal_i == 1, "Balanced", "Sparse tails"), levels = c("Balanced", "Sparse tails"))
  
  label_positions <- rbind(label_positions, data.frame(
    ncat = ncat_f, bal = bal_f, label = fl$label, x = fl$x, y = fl$y, stringsAsFactors = FALSE
  ))
  
  # leader lines only for the moved labels (K = 10, categories 5 and 6)
  moved <- fl[fl$k %in% c(5, 6) & ncat_i == 10, ]
  if (nrow(moved) > 0) {
    leader_lines <- rbind(leader_lines, data.frame(ncat = ncat_f, bal = bal_f, x = moved$x_true, xend = moved$x, y = moved$y_true, yend = moved$y - 0.002, stringsAsFactors = F))
  }
}

p1 <- ggplot(dens_df, aes(x = x, y = dens)) +
  geom_area(aes(fill = y_cat), position = "stack") +
  geom_line(color = "white", linewidth = 0.5) +
  geom_vline(data = thr_df, aes(xintercept = threshold), color = "white", linetype = "dashed", linewidth = 0.6) +
  facet_grid(ncat ~ bal) +
  geom_segment(data = leader_lines, aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, color = "grey40", linewidth = 0.3) +
  geom_text(data = label_positions, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 6, color = "black", alpha = 0.8) +
  labs(x = "U", y = "Density", fill = "Y") +
  theme_minimal(base_size = 12) +
  theme(strip.background = element_rect(fill = "lightgrey", colour = NA),
        legend.position = "bottom", legend.key = element_blank(),
        text = element_text(size = 30)) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_fill_met_d(name = "Benedictus", direction = -1)
p1
ggsave("p1.png", plot = p1, path = fold_graphs, width = 45, height = 30, units = "cm")

### Generate the sampling and then based on that we create the rest

Dat <- vector(mode = "list", length = nsim)
pb <- progress_bar$new(total = nrow(param))

for (s in 1:nsim){
  pb$tick()
  rho_ux_s <- param$rho_ux[s]
  phi_s    <- param$phi[s]
  
  # Respondent-domain moments (standard normalization)
  mu1_u <- 0; sigma1_uu <- 1
  mu1_x <- 0; sigma1_xx <- 1
  
  # Nonrespondent-domain proxy mean shift (fixed design choice)
  mu0_x <- param$shift_x[s]
  sigma0_xx <- 1
  
  # Back out mu0_u, sigma0_uu from Eqs. (5)-(6) using the true phi
  g_phi <- (phi_s + (1-phi_s)*rho_ux_s) / (phi_s*rho_ux_s + (1-phi_s))
  mu0_u <- mu1_u + g_phi*(mu0_x - mu1_x)
  sigma0_uu <- 1 + (g_phi^2)*(sigma0_xx - sigma1_xx)
  
  # Draw S first, then U,X within each stratum
  S <- rbinom(N, 1, pi)
  n1 <- sum(S == 1); n0 <- N - n1
  
  R <- matrix(c(1, rho_ux_s, rho_ux_s, 1), nrow = 2)
  
  d1 <- MASS::mvrnorm(n1, mu = c(mu1_u, mu1_x),
                      Sigma = R * outer(c(1, sqrt(sigma1_xx)), c(1, sqrt(sigma1_xx))))
  d0 <- MASS::mvrnorm(n0, mu = c(mu0_u, mu0_x),
                      Sigma = R * outer(c(sqrt(sigma0_uu), sqrt(sigma0_xx)), c(sqrt(sigma0_uu), sqrt(sigma0_xx))))
  
  dat <- rbind(data.frame(U = d1[,1], X = d1[,2], S = 1),data.frame(U = d0[,1], X = d0[,2], S = 0))
  
  # Cut U into the ordinal target using this scenario's thresholds
  thr <- cuts$thresholds[[s]]
  dat$Y <- factor(cut(dat$U, breaks = c(-Inf, thr, Inf), labels = 1:param$ncat[s]),
                  levels = 1:param$ncat[s])
  
  Dat[[s]] <- dat
}

#### - (III) MUBP function - ####

mubp_k <- function(data, summary = c("arithmetic", "closed_geometric")){
  # 1. Split the data into sample and out of sample
  s_data <- data[data$S == 1,]
  ns_data <- data[data$S == 0,]
  
  # 2. Fit ordinal response model to obtain the cutpoints  
  reg0 <- polr(Y ~ 1, data = s_data, method = "probit")
  tau_hat <- reg0$zeta
  all_taus <- c(-Inf, tau_hat, Inf)
  K <- length(reg0$zeta) + 1
  
  # 3. Proxy moments 
  mu1_x <- mean(s_data$X)
  sigma1_xx <- var(s_data$X)
  mu0_x <- mean(ns_data$X)
  sigma0_xx <- var(ns_data$X)
  rho_ux <- polyserial(s_data$X, s_data$Y)
  pi_sample <- mean(data$S == 1)
  
  # 4. Monte Carlo integration over phi ~ Uniform(0,1)
  M <- 5000
  phi <- runif(M)
  
  g_phi <- ((phi + (1 - phi) * rho_ux) / (phi * rho_ux + (1 - phi)))
  mu0_u <- g_phi * ((mu0_x - mu1_x) / sqrt(sigma1_xx))
  sigma0_uu <- 1 + (g_phi^2) * ((sigma0_xx - sigma1_xx) / sigma1_xx)
  
  # 5. Discard the invalid cases
  valid <- sigma0_uu > 0
  m_valid <- sum(valid)
  m_dropped <- M - m_valid
  
  if (m_valid == 0) {
    out <- rep(NA_real_, K)
    attr(out, "m_dropped") <- M
    return(out)
  }
  
  g_phi <- g_phi[valid]; mu0_u <- mu0_u[valid]; sigma0_uu <- sigma0_uu[valid]
  
  # 6. Evaluate the thresholds
  z1_mat <- outer(rep(1, m_valid), all_taus)
  z0_mat <- outer(1 / sqrt(sigma0_uu), all_taus) - outer(mu0_u / sqrt(sigma0_uu), rep(1, length(all_taus)))
  p1_mat <- pnorm(z1_mat); p0_mat <- pnorm(z0_mat)
  v_k_1_mat <- p1_mat[, -1, drop = F] - p1_mat[, -ncol(p1_mat), drop = F]
  v_k_0_mat <- p0_mat[, -1, drop = F] - p0_mat[, -ncol(p0_mat), drop = F]
  v_k_pop_mat <- pi_sample * v_k_1_mat + (1 - pi_sample) * v_k_0_mat
  
  # 7. Summarize the draws of the population proportions
  v_k_1 <- v_k_1_mat[1,]
  if (summary == "arithmetic"){
    v_k_pop <- colMeans(v_k_pop_mat)
  } else{
    log_v <- log(pmax(v_k_pop_mat, .Machine$double.xmin)) # Guard against log(0)
    g <- exp(colMeans(log_v))
    v_k_pop <- g/sum(g)
  }
  
  # 8. MUBP
  out <- v_k_1 - v_k_pop
  names(out) <- paste0("C", 1:K)
  attr(out, "m_dropped") <- m_dropped
  attr(out, "summary") <- summary
  return(out)
}

#### - (IV) Estimation - ####

Results <- vector(mode = "list", length = nsim)

pb <- progress_bar$new(total = nsim)
for (s in 1:nsim){
  pb$tick()
  
  t_dat <- Dat[[s]]
  t_K <- param$ncat[[s]]
  
  # 1. Ensure that the target variable is a factor
  t_dat$Y <- factor(t_dat$Y, levels = 1:t_K)
  
  # 2. Store the true values
  v1_tab   <- table(factor(t_dat$Y[t_dat$S == 1], levels = 1:t_K))
  vpop_tab <- table(factor(t_dat$Y,             levels = 1:t_K))
  v1_true   <- as.numeric(prop.table(v1_tab))   
  vpop_true <- as.numeric(prop.table(vpop_tab)) 
  mubp_true <- v1_true - vpop_true              
  
  # 3. Estimate the MUBP
  mubp_est <- mubp_k(t_dat, summary = param$summary[[s]])
  summary <- as.character(param$summary[[s]])
  m_dropped <- attr(mubp_est, "m_dropped")
  
  # 4. Store all the information
  Results[[s]] <- data.frame(sim_id = s, d = param$D[s], rho_ux=param$rho_ux[s], phi = param$phi[s], ncat = param$ncat[s], bal = param$bal[s],
                             category = paste0("C", 1:t_K), v1_true   = as.numeric(v1_true), vpop_true = as.numeric(vpop_true), 
                             mubp_true = mubp_true, mubp_est  = as.numeric(mubp_est), mubp_err  = as.numeric(mubp_est) - mubp_true, M_drop = m_dropped, shift_x = param$shift_x[s], Sum = summary)
  
}

sim_results <- bind_rows(Results)
write.table(sim_results, file = paste0(fold_data, "/Simulated_Data_N", N, "_D", D, ".txt"))

sim_results <- read.table("Simulated_Data_N10000_D200.txt")

#### - (V) Graphs and Tables - ####

summary_agg <- sim_results %>%
  group_by(sim_id, d, rho_ux, phi, ncat, bal, shift_x) %>%
  summarise(rmse = sqrt(mean(mubp_err^2)), mae  = mean(abs(mubp_err)),.groups = "drop") %>%
  group_by(rho_ux, phi, ncat, bal, shift_x) %>%
  summarise(mean_rmse = mean(rmse), sd_rmse = sd(rmse), mean_mae  = mean(mae),  sd_mae  = sd(mae), .groups = "drop") %>%
  mutate(
    ncat_f = factor(paste0(ncat), levels = c(3, 5, 10)),
    bal_f  = factor(ifelse(bal == 1, "Balanced", "Sparse tails"), levels = c("Balanced", "Sparse tails")),
    phi_f  = factor(paste0("phi == ", phi)),
    rho_f  = factor(paste0("rho[UX] == ", rho_ux)),
    shift_f   = factor(shift_x), group_bp = interaction(bal_f, shift_x, sep = " | "))

p2 <- ggplot(summary_agg, aes(x = ncat_f, y = mean_rmse, color = shift_f, shape = bal_f, group = group_bp)) +
  #geom_line(aes(linetype = pi_f), linewidth = 0.7, position = position_dodge(width = 0.3)) +
  geom_point(size = 2.8, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = pmax(0, mean_rmse - sd_rmse), ymax = mean_rmse + sd_rmse),
                width = 0.15, position = position_dodge(width = 0.3), alpha = 0.5) +
  facet_grid(phi_f ~ rho_f, labeller = label_parsed) +
  scale_color_paletteer_d("nationalparkcolors::Acadia") +
  labs(x = "K", y = "RMSE", color = expression(Delta[x]), shape = "Category Distribution") +
  theme_minimal(base_size = 13) +
  theme(text = element_text(size = 20), strip.background = element_rect(fill = "grey90", color = NA),
        strip.text = element_text(face = "bold", size = 15), legend.position = "bottom", legend.box = "vertical",
        panel.grid.minor = element_blank())
p2
ggsave("p2.png", plot = p2, path = fold_graphs, width = 20, height = 20, units = "cm")

summary_cat <- sim_results %>%
  group_by(rho_ux, phi, ncat, bal, category, shift_x) %>%
  summarise(mean_err = mean(mubp_err), rmse = sqrt(mean(mubp_err^2)),mae = mean(abs(mubp_err)), .groups = "drop") %>%
  mutate(ncat_f = factor(paste0("K = ", ncat), levels = c("K = 3", "K = 5", "K = 10")), 
         bal_f  = factor(ifelse(bal == 1, "Balanced", "Sparse tails"), levels = c("Balanced", "Sparse tails")),
         phi_f  = factor(paste0("phi == ", phi)), rho_f  = factor(paste0("rho[UX] == ", rho_ux)), 
         shift_f   = factor(shift_x), cat_num = as.numeric(gsub("C", "", category)))

plot_by_K <- function(K_val) {
  df_k <- summary_cat %>% 
    filter(ncat == K_val)
  ggplot(df_k, aes(x = factor(cat_num), y = rmse, color = shift_f, shape = bal_f, group = interaction(shift_f, bal_f))) +
    geom_line(aes(linetype = shift_f), linewidth = 0.7, position = position_dodge(width = 0.25)) +
    geom_point(size = 2.3, position = position_dodge(width = 0.25)) +
    facet_grid(phi_f ~ rho_f, labeller = label_parsed) +  scale_color_paletteer_d("nationalparkcolors::Acadia") +
    labs(x = "Category", y = "RMSE", color = expression(Delta[x]), shape = "Category Distribution", linetype =expression(Delta[x])) +
    theme_minimal(base_size = 12) + theme(text = element_text(size = 20), strip.text = element_text(face = "bold", size = 15),
                                          strip.background = element_rect(fill = "grey90", color = NA), legend.position = "bottom",
                                          legend.box = "vertical", panel.grid.minor = element_blank())
}

p3  <- plot_by_K(3); p4  <- plot_by_K(5); p5 <- plot_by_K(10)
p3; p4; p5
ggsave("p3.png", plot = p3, path = fold_graphs, width = 20, height = 20, units = "cm")
ggsave("p4.png", plot = p4, path = fold_graphs, width = 20, height = 20, units = "cm")
ggsave("p5.png", plot = p5, path = fold_graphs, width = 20, height = 20, units = "cm")



summary_decomp <- sim_results %>%
  group_by(rho_ux, phi, ncat, bal, shift_x, category) %>%
  summarise(
    bias = mean(mubp_est) - mean(mubp_true),
    variance = var(mubp_est),
    mse = mean(mubp_err^2),
    .groups = "drop"
  ) %>%
  mutate(
    bias_sq = bias^2,
    check_diff = mse - (bias_sq + variance)   # sanity check, should be ~0
  ) %>%
  pivot_longer(cols = c(bias_sq, variance), names_to = "component", values_to = "value") %>%
  mutate(
    component = factor(component, levels = c("variance", "bias_sq"),
                       labels = c("Variance", "Bias²")),
    cat_num = as.numeric(gsub("C", "", category)),
    ncat_f  = factor(paste0("K = ", ncat), levels = c("K = 3", "K = 5", "K = 10")),
    bal_f   = factor(ifelse(bal == 1, "Balanced", "Sparse tails"), levels = c("Balanced", "Sparse tails")),
    phi_f   = factor(paste0("phi == ", phi)),
    rho_f   = factor(paste0("rho[UX] == ", rho_ux)),
    shift_f = factor(shift_x)
  )

summary_decomp_pct <- summary_decomp %>%
  mutate(shift_f = factor(paste0("Delta[x] == ", shift_x))) %>%   # relabeled for parsing
  group_by(rho_ux, phi, ncat, bal, shift_x, category, ncat_f, bal_f, phi_f, rho_f, shift_f, cat_num) %>%
  mutate(pct = value / sum(value) * 100) %>%
  ungroup()

plot_decomp_pct <- function(K_val){
  df_k <- summary_decomp_pct %>% filter(ncat == K_val)
  
  ggplot(df_k, aes(x = factor(cat_num), y = pct, fill = component)) +
    geom_col(position = "stack") +
    facet_nested(rho_f + shift_f ~ bal_f + phi_f,
                 labeller = labeller(phi_f = label_parsed, rho_f = label_parsed, shift_f = label_parsed)) +
    scale_fill_manual(values = c("Variance" = "#FED789FF", "Bias²" ="#023743FF")) +
    labs(x = "Category", y = "% of MSE", fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(text = element_text(size = 20),
          strip.background = element_rect(fill = "grey90", color = "white"),
          strip.text = element_text(face = "bold", size = 20),
          legend.position = "bottom",
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0.3, "lines"))
}

#c("#FED789FF","#023743FF", "#72874EFF", "#476F84FF", "#A4BED5FF")
#paquin_cols <-MetBrewer::MetPalettes$Paquin[[1]]

p_decomp_3  <- plot_decomp_pct(3)
p_decomp_5  <- plot_decomp_pct(5)
p_decomp_10 <- plot_decomp_pct(10)

p_decomp_3; p_decomp_5; p_decomp_10
# Manually save with 1200x1500
ggsave("p6.png", plot = p_decomp_3, path = fold_graphs, width = 30, height = 50, units = "cm")
ggsave("p7.png", plot = p_decomp_5, path = fold_graphs, width = 30, height = 50, units = "cm")
ggsave("p8.png", plot = p_decomp_10, path = fold_graphs, width = 30, height = 50, units = "cm")


##### Relative plots

summary_cat <- sim_results %>%
  group_by(rho_ux, phi, ncat, bal, category, shift_x) %>%
  summarise(mean_err = mean(mubp_err),
            rmse     = sqrt(mean(mubp_err^2)),
            mae      = mean(abs(mubp_err)),
            v_true   = mean(vpop_true),          # true population proportion of the category
            .groups  = "drop") %>%
  mutate(rrmse    = rmse / v_true,            # relative RMSE
         rel_bias = mean_err / v_true,        # relative bias (keeps the sign)
         ncat_f  = factor(paste0("K = ", ncat), levels = c("K = 3", "K = 5", "K = 10")),
         bal_f   = factor(ifelse(bal == 1, "Balanced", "Sparse tails"), levels = c("Balanced", "Sparse tails")),
         phi_f   = factor(paste0("phi == ", phi)),
         rho_f   = factor(paste0("rho[UX] == ", rho_ux)),
         shift_f = factor(shift_x),
         cat_num = as.numeric(gsub("C", "", category)))

plot_by_K <- function(K_val, metric = "rmse", y_lab = "RMSE") {
  df_k <- summary_cat %>% filter(ncat == K_val)
  ggplot(df_k, aes(x = factor(cat_num), y = .data[[metric]], color = shift_f, shape = bal_f,
                   group = interaction(shift_f, bal_f))) +
    geom_line(aes(linetype = shift_f), linewidth = 0.7, position = position_dodge(width = 0.25)) +
    geom_point(size = 2.3, position = position_dodge(width = 0.25)) +
    facet_grid(phi_f ~ rho_f, labeller = label_parsed) +
    scale_color_paletteer_d("nationalparkcolors::Acadia") +
    labs(x = "Category", y = y_lab, color = expression(Delta[x]),
         shape = "Category Distribution", linetype = expression(Delta[x])) +
    theme_minimal(base_size = 12) +
    theme(text = element_text(size = 20),
          strip.text = element_text(face = "bold", size = 15),
          strip.background = element_rect(fill = "grey90", color = NA),
          legend.position = "bottom", legend.box = "vertical",
          panel.grid.minor = element_blank())
}

# Separate plots, same layout as before
p3  <- plot_by_K(3);  p4  <- plot_by_K(5);  p5  <- plot_by_K(10)
p3r <- plot_by_K(3,  "rrmse", "Relative RMSE")
p4r <- plot_by_K(5,  "rrmse", "Relative RMSE")
p5r <- plot_by_K(10, "rrmse", "Relative RMSE")
ggsave("p3r.png", plot = p3r, path = fold_graphs, width = 20, height = 20, units = "cm")
ggsave("p4r.png", plot = p4r, path = fold_graphs, width = 20, height = 20, units = "cm")
ggsave("p5r.png", plot = p5r, path = fold_graphs, width = 20, height = 20, units = "cm")
