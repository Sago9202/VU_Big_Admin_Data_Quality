## --------------------------------------------------------------- ##
## Departure from Ignorable Sample Selection with Ordinal Outcomes ##
## Script 1 - Simulation                                           ##
## Written by: Santiago Gómez-Echeverry                            ##
## Last modified: 15/10/2026                                       ##
## --------------------------------------------------------------- ##

#### - (I) Working space and packages - ####

rm(list = ls())
packs <- c('MASS', 'dplyr', 'tidyr', 'purrr', 'tibble', 'ggplot2', 'MetBrewer', 'progress', 'polycor', 'paletteer', 'ggh4x', 'patchwork')
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

setwd(fold_graphs)

## 1. Aggregate plots

p2a <- sim_results %>%
  dplyr::filter(Sum == "arithmetic") |> 
  group_by(sim_id, d, rho_ux, phi, ncat, bal, shift_x) %>%
  summarise(rmse = sqrt(mean(mubp_err^2)), mae  = mean(abs(mubp_err)),.groups = "drop") %>%
  group_by(rho_ux, phi, ncat, bal, shift_x) %>%
  summarise(mean_rmse = mean(rmse), sd_rmse = sd(rmse), mean_mae  = mean(mae),  sd_mae  = sd(mae), .groups = "drop") %>%
  mutate(ncat_f = factor(paste0(ncat), levels = c(3, 5, 10)), bal_f  = factor(ifelse(bal == 1, "Balanced", "Sparse tails"), levels = c("Balanced", "Sparse tails")),
         phi_f  = factor(paste0("phi == ", phi)), rho_f  = factor(paste0("rho[ux] == ", rho_ux)), shift_f   = factor(shift_x), group_bp = interaction(bal_f, shift_x, sep = " | ")) |> 
  ggplot(aes(x = ncat_f, y = mean_rmse, color = shift_f, shape = bal_f, group = group_bp)) +
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
ggsave("p2a.png", plot = p2a, path = fold_graphs, width = 20, height = 20, units = "cm")
p2b <- sim_results %>%
  dplyr::filter(Sum == "closed_geometric") |> 
  group_by(sim_id, d, rho_ux, phi, ncat, bal, shift_x) %>%
  summarise(rmse = sqrt(mean(mubp_err^2)), mae  = mean(abs(mubp_err)),.groups = "drop") %>%
  group_by(rho_ux, phi, ncat, bal, shift_x) %>%
  summarise(mean_rmse = mean(rmse), sd_rmse = sd(rmse), mean_mae  = mean(mae),  sd_mae  = sd(mae), .groups = "drop") %>%
  mutate(ncat_f = factor(paste0(ncat), levels = c(3, 5, 10)), bal_f  = factor(ifelse(bal == 1, "Balanced", "Sparse tails"), levels = c("Balanced", "Sparse tails")),
         phi_f  = factor(paste0("phi == ", phi)), rho_f  = factor(paste0("rho[ux] == ", rho_ux)), shift_f   = factor(shift_x), group_bp = interaction(bal_f, shift_x, sep = " | ")) |> 
  ggplot(aes(x = ncat_f, y = mean_rmse, color = shift_f, shape = bal_f, group = group_bp)) +
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
ggsave("p2b.png", plot = p2b, path = fold_graphs, width = 20, height = 20, units = "cm")


## 2. RMSE and RRMSE per category

a_cat <- sim_results |> 
  dplyr::filter(Sum == "arithmetic") |> 
  group_by(rho_ux, phi, ncat, bal, category, shift_x) |> 
  summarise(mean_err = mean(mubp_err), rmse     = sqrt(mean(mubp_err^2)), mae = mean(abs(mubp_err)), v_true = mean(vpop_true), .groups  = "drop") |> 
  mutate(rrmse    = rmse / v_true, rel_bias = mean_err / v_true, ncat_f  = factor(paste0("K = ", ncat), levels = c("K = 3", "K = 5", "K = 10")),
         bal_f   = factor(ifelse(bal == 1, "Balanced", "Sparse tails"), levels = c("Balanced", "Sparse tails")), phi_f   = factor(paste0("phi == ", phi)),
         rho_f   = factor(paste0("rho[ux] == ", rho_ux)), shift_f = factor(shift_x), cat_num = as.numeric(gsub("C", "", category)))

g_cat <- sim_results |> 
  dplyr::filter(Sum == "closed_geometric") |> 
  group_by(rho_ux, phi, ncat, bal, category, shift_x) |> 
  summarise(mean_err = mean(mubp_err), rmse     = sqrt(mean(mubp_err^2)), mae = mean(abs(mubp_err)), v_true = mean(vpop_true), .groups  = "drop") |> 
  mutate(rrmse    = rmse / v_true, rel_bias = mean_err / v_true, ncat_f  = factor(paste0("K = ", ncat), levels = c("K = 3", "K = 5", "K = 10")),
         bal_f   = factor(ifelse(bal == 1, "Balanced", "Sparse tails"), levels = c("Balanced", "Sparse tails")), phi_f   = factor(paste0("phi == ", phi)),
         rho_f   = factor(paste0("rho[ux] == ", rho_ux)), shift_f = factor(shift_x), cat_num = as.numeric(gsub("C", "", category)))


panel_plot <- function(df, metric, y_lab, y_lim = NULL) {
  ggplot(df, aes(x = cat_num, y = .data[[metric]], color = shift_f, shape = bal_f, linetype = bal_f,
                 group = interaction(shift_f, bal_f))) +
    geom_line(linewidth = 0.45) +
    geom_point(size = 1.3) +
    facet_grid(phi_f ~ rho_f, labeller = label_parsed) +
    scale_x_continuous(breaks = function(l) seq(1, floor(l[2]))) +
    coord_cartesian(ylim = y_lim) +
    scale_color_paletteer_d("nationalparkcolors::Acadia") +
    scale_linetype_manual(values = c("Balanced" = "solid", "Sparse tails" = "22")) +
    labs(x = "Category", y = y_lab, color = expression(Delta[x]),
         shape = "Category distribution", linetype = "Category distribution") +
    theme_bw(base_size = 9) +
    theme(strip.background = element_rect(fill = "grey92", color = NA),
          strip.text = element_text(size = 8),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")
}

fig_by_K <- function(K_val) {
  a_df <- filter(a_cat, ncat == K_val)
  g_df <- filter(g_cat, ncat == K_val)
  
  # Keep the same range per row
  lim_rmse  <- c(0, max(a_df$rmse,  g_df$rmse,  na.rm = TRUE))
  lim_rrmse <- c(0, max(a_df$rrmse, g_df$rrmse, na.rm = TRUE))
  
  (panel_plot(a_df, "rmse",  "RMSE", lim_rmse)  + ggtitle("(a) Arithmetic mean - Absolute error") | panel_plot(g_df, "rmse",  "RMSE", lim_rmse)  + ggtitle("(b) Normalized geometric mean - Absolute error")) /
  (panel_plot(a_df, "rrmse", "Relative RMSE", lim_rrmse) + ggtitle("(c) Arithmetic mean - Relative error") | panel_plot(g_df, "rrmse", "Relative RMSE", lim_rrmse) + ggtitle("(d) Normalized geometric mean - Relative error")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", plot.title = element_text(size = 9, face = "bold"))
}

for (K in c(3, 5, 10)){
  ggsave(sprintf("fig_rmse_K%d.pdf", K), fig_by_K(K), width = 12.6, height = 8.2, device = cairo_pdf)
}

## 3. Bias and Variance decomposition

a_decomp <- sim_results |> 
  dplyr::filter(Sum == "arithmetic") |>
  group_by(rho_ux, phi, ncat, bal, shift_x, category) |> 
  summarise(bias = mean(mubp_est) - mean(mubp_true), variance = var(mubp_est), mse = mean(mubp_err^2), .groups = "drop") |> 
  mutate(bias_sq = bias^2, check_diff = mse - (bias_sq + variance)) |> 
  pivot_longer(cols = c(bias_sq, variance), names_to = "component", values_to = "value") |> 
  mutate(component = factor(component, levels = c("variance", "bias_sq"), labels = c("Variance", "Bias²")), cat_num = as.numeric(gsub("C", "", category)), 
         ncat_f  = factor(paste0("K = ", ncat), levels = c("K = 3", "K = 5", "K = 10")),bal_f = factor(ifelse(bal == 1, "Balanced", "Sparse tails"), levels = c("Balanced", "Sparse tails")),
         phi_f   = factor(paste0("phi == ", phi)), rho_f   = factor(paste0("rho[ux] == ", rho_ux)),  shift_f = factor(shift_x))

g_decomp <- sim_results |> 
  dplyr::filter(Sum == "closed_geometric") |>
  group_by(rho_ux, phi, ncat, bal, shift_x, category) |> 
  summarise(bias = mean(mubp_est) - mean(mubp_true), variance = var(mubp_est), mse = mean(mubp_err^2), .groups = "drop") |> 
  mutate(bias_sq = bias^2, check_diff = mse - (bias_sq + variance)) |> 
  pivot_longer(cols = c(bias_sq, variance), names_to = "component", values_to = "value") |> 
  mutate(component = factor(component, levels = c("variance", "bias_sq"), labels = c("Variance", "Bias²")), cat_num = as.numeric(gsub("C", "", category)), 
         ncat_f  = factor(paste0("K = ", ncat), levels = c("K = 3", "K = 5", "K = 10")),bal_f = factor(ifelse(bal == 1, "Balanced", "Sparse tails"), levels = c("Balanced", "Sparse tails")),
         phi_f   = factor(paste0("phi == ", phi)), rho_f   = factor(paste0("rho[ux] == ", rho_ux)),  shift_f = factor(shift_x))


prep_decomp <- function(summary_decomp) {
  summary_decomp |> 
    select(rho_ux, phi, ncat, bal, shift_x, cat_num, component, value) |> 
    pivot_wider(names_from = component, values_from = value) |> 
    mutate(bias_share = `Bias²` / (`Bias²` + Variance), bal_f   = factor(ifelse(bal == 1, "Balanced", "Sparse tails"), levels = c("Balanced", "Sparse tails")),
           rho_f   = factor(paste0("rho[ux] == ", rho_ux)), shift_f = factor(paste0("Delta[x] == ", shift_x), levels = paste0("Delta[x] == ", c(-0.5, 0, 0.5))),
           phi_f   = factor(paste0("phi == ", phi), levels = rev(paste0("phi == ", c(0.2, 0.5, 0.8)))))
}

plot_bias_heat <- function(df, K_val, title) {
  df <- filter(df, ncat == K_val)
  x_breaks <- if (K_val > 5) c(1, 5, 10) else seq_len(K_val)
  ggplot(df, aes(x = factor(cat_num), y = phi_f, fill = bias_share)) +
    geom_tile(color = "white", linewidth = 0.3) +
    facet_nested(shift_f ~ bal_f + rho_f, switch = "y",labeller = labeller(rho_f = label_parsed, shift_f = label_parsed)) +
    scale_fill_gradientn(colours = c("#FED789", "#C9B97A", "#72874E", "#023743"), limits = c(0, 1), labels = scales::percent,
                         name = expression("Share of MSE due to bias"^2)) +
    scale_y_discrete(labels = function(x) parse(text = x)) +
    scale_x_discrete(breaks = x_breaks) +
    labs(x = "Category", y = NULL, title = title) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), panel.spacing.x = unit(0.25, "lines"), panel.spacing.y = unit(0.6, "lines"),
          strip.background = element_rect(fill = "grey92", color = "white"), strip.text = element_text(size = 8),
          strip.placement = "outside", legend.position = "bottom", legend.key.width = unit(1.6, "cm"), legend.title = element_text(vjust = 0.8))
}

fig_heat_by_K <- function(K_val, a_decomp, g_decomp) {
  a_df <- prep_decomp(a_decomp)
  g_df <- prep_decomp(g_decomp)
  (plot_bias_heat(a_df, K_val, "(a) Arithmetic mean") / plot_bias_heat(g_df, K_val, "(b) Normalized geometric mean")) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom", plot.title = element_text(size = 9, face = "bold"))
}

for (K in c(3, 5, 10)){
  ggsave(sprintf("fig_bias_heat_K%d.pdf", K), fig_heat_by_K(K, a_decomp, g_decomp), width = 6.3, height = 7, device = cairo_pdf)
}

## 4. Explanation of the tails

rho <- 0.2; dx <- 0.5; K <- 10; pi_s <- 0.05
h_fun     <- function(phi) (phi + (1 - phi) * rho) / (phi * rho + (1 - phi))
h_assumed <- (1 + rho) * log(1 / rho) / (1 - rho) - 1
m_true <- h_fun(0.5) * dx; m_ass <- h_assumed * dx
f_pop <- function(u, m) pi_s * dnorm(u) + (1 - pi_s) * dnorm(u, m)
F_pop <- function(u, m) pi_s * pnorm(u) + (1 - pi_s) * pnorm(u, m)
get_probs <- function(bal, gamma = 2) {
  if (bal) return(rep(1 / K, K)); k <- 1:K; w <- (k * (K + 1 - k))^gamma; w / sum(w) }
dists <- c("Balanced" = TRUE, "Sparse tails" = FALSE)

u <- seq(-3.5, 4.5, length.out = 1600)

# Categories, thresholds and differences
cats <- bind_rows(lapply(names(dists), function(nm) {
  brk <- c(-Inf, qnorm(cumsum(get_probs(dists[[nm]]))[-K]), Inf)
  data.frame(dist = nm, k = 1:K, lo = brk[-(K + 1)], hi = brk[-1])
})) %>%
  mutate(dist = factor(dist, levels = names(dists)),
         diff = 100 * ((F_pop(hi, m_ass) - F_pop(lo, m_ass)) - (F_pop(hi, m_true) - F_pop(lo, m_true))),
         x = case_when(k == 1 ~ hi - 0.55, k == K ~ lo + 0.65, TRUE ~ (lo + hi) / 2),
         y = ifelse(k %% 2 == 0, 0.43, 0.405))

# Colours: within each panel, rank the under- and overestimated categories by size
# and spread them over the pink and blue halves of Benedictus (near-white middle dropped)
bene  <- met.brewer("Benedictus")
pinks <- rev(bene[1:5])   # light -> dark pink
blues <- bene[9:13]       # light -> dark blue

cats <- cats %>%
  mutate(under = diff <= 0,
         mag   = round(abs(diff), 1)) %>%        # same rounding as the printed labels
  group_by(under) %>%                            # rank across both panels together
  mutate(r = dense_rank(mag),
         fillcol = if (first(under)) colorRampPalette(pinks)(max(r))[r]
         else              colorRampPalette(blues)(max(r))[r]) %>%
  ungroup()

thr <- cats %>%
  filter(is.finite(hi)) %>%
  select(dist, tau = hi)

# One ribbon piece per category, cut exactly at the thresholds
dens <- cats %>%
  select(dist, k, lo, hi, fillcol) %>%
  group_by(dist, k, fillcol) %>%
  reframe(u = sort(unique(c(lo[is.finite(lo)],
                            .env$u[.env$u > lo & .env$u < hi],
                            hi[is.finite(hi)])))) %>%
  mutate(ass = f_pop(u, m_ass), true = f_pop(u, m_true))

curves <- data.frame(u = u, ass = f_pop(u, m_ass), true = f_pop(u, m_true)) %>%
  crossing(dist = factor(names(dists), levels = names(dists))) %>%
  pivot_longer(c(ass, true), names_to = "curve", values_to = "d")

p <- ggplot() +
  geom_ribbon(data = dens,
              aes(u, ymin = pmin(ass, true), ymax = pmax(ass, true),
                  fill = fillcol, group = k),
              alpha = 0.9) +
  geom_vline(data = thr, aes(xintercept = tau), colour = "grey55", linewidth = 0.3) +
  geom_line(data = curves, aes(u, d, linetype = curve), linewidth = 0.6, colour = "grey15") +
  geom_text(data = cats, aes(x, 0.465, label = k), size = 4, colour = "grey30") +
  geom_text(data = cats, aes(x, y, label = sprintf("%+.1f", diff), colour = diff <= 0),
            size = 4, show.legend = FALSE) +
  facet_wrap(~ dist) +
  scale_fill_identity(
    guide  = "legend",
    breaks = c(cats$fillcol[which.min(cats$diff)],   
               cats$fillcol[which.max(cats$diff)]),  
    labels = c("Underestimation", "Overestimation")
  ) +
  scale_colour_manual(values = c(`TRUE` = "#9a133d", `FALSE` = "#1a318b")) +
  scale_linetype_manual(values = c(ass = "dashed", true = "solid"),
                        labels = c("Assumed by the MUBP",
                                   expression("True (" * phi1 == 0.8 * ")"))) +
  scale_y_continuous(limits = c(0, 0.485), expand = c(0, 0)) +
  labs(x = "Latent variable U in the population", y = "Density",
       fill = NULL, linetype = NULL) +
  theme_bw(base_size = 20) +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        strip.background = element_rect(fill = "grey92", colour = NA))
p
ggsave("fig_area.png", p, width = 12, height = 6, dpi = 300)

