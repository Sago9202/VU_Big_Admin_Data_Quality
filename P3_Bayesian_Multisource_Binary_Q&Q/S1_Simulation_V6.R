## ------------------------------------------------- ##
## Bayesian Inference for proportions with ME and SE ##
## Script 1: Simulation Study                        ##
## Written by: Santiago Gómez-Echeverry              ##
## Last update: 17/09/2026                           ##
## ------------------------------------------------- ##

#### - (I) Working space and packages - ####
rm(list = ls())
packs <- c('tidyr', 'dplyr', 'rstan', 'posterior', 'bayesplot', 'reshape', 'ggplot2', 'priorsense', 'MASS', 'viridis', 'patchwork', 'kableExtra')
ipacks <- packs %in% rownames(installed.packages())
if(any(ipacks == F)){install.packages(packs[!ipacks])}
invisible(lapply(packs, library, character.only = T))
sys <- Sys.info()
fold_data <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/3 - Dichotomous Multisource TE/3 - Data")
fold_graphs <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/3 - Dichotomous Multisource TE/4 - Graphs & Tables")

# Stan specifications
rstan_options(auto_write = T)
options(mc.cores = parallel::detectCores())

#### - (II) Data generation  - #### 

set.seed(42)
D <- 5
N <- 200    # No. of observations
pi_X <- 0.1  # Latent target proportion
K <- 3       # No. of items
q_idx <- 2   # Item from the non-prob sample
seeds <- sample.int(1e6, D)

# 1. Condition generation

# Define scenario levels
scenarios <- expand.grid(me_lev =c("L", "M", "H"), se_lev = c("L","M", "H"), stringsAsFactors = F)
sim_data <- list()

# Function to compute delta
compute_deltaR <- function(eta, theta_q_Se, theta_q_Sp, N = 20000, intercept = -1, betaZ = 0.5, mu0 = -0.5, mu1 = 0.5, sigma0 = 1, sigma1 = 1) {
  X <- rbinom(N, 1, 0.4)
  Yq <- ifelse(X == 1, rbinom(N, 1, theta_q_Se), rbinom(N, 1, 1 - theta_q_Sp))
  Z <- ifelse(X == 1, rnorm(N, mean = mu1, sd = sigma1), rnorm(N, mean = mu0, sd = sigma0))
  pR <- pnorm(intercept + betaZ * Z + eta * X)
  Rq <- rbinom(N, 1, pR)
  return(mean(Rq[Yq==1]) - mean(Rq[Yq==0]))
}

# Calibration to find delta*
find_eta_for_deltaR <- function(target_delta, theta_q_Se, theta_q_Sp,intercept = -1, betaZ = 0.5,lower = -10, upper = 10, tol = 1e-3, mu0 = -0.5, mu1 = 0.5, sigma0 = 1, sigma1 = 1) {
  # Wrapper for uniroot
  f <- function(eta) compute_deltaR(eta, theta_q_Se, theta_q_Sp, intercept = intercept, betaZ = betaZ,
                                    mu0 = mu0, mu1 = mu1, sigma0 = sigma0, sigma1 = sigma1) - target_delta
  
  # Try uniroot; if bracketing fails, adjust delta
  safe_uniroot <- function(f, lower, upper) {
    left <- f(lower); right <- f(upper)
    if (is.na(left) | is.na(right) | left*right > 0) {
      # Not bracketed, then reduce delta automatically
      warning("Target delta not achievable. Adjusting to closest possible value.")
      # Find approximate achievable delta using a large eta grid
      eta_grid <- seq(-10, 10, length.out = 200)
      deltas <- sapply(eta_grid, function(e) compute_deltaR(e, theta_q_Se, theta_q_Sp))
      closest <- eta_grid[which.min(abs(deltas - target_delta))]
      return(closest)
    }
    # Otherwise safe
    return(uniroot(f, lower = lower, upper = upper, tol = tol)$root)
  }
  eta_root <- safe_uniroot(f, lower, upper)
  return(eta_root)
}

# Error conditions
ME_targets <- list(Low  = list(Se = 0.56, Sp = 0.55), Mod  = list(Se = 0.90, Sp = 0.85), High = list(Se = 0.85, Sp = 0.85))
SE_targets <- list(Low  = 0.03, Mod = 0.06, High = 0.10)
scenarios <- expand.grid(me_lev = c("Low", "Mod", "High"), se_lev = c("Low", "Mod", "High"), stringsAsFactors = F)

calibrated <- data.frame(Scenario=character(), Se=numeric(), Sp=numeric(), Eta=numeric(), Realized_DeltaR=numeric(), stringsAsFactors=F)

for (s in 1:nrow(scenarios)) {
  me <- scenarios$me_lev[s]; se <- scenarios$se_lev[s]
  Se <- ME_targets[[me]]$Se
  Sp <- ME_targets[[me]]$Sp
  target_delta <- SE_targets[[se]]
  eta_cal <- find_eta_for_deltaR(target_delta, Se, Sp)
  realized_deltaR <- compute_deltaR(eta_cal, Se, Sp)
  
  calibrated[s, ] <- list(Scenario = paste0("ME_", me, "_SE_", se), Se = Se, Sp = Sp, Eta = eta_cal, Realized_DeltaR = realized_deltaR)
}

print(calibrated) # Check the calibration

# Let's export the calibrated values for the table in the paper
calibrated[,1] <- 1:nrow(scenarios)
calibrated[,4:5] <- round(calibrated[,4:5],3)
Tab_Cal <- kable(calibrated, format = "latex")
Tab_Cal

# 2. Data generation

mu0_true <- -0.5
mu1_true <- 0.5
sigma0_true <- 1.0
sigma1_true <- 1.0

for (s in 1:nrow(scenarios)){
  # True (Population) values 
  theta_vals <- matrix(c(0.15, (1-calibrated$Sp[s]), 0.20, 0.85, calibrated$Se[s], 0.80), ncol = 2, byrow = F)
  eta2_true <- calibrated$Eta[s] # Selection error level
  Se_q <- theta_vals[q_idx, 2]; Sp_q <- calibrated$Sp[s]; J_q <- Se_q + Sp_q -1 
  intercept <- -1.0; betaZ <- 0.5 
  
  pR_X1 <- pnorm((intercept + betaZ*mu1_true + eta2_true) / sqrt(1 + betaZ^2 * sigma1_true^2))
  pR_X0 <- pnorm((intercept + betaZ*mu0_true) / sqrt(1 + betaZ^2 * sigma0_true^2))
  fq_t <- pR_X1*pi_X + pR_X0*(1-pi_X)
  pX_q_t <- pi_X*pR_X1 / fq_t
  
  pX1_Yq1 <- pi_X*Se_q / (pi_X*Se_q + (1-pi_X)*(1-Sp_q))
  pX1_Yq0 <- pi_X*(1-Se_q) / (pi_X*(1-Se_q) + (1-pi_X)*Sp_q)
  Delta_RY_t <- (pX1_Yq1 - pX1_Yq0) * (pR_X1 - pR_X0)  
  pYq_t <- pi_X*Se_q + (1-pi_X)*(1-Sp_q)                   
  sigma2_Yq_t <- pYq_t * (1 - pYq_t)
  
  EM_t <- 1 - Sp_q - pX_q_t*(1 - J_q)             
  ES_t <- Delta_RY_t*(1/fq_t)*sigma2_Yq_t*(1/J_q)  
  EC_t <- EM_t + ES_t                                 
  
  sim_data[[s]] <- vector("list", D)
  
  for (d in 1:D){
    set.seed(seeds[d])
    X <- rbinom(N, 1, pi_X)               # Latent target variable
    Y <- matrix(NA, nrow = N, ncol = K)   # Observed target variable
    for (k in 1:K){
      Y[,k] <- ifelse(X==1, rbinom(N, 1, theta_vals[k,2]), rbinom(N, 1, theta_vals[k,1]))
    }
    Z <- ifelse(X == 1, rnorm(N, mean = mu1_true, sd = sigma1_true), rnorm(N, mean = mu0_true, sd = sigma0_true))
    # MNAR process (Y2)
    pR <- pnorm(-1.0+0.5*Z+eta2_true*X)
    Rq <- rbinom(N, 1, pR)
    fq <- mean(Rq)
    Y_obs <- Y
    p_miss <- 0.3
    Y_obs[Rq==0, q_idx] <- NA
    
    # MAR process (Y2^{-})
    for (i in setdiff(1:K, q_idx)){
      miss_k <- rbinom(N, 1, p_miss)
      Y_obs[miss_k ==1, i] <- NA
    }
    
    # We need a placeholder for the NA's since Stan doesn't accept missing values
    obs <- ifelse(is.na(Y_obs), 0, 1)
    Y_obs[is.na(Y_obs)] <- 0  
    sim_data[[s]][[d]] <- list(stan_data = list(N = N, K = K, C = 2,  Y = Y_obs,  R = Rq, q_idx = q_idx, fq = fq, Z = Z, obs = obs),
                          true_values = list(EM = EM_t, ES = ES_t, EC = EC_t, eta2 = eta2_true, theta = theta_vals),
                          label = paste0("ME_", scenarios$me_lev[s], "_SE_", scenarios$se_lev[s], "_d_", d)) 
  }
}

#### - (III) Stan model and estimation  - ####

# 1. Stan model
stan_code <- " 
data { 
  int<lower=1> N; 
  int<lower=1> K; 
  int<lower=1> C; 
  int<lower=0, upper=1> Y[N, K]; 
  int<lower=0, upper=1> R[N]; 
  int<lower=1, upper=K> q_idx; 
  real<lower=0, upper=1> fq; 
  vector[N] Z; 
  int<lower=0, upper=1> obs[N, K]; 
} 

parameters { 
  simplex[C] pi; 
  vector[K] beta0; 
  vector<lower=0>[K] beta1;  
  real eta0; 
  real eta1; 
  real eta2;            
  vector[C] mu_z;
  vector<lower=0>[C] sigma_z;
} 

transformed parameters { 
  matrix[K, 2] theta;  
  for (k in 1:K) { 
    // Class 1 (X=0) success probability 
    theta[k, 1] = inv_logit(beta0[k]); 
    // Class 2 (X=1) success probability, forced higher by beta1 > 0 
    theta[k, 2] = inv_logit(beta0[k] + beta1[k]); 
  } 
} 

model { 
  // 1. Priors: 
  {{PRIORS}} 
  
  // Pre-calculate linear predictors for selection to avoid redundant math 
  vector[N] mu_base = eta0 + eta1 * Z; 
  vector[N] mu_y1 = mu_base + eta2; 
  
  // 2. Likelihood 
  for (i in 1:N) { 
    vector[C] log_theta_classes; 
    for (c in 1:C) { 
      real lp = log(pi[c]); 
      
      lp += normal_lpdf(Z[i] | mu_z[c], sigma_z[c]);
      
      // Measurement model: Only use items where obs == 1 
      for (k in 1:K) { 
        if (obs[i, k] == 1)  
          lp += bernoulli_lpmf(Y[i, k] | theta[k, c]); 
      } 
      
      real mu_sel = (c == 1) ? mu_base[i] : mu_y1[i]; 
      if (R[i] == 1) { 
        lp += std_normal_lcdf(mu_sel); 
      } else { 
        lp += std_normal_lccdf(mu_sel); 
      } 
        
      log_theta_classes[c] = lp; 
    } 
    // Increment total log-probability with marginalized mixture likelihood 
    target += log_sum_exp(log_theta_classes); 
  } 
} 

generated quantities { 
  vector[N] log_lik;
  real Se_q = theta[q_idx, 2]; 
  real Sp_q = 1.0 - theta[q_idx, 1]; 
  real J_q = Se_q + Sp_q - 1.0; 
  
  vector[N] p_x1_given_y; 
  { 
    vector[N] mu_base = eta0 + eta1 * Z; 
    vector[N] mu_y1 = mu_base + eta2; 
    
    for (i in 1:N) { 
      vector[2] logp; 
      for (c in 1:2) { 
        logp[c] = log(pi[c]); 
        
        logp[c] += normal_lpdf(Z[i] | mu_z[c], sigma_z[c]);
        
        for (k in 1:K) { 
          if (obs[i, k] == 1) logp[c] += bernoulli_lpmf(Y[i, k] | theta[k, c]); 
        } 
        
        real mu_sel = (c == 1) ? mu_base[i] : mu_y1[i]; 
        if (R[i] == 1) { 
          logp[c] += std_normal_lcdf(mu_sel); 
        } else { 
          logp[c] += std_normal_lccdf(mu_sel); 
        } 
      } 
      log_lik[i] = log_sum_exp(logp);
      p_x1_given_y[i] = softmax(logp)[2]; 
    } 
  } 
  
  real piX_given_q = 0.0; 
  int n_q = sum(R); 
  
  if (n_q > 0) { 
    for (i in 1:N) { 
      if (R[i] == 1) { 
        piX_given_q += p_x1_given_y[i]; 
      } 
    } 
    piX_given_q /= n_q; 
  } 
  
  real sigma2_X = pi[2] * (1.0 - pi[2]); 
  real EM = 1.0 - Sp_q - piX_given_q * (1.0 - J_q);  
  real Delta_R_X = Phi((eta0 + eta1 * mu_z[2] + eta2) / sqrt(1.0 + square(eta1) * square(sigma_z[2]))) - 
                   Phi((eta0 + eta1 * mu_z[1]) / sqrt(1.0 + square(eta1) * square(sigma_z[1]))); 
  real ES = (Delta_R_X * sigma2_X) / fq; 
  real EC = EM + ES; 
} 
"  
# We will consider three different prior for sensitivity analyses
prior_blocks <- list(
  weak = "
  pi      ~ dirichlet(rep_vector(0.5, C));
  beta0   ~ normal(-1.4, 1.5);
  beta1   ~ normal(3.0, 1.5);
  eta0    ~ normal(-1, 1);
  eta1    ~ normal(0.5, 1);
  eta2    ~ normal(0, 1.5);
  mu_z    ~ normal(0, 1);
  sigma_z ~ exponential(1);
  "
  # 
  # base = "
  # pi      ~ dirichlet(rep_vector(1.0, C));
  # beta0   ~ normal(-1.4, 1);
  # beta1   ~ normal(3.0, 1);
  # eta0    ~ normal(-1, 0.75);
  # eta1    ~ normal(0.5, 0.75);
  # eta2    ~ normal(0, 1);
  # mu_z    ~ normal(0, 1);
  # sigma_z ~ exponential(1);
  # "
  #,
  # 
  # strong = "
  # pi      ~ dirichlet(rep_vector(2.0, C));
  # beta0   ~ normal(-1.5, 0.5);
  # beta1   ~ normal(3.5, 0.7);
  # eta0    ~ normal(-1, 0.5);
  # eta1    ~ normal(0.5, 0.5);
  # eta2    ~ normal(0, 0.5);
  # mu_z[1] ~ normal(-0.5, 0.5);
  # mu_z[2] ~ normal(0.5, 0.5);
  # sigma_z ~ exponential(2);
  # "
)
# 2. Estimation

diag_pars <- c("pi", "beta0", "beta1", "eta0", "eta1", "eta2", "mu_z", "sigma_z",
               "EM", "ES", "EC", "Se_q", "Sp_q", "J_q")
report_pars <- c("EM", "ES", "EC", "Se_q", "Sp_q", "J_q")  # subset we attach true values to

max_treedepth_ctrl <- 10  # default max_treedepth; update if you pass a different `control` to sampling()
n_chains <- 3; n_iter <- 4000; n_warmup <- 1000
n_postwarmup_draws <- n_chains * (n_iter - n_warmup)

results_list <- list()       # Posterior summaries for EM/ES/EC/Se_q/Sp_q/J_q, with true values attached
diagnostics_list <- list()   # One row per fit: divergences, treedepth, E-BFMI, worst Rhat/ESS across ALL diag_pars
checkpoint_folder <- file.path(fold_data, "checkpoints")
if(!dir.exists(checkpoint_folder)) dir.create(checkpoint_folder, recursive = T)

for(p in names(prior_blocks)) {
  cat("\nCompiling model for prior:", p, "\n")
  this_code <- gsub("{{PRIORS}}", prior_blocks[[p]], stan_code, fixed = T)
  stan_file <- file.path(tempdir(), paste0("model_", p, ".stan"))
  writeLines(this_code, stan_file)
  model <- stan_model(file = stan_file)
  
  for(s in 1:length(sim_data)) {
    
    for(r in 1:D) {
      curr_data <- sim_data[[s]][[r]]
      cat("  Estimating Scenario:", curr_data$label, "\n")
      fit <- sampling(model, data = curr_data$stan_data, chains = n_chains, cores = n_chains,
                      iter = n_iter, warmup = n_warmup, seed = 42, refresh = 500)
      
      # --- Parameter-level summary (mean, sd, quantiles, Rhat, bulk/tail ESS) via posterior ---
      draws_arr <- as.array(fit, pars = diag_pars)
      param_diag <- posterior::summarise_draws(draws_arr)
      
      # Keep only the error-decomposition parameters for the results table, with true values attached
      fit_summary <- as.data.frame(param_diag[param_diag$variable %in% report_pars, ])
      fit_summary$prior <- p
      fit_summary$scenario <- curr_data$label
      fit_summary$draw <- r
      fit_summary$true_EM <- curr_data$true_values$EM
      fit_summary$true_ES <- curr_data$true_values$ES
      fit_summary$true_EC <- curr_data$true_values$EC
      results_list[[length(results_list) + 1]] <- fit_summary
      
      # --- Fit-level HMC diagnostics: divergences, treedepth saturation, E-BFMI ---
      sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
      n_divergent <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
      n_max_treedepth <- sum(sapply(sp, function(x) sum(x[, "treedepth__"] >= max_treedepth_ctrl)))
      ebfmi_chain <- sapply(sp, function(x) {
        e <- x[, "energy__"]
        mean(diff(e)^2) / var(e)
      })
      min_ebfmi <- min(ebfmi_chain)
      
      diag_row <- data.frame(
        prior = p, scenario = curr_data$label, draw = r,
        n_divergent = n_divergent,
        pct_divergent = n_divergent / n_postwarmup_draws,
        n_max_treedepth = n_max_treedepth,
        min_ebfmi = min_ebfmi,
        max_rhat = max(param_diag$rhat, na.rm = T),
        min_ess_bulk = min(param_diag$ess_bulk, na.rm = T),
        min_ess_tail = min(param_diag$ess_tail, na.rm = T)
      )
      diagnostics_list[[length(diagnostics_list) + 1]] <- diag_row
      
      # Flag problems immediately in the console so a long run doesn't hide issues
      if(n_divergent > 0 || diag_row$max_rhat > 1.01 || min_ebfmi < 0.2 || diag_row$min_ess_bulk < 400) {
        cat("    [!] Convergence flag -- divergent:", n_divergent,
            "| max_treedepth hits:", n_max_treedepth,
            "| max Rhat:", round(diag_row$max_rhat, 3),
            "| min E-BFMI:", round(min_ebfmi, 3),
            "| min ESS_bulk:", round(diag_row$min_ess_bulk, 0), "\n")
      }
      
      rm(fit, draws_arr, sp)
    }
    
    # Checkpoint after each scenario is fully estimated for this prior
    saveRDS(dplyr::bind_rows(results_list), file.path(checkpoint_folder, paste0("results_checkpoint_", p, "_scenario", s, ".rds")))
    saveRDS(dplyr::bind_rows(diagnostics_list), file.path(checkpoint_folder, paste0("diagnostics_checkpoint_", p, "_scenario", s, ".rds")))
  }
}

results_all <- dplyr::bind_rows(results_list)
diagnostics_all <- dplyr::bind_rows(diagnostics_list)
saveRDS(results_all, file.path(fold_data, "simulation_results_all_draws.rds"))
saveRDS(diagnostics_all, file.path(fold_data, "simulation_diagnostics_all_draws.rds"))

# --- Convergence summary table, aggregated across draws within each prior x scenario ---
# Meant to be reported alongside the main results (analogous to Table 2/App. E in the paper)
convergence_summary <- diagnostics_all %>%
  group_by(prior, scenario) %>%
  summarise(
    n_fits = n(),
    pct_fits_with_divergences = mean(n_divergent > 0),
    mean_pct_divergent_iters = mean(pct_divergent),
    pct_fits_maxtreedepth = mean(n_max_treedepth > 0),
    pct_fits_rhat_over_1.01 = mean(max_rhat > 1.01),
    mean_max_rhat = mean(max_rhat),
    pct_fits_low_ebfmi = mean(min_ebfmi < 0.2),
    mean_min_ess_bulk = mean(min_ess_bulk),
    mean_min_ess_tail = mean(min_ess_tail),
    .groups = "drop"
  )

print(convergence_summary)
Tab_Conv <- kable(convergence_summary, format = "latex", digits = 3)
Tab_Conv

# Store the results! -> This is key since they take a while to run
setwd(fold_data)
saveRDS(fitted_models, file = "Fitted_Models_120626.rds")

# In case we don't want to run, but just load the results instead
fitted_models <- readRDS("Fitted_Models_120626.rds")

#### - (IV) Graphs and Tables - ####

### - Figures Appendix D: Convergence - ###

color_scheme_set("viridis")
fit_base <- fitted_models[[2]]
pars_to_plot <- c("pi[1]", "pi[2]", "beta0[1]", "beta0[2]", "beta0[3]", "beta1[1]", "beta1[2]","beta1[3]", 
                  "eta0", "eta1", "eta2", "mu_z[1]", "mu_z[2]", "sigma_z[1]", "sigma_z[2]")

for (i in 1:length(fit_base)){
  t_plots <- mcmc_trace(fit_base[[i]], pars = pars_to_plot,facet_args = list(nrow = 3)) +
    theme_minimal() + theme(strip.text = element_text(size = 18)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  d_plots <- mcmc_dens_overlay(fit_base[[i]], pars = pars_to_plot, facet_args = list(nrow = 3)) +
    theme_minimal() + theme(strip.text = element_text(size = 18)) +
    theme(legend.position = "right")
  combined <- t_plots / d_plots + plot_layout(guides = "collect")
  ggsave(paste0("Fig_conv_S", i, ".png"), plot = combined, path = fold_graphs, width = 25, height = 20, units = "cm")
  rm(t_plots, d_plots, combined)
}

### - Figure 1  - ###

# For this plot we need a lot of prep. We need to get a the posteriors and made several precalculations

fit_base <- fitted_models[[2]]
all_post <- lapply(seq_along(fit_base), function(i) {
  fit <- fit_base[[i]]
  draws <- as_draws_df(fit, variable = c("EM", "ES"))
  draws$Scenario <- sim_data[[i]]$label
  draws$EM_true <- sim_data[[i]]$true_values$EM
  draws$ES_true <- sim_data[[i]]$true_values$ES
  draws
})
df_post <- bind_rows(all_post)

# Compute local density for coloring
comp_dens <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x, y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  dens$z[cbind(ix, iy)]
}

df_post <- df_post %>%
  group_by(Scenario) %>%
  mutate(density = comp_dens(EM, ES), density = (density - min(density)) / (max(density) - min(density))) %>%
  ungroup()

# Posterior EC
df_post <- df_post %>%
  mutate(EC = EM + ES)

# Posterior mean point per scenario
df_ec_post <- df_post %>%
  group_by(Scenario) %>%
  summarise(EM_bar = median(EM), ES_bar = median(ES), .groups = "drop")

# Compute EC line per scenario
df_ec <- df_post %>%
  group_by(Scenario) %>%
  summarise(EC = unique(EM_true + ES_true), .groups = "drop") %>%
  rowwise() %>%
  mutate(EM_seq = list(seq(min(df_post$EM[df_post$Scenario==Scenario]), max(df_post$EM[df_post$Scenario==Scenario]), length.out = 100)), ES_seq = list(EC - EM_seq)) %>%
  unnest(cols = c(EM_seq, ES_seq))

# Distance to true EC line: ES = EC_true - EM
df_post <- df_post %>%
  group_by(Scenario) %>%
  mutate(EC_true = unique(EM_true + ES_true),EC_dist = abs(ES - (EC_true - EM))) %>%
  ungroup()

corr_df <- do.call(rbind,by(df_post, df_post$Scenario, function(d) {data.frame(Scenario = unique(d$Scenario), rho = cor(d$EM, d$ES))}))
corr_df$label <- sprintf("rho == %.2f", corr_df$rho)
pos_df <- data.frame(Scenario = unique(df_post$Scenario), x = min(df_post$EM), y = min(df_post$ES))
corr_df <- merge(corr_df, pos_df, by = "Scenario")

facet_labels <- c("ME_Low_SE_Low" = 'Se :0.95~\"/\"~Delta^"*" :0.06', "ME_Low_SE_Mod" = 'Se :0.95~\"/\"~Delta^"*" :\"0.09\"', "ME_Low_SE_High" = 'Se :0.95~\"/\"~Delta^"*" :0.12', 
                  "ME_Mod_SE_Low" = 'Se :0.90~\"/\"~Delta^"*" :0.06', "ME_Mod_SE_Mod" = 'Se :0.90~\"/\"~Delta^"*" :\"0.09\"', "ME_Mod_SE_High" = 'Se :0.90~\"/\"~Delta^"*" :0.12', 
                  "ME_High_SE_Low" = 'Se :0.85~\"/\"~Delta^"*" :0.06', "ME_High_SE_Mod" = 'Se :0.85~\"/\"~Delta^"*" :\"0.09\"', "ME_High_SE_High" = 'Se :0.85~\"/\"~Delta^"*" :0.12')
# Get the error settings' order
ordered_scenarios <- c("ME_Low_SE_Low",  "ME_Low_SE_Mod",  "ME_Low_SE_High","ME_Mod_SE_Low", "ME_Mod_SE_Mod", "ME_Mod_SE_High", "ME_High_SE_Low", "ME_High_SE_Mod", "ME_High_SE_High")
df_post$Scenario <- factor(df_post$Scenario, levels = ordered_scenarios)

# Apply to your auxiliary dataframes so they match the facets
df_ec$Scenario <- factor(df_ec$Scenario, levels = ordered_scenarios)
df_ec_post$Scenario <- factor(df_ec_post$Scenario, levels = ordered_scenarios)
corr_df$Scenario <- factor(corr_df$Scenario, levels = ordered_scenarios)

# Contour Plot
p1 <- ggplot(df_post, aes(x = EM, y = ES)) +
  geom_point(aes(color = density), alpha = 0.6, size = 1.2) +
  scale_color_viridis(option = "viridis", direction = 1, name = "Posterior density") +
  geom_point(aes(x = EM_true, y = ES_true, shape = "True"),color = "red3",size = 3) +
  geom_line(data = df_ec,aes(x = EM_seq, y = ES_seq,linetype = "Composite error identified set"), 
            color = "red3",linewidth = 0.8) +
  geom_point(data = df_ec_post,aes(x = EM_bar, y = ES_bar, shape = "Posterior median"), color = "gray17", size = 3,stroke = 1.2) +
  facet_wrap(~ Scenario,labeller = as_labeller(facet_labels, label_parsed)) +
  labs(x = "Measurement error", y = "Selection error") +
  scale_shape_manual(name = "Composite error", values = c( "True" = 17,"Posterior median" = 4)) +
  scale_linetype_manual(name = "Lines",values = c("Composite error identified set" = "dashed")) +
  guides(shape = guide_legend(override.aes = list(color = c("gray17", "red3")))) +
  theme_minimal(base_size = 15) +
  geom_text(data = corr_df, aes(x = x, y = y, label = label), parse = T, hjust = 0, vjust = 0, size = 6, color = "black") +
  theme(strip.text = element_text(face = "bold", size = 30), panel.grid = element_line(color = "gray85"), 
        legend.title = element_text(size = 18), legend.text = element_text(size = 15), axis.title = element_text(size = 30),
        axis.text=element_text(size= 25))
p1
ggsave(paste0("Fig_1.png"), plot = p1, path = fold_graphs, width = 50, height = 30, units = "cm")

### - Figure 2  - ###

# Get the posteriors of EC
ec_summ <- df_post %>%
  group_by(Scenario) %>%
  summarise(EC_true = unique(EM_true + ES_true), EC_med  = median(EC), EC_lo95 = quantile(EC, 0.025),
            EC_hi95 = quantile(EC, 0.975),EC_lo80 = quantile(EC, 0.10),EC_hi80 = quantile(EC, 0.90),.groups = "drop")

xlims <- df_post %>%
  group_by(Scenario) %>%
  summarise(xmin = quantile(EC, 0.005), xmax = quantile(EC, 0.995),.groups = "drop")

df_post <- df_post %>%
  left_join(xlims, by = "Scenario")

p2 <- ggplot(df_post, aes(x = EC)) +
  geom_density(fill = "#20A387FF", color = NA, alpha = 0.45, adjust = 1) +
  geom_segment(data = ec_summ, aes(x = EC_lo95, xend = EC_hi95, y = 0, yend = 0), linewidth = 1.1, color = "gray40") +
  geom_segment(data = ec_summ, aes(x = EC_lo80, xend = EC_hi80, y = 0, yend = 0), linewidth = 3, color = "#238A8DFF") +
  geom_vline(data = ec_summ, aes(xintercept = EC_med, color = "Posterior Median"), linewidth = 1) +
  geom_vline(data = ec_summ, aes(xintercept = EC_true, color = "True Value"), linetype = "dashed", linewidth = 1.2) +
  scale_color_manual(name = "Error", values = c("Posterior Median" = "black", "True Value" = "firebrick4"),
                     guide = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
  facet_wrap(~ Scenario, labeller = as_labeller(facet_labels, label_parsed), scales = "free_x") +
  coord_cartesian(xlim = c(min(df_post$xmin), max(df_post$xmax))) +
  labs(x = expression("Composite error"), y = "Posterior density") +
  theme_minimal(base_size = 15) +
  theme(axis.title = element_text(size = 30), strip.text = element_text(face = "bold", size = 30), panel.grid = element_blank(),
        axis.text=element_text(size= 25))
p2
ggsave(paste0("Fig_2.png"), plot = p2, path = fold_graphs, width = 50, height = 30, units = "cm")

### - Table 1 - ###
pars_interest <- c("Se_q","eta2","EM","ES","EC")
results <- list()
indx <- 1
for(p in names(fitted_models)){
  for(s in seq_along(fitted_models[[p]])){
    fit <- fitted_models[[p]][[s]]
    summ <- summary(fit, pars = pars_interest)$summary
    df <- data.frame(Scenario = sim_data[[s]]$label, Prior = p, Parameter = pars_interest, Mean = summ[pars_interest,"mean"], 
                     LCI = summ[pars_interest,"2.5%"], UCI = summ[pars_interest,"97.5%"])
    results[[indx]] <- df
    indx <- indx+1
  }
}
results <- bind_rows(results)
truth_list <- list()
for(s in seq_along(sim_data)){
  t_vals <- sim_data[[s]]$true_values
  theta <- t_vals$theta
  truth_df <- data.frame(Scenario = sim_data[[s]]$label, Parameter = pars_interest, Truth = c(theta[2,2], t_vals$eta2, t_vals$EM, t_vals$ES, t_vals$EC))
  truth_list[[s]] <- truth_df
}

truth <- bind_rows(truth_list)
results <- left_join(results, truth, by=c("Scenario","Parameter"))

table_wide <- results %>%
  pivot_wider(names_from = Prior, values_from = c(Mean,LCI,UCI))
table_wide$Estimate_weak <- sprintf("%.3f [%.3f, %.3f]", table_wide$Mean_weak, table_wide$LCI_weak, table_wide$UCI_weak)
table_wide$Estimate_base <- sprintf("%.3f [%.3f, %.3f]", table_wide$Mean_base, table_wide$LCI_base, table_wide$UCI_base)
table_wide$Estimate_informative <- sprintf("%.3f [%.3f, %.3f]", table_wide$Mean_strong, table_wide$LCI_strong, table_wide$UCI_strong)
table_wide$Truth <- sprintf("%.3f", table_wide$Truth)

table_wide <- table_wide %>%
  dplyr::select(!dplyr::starts_with("Mean") & !dplyr::starts_with("LCI") & !dplyr::starts_with("UCI"))

# Alternative with more control over formatting
kable(table_wide, format = "latex", booktabs = T, longtable = F,  linesep = "", digits = 3, col.names = c("Scenario", "Parameter", "True value", "Weakly Inf.", "Base", "Strongly Inf."), 
      align = c("l", "l", "c", "c", "c", "c"), caption = "Posterior estimates under different prior specifications", label = "tab:posterior_estimates") %>%
  add_header_above(c(" " = 3, "Posterior mean [95% credible interval]" = 3), bold = T, line = T) %>%
  collapse_rows(columns = 1, latex_hline = "major",  valign = "middle") %>%
  kable_styling(latex_options = c("hold_position", "scale_down"), font_size = 10, bootstrap_options = c("condensed")) %>%
  row_spec(0, bold = T) %>%
  row_spec(0, extra_latex_after = "\\midrule") %>%
  column_spec(2, italic = TRUE) %>%
  column_spec(4:6, width = "3cm")

### - Figure C1: ME and SE relationship - ###

## Effect of J_q on Em
posterior_pi <- function(p, Se, Sp){
  (p * Se) / (p * Se + (1 - p) * (1 - Sp))
}

# Partial derivative wrt Se
partial_pi <- function(p, Se, Sp){
  D <- p * Se + (1 - p) * (1 - Sp)
  p * (1 - p) * (1 - Sp) / D^2
}

# Grid of parameters
p_vals  <- c(0.2, 0.4, 0.6)
Sp_vals <- c(0.6, 0.8, 0.9)
Se_vals <- seq(0.3, 0.95, by = 0.01)

# Build full simulation grid
df <- expand.grid(Se = Se_vals,p  = p_vals,Sp = Sp_vals) %>%
  mutate(J = Se + Sp - 1, direct = posterior_pi(p, Se, Sp), dpi_dJ = partial_pi(p, Se, Sp),indirect = (1 - J) * dpi_dJ, total_deriv = direct - indirect) %>%
  pivot_longer(cols = c(direct, indirect, total_deriv), names_to = "component", values_to = "value") %>%
  mutate(component = factor(component,levels = c("direct", "indirect", "total_deriv"),labels = c("Direct effect", "Indirect effect", "Total derivative")),
         p = factor(p, levels = p_vals,labels = paste0("p = ", p_vals)),Sp = factor(Sp, levels = Sp_vals, labels = paste0("Sp = ", Sp_vals)))

pC1 <- ggplot(df, aes(x = J, y = value, color = component, linetype = component)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
  facet_grid(p ~ Sp) +
  scale_color_manual(values = c("Direct effect" = "red4","Indirect effect" = "blue4","Total derivative" = "black")) +
  scale_linetype_manual(values = c("Direct effect" = "dashed","Indirect effect" = "dotted","Total derivative" = "solid")) +
  labs(x = expression(J[q]),y = expression(paste(partialdiff, epsilon[M], "/", partialdiff, J[q])),color = "Component",linetype = "Component") +
  theme_minimal(base_size = 14) + theme(text = element_text(size =20), legend.position = "bottom",strip.text = element_text(face = "bold"))
pC1
ggsave(paste0("Fig_C1.png"), plot = pC1, path = fold_graphs, width = 30, height = 30, units = "cm")

#############


cov <- results_all %>%
  filter(variable %in% c("EM", "ES", "EC")) %>%
  mutate(
    true = case_when(variable == "EM" ~ true_EM,
                     variable == "ES" ~ true_ES,
                     variable == "EC" ~ true_EC),
    covered  = true >= q5 & true <= q95,
    variable = factor(variable, levels = c("EM", "ES", "EC"))
  ) %>%
  group_by(scenario, variable) %>%
  summarise(n = n(), k = sum(covered), .groups = "drop") %>%
  mutate(coverage = k / n,
         # exact (Clopper–Pearson) 95% interval for the coverage itself
         lo = qbeta(0.025, k, n - k + 1),
         hi = qbeta(0.975, k + 1, n - k))

cov_table <- cov %>%
  mutate(cell = sprintf("%.2f [%.2f, %.2f]", coverage, lo, hi)) %>%
  dplyr::select(scenario, variable, cell) %>%
  pivot_wider(names_from = variable, values_from = cell)

strip_lab <- cov %>%
  arrange(scenario, variable) %>%
  group_by(scenario) %>%
  summarise(lab = paste0(pretty_scen(first(scenario)), "\n",
                         "90% CrI coverage: ",
                         paste(sprintf("%s %d/%d", variable, k, n), collapse = " · ")),
            .groups = "drop")
strip_lab <- setNames(strip_lab$lab, strip_lab$scenario)

ggplot(wide, aes(EM, ES)) +
  geom_abline(data = truth, aes(intercept = true_EC, slope = -1),
              colour = cols["True value"], linetype = "22", linewidth = 0.4) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", bins = 4) +
  geom_point(colour = est_col, size = 0.7, alpha = 0.35) +
  geom_point(data = marks, aes(colour = type, shape = type),
             size = 2.8, stroke = 1.2) +
  geom_xsidedensity(fill = fill_col, colour = NA, alpha = 0.7) +
  geom_ysidedensity(fill = fill_col, colour = NA, alpha = 0.7) +
  geom_xsidevline(data = marks, aes(xintercept = EM, colour = type),
                  linewidth = 0.6, show.legend = FALSE) +
  geom_ysidehline(data = marks, aes(yintercept = ES, colour = type),
                  linewidth = 0.6, show.legend = FALSE) +
  facet_wrap(~ scenario, labeller = as_labeller(strip_lab)) +
  scale_fill_gradient(low = "#F1F5FA", high = "#8DB0D3", guide = "none") +
  scale_colour_manual(values = cols, name = NULL) +
  scale_shape_manual(values = shapes, name = NULL) +
  scale_xsidey_continuous(breaks = NULL) +
  scale_ysidex_continuous(breaks = NULL) +
  labs(x = "EM", y = "ES",
       caption = "Dashed line: true EC (EM + ES)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor        = element_blank(),
    panel.grid.major        = element_line(colour = "grey95"),
    panel.border            = element_rect(colour = "grey85", fill = NA),
    strip.text              = element_text(hjust = 0, size = 9, lineheight = 1.1),
    legend.position         = "top",
    legend.justification    = "left",
    plot.caption            = element_text(colour = "grey50", size = 8.5),
    ggside.panel.scale      = 0.15,
    ggside.panel.border     = element_blank(),
    ggside.panel.grid       = element_blank(),
    ggside.panel.background = element_blank()
  )
