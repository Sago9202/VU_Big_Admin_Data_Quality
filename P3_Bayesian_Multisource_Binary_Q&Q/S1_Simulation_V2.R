## ---------------------------------------- ##
## Modeling Measurement and Selection Error ##
## Script 1 - Simulation 1                  ##
## Written by: Santiago Gómez-Echeverry     ##
## Last modified: 08/03/2026                ##
## ---------------------------------------- ##

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
N <- 200     # No. of observations
pi_X <- 0.4  # Latent target proportion
K <- 3       # No. of items
q_idx <- 2   # Item from the non-prob sample
# Define scenario levels
scenarios <- expand.grid(me_lev =c("Low", "High"), se_lev = c("Low", "High"), stringsAsFactors = F)
sim_data <- list()

compute_deltaR <- function(eta, theta_q_Se, theta_q_Sp, N = 20000, intercept = -1, betaZ = 0.5) {
  X <- rbinom(N, 1, 0.4)
  Yq <- ifelse(X == 1, rbinom(N, 1, theta_q_Se), rbinom(N, 1, 1 - theta_q_Sp))
  Z <- rnorm(N)
  pR <- plogis(intercept + betaZ * Z + eta * Yq)
  Rq <- rbinom(N, 1, pR)
  return(mean(Rq[X==1]) - mean(Rq[X==0]))
}

# Find eta robustly
find_eta_for_deltaR <- function(target_delta, theta_q_Se, theta_q_Sp,intercept = -1, betaZ = 0.5,lower = -10, upper = 10, tol = 1e-3) {
  # Wrapper for uniroot
  f <- function(eta) compute_deltaR(eta, theta_q_Se, theta_q_Sp, 
                                    intercept = intercept, betaZ = betaZ) - target_delta
  
  # Try uniroot; if bracketing fails, adjust delta
  safe_uniroot <- function(f, lower, upper) {
    left <- f(lower); right <- f(upper)
    if (is.na(left) | is.na(right) | left*right > 0) {
      # Not bracketed → reduce delta automatically
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

# Example: targets
ME_targets <- list(Low  = list(Se = 0.95, Sp = 0.85),
                   High = list(Se = 0.85, Sp = 0.85))
SE_targets <- list(Low  = 0.05, High = 0.15)

scenarios <- expand.grid(me_lev = c("Low", "High"),
                         se_lev = c("Low", "High"),
                         stringsAsFactors = FALSE)

calibrated <- data.frame(Scenario=character(),
                         Se=numeric(), Sp=numeric(), Eta=numeric(), Realized_DeltaR=numeric(),
                         stringsAsFactors=FALSE)

for (s in 1:nrow(scenarios)) {
  me <- scenarios$me_lev[s]; se <- scenarios$se_lev[s]
  Se <- ME_targets[[me]]$Se
  Sp <- ME_targets[[me]]$Sp
  target_delta <- SE_targets[[se]]
  eta_cal <- find_eta_for_deltaR(target_delta, Se, Sp)
  realized_deltaR <- compute_deltaR(eta_cal, Se, Sp)
  calibrated[s, ] <- list(Scenario = paste0("Meas_", me, "_Sel_", se), Se = Se, Sp = Sp, 
                          Eta = eta_cal, Realized_DeltaR = realized_deltaR)
}

print(calibrated)
calibrated[,1] <- 1:nrow(scenarios)
calibrated[,4:5] <- round(calibrated[,4:5],3)
Tab_Cal <- kable(calibrated, format = "latex")
Tab_Cal

for (s in 1:nrow(scenarios)){
  theta_vals <- matrix(c(0.15, (1-calibrated$Sp[s]), 0.20, 0.85, calibrated$Se[s], 0.80), ncol = 2, byrow = F)
  
  # Selection error level
  eta2_true <- calibrated$Eta[s]
  
  # Remainder of the data generation
  X <- rbinom(N, 1, pi_X)               # Latent target variable
  Y <- matrix(NA, nrow = N, ncol = K)   # Observed target variable
  for (k in 1:K){
    Y[,k] <- ifelse(X==1, rbinom(N, 1, theta_vals[k,2]), rbinom(N, 1, theta_vals[k,1]))
  }
  #Z <- rnorm(N, 0, 1)  # Should we adjust this variable so that we have a higher correlation with X?
  rho_XZ <- 0.7
  X_std <- (X-mean(X))/sd(X)
  Z <- rho_XZ*X_std+sqrt(1-rho_XZ^2)*rnorm(N)
  pR <- plogis(-1.0+0.5*Z+eta2_true*Y[,q_idx])
  Rq <- rbinom(N, 1, pR)
  fq <- mean(Rq)
  
  Y_obs <- Y
  p_miss <- 0.3
  Y_obs[Rq==0, q_idx] <- NA
  
  for (i in setdiff(1:K, q_idx)){
    miss_k <- rbinom(N, 1, p_miss)
    Y_obs[miss_k ==1, i] <- NA
  }
  
  obs <- ifelse(is.na(Y_obs), 0, 1)
  Y_obs[is.na(Y_obs)] <- 0  # placeholder
  
  # Let's define the true values for later comparison
  Se_q <- theta_vals[q_idx, 2]; Sp_q <- calibrated$Sp[s]; J_q <- Se_q + Sp_q - 1
  Yq <- Y[, q_idx]
  p1 <- mean(Rq[X==1]); p0 <- mean(Rq[X==0])
  Delta_RY <- p1 - p0
  piX_given_q <- mean(X[Rq==1])
  fq <- mean(Rq)
  piYq <- mean(Yq[Rq==1])
  sigma2_Yq <- piYq*(1-piYq)
  EM_t <- 1 - Sp_q - piX_given_q*(1- J_q)
  ES_t <- Delta_RY*(1/fq)*sigma2_Yq*(1/J_q)
  EC_t <- EM_t + ES_t
  sim_data[[s]] <- list(stan_data = list(N = N, K = K, C = 2,  Y = Y_obs,  R = Rq, q_idx = q_idx, fq = fq, Z = Z, obs = obs),
                        true_values = list(EM = EM_t, ES = ES_t, EC = EC_t, eta2 = eta2_true, theta = theta_vals),
                        label = paste0("Meas_", scenarios$me_lev[s], "_Sel_", scenarios$se_lev[s]))
}

sim_data[[1]]$true_values$EM; sim_data[[2]]$true_values$EM;sim_data[[3]]$true_values$EM;sim_data[[4]]$true_values$EM 
sim_data[[1]]$true_values$ES; sim_data[[2]]$true_values$ES;sim_data[[3]]$true_values$ES;sim_data[[4]]$true_values$ES 


#### - (III) Stan model and estimation  - ####

stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> C;
  int<lower=0, upper=1> Y[N,K];
  int<lower=0, upper=1> R[N];
  int<lower=1, upper=K> q_idx;
  real<lower=0, upper=1> fq;
  vector[N] Z;
  int<lower=0, upper=1> obs[N, K];
}

parameters {
  simplex[C] pi;
  vector[K] alpha;
  real<lower=0> gamma[K];
  real eta0;
  real eta1;
  real eta2;         
}

transformed parameters {
  real<lower=0, upper=1> theta[K,2];
  for (k in 1:K) {
    theta[k,1] = inv_logit(alpha[k]);
    theta[k,2] = inv_logit(alpha[k] + gamma[k]);
  }

  real Se_q = theta[q_idx, 2];
  real Sp_q = 1.0 - theta[q_idx, 1];
  real J_q = Se_q + Sp_q - 1.0;
}

model {
  // 1. Priors:
  {{PRIORS}}
  
  // 2. Likelihood
  for (i in 1:N) {
    vector[C] logp;
    for (c in 1:C) {
      real lp = log(pi[c]);
      
      // Measurement model: Only use observed data
      for (k in 1:K) {
        if (obs[i,k] == 1) 
          lp += bernoulli_lpmf(Y[i,k] | theta[k,c]);
      }

      // Selection model (MNAR Logic)
      if (R[i] == 1) {
        lp += bernoulli_lpmf(1 | Phi(eta0 + eta1 * Z[i] + eta2 * Y[i, q_idx]));
      } else {
        // Marginalize over the missing Y[q_idx]
        real logpY1 = bernoulli_lpmf(1 | theta[q_idx, c]);
        real logpY0 = bernoulli_lpmf(0 | theta[q_idx, c]);
        real logpR0_y1 = bernoulli_lpmf(0 | Phi(eta0 + eta1 * Z[i] + eta2 * 1));
        real logpR0_y0 = bernoulli_lpmf(0 | Phi(eta0 + eta1 * Z[i] + eta2 * 0));
        lp += log_sum_exp(logpY1 + logpR0_y1, logpY0 + logpR0_y0);
      }
      logp[c] = lp;
    }
    target += log_sum_exp(logp);
  }
}

generated quantities {
  vector[N] log_lik;
  real p_x1_given_y[N];
  
  for (i in 1:N) {
    vector[C] logp;
    for (c in 1:C) {
      real lp = log(pi[c]);
      // Measurement part (Fixed: checking 'obs' flag)
      for (k in 1:K) {
        if (obs[i,k] == 1) 
          lp += bernoulli_lpmf(Y[i, k] | theta[k, c]);
      }
      // Selection part
      if (R[i] == 1) {
        lp += bernoulli_lpmf(1 | Phi(eta0 + eta1 * Z[i] + eta2 * Y[i, q_idx]));
      } else {
        real logpY1 = bernoulli_lpmf(1 | theta[q_idx, c]);
        real logpY0 = bernoulli_lpmf(0 | theta[q_idx, c]);
        real logpR0_y1 = bernoulli_lpmf(0 | Phi(eta0 + eta1 * Z[i] + eta2 * 1));
        real logpR0_y0 = bernoulli_lpmf(0 | Phi(eta0 + eta1 * Z[i] + eta2 * 0));
        lp += log_sum_exp(logpY1 + logpR0_y1, logpY0 + logpR0_y0);
      }
      logp[c] = lp;
    }
    log_lik[i] = log_sum_exp(logp);
    p_x1_given_y[i] = softmax(logp)[2]; // Probability of Class 2
  }

  // 3. Bias Terms Calculation
  real piX_given_q = 0.0;
  real piY_given_q = 0.0;
  int n_q = sum(R);

  if (n_q > 0) {
    for (i in 1:N) {
      if (R[i] == 1) {
        piX_given_q += p_x1_given_y[i];
        piY_given_q += Y[i, q_idx];
      }
    }
    piX_given_q /= n_q;
    piY_given_q /= n_q;
  }

  real sigma2_Yq = piY_given_q * (1.0 - piY_given_q);
  real EM = 1.0 - Sp_q - piX_given_q * (1.0 - J_q);

  real sum_delta = 0.0;
  for (i in 1:N) {
    sum_delta += (Phi(eta0 + eta1 * Z[i] + eta2 * 1) - 
                  Phi(eta0 + eta1 * Z[i] + eta2 * 0));
  }
  real Delta_R_Y = sum_delta / N;
  real ES = (Delta_R_Y * sigma2_Yq) / (fq * J_q);
  real EC = EM + ES;
} 
" 

#stan_file <- file.path(tempdir(), "stan_code.stan")
#writeLines(stan_code, con = stan_file)
#cat("Stan model written to:", stan_file, "\n")

#stan_model_rstan <- stan_model(file = stan_file)
#fit_rstan <- sampling(stan_model_rstan, data = data_list, seed = 42, chains = 4, iter = 5000, warmup = 1000)

prior_blocks <- list(
  # weak = "
  # pi    ~ dirichlet(rep_vector(0.5, C));
  # alpha ~ normal(0, 1.2);
  # gamma ~ normal(0, 1.2);
  # eta0  ~ normal(0, 1.2);
  # eta1  ~ normal(0, 1.2);
  # eta2  ~ normal(0, 1.2);
  # ",
  # 
  baseline = "
  pi    ~ dirichlet(rep_vector(1.0, C));
  alpha ~ normal(0, 1);
  gamma ~ normal(0.5, 1);
  eta0  ~ normal(0, 1);
  eta1  ~ normal(0, 1);
  eta2  ~ normal(0, 1);
  ",
  
  # informative = "
  # pi    ~ dirichlet(rep_vector(3.0, C));
  # alpha ~ normal(0, 0.5);
  # gamma ~ normal(1, 0.5);
  # eta0  ~ normal(0, 0.5);
  # eta1  ~ normal(0, 0.5);
  # eta2  ~ normal(0.5, 0.5);
  # "
)


fitted_models <- list()

for(p in names(prior_blocks)) {
  
  cat("\nCompiling model for prior:", p, "\n")
  
  this_code <- gsub("{{PRIORS}}", prior_blocks[[p]], stan_code, fixed = TRUE)
  
  stan_file <- file.path(tempdir(), paste0("model_", p, ".stan"))
  writeLines(this_code, stan_file)
  model <- stan_model(file = stan_file)
  
  fitted_models[[p]] <- list()
  
  for(s in 1:length(sim_data)) {
    
    curr_data <- sim_data[[s]]
    cat("  Estimating Scenario:", curr_data$label, "\n")
    
    fitted_models[[p]][[s]] <- sampling(
      model,
      data = curr_data$stan_data,
      chains = 4,
      iter = 4000,
      warmup = 1000,
      seed = 42,
      refresh = 500
    )
  }
}


#### - (IV) Results  - ####

## Convergence
# Before jumping into analyzing the results we should be certain that the model is converging. The easiest way to evaluate this is check the traceplots and the density plots.
# Another option is to also report the R_hat measures

fit_baseline <- fitted_models[[2]]

pars_to_plot <- c("pi[1]", "pi[2]", "alpha[1]", "alpha[2]", "gamma[1]", "gamma[2]", "eta0", "eta1")

for (i in 1:length(fit_baseline)){
 t_plots <- mcmc_trace(fit_baseline[[i]], pars = pars_to_plot,facet_args = list(nrow = 2)) +
    theme_minimal() + theme(strip.text = element_text(size = 18)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  d_plots <- mcmc_dens_overlay(fit_baseline[[i]], pars = pars_to_plot, facet_args = list(nrow = 2)) +
    theme_minimal() + theme(strip.text = element_text(size = 18)) +
    theme(legend.position = "right")
  combined <- t_plots / d_plots + plot_layout(guides = "collect")
  ggsave(paste0("Fig_conv", i, ".png"), plot = combined, path = fold_graphs, width = 25, height = 20, units = "cm")
  rm(t_plots, d_plots, combined)
}

## Contour plots

# Get the posterior draws along with the true values
all_post <- lapply(seq_along(fit_baseline), function(i) {
  fit <- fit_baseline[[i]]
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
  mutate(density = comp_dens(EM, ES),
         density = (density - min(density)) / (max(density) - min(density))) %>%
  ungroup()

# Posterior EC
df_post <- df_post %>%
  mutate(EC = EM + ES)

# Posterior mean point per scenario
df_ec_post <- df_post %>%
  group_by(Scenario) %>%
  summarise(
    EM_bar = median(EM),
    ES_bar = median(ES),
    .groups = "drop")


# Compute EC line per scenario
df_ec <- df_post %>%
  group_by(Scenario) %>%
  summarise(EC = unique(EM_true + ES_true), .groups = "drop") %>%
  rowwise() %>%
  mutate(EM_seq = list(seq(min(df_post$EM[df_post$Scenario==Scenario]),
                           max(df_post$EM[df_post$Scenario==Scenario]), length.out = 100)),
         ES_seq = list(EC - EM_seq)) %>%
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

facet_labels <- c(
  "Meas_Low_Sel_Low"   = 'Se :0.95~\"/\"~Delta^"*" :0.05',
  "Meas_Low_Sel_High"  = 'Se :0.95~\"/\"~Delta^"*" :0.15',
  "Meas_High_Sel_Low"  = 'Se :0.85~\"/\"~Delta^"*" :0.05',
  "Meas_High_Sel_High" = 'Se :0.85~\"/\"~Delta^"*" :0.15')
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
  geom_text(data = corr_df, aes(x = x, y = y, label = label), parse = TRUE, hjust = 0, vjust = 0, size = 6, color = "black") +
  theme(strip.text = element_text(face = "bold", size = 20),
        panel.grid = element_line(color = "gray85"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15))
p1
ggsave(paste0("Fig_1.png"), plot = p1, path = fold_graphs, width = 35, height = 20, units = "cm")

## Summary of the results

# 1. Initialize an empty list to collect results
results_list <- list()

for (i in 1:length(fit_baseline)) {
  
  # A. Get the Stan summary for just the bias components
  # This avoids the N-length variables (y_rep, log_lik)
  fit_sum <- summary(fit_baseline[[i]], pars = c("EM", "ES", "EC"))$summary
  fit_df  <- as.data.frame(fit_sum)
  
  # B. Add the parameter names as a column
  fit_df$parameter <- rownames(fit_df)
  
  # C. Match with the True values from your sim_data list
  trues <- sim_data[[i]]$true_values
  # We map them: EM -> EM_t, ES -> ES_t, EC -> EC_t
  fit_df$true_val <- c(trues$EM, trues$ES, trues$EC)
  
  # D. Add metadata and calculate Bias
  fit_df$scenario <- sim_data[[i]]$label
  fit_df$bias     <- fit_df$mean - fit_df$true_val
  
  # E. Add a "Coverage" flag: Does the 95% interval contain the truth?
  fit_df$covered  <- ifelse(fit_df$true_val >= fit_df$`2.5%` & 
                              fit_df$true_val <= fit_df$`97.5%`, 1, 0)
  
  # Store this scenario's results
  results_list[[i]] <- fit_df
}

# 2. Combine all scenarios into one master data frame
final_results_table <- do.call(rbind, results_list)

# 3. View the top of your results
head(final_results_table)


# Focus on the 'High/High' scenario (usually index 4)
worst_case <- fitted_models[[4]]



summ <- summarise_draws(
  worst_case,
  variables = c("pi[1]","pi[2]","alpha[1]","alpha[2]",
                "gamma[1]","gamma[2]","eta0","eta1","eta2")
)

summ[, c("variable","rhat","ess_bulk","ess_tail")]
# Trace plots for the selection parameters
p3 <- mcmc_trace(worst_case, pars = c("pi[1]","pi[2]","alpha[1]","alpha[2]","gamma[1]","gamma[2]","eta0", "eta1", "eta2"))
p3

traceplot(worst_case, pars = c("pi[1]","pi[2]","alpha[1]","alpha[2]","gamma[1]","gamma[2]","eta0", "eta1", "eta2"))

ggsave(paste0("Fig_3.png"), plot = p3, path = fold_graphs, width = 20, height = 20, units = "cm")

# Posterior density vs True Values
mcmc_areas(worst_case, pars =  c("pi[1]","pi[2]","alpha[1]","alpha[2]","gamma[1]","gamma[2]","eta0", "eta1", "eta2"), prob = 0.95) +
  geom_vline(data = data.frame(parameter =  c("pi[1]","pi[2]","alpha[1]","alpha[2]","gamma[1]","gamma[2]","eta0", "eta1", "eta2"), 
                               vals = unlist(sim_data[[4]]$true_values[1:3])),
             aes(xintercept = vals), linetype = "dashed", color = c("red","blue", "green") )

# Prior sensitivity analysis
sensitivity_tab <- powerscale_sensitivity(fitted_models[[4]])
print(sensitivity_tab)


param_names <- names(rstan::extract(fitted_models[[4]]))
# Filter out generated quantities so you only pass model parameters:
param_names <- setdiff(param_names, c("log_lik","lp__", "lprior", "p_x1_given_y", "y_rep_obs","y_rep_unobs", "piX_given_q", "piY_given_q",
                                      "EM","ES","EC", "eps_J", "n_q", "sum_delta"))

powerscale_plot_quantities(fitted_models[[4]],log_lik_name = "log_lik",log_prior_name = "lprior",variables = param_names)




####
ec_summ <- df_post %>%
  group_by(Scenario) %>%
  summarise(EC_true = unique(EM_true + ES_true), EC_med  = median(EC), EC_lo95 = quantile(EC, 0.025),
            EC_hi95 = quantile(EC, 0.975),EC_lo80 = quantile(EC, 0.10),EC_hi80 = quantile(EC, 0.90),.groups = "drop")


xlims <- df_post %>%
  group_by(Scenario) %>%
  summarise(
    xmin = quantile(EC, 0.005),
    xmax = quantile(EC, 0.995),
    .groups = "drop"
  )

df_post <- df_post %>%
  left_join(xlims, by = "Scenario")

p_ec <- ggplot(df_post, aes(x = EC)) +
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
  theme(axis.title = element_text(size = 20), 
        strip.text = element_text(face = "bold", size = 20), 
        panel.grid = element_blank(),)

p_ec
ggsave(paste0("Fig_2.png"), plot = p_ec, path = fold_graphs, width = 35, height = 20, units = "cm")

# Parameters' table

pars_interest <- c("Se_q","Sp_q","eta2","EM","ES","EC")
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
  true_vals <- sim_data[[s]]$true_values
  theta <- true_vals$theta
  truth_df <- data.frame(
    Scenario = sim_data[[s]]$label,
    Parameter = pars_interest,
    Truth = c(
      theta[2,2],          # Se_q
      1-theta[2,1],        # Sp_q
      true_vals$eta2,
      true_vals$EM,
      true_vals$ES,
      true_vals$EC
    )
  )
  
  truth_list[[s]] <- truth_df
}

truth <- bind_rows(truth_list)
results <- left_join(results, truth, by=c("Scenario","Parameter"))

table_wide <- results %>%
  pivot_wider(names_from = Prior, values_from = c(Mean,LCI,UCI))
table_wide$Estimate_baseline <- sprintf("%.3f [%.3f, %.3f]", table_wide$Mean_baseline, table_wide$LCI_baseline, table_wide$UCI_baseline)
table_wide$Estimate_weak <- sprintf("%.3f [%.3f, %.3f]", table_wide$Mean_weak, table_wide$LCI_weak, table_wide$UCI_weak)
table_wide$Estimate_informative <- sprintf("%.3f [%.3f, %.3f]", table_wide$Mean_informative, table_wide$LCI_informative, table_wide$UCI_informative)
table_wide$Truth <- sprintf("%.3f", table_wide$Truth)

table_wide <- table_wide %>%
  dplyr::select(!dplyr::starts_with("Mean") & !dplyr::starts_with("LCI") & !dplyr::starts_with("UCI"))

kable(table_wide,
      format="latex",
      booktabs=T,
      col.names=c("Scenario","Parameter","True value",
                  "Baseline prior","Weak prior","Informative prior")
) %>%
  add_header_above(c(
    " " = 5,
    "Posterior mean [95\\% credible interval]" = 3
  ))


#### - (IV) Measurement and selection error relationship  - ####

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

# Plot
pA1 <- ggplot(df, aes(x = J, y = value, color = component, linetype = component)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
  facet_grid(p ~ Sp) +
  scale_color_manual(values = c("Direct effect" = "red4","Indirect effect" = "blue4","Total derivative" = "black")) +
  scale_linetype_manual(values = c("Direct effect" = "dashed","Indirect effect" = "dotted","Total derivative" = "solid")) +
  labs(x = expression(J[q]),y = expression(paste(partialdiff, epsilon[M], "/", partialdiff, J[q])),color = "Component",linetype = "Component") +
  theme_minimal(base_size = 14) + theme(text = element_text(size =20), legend.position = "bottom",strip.text = element_text(face = "bold"))
pA1
ggsave(paste0("Fig_A1.png"), plot = pA1, path = fold_graphs, width = 30, height = 30, units = "cm")




