## ------------------------------------------------- ##
## Bayesian Inference for proportions with ME and SE ##
## Script 2: Empirical illustration                  ##
## Written by: Santiago Gómez-Echeverry              ##
## Last update: 05/05/2026                           ##
## ------------------------------------------------- ##

#### - (I) Working space and packages - ####
rm(list = ls())
packs <- c('tidyr', 'dplyr', 'rstan', 'posterior', 'ggplot2', 'ipumsr', 'viridis', 'kableExtra', 'purrr', 'broom.mixed', 'loo', 'bayesplot', 'patchwork')
ipacks <- packs %in% rownames(installed.packages())
if(any(ipacks == F)){install.packages(packs[!ipacks])}
invisible(lapply(packs, library, character.only = T))
sys <- Sys.info()
options(scipen = 999)
fold_data <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/3 - Dichotomous Multisource TE/3 - Data")
fold_graphs <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/3 - Dichotomous Multisource TE/4 - Graphs & Tables")

#### - (II) Data arrangement  - ####

# For this empirical illustration we will use IPUMS data via an API
#set_ipums_api_key("59cba10d8a5da536fc06b59d31e08525f7c34a8d8104321d30caa4a8", save = TRUE)

# We will use information from three sources: (i) the Current Population Survey (CPS), (ii) the American Time Use Survey (ATUS), and the (iii) American Community Survey (ACS)
# To use this data we first need to make the request and then extract the data
cps_req <- define_extract_micro(collection = "cps", description = "CPS Anchors", samples = "cps2023_03s", variables = c("CPSIDP", "INCTOT", "HHINCOME", "INCWAGE", "WKSWORK1", "HOURWAGE", "AGE"))
atus_req <- define_extract_micro(collection = "atus", description = "ATUS Target", samples = "at2023", variables = c("CPSIDP", "EARNWEEK"))
#acs_req <- define_extract_micro(collection = "usa", description = "ACS Benchmark", "us2023a", c("INCWAGE")) # This will be our benchmark

# Loading the data
cps_data  <- read_ipums_micro(download_extract(wait_for_extract(submit_extract(cps_req))))
atus_data <- read_ipums_micro(download_extract(wait_for_extract(submit_extract(atus_req))))
#acs_data  <- read_ipums_micro(download_extract(wait_for_extract(submit_extract(acs_req))))

# We will use the CPS-ATUS link data to create our indicators 
dat_emp <- inner_join(cps_data, atus_data, by = "CPSIDP") %>% 
  mutate(X = ifelse(INCWAGE < 40000, 0, ifelse((INCWAGE == 99999998 | INCWAGE == 99999999), NA, 1))) %>%               # Ground truth
  mutate(Y1 = ifelse(HHINCOME < 80000, 0, ifelse(HHINCOME == 99999999, NA, 1))) %>%                                    # CPS Proxy 1
  mutate(Y2 = ifelse(WKSWORK1*(HOURWAGE*40) < 40000, 0, ifelse((HOURWAGE == 99.99 | HOURWAGE == 999.99), NA, 1))) %>%  # CPS Proxy 2
  mutate(Y3 = ifelse((EARNWEEK * 52) < 40000, 0, ifelse(EARNWEEK == 99999.99, NA, 1))) %>%                             # ATUS target 
  mutate(R3 = ifelse(is.na(Y3), 0, 1)) %>% 
  mutate(Z = as.vector(scale(AGE))) %>% 
  filter(!is.na(X), !is.na(Y1), !is.na(Y2), !is.na(Z))   # Filtering: Ensure the observed 'Truth' and CPS indicators are present

# Check the number of missing values and the Se and Sp so that we report that
sum(!is.na(dat_emp$Y3))


# Now we need to adjust the data format so that we can pass it through Stan
N <- nrow(dat_emp)
K <- 3
C <- 2
q_idx <- 3
fq <- mean(dat_emp$R3)
Y_mat <- as.matrix(dat_emp[, c("Y1", "Y2", "Y3")])
Y_mat[is.na(Y_mat)] <- 0                           # Masked by 'obs'; Stan doesn't allow NA's!
obs_mat <- matrix(1, nrow = N, ncol = K)
obs_mat[, q_idx] <- dat_emp$R3
mu_z_low <- mean(dat_emp$Z[dat_emp$Y1 == 0])
mu_z_high <- mean(dat_emp$Z[dat_emp$Y1 == 1])
sd_z_low  <- sd(dat_emp$Z[dat_emp$Y1 == 0])
sd_z_high  <- sd(dat_emp$Z[dat_emp$Y1 == 1])

## True values
true_pi_pop <- mean(dat_emp$X)                # pi_X
true_pi_sam <- mean(dat_emp$X[dat_emp$R3==1]) # pi_XR
obs_pí <- mean(dat_emp$Y3, na.rm = T)         # pi_YR

true_em <- obs_pí - true_pi_sam       # Selection error
true_es <- true_pi_sam - true_pi_pop  # Measurement error
true_ec <- true_em + true_es          # Composite error
true_Se_q <- mean(dat_emp$Y3[dat_emp$X == 1], na.rm = T)
true_Sp_q <- mean(dat_emp$Y3[dat_emp$X == 0]==0, na.rm = T)
truth_vector <- c(true_pi_pop, true_Se_q, true_Sp_q,true_em, true_es, true_ec)

# Get the true etas
probit_mod <- glm(R3 ~ Z + X, data = dat_emp, family = binomial(link = "probit"))
true_eta0 <- coef(probit_mod)["(Intercept)"]
true_eta1 <- coef(probit_mod)["Z"]
true_eta2 <- coef(probit_mod)["X"]

# Get the true betas
beta0_prior_mean <- numeric(3)
beta1_prior_mean <- numeric(3)
prox <- paste0("Y",1:3)
for(k in 1:3){
  prox_nam <- prox[k]
  se_k <- mean(dat_emp[[prox_nam]][dat_emp$X==1], na.rm = T)
  sp_k <- mean(dat_emp[[prox_nam]][dat_emp$X==0]==0, na.rm = T)
  beta0_prior_mean[k] <- qlogis(1-sp_k)
  beta1_prior_mean[k] <- qlogis(se_k) - qlogis(1-sp_k)
}


# Stan data
stan_data <- list(N = N, K = K, C = C, Y = Y_mat, R = dat_emp$R3, q_idx = q_idx, fq = fq, Z = dat_emp$Z, obs = obs_mat, mu_z_prior = c(mu_z_low, mu_z_high),
                  sd_z_prior = c(sd_z_low, sd_z_high), beta0_prior_mean = beta0_prior_mean,
                  beta1_prior_mean = beta1_prior_mean)

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
  vector[C] mu_z_prior;
  vector[C] sd_z_prior;
  vector[K] beta0_prior_mean;
  vector[K] beta1_prior_mean;
} 

parameters { 
  simplex[C] pi; 
  vector[K] beta0; 
  vector<lower=0>[K] beta1;  
  real eta0; 
  real eta1; 
  real eta2;            
  ordered[C] mu_z;
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
prior_blocks <- list(
  # 1. BASE
  base = "
  pi[2] ~ beta(2, 5);
  beta0 ~ normal(0, 1);
  beta1 ~ normal(0, 1);
  eta0  ~ normal(0, 1);
  eta1  ~ normal(0.5, 1);
  eta2  ~ normal(0, 1);
  mu_z  ~ normal(0, 1);          
  sigma_z ~ exponential(1); 
  ", 
  # 2. SELECTION-CALIBRATED
  is = "
  pi[2] ~ beta(2, 5);
  beta0 ~ normal(0, 1);
  beta1 ~ normal(0, 1);
  eta0  ~ normal(-1.10, 0.1);
  eta1  ~ normal(-0.56, 0.1);
  eta2  ~ normal(3, 0.1);
  mu_z  ~ normal(0, 1);          
  sigma_z ~ exponential(1); 
  ",

  # 3. MEASUREMENT-CALIBRATED
  im = "
  pi[2] ~ beta(2, 5);
  beta0 ~ normal(beta0_prior_mean, 0.5);
  beta1 ~ normal(beta1_prior_mean, 0.5);
  eta0  ~ normal(0, 1);
  eta1  ~ normal(0.5, 1);
  eta2  ~ normal(0, 1);
  mu_z  ~ normal(0, 1);          
  sigma_z ~ exponential(1); 
  ",
  # 4. Selection- and Measurement-informed prior
  ism = "
  pi[2] ~ beta(2, 5);
  beta0 ~ normal(beta0_prior_mean, 0.5);
  beta1 ~ normal(beta1_prior_mean, 0.5);
  eta0  ~ normal(-1.10, 0.1);
  eta1  ~ normal(-0.56, 0.1);
  eta2  ~ normal(3, 0.1);
  mu_z  ~ normal(0, 1);          
  sigma_z ~ exponential(1); 
  "
)



# 2. Estimation
init_fun <- function() {
  list(pi = c(0.7, 0.3),
    # Start the chains near the anchored -0.96 and 2.79, with slight jitter
    beta0 = rep(-0.96, K) + rnorm(K, mean = 0, sd = 0.05),
    beta1 = rep(2.79, K) + rnorm(K, mean = 0, sd = 0.05), 
    eta0 = -0.7 + rnorm(1, 0, 0.05),
    eta1 = 0.4 + rnorm(1, 0, 0.05),
    eta2 = 0.3 + rnorm(1, 0, 0.05),
    mu_z = sort(stan_data$mu_z_prior),
    sigma_z = rep(1.0, C))
}

priors_nam <- names(prior_blocks)
fitted_models <- list()
for (p in priors_nam){
  cat("\nCompiling model for prior:", p, "\n")
  full_mod <- gsub("{{PRIORS}}", prior_blocks[[p]], stan_code, fixed = T)
  stan_file <- file.path(tempdir(), paste0("model_", p, ".stan"))
  writeLines(full_mod, stan_file)
  model <- stan_model(file = stan_file)
  fitted_models[[p]] <- sampling(model, data = stan_data, init = init_fun, chains = 4, cores = 3, iter = 4000, warmup = 1000, seed = 42,
                                 control = list(adapt_delta = 0.95, max_treedepth = 12), refresh = 1000)
}

#### - (IV) Graphs & Tables  - ####

### - Figures Appendix D: Convergence - ###

color_scheme_set("viridis")
pars_to_plot <- c("pi[1]", "pi[2]", "beta0[1]", "beta0[2]", "beta0[3]", "beta1[1]", "beta1[2]","beta1[3]", 
                  "eta0", "eta1", "eta2", "mu_z[1]", "mu_z[2]", "sigma_z[1]", "sigma_z[2]")

for (i in 1:length(fitted_models)){
  t_plots <- mcmc_trace(fitted_models[[1]], pars = pars_to_plot,facet_args = list(nrow = 3)) +
    theme_minimal() + theme(strip.text = element_text(size = 18)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  d_plots <- mcmc_dens_overlay(fitted_models[[1]], pars = pars_to_plot, facet_args = list(nrow = 3)) +
    theme_minimal() + theme(strip.text = element_text(size = 18)) +
    theme(legend.position = "right")
  combined <- t_plots / d_plots + plot_layout(guides = "collect")
  ggsave(paste0("Fig_conv_E", i, ".png"), plot = combined, path = fold_graphs, width = 25, height = 20, units = "cm")
  rm(t_plots, d_plots, combined)
}

### - Table 3 - ###

res_tables <- list()
for (i in 1:length(priors_nam)){
  res_tables[[i]] <- as.data.frame(summary(fitted_models[[i]], pars = c("pi[2]", "Se_q", "Sp_q", "EM", "ES", "EC"), probs = c(0.025, 0.50, 0.975))$summary)
  res_tables[[i]]$Parameter <- rownames(res_tables[[i]])
  res_tables[[i]]$Prior <- priors_nam[[i]]
}

final_res <- do.call(rbind.data.frame, res_tables)
colnames(final_res)[colnames(final_res) == "2.5%"] <- "P2.5"
colnames(final_res)[colnames(final_res) == "50%"] <- "P50"
colnames(final_res)[colnames(final_res) == "97.5%"] <- "P97.5"
final_res$Truth <- rep(truth_vector, times = length(priors_nam))

# 3. Clean up Parameter names and round
table_data <- final_res %>%
  mutate(Parameter = case_when(
    Parameter == "pi2" ~ "$\\pi_{2}$ ",
    Parameter == "Se_q"  ~ "$Se_q$",
    Parameter == "Sp_q"  ~ "$Sp_q$",
    Parameter == "EM"    ~ "$E_M$",
    Parameter == "ES"    ~ "$E_S$",
    Parameter == "EC"    ~ "$E_C$",
    TRUE ~ Parameter)) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>% 
  mutate(Rhat = round(Rhat, 4))

cols_to_show <- c("Parameter", "Truth", "mean", "sd", "P2.5", "P50", "P97.5", "n_eff", "Rhat")

table_data[, cols_to_show] %>%
  kable(format = "latex",booktabs = T, escape = F, longtable = F, linesep = "", row.names = F, 
        caption = "Posterior Estimates vs. Benchmark", col.names = c("Parameter", "Benchmark", "Mean", "SD", "P2.5","P50", "P97.5", "$n_{eff}$", "$\\hat{R}$")) %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), font_size = 10) %>%
  pack_rows(priors_nam[1], 1, 6, latex_gap_space = "0em") %>%
  pack_rows(priors_nam[2], 7, 12, latex_gap_space = "0em") %>%
  pack_rows(priors_nam[3], 13, 18, latex_gap_space = "0em") %>% 
  pack_rows(priors_nam[4], 19, 24, latex_gap_space = "0em")


### - Table 4 - ###

# 1. Extract the log_lik from each fitted model and compute LOO
loo_res <- list()
for (p in names(fitted_models)) {
  log_lik_matrix <- extract_log_lik(fitted_models[[p]], parameter_name = "log_lik")   # Extract log-likelihood matrix [iterations x observations]
  loo_res[[p]] <- loo(log_lik_matrix)  # Calculate PSIS-LOO
  cat("\nLOO Diagnostics for", p, ":\n")
  print(loo_res[[p]])
}

# 2. Compare all models against each other
loo_comp <- loo_compare(loo_res)
m_labs <- c("base" = "$\\mathcal{M}_{0'}$", "is" = "$\\mathcal{M}_{S}$", "im" = "$\\mathcal{M}_{M}$", "ism" = "$\\mathcal{M}_{C}$")
table_loo <- as.data.frame(loo_comp) %>%
  mutate(Model = m_labs[model]) %>%
  mutate(across(where(is.numeric), ~ round(., digits = 3))) %>%
  select(Model, elpd_loo, se_elpd_loo, elpd_diff, se_diff, p_worse, diag_diff) %>%
  kbl(format = "latex", booktabs = TRUE, escape = F, col.names = c("Prior", "elpd", "$\\mathrm{SE}_{elpd}$","$\\Delta$ ELPD", "$\\mathrm{SE}_{\\Delta}$", "p-worse", "diag_{\\Delta}"),
    caption = "Model Comparison via PSIS-LOO Cross-Validation", label = "4", align = "lcccc") %>%
  kable_styling(latex_options = c("hold_position"))
table_loo



