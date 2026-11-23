## ------------------------------------------------- ##
## Bayesian Inference for proportions with ME and SE ##
## Script 2: Empirical illustration                  ##
## Written by: Santiago Gómez-Echeverry              ##
## Last update: 05/05/2026                           ##
## ------------------------------------------------- ##

#### - (I) Working space and packages - ####
rm(list = ls())
packs <- c('tidyr', 'dplyr', 'rstan', 'posterior', 'ggplot2', 'ipumsr', 'viridis', 'kableExtra', 'purrr', 'broom.mixed', 'loo')
ipacks <- packs %in% rownames(installed.packages())
if(any(ipacks == F)){install.packages(packs[!ipacks])}
invisible(lapply(packs, library, character.only = T))
sys <- Sys.info()
options(scipen = 999)
fold_data <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/3 - Dichotomous Multisource TE/3 - Data")
fold_graphs <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/3 - Dichotomous Multisource TE/4 - Graphs & Tables")

#### - (II) Data arrangement  - ####

# For this empirical illustration we will use IPUMS data via an API
#set_ipums_api_key("", save = TRUE)

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
  mutate(X = ifelse(INCWAGE < 60000, 0, ifelse((INCWAGE == 99999998 | INCWAGE == 99999999), NA, 1))) %>%               # Ground truth
  mutate(Y1 = ifelse(HHINCOME < 100000, 0, ifelse(HHINCOME == 99999999, NA, 1))) %>%                                   # CPS Proxy 1
  mutate(Y2 = ifelse(WKSWORK1*(HOURWAGE*40) < 60000, 0, ifelse((HOURWAGE == 99.99 | HOURWAGE == 999.99), NA, 1))) %>%  # CPS Proxy 2
  mutate(Y3 = ifelse((EARNWEEK * 52) < 60000, 0, ifelse(EARNWEEK == 99999.99, NA, 1))) %>%                             # ATUS target 
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

stan_data <- list(N = N, K = K, C = C, Y = Y_mat, R = dat_emp$R3, q_idx = q_idx, fq = fq, Z = dat_emp$Z, obs = obs_mat)

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

#### - (III) Stan model and estimation  - ####

# 1. Stan model
stan_code <- "
data { 
  int<lower=1> N; 
  int<lower=1> K; 
  int<lower=1> C; 
  int<lower=0, upper=1> Y[N, K]; 
  int<lower=0, upper=1> R[N]; 
  int<lower=1> q_idx; 
  real<lower=0, upper=1> fq; 
  vector[N] Z; 
  int<lower=0, upper=1> obs[N, K]; 
} 

parameters { 
  simplex[C] pi; 
  vector[K] alpha; 
  vector<lower=0>[K] gamma;  
  real eta0; 
  real eta1; 
  real<lower=0> eta2;            
} 

transformed parameters { 
  matrix[K, 2] theta;  
  for (k in 1:K) { 
    theta[k, 1] = inv_logit(alpha[k]); 
    theta[k, 2] = inv_logit(alpha[k] + gamma[k]); 
  } 
} 

model { 
  {{PRIORS}}

  vector[N] mu_base = eta0 + eta1 * Z; 
  vector[N] mu_y1 = mu_base + eta2; 
  
  for (i in 1:N) { 
    vector[C] log_theta_classes; 
    for (c in 1:C) { 
      real lp = log(pi[c]); 
      
      for (k in 1:K) { 
        if (obs[i, k] == 1)  
          lp += bernoulli_lpmf(Y[i, k] | theta[k, c]); 
      } 
      
      if (R[i] == 1) { 
        lp += std_normal_lcdf(mu_base[i] + eta2 * Y[i, q_idx]); 
      } else { 
        real lp_y1 = bernoulli_lpmf(1 | theta[q_idx, c]) + std_normal_lccdf(mu_y1[i]); 
        real lp_y0 = bernoulli_lpmf(0 | theta[q_idx, c]) + std_normal_lccdf(mu_base[i]); 
        lp += log_sum_exp(lp_y1, lp_y0); 
      } 
      log_theta_classes[c] = lp; 
    } 
    target += log_sum_exp(log_theta_classes); 
  } 
} 

generated quantities { 
  vector[N] log_lik;
  real Se_q = theta[q_idx, 2]; 
  real Sp_q = 1.0 - theta[q_idx, 1]; 
  real J_q = Se_q + Sp_q - 1.0; 
  
  vector[N] p_x1_given_y; 
  vector[N] mu_base_g = eta0 + eta1 * Z; 
  vector[N] mu_y1_g = mu_base_g + eta2; 

  for (i in 1:N) { 
    vector[2] logp; 
    for (c in 1:2) { 
      logp[c] = log(pi[c]); 
      for (k in 1:K) { 
        if (obs[i, k] == 1) 
          logp[c] += bernoulli_lpmf(Y[i, k] | theta[k, c]); 
      } 
      if (R[i] == 1) { 
        logp[c] += std_normal_lcdf(mu_base_g[i] + eta2 * Y[i, q_idx]); 
      } else { 
        real lp_y1 = bernoulli_lpmf(1 | theta[q_idx, c]) + std_normal_lccdf(mu_y1_g[i]); 
        real lp_y0 = bernoulli_lpmf(0 | theta[q_idx, c]) + std_normal_lccdf(mu_base_g[i]); 
        logp[c] += log_sum_exp(lp_y1, lp_y0); 
      } 
    } 
    log_lik[i] = log_sum_exp(logp);
    p_x1_given_y[i] = softmax(logp)[2]; 
  } 

  real piX_given_q = 0.0; 
  int n_q = sum(R); 
  if (n_q > 0) { 
    for (i in 1:N) 
      if (R[i] == 1) 
        piX_given_q += p_x1_given_y[i]; 
    piX_given_q /= n_q; 
  } 

  real piY_pop   = pi[2] * Se_q + pi[1] * (1.0 - Sp_q);
  real sigma2_Yq = piY_pop * (1.0 - piY_pop);
  real EM = 1.0 - Sp_q - piX_given_q * (1.0 - J_q); 
  real Delta_R_Y = Phi((eta0 + eta2) / sqrt(1.0 + eta1^2)) - Phi(eta0 / sqrt(1.0 + eta1^2)); 
  real ES = (Delta_R_Y * sigma2_Yq) / (fq * J_q); 
  real EC = EM + ES; 
} 
"
prior_blocks <- list(
  # 1. BASE: Keep as is - it's your control group that shows the bias.
  base = "
  pi[2] ~ beta(2, 5);
  alpha ~ normal(0, 1);
  gamma ~ normal(0.5, 1);
  eta0  ~ normal(0, 1);
  eta1  ~ normal(0.5, 1);
  eta2  ~ normal(0, 1);
  ",
  
  # 2. SELECTION-CALIBRATED: Aiming for ES = 0.149
  is = "
  pi[2] ~ beta(2, 5);
  alpha ~ normal(0, 1);
  gamma ~ normal(0.5, 1);
  eta0  ~ normal(-0.7, 0.1);
  eta1  ~ normal(0.4, 0.1); 
  eta2  ~ normal(0.3, 0.1); 
  ",
  
  # 3. MEASUREMENT-CALIBRATED: Aiming for EM = 0.092
  im = "
  pi[2] ~ beta(2, 5);
  alpha ~ normal(-1.73, 0.1); 
  gamma ~ normal(3.34, 0.1);
  eta0  ~ normal(0, 1);
  eta1  ~ normal(0.5, 1);
  eta2  ~ normal(0, 1);
  ",
  # 4. Selection- and Measurement-informed prior
  ism = "
  pi[2] ~ beta(2, 5);
  alpha ~ normal(-1.73, 0.1); 
  gamma ~ normal(3.34, 0.1);
  eta0  ~ normal(-0.7, 0.1);
  eta1  ~ normal(0.4, 0.1); 
  eta2  ~ normal(0.3, 0.1); 
  "
)

# 2. Estimation

priors_nam <- names(prior_blocks)
fitted_models <- list()
for (p in priors_nam){
  cat("\nCompiling model for prior:", p, "\n")
  full_mod <- gsub("{{PRIORS}}", prior_blocks[[p]], stan_code, fixed = T)
  stan_file <- file.path(tempdir(), paste0("model_", p, ".stan"))
  writeLines(full_mod, stan_file)
  model <- stan_model(file = stan_file)
  fitted_models[[p]] <- sampling(model, data = stan_data, chains = 3, cores = 3, iter = 4000, warmup = 1000, seed = 42, refresh = 5000)
}

#### - (IV) Graphs & Tables  - ####

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

cols_to_show <- c("Parameter", "Benchmark", "mean", "sd", "P2.5", "P50", "P97.5", "n_eff", "Rhat")

table_data[, cols_to_show] %>%
  kable(format = "latex",booktabs = T, escape = F, longtable = F, linesep = "", row.names = F, 
        caption = "Posterior Estimates vs. Benchmark", col.names = c("Parameter", "Benchmark", "Mean", "SD", "P2.5","P50", "P97.5", "$n_{eff}$", "$\\hat{R}$")) %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), font_size = 10) %>%
  pack_rows(priors_nam[1], 1, 6, latex_gap_space = "0em") %>%
  pack_rows(priors_nam[2], 7, 12, latex_gap_space = "0em") %>%
  pack_rows(priors_nam[3], 13, 18, latex_gap_space = "0em") %>% 
  pack_rows(priors_nam[4], 19, 24, latex_gap_space = "0em")


## LOO

# 1. Extract the log_lik from each fitted model and compute LOO
loo_res <- list()
for (p in names(fitted_models)) {
  log_lik_matrix <- extract_log_lik(fitted_models[[p]], parameter_name = "log_lik")   # Extract log-likelihood matrix [iterations x observations]
  loo_res[[p]] <- loo(log_lik_matrix)  # Calculate PSIS-LOO
  cat("\nLOO Diagnostics for", p, ":\n")
  print(loo_results[[p]])
}

# 2. Compare all models against each other
loo_comp <- loo_compare(loo_results)
m_labs <- c("base" = "$\\mathcal{M}_{0}$", "is" = "$\\mathcal{M}_{S}$", "im" = "$\\mathcal{M}_{M}$", "ism" = "$\\mathcal{M}_{C}$")
table_loo <- as.data.frame(loo_comp) %>%
  mutate(Model = m_labs[rownames(.)]) %>%
  select(Model, elpd_loo, se_elpd_loo, elpd_diff, se_diff) %>%
  kbl(format = "latex", booktabs = TRUE, escape = F, col.names = c("Prior", "ELPD", "$\\mathrm{SE}_{ELPD}$","$\\Delta$ ELPD", "$\\mathrm{SE}_{\\Delta}$"),
    caption = "Model Comparison via PSIS-LOO Cross-Validation", label = "4", align = "lcccc") %>%
  kable_styling(latex_options = c("hold_position"))



