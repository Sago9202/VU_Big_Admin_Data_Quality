## ------------------------------------------------- ##
## Bayesian Inference for proportions with ME and SE ##
## Script 1: Simulation Study                        ##
## Written by: Santiago Gómez-Echeverry              ##
## Last update: 19/10/2026                           ##
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

#### - (II) Data generation  - ####

## ---------------------------------------------------------------------- ##
## Simulation: selection, measurement and joint models                   ##
## ---------------------------------------------------------------------- ##

D <- 200                    # Replications
K <- 3                      # No. of items
q_idx <- 2                  # Item from the non-prob sample
seeds <- sample.int(1e6, D) # Seeds for random sampling
sigma0 <- 1.0
sigma1 <- 1.0
eta1 <- 0.5                 # eta0 is solved so that Pr(R = 1) = fq
fq <- 0.2

scenarios <- expand.grid(Delta_RY = c(0.03, 0.06, 0.09), J = c(0.3, 0.5, 0.8),
                         N = c(200, 500, 1000), pi = c(0.2, 0.4, 0.6),
                         delta_Z = c(0.6, 0.8, 1))

## find_cell(): every element of one condition, with fq fixed
##   Z | X = x ~ N(mu_x, sigma_x^2),  mu0 = -delta_Z/2,  mu1 = delta_Z/2
##   Pr(R=1 | X=x) = Phi((eta0 + eta1 mu_x + eta2 x) / s_x),  s_x = sqrt(1 + eta1^2 sigma_x^2)
##   Delta_RY = C * Delta_RX,  C = J sigma2_X / sigma2_Y   (Se = Sp = (1 + J)/2)
find_cell <- function(Delta_RY, J, N, pi, delta_Z) {
  # Indicator q
  Se <- (1 + J) / 2; Sp <- Se
  pi_Y <- pi * Se + (1 - pi) * (1 - Sp)
  C <- J * pi * (1 - pi) / (pi_Y * (1 - pi_Y))
  
  # Bounds (from 0 < p0, p1 < 1 with fq fixed)
  RX_min <- max(-fq / (1 - pi), -(1 - fq) / pi)
  RX_max <- min( fq / pi,        (1 - fq) / (1 - pi))
  RY_min <- C * RX_min
  RY_max <- C * RX_max
  feasible <- Delta_RY > RY_min && Delta_RY < RY_max
  
  # Selection probabilities
  Delta_RX <- Delta_RY / C
  p0 <- fq - pi * Delta_RX                  # Pr(R=1 | X=0)
  p1 <- p0 + Delta_RX                       # Pr(R=1 | X=1)
  
  # Selection parameters
  mu0 <- -delta_Z / 2; mu1 <- delta_Z / 2
  s0 <- sqrt(1 + eta1^2 * sigma0^2); s1 <- sqrt(1 + eta1^2 * sigma1^2)
  eta0 <- if (feasible) s0 * qnorm(p0) - eta1 * mu0 else NA
  a1   <- (eta0 + eta1 * mu1) / s1
  eta2 <- if (feasible) s1 * (qnorm(p1) - a1) else NA      # eta2 = s1 [Phi^{-1}(Delta_RX + p0) - a1]
  
  # True errors
  pi_XR <- pi * p1 / fq
  pi_YR <- pi_XR * Se + (1 - pi_XR) * (1 - Sp)
  
  data.frame(Delta_RY, J, N, pi, delta_Z, Se, Sp, mu0, mu1,
             eta0, eta1, eta2, fq, n_q = N * fq, p0, p1, Delta_RX,
             RX_min, RX_max, RY_min, RY_max, feasible,
             E_S = pi_XR - pi, E_M = pi_YR - pi_XR, E_C = pi_YR - pi)
}

## Apply to every condition
conditions <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i)
  do.call(find_cell, as.list(scenarios[i, ]))))
conditions <- cbind(condition = seq_len(nrow(conditions)), conditions)

options(width = 250)
print(conditions, digits = 3, row.names = FALSE)
cat("\nFeasible conditions:", sum(conditions$feasible), "of", nrow(conditions), "\n")

## ---------------------------------------------------------------------- ##
## Data simulation                                                        ##
## ---------------------------------------------------------------------- ##
sim_data <- vector("list", nrow(conditions))

for (s in 1:nrow(conditions)) {
  cond <- conditions[s, ]
  
  # True (population) values for this condition
  N      <- cond$N
  pi_X   <- cond$pi
  mu0_true <- cond$mu0; mu1_true <- cond$mu1
  sigma0_true <- sigma0; sigma1_true <- sigma1
  eta0_true <- cond$eta0; eta2_true <- cond$eta2
  # theta_vals: column 1 = Pr(Y_k = 1 | X = 0) = 1 - Sp_k, column 2 = Pr(Y_k = 1 | X = 1) = Se_k
  theta_vals <- matrix(c(0.15, 1 - cond$Sp, 0.20,
                         0.85, cond$Se,     0.80), ncol = 2, byrow = FALSE)
  
  EM_t <- cond$E_M
  ES_t <- cond$E_S
  EC_t <- cond$E_C
  
  # True parameter values on the scale of the Stan parameters. They are the centres of the
  # strong prior ("strong and correct"). Only the joint model receives them, and only the
  # strong prior uses them.
  prior_m <- list(pi_m    = c(1 - pi_X, pi_X),                                   # (Pr(X=0), Pr(X=1))
                  beta0_m = qlogis(theta_vals[, 1]),                             # logit(1 - Sp_k)
                  beta1_m = qlogis(theta_vals[, 2]) - qlogis(theta_vals[, 1]),   # logit(Se_k) - logit(1 - Sp_k)
                  eta0_m  = eta0_true,
                  eta1_m  = eta1,
                  eta2_m  = eta2_true,
                  mu_m    = c(mu0_true, mu1_true),
                  sigma_m = sigma0_true)
  
  sim_data[[s]] <- vector("list", D)
  
  for (d in 1:D) {
    set.seed(seeds[d])
    X <- rbinom(N, 1, pi_X)               # Latent target variable
    Y <- matrix(NA, nrow = N, ncol = K)   # Observed target variable
    for (k in 1:K) {
      Y[, k] <- ifelse(X == 1, rbinom(N, 1, theta_vals[k, 2]), rbinom(N, 1, theta_vals[k, 1]))
    }
    Z <- ifelse(X == 1, rnorm(N, mean = mu1_true, sd = sigma1_true),
                rnorm(N, mean = mu0_true, sd = sigma0_true))
    
    # MNAR process (Y2)
    pR <- pnorm(eta0_true + eta1 * Z + eta2_true * X)
    Rq <- rbinom(N, 1, pR)
    fq_obs <- mean(Rq)                    # realized fraction (not 'fq': that is the design value)
    Y_obs <- Y
    Y_obs[Rq == 0, q_idx] <- NA
    
    # MAR process (Y2^{-})
    p_miss <- 0.3
    for (i in setdiff(1:K, q_idx)) {
      miss_k <- rbinom(N, 1, p_miss)
      Y_obs[miss_k == 1, i] <- NA
    }
    
    # Placeholder for the NA's since Stan doesn't accept missing values
    obs <- ifelse(is.na(Y_obs), 0, 1)
    Y_obs[is.na(Y_obs)] <- 0
    
    sim_data[[s]][[d]] <- list(
      stan_data   = list(N = N, K = K, C = 2, Y = Y_obs, R = Rq, q_idx = q_idx,
                         fq = fq_obs, Z = Z, obs = obs,
                         prior_m = prior_m),       # centres of the strong prior (joint model only)
      true_values = list(EM = EM_t, ES = ES_t, EC = EC_t, eta0 = eta0_true, eta2 = eta2_true,
                         theta = theta_vals, pi = pi_X, fq = fq, Delta_RY = cond$Delta_RY),
      label = sprintf("C%03d_DRY%.2f_J%.1f_N%d_pi%.1f_dZ%.1f_d%03d", s, cond$Delta_RY,
                      cond$J, N, pi_X, cond$delta_Z, d))
  }
}

## ---------------------------------------------------------------------- ##
## Model estimation                                                       ##
## ---------------------------------------------------------------------- ##
library(cmdstanr); library(posterior); library(dplyr)
# cmdstanr compiles with CmdStan's own toolchain (-O3) and runs each chain as a separate
# CmdStan process. One-time setup:
#   install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
#   cmdstanr::install_cmdstan()
if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))
  stop("CmdStan not found: run cmdstanr::install_cmdstan() once, ",
       "or cmdstanr::set_cmdstan_path('<folder>') if it is installed elsewhere.")
# fold_data: folder for checkpoints and results (defined in your workspace header)

# 1. Stan models

# 1.1 Selection model only: Y_q taken as X. Only Z, R and Y_q (for R = 1) enter. Vectorised; Phi_approx.
# Z is binned (see z_bin_width): units with the same (R, Y_q, Z-bin) are collapsed into one
# pattern with weight w, so the likelihood runs over the patterns instead of the N units.
stan_code_sel <- "
functions {
  // log(Phi_approx(x)) = log_inv_logit(0.07056 x^3 + 1.5976 x), vectorised (|Phi_approx - Phi| < 1.4e-4)
  vector log_Phi_approx(vector x) {
    return log_inv_logit(0.07056 * x .* x .* x + 1.5976 * x);
  }
}

data {
  int<lower=0> P1;                      // patterns with R = 1: (Z-bin, Y_q)
  int<lower=0> P0;                      // patterns with R = 0: Z-bin
  int<lower=1> C;
  vector[P1] Z1;                        // bin midpoint
  array[P1] int<lower=0, upper=1> Y1;   // Y_q for R = 1, taken as X
  vector<lower=0>[P1] w1;               // no. of units with the pattern
  vector[P0] Z0;
  vector<lower=0>[P0] w0;
  real<lower=0, upper=1> fq;
}

transformed data {
  vector[P1] x1 = to_vector(Y1);
  real n_x1 = dot_product(w1, x1);
  real n_x0 = sum(w1) - n_x1;
}

parameters {
  simplex[C] pi;
  real eta0;
  real eta1;
  real eta2;
  vector[C] mu_z;
  real<lower=0> sigma_z;              // common SD of Z in both classes
}

model {
  // 1. Priors:
  {{PRIORS}}

  // 2a. Likelihood for R = 1: X observed, no mixture
  target += n_x0 * log(pi[1]) + n_x1 * log(pi[2]);
  target += dot_product(w1, -0.918938533204672742 - log(sigma_z)
                            - 0.5 * square((Z1 - mu_z[1] - (mu_z[2] - mu_z[1]) * x1) / sigma_z));
  target += dot_product(w1, log_Phi_approx(eta0 + eta1 * Z1 + eta2 * x1));

  // 2b. Likelihood for R = 0: X latent, marginalised over the two classes
  {
    vector[P0] s0 = eta0 + eta1 * Z0;
    vector[P0] la = log(pi[1]) - 0.918938533204672742 - log(sigma_z)
                    - 0.5 * square((Z0 - mu_z[1]) / sigma_z);
    vector[P0] lb = log(pi[2]) - 0.918938533204672742 - log(sigma_z)
                    - 0.5 * square((Z0 - mu_z[2]) / sigma_z);
    la += log_Phi_approx(-s0);             // log(1 - Phi_approx(s0))
    lb += log_Phi_approx(-(s0 + eta2));
    target += dot_product(w0, la + log1p_exp(lb - la));
  }
}

generated quantities {
  real pi_X = pi[2];
  real sigma2_X = pi[2] * (1.0 - pi[2]);
  real Delta_R_X = Phi((eta0 + eta1 * mu_z[2] + eta2) / sqrt(1.0 + square(eta1) * square(sigma_z))) -
                   Phi((eta0 + eta1 * mu_z[1]) / sqrt(1.0 + square(eta1) * square(sigma_z)));
  real ES = (Delta_R_X * sigma2_X) / fq;
}
"

# 1.2 Measurement model only: latent class model on the K items.
# Units are collapsed into unique response patterns (incl. missingness and R),
# so the likelihood loops over P patterns (at most 3^K * 2) instead of N units.
stan_code_meas <- "
data {
  int<lower=1> P;                       // number of unique response patterns
  int<lower=1> K;
  int<lower=1> C;
  array[P, K] int<lower=0, upper=1> Y;        // pattern responses (0 placeholder if missing)
  array[P, K] int<lower=0, upper=1> obs;      // pattern observation indicators
  vector<lower=0>[P] w;                 // no. of units with each pattern
  vector<lower=0>[P] wq;                // no. of units with each pattern that are in q (R = 1)
  int<lower=1, upper=K> q_idx;
}

parameters {
  simplex[C] pi;
  vector[K] beta0;
  vector<lower=0>[K] beta1;
}

transformed parameters {
  matrix[K, 2] theta;
  for (k in 1:K) {
    theta[k, 1] = inv_logit(beta0[k]);
    theta[k, 2] = inv_logit(beta0[k] + beta1[k]);
  }
}

model {
  // 1. Priors:
  {{PRIORS}}

  // 2. Likelihood (weighted by pattern frequency)
  for (p in 1:P) {
    vector[C] lp = log(pi);
    for (c in 1:C) {
      for (k in 1:K) {
        if (obs[p, k] == 1) lp[c] += bernoulli_lpmf(Y[p, k] | theta[k, c]);
      }
    }
    target += w[p] * log_sum_exp(lp);
  }
}

generated quantities {
  real pi_X = pi[2];
  real Se_q = theta[q_idx, 2];
  real Sp_q = 1.0 - theta[q_idx, 1];
  real J_q = Se_q + Sp_q - 1.0;

  // Average posterior P(X = 1 | Y) among units in the non-prob sample
  real piX_given_q = 0.0;
  for (p in 1:P) {
    if (wq[p] > 0) {
      vector[C] lp = log(pi);
      for (c in 1:C) {
        for (k in 1:K) {
          if (obs[p, k] == 1) lp[c] += bernoulli_lpmf(Y[p, k] | theta[k, c]);
        }
      }
      piX_given_q += wq[p] * softmax(lp)[2];
    }
  }
  piX_given_q /= sum(wq);

  real EM = 1.0 - Sp_q - piX_given_q * (1.0 - J_q);
}
"

# 1.3 Joint model (ME + SE). Vectorised; probit selection evaluated with Phi_approx
# Z is binned (see z_bin_width): units with the same (Y, obs, R, Z-bin) are collapsed into one
# pattern with weight w, so the likelihood runs over the patterns instead of the N units.
# The data block also holds the centres of the strong prior (*_m): the true values of the
# condition. The weak and base priors do not use them.
stan_code_joint <- "
functions {
  // log(Phi_approx(x)) = log_inv_logit(0.07056 x^3 + 1.5976 x), vectorised (|Phi_approx - Phi| < 1.4e-4)
  vector log_Phi_approx(vector x) {
    return log_inv_logit(0.07056 * x .* x .* x + 1.5976 * x);
  }

  // Log-likelihood contribution of every unit under each class (N x 2)
  // sgn = +1 if R = 1 and -1 if R = 0, so that log Pr(R | class) = log_Phi_approx(sgn .* index)
  matrix class_loglik(vector Z, matrix A, matrix B, vector sgn, vector pi, matrix theta,
                      real eta0, real eta1, real eta2, vector mu_z, real sigma_z) {
    int N = rows(Z);
    matrix[N, 2] lp;
    vector[N] s = eta0 + eta1 * Z;
    for (c in 1:2) {
      // class weight + auxiliary Z + measurement model (A: observed & Y = 1, B: observed & Y = 0)
      lp[, c] = log(pi[c]) - 0.918938533204672742 - log(sigma_z)
                - 0.5 * square((Z - mu_z[c]) / sigma_z)
                + A * log(theta[, c]) + B * log1m(theta[, c]);
    }
    // selection model (probit, via Phi_approx)
    lp[, 1] += log_Phi_approx(sgn .* s);
    lp[, 2] += log_Phi_approx(sgn .* (s + eta2));
    return lp;
  }
}

data {
  int<lower=1> N;                       // number of patterns
  int<lower=1> K;
  int<lower=1> C;
  array[N, K] int<lower=0, upper=1> Y;
  array[N] int<lower=0, upper=1> R;
  int<lower=1, upper=K> q_idx;
  real<lower=0, upper=1> fq;
  vector[N] Z;
  array[N, K] int<lower=0, upper=1> obs;
  vector<lower=0>[N] w;                 // no. of units with the pattern
  // Centres of the strong prior: true values of the condition (unused by weak and base)
  simplex[C] pi_m;                      // (Pr(X = 0), Pr(X = 1))
  vector[K] beta0_m;                    // logit(1 - Sp_k)
  vector<lower=0>[K] beta1_m;           // logit(Se_k) - logit(1 - Sp_k)
  real eta0_m;
  real eta1_m;
  real eta2_m;
  vector[C] mu_m;                       // (mu_0, mu_1)
  real<lower=0> sigma_m;
}

transformed data {
  matrix[N, K] A;      // 1 if item observed and Y = 1
  matrix[N, K] B;      // 1 if item observed and Y = 0
  vector[N] Rv;        // R as a vector
  vector[N] sgn;       // +1 if R = 1, -1 if R = 0
  for (i in 1:N) {
    Rv[i] = R[i];
    sgn[i] = 2 * R[i] - 1;
    for (k in 1:K) {
      A[i, k] = obs[i, k] * Y[i, k];
      B[i, k] = obs[i, k] * (1 - Y[i, k]);
    }
  }
}

parameters {
  simplex[C] pi;
  vector[K] beta0;
  vector<lower=0>[K] beta1;
  real eta0;
  real eta1;
  real eta2;
  vector[C] mu_z;
  real<lower=0> sigma_z;              // common SD of Z in both classes
}

transformed parameters {
  matrix[K, 2] theta;
  for (k in 1:K) {
    theta[k, 1] = inv_logit(beta0[k]);
    theta[k, 2] = inv_logit(beta0[k] + beta1[k]);
  }
}

model {
  // 1. Priors:
  {{PRIORS}}

  // 2. Likelihood: log(exp(lp1) + exp(lp2)) for every unit, vectorised
  {
    matrix[N, 2] lp = class_loglik(Z, A, B, sgn, pi, theta, eta0, eta1, eta2, mu_z, sigma_z);
    target += dot_product(w, lp[, 1] + log1p_exp(lp[, 2] - lp[, 1]));
  }
}

generated quantities {
  real pi_X = pi[2];
  real Se_q = theta[q_idx, 2];
  real Sp_q = 1.0 - theta[q_idx, 1];
  real J_q = Se_q + Sp_q - 1.0;
  real piX_given_q;
  {
    matrix[N, 2] lp = class_loglik(Z, A, B, sgn, pi, theta, eta0, eta1, eta2, mu_z, sigma_z);
    vector[N] p_x1_given_y = inv_logit(lp[, 2] - lp[, 1]);   // Pr(X = 1 | data) per unit
    piX_given_q = dot_product(p_x1_given_y, w .* Rv) / dot_product(w, Rv);
  }
  real sigma2_X = pi[2] * (1.0 - pi[2]);
  real EM = 1.0 - Sp_q - piX_given_q * (1.0 - J_q);
  real Delta_R_X = Phi((eta0 + eta1 * mu_z[2] + eta2) / sqrt(1.0 + square(eta1) * square(sigma_z))) -
                   Phi((eta0 + eta1 * mu_z[1]) / sqrt(1.0 + square(eta1) * square(sigma_z)));
  real ES = (Delta_R_X * sigma2_X) / fq;
  real EC = EM + ES;
}
"

# 2. Data builders: each model receives only the data it uses

# Z binning: Z is replaced by the midpoint of its bin (bins of width z_bin_width on a fixed
# grid), then identical units are collapsed with a count. Rounding adds z_bin_width^2 / 12
# to the variance of Z (0.0002 for 0.05, vs sigma_z^2 = 1). z_bin_width <- NULL: no binning.
z_bin_width <- 0.1
bin_z <- function(z, h = z_bin_width) if (is.null(h)) z else (floor(z / h) + 0.5) * h

make_data_sel <- function(sd) {
  r1 <- sd$R == 1
  zb <- bin_z(sd$Z)
  a1 <- aggregate(list(w = rep(1, sum(r1))),
                  list(Z = zb[r1], Y = as.integer(sd$Y[r1, sd$q_idx])), length)
  a0 <- aggregate(list(w = rep(1, sum(!r1))), list(Z = zb[!r1]), length)
  list(P1 = nrow(a1), P0 = nrow(a0), C = sd$C,
       Z1 = a1$Z, Y1 = as.integer(a1$Y), w1 = as.numeric(a1$w),
       Z0 = a0$Z, w0 = as.numeric(a0$w), fq = sd$fq)
}

# 2b. Starting values for the selection model, computed from each dataset
# (instead of Stan's random inits), with a small chain-specific shift so the
# chains do not start at the same point.
make_init_sel <- function(d, chain) {
  j  <- c(-0.1, 0, 0.1)[(chain - 1) %% 3 + 1]
  Z1 <- rep(d$Z1, d$w1); Y1 <- rep(d$Y1, d$w1); Z0 <- rep(d$Z0, d$w0)  # back to unit level
  z  <- c(Z1, Z0)
  R  <- c(rep(1, length(Z1)), rep(0, length(Z0)))
  p  <- min(max(mean(Y1), 0.05), 0.95)                       # share of Y_q = 1 among R = 1
  m0 <- if (any(Y1 == 0)) mean(Z1[Y1 == 0]) else mean(z)     # Z mean by Y_q among R = 1
  m1 <- if (any(Y1 == 1)) mean(Z1[Y1 == 1]) else mean(z)
  g  <- unname(suppressWarnings(glm(R ~ z, family = binomial(link = "probit")))$coefficients)
  list(pi      = c(1 - p - j / 2, p + j / 2),
       eta0    = g[1] + j,
       eta1    = g[2] + j,
       eta2    = j,
       mu_z    = c(m0, m1) + j,
       sigma_z = sd(z) * (1 + j))
}

make_data_meas <- function(sd) {
  Ycode <- ifelse(sd$obs == 1, sd$Y, -1)          # -1 marks a missing item
  pat <- as.data.frame(Ycode)
  pat$R <- sd$R
  agg <- pat %>% group_by(across(everything())) %>% summarise(w = n(), .groups = "drop")
  Yp <- as.matrix(agg[, 1:sd$K])
  obsp <- (Yp != -1) * 1L
  Yp[Yp == -1] <- 0
  storage.mode(Yp) <- "integer"; storage.mode(obsp) <- "integer"
  list(P = nrow(agg), K = sd$K, C = sd$C, Y = Yp, obs = obsp,
       w = as.numeric(agg$w), wq = as.numeric(agg$w * agg$R), q_idx = sd$q_idx)
}

make_data_joint <- function(sd) {
  u <- data.frame(Y = sd$Y, obs = sd$obs, R = sd$R, Z = bin_z(sd$Z))
  a <- aggregate(list(w = rep(1, nrow(u))), u, length)
  Yp <- as.matrix(a[, paste0("Y.", 1:sd$K)]);   storage.mode(Yp) <- "integer"
  op <- as.matrix(a[, paste0("obs.", 1:sd$K)]); storage.mode(op) <- "integer"
  dimnames(Yp) <- dimnames(op) <- NULL
  c(list(N = nrow(a), K = sd$K, C = sd$C, Y = Yp, R = as.integer(a$R), q_idx = sd$q_idx,
         fq = sd$fq, Z = a$Z, obs = op, w = as.numeric(a$w)),
    sd$prior_m)                                    # centres of the strong prior
}

# 3. Priors
# Split into blocks so each model only receives the priors for its own parameters.
# The selection-only and measurement-only models use the base prior only; the joint
# model is fitted with all three (see 'priors' in model_specs).
prior_components <- list(
  
  # Weak: flat on the probability scale.
  #  - pi: Dirichlet(1, 1) = uniform on Pr(X = 1).
  #  - eta: N(0, 1) on the probit scale makes the implied selection probabilities roughly
  #    uniform, and centring at 0 encodes neither the sampling rate nor a direction of selection.
  #  - beta: SD 1.5 on the logit scale spreads 1 - Sp and Se widely without piling mass at 0 or 1.
  weak = list(
    pi   = "pi      ~ dirichlet(rep_vector(1.0, C));",
    meas = "beta0   ~ normal(-1.4, 1.5);
  beta1   ~ normal(3.0, 1.5);",
    sel  = "eta0    ~ normal(0, 1);
  eta1    ~ normal(0, 1);
  eta2    ~ normal(0, 1);",
    z    = "mu_z    ~ normal(0, 1);
  sigma_z ~ exponential(1);"
  ),
  
  # Base: weakly informative and roughly centred. eta0 ~ N(-1, 0.75) encodes a sampling
  # rate of about 20%; beta0 ~ N(-1.4, 1) and beta1 ~ N(3, 1) describe items of moderate
  # quality (1 - Sp about 0.2, Se about 0.8).
  base = list(
    pi   = "pi      ~ dirichlet(rep_vector(1.0, C));",
    meas = "beta0   ~ normal(-1.4, 1);
  beta1   ~ normal(3.0, 1);",
    sel  = "eta0    ~ normal(-1, 0.75);
  eta1    ~ normal(0.5, 0.75);
  eta2    ~ normal(0, 1);",
    z    = "mu_z    ~ normal(0, 1);
  sigma_z ~ exponential(1);"
  ),
  
  # Strong: accurate and informative. Every prior is centred on the true value of the
  # condition (passed as data: *_m), with a small SD:
  #  - pi: Dirichlet(20 * pi_true), worth 20 prior observations (SD 0.09 at pi = 0.2,
  #    0.11 at pi = 0.4 and 0.6).
  #  - beta, eta, mu_z: SD 0.5 on the logit / probit / Z scale.
  #  - sigma_z: SD 0.25 around the true value 1.
  strong = list(
    pi   = "pi      ~ dirichlet(20 * pi_m);",
    meas = "beta0   ~ normal(beta0_m, 0.5);
  beta1   ~ normal(beta1_m, 0.5);",
    sel  = "eta0    ~ normal(eta0_m, 0.5);
  eta1    ~ normal(eta1_m, 0.5);
  eta2    ~ normal(eta2_m, 0.5);",
    z    = "mu_z    ~ normal(mu_m, 0.5);
  sigma_z ~ normal(sigma_m, 0.25);"
  )
)

# 4. Model specifications (run in this order: selection, measurement, joint;
#    each model goes through all conditions before the next one starts).
#    'priors': the prior specifications each model is fitted with. Only the joint model
#    is compared across priors; the single-source models use the base prior.
model_specs <- list(
  selection = list(
    code = stan_code_sel, make_data = make_data_sel, make_init = make_init_sel,
    # Caps the work per iteration at 2^8 - 1 = 255 leapfrog steps (default 10 -> 1023).
    # Affects efficiency (ESS), not the target distribution; see n_max_treedepth.
    max_treedepth = 8,
    prior_parts = c("pi", "sel", "z"),
    priors   = "base",
    cond_ids = conditions$condition,             # uses Z, R and Y_q: all five factors matter
    pooled   = character(0),
    diag_pars = c("pi", "eta0", "eta1", "eta2", "mu_z", "sigma_z", "ES", "pi_X"),
    report_pars = c("ES", "pi_X")
  ),
  measurement = list(
    code = stan_code_meas, make_data = make_data_meas,
    prior_parts = c("pi", "meas"),
    priors   = "base",
    cond_ids = conditions$condition,             # all conditions
    pooled   = character(0),
    diag_pars = c("pi", "beta0", "beta1", "EM", "Se_q", "Sp_q", "J_q", "pi_X"),
    report_pars = c("EM", "Se_q", "Sp_q", "J_q", "pi_X")
  ),
  joint = list(
    code = stan_code_joint, make_data = make_data_joint,
    # Caps the work per iteration at 2^8 - 1 = 255 leapfrog steps (default 10 -> 1023).
    # Affects efficiency (ESS), not the target distribution; see n_max_treedepth.
    max_treedepth = 8,
    prior_parts = c("pi", "meas", "sel", "z"),
    priors   = c("base", "weak", "strong"),      # base first, then the prior-sensitivity runs
    cond_ids = conditions$condition,           # all five factors matter
    pooled   = character(0),
    diag_pars = c("pi", "beta0", "beta1", "eta0", "eta1", "eta2", "mu_z", "sigma_z",
                  "EM", "ES", "EC", "Se_q", "Sp_q", "J_q", "pi_X"),
    report_pars = c("EM", "ES", "EC", "Se_q", "Sp_q", "J_q", "pi_X")
  )
)

# 5. Estimation

max_treedepth_default <- 10  # used for models without their own max_treedepth in model_specs
n_chains <- 3; n_iter <- 2000; n_warmup <- 1000   # 3 x 1000 = 3000 post-warmup draws per fit
n_warmup_short <- 200   # warm-up for replications 2..D, started from the pilot's adaptation
sampler_metric <- "diag_e"    # Stan default (dense_e gave no speed-up)
n_workers <- max(n_chains, parallel::detectCores() - 1)  # R workers; each runs one chain at a time
n_postwarmup_draws <- n_chains * (n_iter - n_warmup)

cond_vars <- c("Delta_RY", "J", "N", "pi", "delta_Z")   # design factors carried into the results

results_list <- list()       # Posterior summaries for each model's report_pars, with true values attached
diagnostics_list <- list()   # One row per fit: divergences, treedepth, E-BFMI, worst Rhat/ESS across diag_pars
checkpoint_folder <- file.path(fold_data, "checkpoints")
if(!dir.exists(checkpoint_folder)) dir.create(checkpoint_folder, recursive = T)

# Persistent cluster of n_workers R processes, started once.
# All chains of all replications of a condition are queued together, so every core
# stays busy (not only n_chains of them). Worker errors are returned and printed;
# worker output also goes to stan_workers.log.
cl <- parallel::makeCluster(n_workers, outfile = file.path(checkpoint_folder, "stan_workers.log"))
invisible(parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(cmdstanr))))
# Workers are separate R sessions: give them the same CmdStan installation as this session
invisible(parallel::clusterCall(cl, cmdstanr::set_cmdstan_path, cmdstanr::cmdstan_path()))


# Progress line for one condition. Time left = mean time per finished chain x chains
# still to start / workers, plus the time the running chains still need.
show_progress <- function(done, total, t0, chain_secs = numeric(0), n_workers = 1) {
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  left <- if (length(chain_secs) > 0) {
    avg <- mean(chain_secs) / 60                                  # minutes per chain
    max(avg * ceiling((total - done) / n_workers), 0)
  } else NA
  w <- 30; k <- round(w * done / total)
  cat(sprintf("\r    [%s%s] %3.0f%% | %d/%d chains | %.1f min elapsed | %s min per chain | ~%s min left   ",
              strrep("=", k), strrep(" ", w - k), 100 * done / total, done, total, el,
              if (length(chain_secs)) sprintf("%.1f", mean(chain_secs) / 60) else "?",
              if (is.na(left)) "?" else sprintf("%.0f", left)))
  if (done == total) cat("\n")
  utils::flush.console()
}

# Runs one chain on a worker with CmdStan; returns its draws, sampler diagnostics and the
# adapted step size and inverse metric (mass matrix)
run_one_chain <- function(task) {
  out_dir <- file.path(tempdir(), "cmdstan_out")
  dir.create(out_dir, showWarnings = FALSE)
  tryCatch({
    fit <- model$sample(data = task$dat, chains = 1, chain_ids = task$chain, seed = 42,
                        iter_warmup = task$warmup, iter_sampling = task$sampling,
                        init = if (is.null(task$init)) NULL else list(task$init),
                        metric = task$metric, step_size = task$step_size, inv_metric = task$inv_metric,
                        max_treedepth = task$max_treedepth, refresh = 0,
                        show_messages = FALSE, output_dir = out_dir)
    tm  <- tryCatch(fit$time()$chains, error = function(e) NULL)   # seconds CmdStan spent on this chain
    res <- list(draws      = fit$draws(variables = task$pars, format = "draws_array"),
                sdiag      = fit$sampler_diagnostics(format = "draws_array"),
                step_size  = fit$metadata()$step_size_adaptation[1],
                inv_metric = fit$inv_metric(matrix = task$metric == "dense_e")[[1]],   # vector for diag_e, matrix for dense_e
                time_warmup   = if (is.null(tm)) NA_real_ else tm$warmup[1],
                time_sampling = if (is.null(tm)) NA_real_ else tm$sampling[1])
    unlink(fit$output_files())                     # keep the temp folder small
    res
  }, error = function(e) paste("R error:", conditionMessage(e)))
}

# Sends a list of chain tasks to the workers (load-balanced) and updates the progress line
run_tasks <- function(tasks, prog) {
  n <- length(tasks); p <- length(cl); res <- vector("list", n)
  started <- rep(Sys.time(), n)
  for (i in seq_len(min(n, p))) {
    started[i] <- Sys.time()
    parallel:::sendCall(cl[[i]], run_one_chain, list(tasks[[i]]), tag = i)
  }
  for (i in seq_len(n)) {
    d <- parallel:::recvOneResult(cl)
    prog$secs <- c(prog$secs, as.numeric(difftime(Sys.time(), started[d$tag], units = "secs")))
    nxt <- i + min(n, p)
    if (nxt <= n) {
      started[nxt] <- Sys.time()
      parallel:::sendCall(cl[[d$node]], run_one_chain, list(tasks[[nxt]]), tag = nxt)
    }
    res[d$tag] <- list(d$value)
    prog$done <- prog$done + 1
    show_progress(prog$done, prog$total, prog$t0, prog$secs, p)
  }
  res
}

# Fits one model to all replications of a condition.
# 1. Pilot: replication 1 with the full warm-up (n_warmup).
# 2. Replications 2..D: short warm-up (n_warmup_short) that starts from the pilot's adapted
#    step size and inverse metric (averaged over its chains);
#    adaptation stays switched on.
#    Sampling iterations are the same for all replications.
# Returns, per replication, the combined draws and sampler diagnostics (NULL if a chain failed).
run_condition <- function(data_list, make_init = NULL, max_treedepth = max_treedepth_default, pars = NULL) {
  D_ <- length(data_list)
  make_task <- function(r, chain, warmup, step_size = NULL, inv_metric = NULL) {
    dat <- data_list[[r]]
    list(dat = dat, chain = chain, warmup = warmup, sampling = n_iter - n_warmup, metric = sampler_metric,
         max_treedepth = max_treedepth, pars = pars,
         init = if (is.null(make_init)) NULL else make_init(dat, chain),
         step_size = step_size, inv_metric = inv_metric)
  }
  prog <- new.env()
  prog$done <- 0; prog$total <- D_ * n_chains; prog$t0 <- Sys.time(); prog$secs <- numeric(0)
  show_progress(0, prog$total, prog$t0)
  
  # 1. Pilot
  pilot <- run_tasks(lapply(seq_len(n_chains), function(ch) make_task(1, ch, n_warmup)), prog)
  ok_p  <- vapply(pilot, function(x) is.list(x) && !is.null(x$draws), logical(1))
  if (any(ok_p)) {
    ss <- mean(vapply(pilot[ok_p], function(x) x$step_size, numeric(1)))
    im <- Reduce(`+`, lapply(pilot[ok_p], function(x) x$inv_metric)) / sum(ok_p)
    w2 <- n_warmup_short
  } else {                                          # pilot failed: full warm-up for all
    ss <- NULL; im <- NULL; w2 <- n_warmup
  }
  
  # 2. Remaining replications
  rest <- list()
  if (D_ > 1) {
    jobs <- expand.grid(chain = seq_len(n_chains), r = 2:D_)
    rest <- run_tasks(lapply(seq_len(nrow(jobs)), function(j)
      make_task(jobs$r[j], jobs$chain[j], w2, ss, im)), prog)
  }
  res <- c(pilot, rest)
  rep_of <- c(rep(1, n_chains), if (D_ > 1) rep(2:D_, each = n_chains))
  
  lapply(seq_len(D_), function(r) {
    rr <- res[rep_of == r]
    ok <- vapply(rr, function(x) is.list(x) && !is.null(x$draws), logical(1))
    if (!all(ok)) {
      cat("    [!] Replication", r, "- chains failed:", which(!ok), "\n")
      for (x in rr[!ok]) if (is.character(x)) cat("       ", x, "\n")
      return(NULL)
    }
    tw <- vapply(rr, function(x) x$time_warmup,   numeric(1))
    ts <- vapply(rr, function(x) x$time_sampling, numeric(1))
    list(draws = do.call(posterior::bind_draws, c(lapply(rr, `[[`, "draws"), along = "chain")),
         sdiag = do.call(posterior::bind_draws, c(lapply(rr, `[[`, "sdiag"), along = "chain")),
         time  = c(warmup   = sum(tw),             # seconds, summed over the chains
                   sampling = sum(ts),
                   slowest_chain = max(tw + ts)))   # seconds of the slowest chain
  })
}

# Wrapped in tryCatch(finally = ...) so the cluster is shut down even if the run
# stops early (error or interrupt). Assignments inside still go to the workspace.
# Order: selection (base) -> measurement (base) -> joint (base, weak, strong).
tryCatch({
  for(m in names(model_specs)) {
    spec <- model_specs[[m]]
    for(p in spec$priors) {
      max_treedepth_ctrl <- if (is.null(spec$max_treedepth)) max_treedepth_default else spec$max_treedepth
      cat("\nCompiling model:", m, "| prior:", p, "\n")
      this_priors <- paste(unlist(prior_components[[p]][spec$prior_parts]), collapse = "\n  ")
      this_code <- gsub("{{PRIORS}}", this_priors, spec$code, fixed = T)
      stan_dir <- file.path(checkpoint_folder, "stan_models")
      dir.create(stan_dir, showWarnings = FALSE)
      stan_file <- cmdstanr::write_stan_file(this_code, dir = stan_dir, basename = paste0("model_", m, "_", p))
      model <- cmdstan_model(stan_file)              # compiled once; recompiles only if the code changed
      parallel::clusterExport(cl, "model")           # send the compiled model to the workers once
      
      all_pars <- union(spec$diag_pars, spec$report_pars)
      
      for(s in spec$cond_ids) {
        cond_info <- conditions[s, c("condition", cond_vars)]
        cond_info[spec$pooled] <- NA               # factor that does not affect this model's data
        
        cat("  Estimating model:", m, "| prior:", p, "| condition:", s, "(", D, "replications )\n")
        t_start <- Sys.time()
        fits <- run_condition(lapply(sim_data[[s]], function(x) spec$make_data(x$stan_data)),
                              make_init = spec$make_init, max_treedepth = max_treedepth_ctrl,
                              pars = all_pars)
        n_diag_before <- length(diagnostics_list)
        
        for(r in 1:D) {
          curr_data <- sim_data[[s]][[r]]
          tv <- curr_data$true_values
          fit <- fits[[r]]
          
          # Skip (and record) fits without draws, so one failure does not stop the run
          if (is.null(fit)) {
            cat("    [!] No samples for", curr_data$label, "-- skipped\n")
            diagnostics_list[[length(diagnostics_list) + 1]] <- data.frame(
              model = m, prior = p, scenario = curr_data$label, cond_info, draw = r,
              n_divergent = NA, pct_divergent = NA, n_max_treedepth = NA, mean_leapfrog = NA, min_ebfmi = NA,
              max_rhat = NA, min_ess_bulk = NA, min_ess_tail = NA, failed = TRUE,
              pilot = (r == 1), time_warmup_s = NA, time_sampling_s = NA, time_total_s = NA,
              time_slowest_chain_s = NA, row.names = NULL)
            next
          }
          
          # --- Parameter-level summary (mean, sd, quantiles, Rhat, bulk/tail ESS) via posterior ---
          draws_arr <- fit$draws
          param_sum <- posterior::summarise_draws(draws_arr)
          base_name <- sub("\\[.*$", "", param_sum$variable)
          param_diag <- param_sum[base_name %in% spec$diag_pars, ]
          
          # True values for every reported quantity (pi_X from this condition)
          true_vec <- c(EM = tv$EM, ES = tv$ES, EC = tv$EC,
                        Se_q = tv$theta[q_idx, 2], Sp_q = 1 - tv$theta[q_idx, 1],
                        J_q = tv$theta[q_idx, 2] - tv$theta[q_idx, 1], pi_X = tv$pi)
          
          fit_summary <- as.data.frame(param_sum[param_sum$variable %in% spec$report_pars, ])
          fit_summary$model <- m
          fit_summary$prior <- p
          fit_summary$scenario <- curr_data$label
          fit_summary <- cbind(fit_summary, cond_info[rep(1, nrow(fit_summary)), ], row.names = NULL)
          fit_summary$draw <- r
          fit_summary$true_value <- unname(true_vec[fit_summary$variable])
          fit_summary$true_EM <- tv$EM
          fit_summary$true_ES <- tv$ES
          fit_summary$true_EC <- tv$EC
          results_list[[length(results_list) + 1]] <- fit_summary
          
          # --- Fit-level HMC diagnostics: divergences, treedepth saturation, E-BFMI ---
          sp_div <- posterior::extract_variable_matrix(fit$sdiag, "divergent__")   # iterations x chains
          sp_td  <- posterior::extract_variable_matrix(fit$sdiag, "treedepth__")
          sp_en  <- posterior::extract_variable_matrix(fit$sdiag, "energy__")
          sp_lf  <- posterior::extract_variable_matrix(fit$sdiag, "n_leapfrog__")   # steps per iteration
          n_divergent <- sum(sp_div)
          n_max_treedepth <- sum(sp_td >= max_treedepth_ctrl)
          ebfmi_chain <- apply(sp_en, 2, function(e) mean(diff(e)^2) / var(e))
          min_ebfmi <- min(ebfmi_chain)
          
          diag_row <- data.frame(
            model = m, prior = p, scenario = curr_data$label, cond_info, draw = r,
            n_divergent = n_divergent,
            pct_divergent = n_divergent / n_postwarmup_draws,
            n_max_treedepth = n_max_treedepth,
            mean_leapfrog = mean(sp_lf),
            min_ebfmi = min_ebfmi,
            max_rhat = max(param_diag$rhat, na.rm = T),
            min_ess_bulk = min(param_diag$ess_bulk, na.rm = T),
            min_ess_tail = min(param_diag$ess_tail, na.rm = T),
            failed = FALSE,
            # Computation time of this fit (seconds, from CmdStan). Replication 1 is the pilot
            # with the full warm-up; the others reuse its adaptation and have a short warm-up.
            pilot                = (r == 1),
            time_warmup_s        = unname(fit$time["warmup"]),         # summed over the chains
            time_sampling_s      = unname(fit$time["sampling"]),       # summed over the chains
            time_total_s         = unname(fit$time["warmup"] + fit$time["sampling"]),
            time_slowest_chain_s = unname(fit$time["slowest_chain"]),
            row.names = NULL
          )
          diagnostics_list[[length(diagnostics_list) + 1]] <- diag_row
          
          # Flag problems immediately in the console so a long run doesn't hide issues
          if(n_divergent > 0 || diag_row$max_rhat > 1.01 || min_ebfmi < 0.2 || diag_row$min_ess_bulk < 400) {
            cat("    [!]", curr_data$label, "-- divergent:", n_divergent,
                "| max_treedepth hits:", n_max_treedepth,
                "| max Rhat:", round(diag_row$max_rhat, 3),
                "| min E-BFMI:", round(min_ebfmi, 3),
                "| min ESS_bulk:", round(diag_row$min_ess_bulk, 0), "\n")
          }
          
          rm(fit, draws_arr, sp_div, sp_td, sp_en, sp_lf)
        }
        rm(fits)
        
        # One summary line per condition: how many of the D fits were flagged or failed
        d_cond <- dplyr::bind_rows(diagnostics_list[(n_diag_before + 1):length(diagnostics_list)])
        flagged <- with(d_cond, !failed & (n_divergent > 0 | max_rhat > 1.01 | min_ebfmi < 0.2 | min_ess_bulk < 400))
        cond_wall_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
        cat(sprintf("    done: %d fits | %d flagged | %d failed | %.0f steps/iteration | %.1f min\n",
                    nrow(d_cond), sum(flagged), sum(d_cond$failed), mean(d_cond$mean_leapfrog, na.rm = TRUE),
                    cond_wall_min))
        
        # Store the condition's wall-clock time (the minutes printed above) and the number
        # of parallel workers on every fit of this condition, next to the flags
        for (i in (n_diag_before + 1):length(diagnostics_list)) {
          diagnostics_list[[i]]$cond_wall_min <- cond_wall_min
          diagnostics_list[[i]]$n_workers     <- n_workers
        }
        
        # Checkpoint after each condition is fully estimated for this model and prior
        saveRDS(dplyr::bind_rows(results_list), file.path(checkpoint_folder, paste0("results_checkpoint_", m, "_", p, "_condition", s, ".rds")))
        saveRDS(dplyr::bind_rows(diagnostics_list), file.path(checkpoint_folder, paste0("diagnostics_checkpoint_", m, "_", p, "_condition", s, ".rds")))
      }
    }
  }
  
}, finally = parallel::stopCluster(cl))

results_all <- dplyr::bind_rows(results_list)
diagnostics_all <- dplyr::bind_rows(diagnostics_list)
saveRDS(results_all, file.path(fold_data, "results_all_models.rds"))
saveRDS(diagnostics_all, file.path(fold_data, "diagnostics_all_models.rds"))