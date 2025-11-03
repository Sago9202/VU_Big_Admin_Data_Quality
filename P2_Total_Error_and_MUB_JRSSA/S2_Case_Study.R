##########################################
## Total Error Case Study               ##
## Written by: Santiago Gómez-Echeverry ##
## Last modified: 08/10/2025            ##
##########################################

#### - (I) Working space and packages - ####

rm(list = ls())  # Remove everything in the working environment
packs <- c('data.table', 'ggplot2', 'tidyr', 'tidySEM', 'dplyr', 'forcats', 'reshape2', 'lavaan', 'semTools', 'modelr', 'reshape', 'stringr', 'cols4all', 'xtable', 'progress')
ipacks <- packs %in% rownames(installed.packages()) 
if(any(ipacks == F)){install.packages(packs[!ipacks])}
invisible(lapply(packs, library, character.only = T))
sys <- Sys.info()
fold_data <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Italian Data/Total Error")
fold_graphs <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/2 - Continuous Multisource TE/4 - Graphs & Tables")

compute_phi <- function(rho_XS, rho_XZ, rho_ZS) {
  num <- rho_XS - rho_XZ * rho_ZS
  denom <- (rho_XS - rho_XZ * rho_ZS) + (rho_ZS - rho_XZ * rho_XS)
  return(num / denom)
}

#### - (II) Data arrangement  - ####
setwd(fold_data)
full_dat <- fread("TotErr_Rel20.csv") # We will use fread since we have a very large data set
full_N <- nrow(full_dat)   
full_N  # Number of observations is 24,484,611

# As a first exercise, let us use a limited number of observations so that we can perform faster analyses
# To do this, let us first see how many unique values we have

n_prid <- length(unique(full_dat$individual_id))
p_sample <- 0.1 # Percentage of the sample
sample_dat <- full_dat[full_dat$individual_id %in% sample(unique(full_dat$individual_id), n_prid*p_sample, replace = F, prob = NULL),]

# Note that we are sampling based on the id's so, in theory, we should be able to see entire trajectories as in the full data.

# Let us consider age as our target variable first
full_dat$survey_yq <- paste0(full_dat$survey_year, full_dat$survey_quarter)
colnames(full_dat)[colnames(full_dat)=="age31_12_rbiadm"] <- "age_adm"
full_dat <- full_dat %>% 
  filter(age_adm>=18)
# table(full_dat$gender_adm)
# full_dat$gender_adm[full_dat$gender_adm == ""] <- NA # Properly assing the NA's in the gender variable
# table(full_dat$gender_adm)

#  Now let's turn to income

full_dat$inc_adm <- full_dat$pe_Wages_Salaries_adm/12
# full_dat$inc_adm[full_dat$survey_month == 1] <- full_dat$gross_m_wage_aj_adm_1[full_dat$survey_month == 1] # This is wrong. We need to account for the year
# for (i in 2:12){
#  full_dat$inc_adm[full_dat$survey_month == i] <- full_dat[[paste0("gross_m_wage_aj_adm_", i)]][full_dat$survey_month == i]
# }

full_dat$inc_lfs <- full_dat$net_wage
full_dat$inc_lfs[full_dat$inc_lfs==99999 | full_dat$inc_lfs==99998 | full_dat$inc_lfs==99997] <- NA 
full_dat$inc_hh <- full_dat$hh_total_income_adm/12
cor(full_dat$inc_adm, full_dat$inc_lfs, use = "complete")

reg <- lm(inc_adm ~ inc_hh + hhsize_adm + factor(gender_rbiadm) + factor(edu_isced_rbiadm), data = full_dat)
summary(reg)
full_dat <- add_predictions(full_dat, reg, var = "aux")
cor(full_dat$aux[full_dat$survey_year == 2019 & full_dat$survey_month == 1], 
    full_dat$inc_adm[full_dat$survey_year == 2019 & full_dat$survey_month == 1], use = "complete.obs")

Dat <- vector(mode = "list", length = 3)
samps <- c(5, 10, 20)

sample_dat %>%
  summarise(unique_ids = n_distinct(paste(individual_id, survey_month, survey_year, survey_quarter))) %>%
  mutate(all_unique = unique_ids == n())

full_dat <- full_dat %>% 
  mutate(unq_id = paste0(survey_year, survey_month, individual_id))
full_dat$idn <- duplicated(full_dat$unq_id)
table(full_dat$idn)/nrow(full_dat) # 1.5% of the cases are duplicates - Why? 
full_dat <- full_dat[idn == F,]
full_dat$inc_adm_full <- (full_dat$inc_adm - mean(full_dat$inc_adm, na.rm = T))/sd(full_dat$inc_adm, na.rm = T)
full_dat$inc_lfs_full <- (full_dat$inc_lfs - mean(full_dat$inc_lfs, na.rm = T))/sd(full_dat$inc_lfs, na.rm = T)

w_lfs <- glm(inc_lfs_full ~ inc_adm_full, data = full_dat) %>%
  resid() %>%
  var()
w_lfs <- 1/w_lfs

w_adm <- glm(inc_adm_full ~ inc_lfs_full, data = full_dat) %>%
  resid() %>%
  var()
w_adm <- 1/w_adm

full_dat$inc_bench <- (w_adm*full_dat$inc_adm_full + w_lfs*full_dat$inc_lfs_full)/(w_adm + w_lfs)

full_data_long <- full_dat %>%
  dplyr::select(individual_id, survey_yq, inc_lfs, inc_adm, inc_lfs_full, inc_adm_full, inc_bench, age_lfs, age_adm, aux, sample5, sample10, sample20) %>% 
  pivot_wider(names_from = survey_yq, values_from = c(inc_lfs, inc_adm, inc_lfs_full, inc_adm_full, inc_bench, age_lfs, age_adm, aux)) %>% 
  mutate(across(
    starts_with("inc_lfs_") | starts_with("inc_adm_") |starts_with("inc_bench_") | starts_with("aux_") | starts_with("age_lfs_") | starts_with("age_adm_"),  # select all income columns
    ~ (.-mean(., na.rm = TRUE)) / sd(., na.rm = TRUE),    # standardize each column
    .names = "{.col}_std"                                 # create new columns with _std suffix
  )) 

mod0 <- 'linc =~ lmbd1*inc_adm_full_20192_std + lmbd2*inc_lfs_full_20192_std + lmbd3*inc_lfs_full_20191_std
         lage =~ dlt1*age_adm_20192_std + dlt2*age_lfs_20192_std + dlt3*age_lfs_20191_std
         
         linc ~~ lage
         # Effects coding constraint
         lmbd1 + lmbd2 + lmbd3 == 3
         dlt1 + dlt2 + dlt3 == 3
        '
fit <- sem(model = mod0, data = full_data_long, meanstructure = T, effect.coding = T)
summary(fit)
x_hat_std <- (lavPredict(fit)[,"linc"]- mean(lavPredict(fit)[,"linc"]))/sd(lavPredict(fit)[,"linc"])
mx_hat <- mean(x_hat_std)
s2x_hat <- var(x_hat_std)
my_hat <- lavPredictY(fit, ynames = "inc_adm_full_20192_std", xnames = c("inc_lfs_full_20192_std", "inc_lfs_full_20191_std")) %>% 
  mean()
s2y_hat <- lavPredictY(fit, ynames = "inc_adm_full_20192_std", xnames = c("inc_lfs_full_20192_std", "inc_lfs_full_20191_std")) %>% 
  var() %>% 
  as.numeric()
rel <- reliability(fit)
gamma_hat <- rel[rownames(rel)=="omega", colnames(rel)=="linc"]
parms <- as.data.frame(parameterEstimates(fit, standardized = T))
lambda <- parms[parms$op == "=~" & parms$rhs == "inc_adm_full_20192_std", colnames(parms)=="std.all"]
Fit <- fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic", "bic"))
full_data_long$inc_bench2_20192 <- mx_hat + gamma_hat*(s2x_hat/s2y_hat)*(full_data_long$inc_adm_full_20192_std - my_hat)
full_data_long$inc_bench2_20192_std <- (full_data_long$inc_bench2_20192 - mean(full_data_long$inc_bench2_20192, na.rm = T))/sd(full_data_long$inc_bench2_20192, na.rm = T)

for (i in 1:length(samps)){
  samp_size <- samps[i]
  t_data <- full_data_long  
  cutoff <- quantile(t_data$inc_adm_20192, probs = 0.70, na.rm = TRUE)
  trows_np <- t_data$inc_adm_20192 < cutoff & !is.na(t_data$inc_adm_20192)        # Identify non-probability sample rows (everyone below cutoff)
  tcols_np <- grepl("inc_adm", names(t_data)) | grepl("age_adm", names(t_data))   # Columns to blank out for the non-probability sample
  t_data[trows_np, tcols_np] <- NA                                                # Assign NA in those columns for the non-probability sample

  # And let's use the probability sample of the given size 
  trows_p <- t_data[[paste0("sample", samp_size)]] == 1
  tcols_p <- grepl("inc_lfs", colnames(t_data))| grepl("age_2", colnames(t_data))
  t_data[trows_p, tcols_p] <- NA
  
  t_data$S <- ifelse(is.na(t_data$inc_adm_20192_std), 0, 1)
  Dat[[i]] <- t_data
}

test_data <- Dat[[1]]
cor_long <- test_data[, c(grep("20192", colnames(test_data)), grep("20191", colnames(test_data)))] %>%
  filter(age_adm_20192>=18) %>% 
  dplyr::select(!ends_with("std")) %>%
  dplyr::select(!starts_with("inc_bench")) %>% 
  dplyr::select(!ends_with("full_20191")) %>% 
  dplyr::select(!ends_with("full_20192")) %>%
  cor(use = "pairwise.complete.obs") %>% 
  melt()
colnames(cor_long) <- c("Var1", "Var2", "Correlation")

corg <- ggplot(cor_long, aes(Var1, Var2, fill = Correlation)) +
  geom_tile() +
  scale_fill_viridis_c(option = "mako", name = "Correlation", limits = c(0, 1), direction = -1) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Corr_IT.png", plot = corg, path = fold_graphs, width = 15, height = 10, units = "cm")

#### - (III)  Estimation  - ####

mod1 <- 'linc =~ lmbd1*inc_adm_20192_std + lmbd2*inc_lfs_20192_std + lmbd3*inc_lfs_20191_std
         lage =~ dlt1*age_adm_20192_std + dlt2*age_lfs_20192_std + dlt3*age_lfs_20191_std
         
         linc ~~ lage
         # Effects coding constraint
         lmbd1 + lmbd2 + lmbd3 == 3
         dlt1 + dlt2 + dlt3 == 3
        '

Res <- matrix(NA, nrow = 9, ncol = 18) # Rows = sample size (3), Columns =
if(!exists("BootRes")) BootRes <- list()
r_methods <- c("monte_carlo", "naive", "adjusted")
set.seed(42) # For reproducibility
for (i in 1:3){
  for (j in 1:3){
    ### - Direct estimation - ###
    t_data <- Dat[[i]]
    S <- t_data$S
    methd <- r_methods[j]
    # 'Actual' errors
    m_err <- mean(full_data_long$inc_adm_20192_std, na.rm = T) - mean(full_data_long$inc_bench_20192_std[S==1], na.rm = T)
    r_err <- mean(full_data_long$inc_bench_20192_std[S==1], na.rm = T) - mean(full_data_long$inc_bench_20192, na.rm = T)
    t_err <- m_err + r_err
  
    # -  Step 1 - Measurement error - #
    fit <- sem(model = mod1, data = t_data, meanstructure = T, effect.coding = T)
    summary(fit)
    x_hat_std <- (lavPredict(fit)[,"linc"]- mean(lavPredict(fit)[,"linc"]))/sd(lavPredict(fit)[,"linc"])
    mx_hat <- mean(x_hat_std)
    s2x_hat <- var(x_hat_std)
    my_hat <- lavPredictY(fit, ynames = "inc_adm_20192_std", xnames = c("inc_lfs_20192_std", "inc_lfs_20191_std")) %>% 
      mean()
    s2y_hat <- lavPredictY(fit, ynames = "inc_adm_20192_std", xnames = c("inc_lfs_20192_std", "inc_lfs_20191_std")) %>% 
      var() %>% 
      as.numeric()
    rel <- reliability(fit)
    gamma_hat <- rel[rownames(rel)=="omega", colnames(rel)=="linc"]
    parms <- as.data.frame(parameterEstimates(fit, standardized = T))
    lambda <- parms[parms$op == "=~" & parms$rhs == "inc_adm_20192_std", colnames(parms)=="std.all"]
    Fit <- fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic", "bic"))
    x_hat <- mx_hat + gamma_hat*(s2x_hat/s2y_hat)*(t_data$inc_adm_20192_std - my_hat)

    m_err_hat <- mean(t_data$inc_adm_20192_std, na.rm = T) - mean(x_hat[S==1], na.rm = T)

    # -  Step 2 - Representation error - #
    Z <- t_data$aux_20192_std
    r_yz <- cor(t_data$inc_adm_20192_std, Z, use = "complete.obs")
    r_zs <- cor(Z, S, use = "complete.obs")
    if (methd == "monte_carlo"){
      M <- 2000
      grid_n    <- 10000
      rho_grid  <- seq(-1, 1, length.out = grid_n)
      psi_grid  <- sapply(rho_grid, function(rho_xs) {compute_phi(rho_xs, r_yz, r_zs)})
      psi_valid <- psi_grid[is.finite(psi_grid) & psi_grid >= 0 & psi_grid <= 1]
      psi_min <- min(psi_valid); psi_max <- max(psi_valid)
      phi_guess <- runif(M, min = psi_min, max = psi_max)
      # MUB
      g_hat <- (phi_guess + (1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess))
      c_obs <- sd(t_data$inc_adm_20192_std, na.rm = T)/sd(Z, na.rm = T)
      r_err_hat <- (g_hat*c_obs*(mean(Z[S == 1], na.rm = T) - mean(Z, na.rm = T))) %>% 
        mean(na.rm = T)
      
      # MUBM  
      g_hat_m <- (phi_guess*sqrt(gamma_hat) + gamma_hat*(1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess)*sqrt(gamma_hat))
      c_obs_m <- (sd(t_data$inc_adm_20192_std, na.rm = T)*sqrt(gamma_hat))/sd(Z, na.rm = T)
      rm_err_hat <- g_hat_m*c_obs_m*(mean(Z[S == 1], na.rm = T) - mean(Z, na.rm = T))
      rm_err_hat <- mean(rm_err_hat, na.rm = T)/lambda
    }
    if (methd == "naive"){
      r_ys <- cor(t_data$inc_lfs_20192_std, S, use = "complete.obs")
      phi_guess <- compute_phi(r_ys, r_yz, r_zs) %>% 
        floor()
      # MUB
      g_hat <- (phi_guess + (1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess))
      c_obs <- sd(t_data$inc_adm_20192_std, na.rm = T)/sd(Z, na.rm = T)
      r_err_hat <- g_hat*c_obs*(mean(Z[S == 1], na.rm = T) - mean(Z, na.rm = T))
    
      # MUBM  
      g_hat_m <- (phi_guess*sqrt(gamma_hat) + gamma_hat*(1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess)*sqrt(gamma_hat))
      c_obs_m <- (sd(t_data$inc_adm_20192_std, na.rm = T)*sqrt(gamma_hat))/sd(Z, na.rm = T)
      rm_err_hat <- g_hat_m*c_obs_m*(mean(Z[S == 1], na.rm = T) - mean(Z, na.rm = T))/lambda
    }
    if (methd == "adjusted") {
      reg_coef <- coef(glm(S~t_data$inc_lfs_20192_std + Z, family = "binomial"))
      phi_guess <- reg_coef[2]/(reg_coef[2] + reg_coef[3])
      # MUB
      g_hat <- (phi_guess + (1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess))
      c_obs <- sd(t_data$inc_adm_20192_std, na.rm = T)/sd(Z, na.rm = T)
      r_err_hat <- g_hat*c_obs*(mean(Z[S == 1], na.rm = T) - mean(Z, na.rm = T))
      # MUBM  
      g_hat_m <- (phi_guess*sqrt(gamma_hat) + gamma_hat*(1 - phi_guess)*r_yz)/(phi_guess*r_yz + (1 - phi_guess)*sqrt(gamma_hat))
      c_obs_m <- (sd(t_data$inc_adm_20192_std, na.rm = T)*sqrt(gamma_hat))/sd(Z, na.rm = T)
      rm_err_hat <- g_hat_m*c_obs_m*(mean(Z[S == 1], na.rm = T) - mean(Z, na.rm = T))/lambda
    }
    t_err_hat <- m_err_hat + rm_err_hat
    
    ### - Bootstrapping - ###
    B <- 300               # no. of bootstrap reps.
    M_boot <- 500          # phi draws inside bootstrap for Monte Carlo 
    grid_n_boot <- 2000    # grid used to find phi-range in bootstrap 
    
    # Vectors for the bootstrap draws
    abs_m <- abs_r <- abs_rm <- abs_t <- boot_m <-  boot_r <- boot_rm <- boot_t <-  boot_gamma <- numeric(B)
    nrows <- nrow(t_data)
    pb <- progress_bar$new(format = paste0("Sample ", i, ", Method ", j, " [:bar] :percent ETA: :eta"), total = B, clear = FALSE, width = 80)
    for(b in 1:B){
      pb$tick() 
      # sample rows with replacement
      ids_b <- sample(seq_len(nrows), size = nrows, replace = T)
      t_b <- t_data[ids_b,]
      Z_b <- t_b$aux_20192_std
      S_b <- t_b$S
      
      # 1) re-fit S1 on bootstrap sample 
      fit_b <- sem(model = mod1, data = t_b, meanstructure = T, effect.coding = T)
      x_hat_std_b <- (lavPredict(fit)[,"linc"]- mean(lavPredict(fit)[,"linc"]))/sd(lavPredict(fit)[,"linc"])
      mx_hat_b <- mean(x_hat_std_b)
      s2x_hat_b <- var(x_hat_std_b)
      my_hat_b <- lavPredictY(fit_b, ynames = "inc_adm_20192_std", xnames = c("inc_lfs_20192_std", "inc_lfs_20191_std")) %>% 
        mean()
      s2y_hat_b <- lavPredictY(fit_b, ynames = "inc_adm_20192_std", xnames = c("inc_lfs_20192_std", "inc_lfs_20191_std")) %>% 
        var() %>% 
        as.numeric()
      rel_b <- reliability(fit_b)
      gamma_hat_b <- rel[rownames(rel)=="omega", colnames(rel)=="linc"]
      parms_b <- as.data.frame(parameterEstimates(fit, standardized = T))
      lambda_b <- parms[parms$op == "=~" & parms$rhs == "inc_adm_20192_std", colnames(parms)=="std.all"]
      x_hat_b <- mx_hat_b + gamma_hat_b*(s2x_hat_b/s2y_hat_b)*(t_b$inc_adm_20192_std - my_hat_b)
      m_err_hat_b <- mean(t_b$inc_adm_20192_std, na.rm = T) - mean(x_hat_b[S==1], na.rm = T)
      
      # reliability gamma and lambda for bootstrap
      rel_b <- tryCatch(reliability(fit_b), error = function(e) NULL)
      gamma_b <- if(!is.null(rel_b)) as.numeric(rel_b["omega", "linc"]) else NA
      parms_b <- as.data.frame(parameterEstimates(fit_b, standardized = T))
      lambda_b <- as.numeric(parms_b[parms_b$op == "=~" & parms_b$rhs == "inc_adm_20192_std", "std.all"])
      # 2) Compute S2 on bootstrap sample
      r_yz_b <- cor(t_b$inc_adm_20192_std, Z_b, use = "complete.obs")
      r_zs_b <- cor(Z_b, S_b, use = "complete.obs")
      
      if(methd == "monte_carlo"){
        rho_grid_b <- seq(-1, 1, length.out = grid_n_boot)
        psi_grid_b <- sapply(rho_grid_b, function(rho_xs) { compute_phi(rho_xs, r_yz_b, r_zs_b) })
        psi_valid_b <- psi_grid_b[is.finite(psi_grid_b) & psi_grid_b >= 0 & psi_grid_b <= 1]
        if(length(psi_valid_b) == 0){
          phi_vec_b <- runif(M_boot, 0, 1)
        } else {
          phi_vec_b <- runif(M_boot, min = min(psi_valid_b), max = max(psi_valid_b))
        }
        # MUB
        g_hat_b <- (phi_vec_b + (1 - phi_vec_b) * r_yz_b)/(phi_vec_b * r_yz_b + (1 - phi_vec_b))
        c_obs_b <- (sd(t_b$inc_adm_20192_std, na.rm = T)) / sd(Z_b, na.rm = T)
        r_err_vec_b <- g_hat_b * c_obs_b * (mean(Z_b[S_b == 1], na.rm = T) - mean(Z_b, na.rm = T))
        r_err_hat_b <- mean(r_err_vec_b, na.rm = T)/lambda_b
        
        # MUBM
        g_hat_m_b <- (phi_vec_b * sqrt(gamma_b) + gamma_b * (1 - phi_vec_b) * r_yz_b)/(phi_vec_b * r_yz_b + (1 - phi_vec_b) * sqrt(gamma_b))
        c_obs_m_b <- (sd(t_b$inc_adm_20192_std, na.rm = T) * sqrt(gamma_b)) / sd(Z_b, na.rm = T)
        rm_err_vec_b <- g_hat_m_b * c_obs_m_b * (mean(Z_b[S_b == 1], na.rm = T) - mean(Z_b, na.rm = T))
        rm_err_hat_b <- mean(rm_err_vec_b, na.rm = T)/lambda_b
      }
      if(methd == "naive"){
        r_ys_b <- suppressWarnings(cor(t_b$inc_lfs_20192_std, S_b, use = "complete.obs"))
        phi_guess_b <- compute_phi(r_ys_b, r_yz_b, r_zs_b)
        # MUB
        g_hat_b <- (phi_guess_b + (1 - phi_guess_b) * r_yz_b)/(phi_guess_b * r_yz_b + (1 - phi_guess_b))
        c_obs_b <- (sd(t_b$inc_adm_20192_std, na.rm = T)) / sd(Z_b, na.rm = T)
        r_err_hat_b <- g_hat_b * c_obs_b * (mean(Z_b[S_b == 1], na.rm = T) - mean(Z_b, na.rm = T))/lambda_b
        
        # MUBM
        g_hat_m_b <- (phi_guess_b * sqrt(gamma_b) + gamma_b * (1 - phi_guess_b) * r_yz_b)/(phi_guess_b * r_yz_b + (1 - phi_guess_b) * sqrt(gamma_b))
        c_obs_m_b <- (sd(t_b$inc_adm_20192_std, na.rm = T) * sqrt(gamma_b)) / sd(Z_b, na.rm = T)
        rm_err_hat_b <- g_hat_m_b * c_obs_m_b * (mean(Z_b[S_b == 1], na.rm = T) - mean(Z_b, na.rm = T))/lambda_b
      }
      
      if(methd == "adjusted"){
        glm_b <- coef(glm(S_b ~ t_b$inc_lfs_20192_std + Z_b, family = "binomial"))
        phi_guess_b <- as.numeric(glm_b[2] / (glm_b[2] + glm_b[3]))
        # MUB
        g_hat_b <- (phi_guess_b + (1 - phi_guess_b) * r_yz_b)/(phi_guess_b * r_yz_b + (1 - phi_guess_b))
        c_obs_b <- (sd(t_b$inc_adm_20192_std, na.rm = T)) / sd(Z_b, na.rm = T)
        r_err_hat_b <- g_hat_b * c_obs_b * (mean(Z_b[S_b == 1], na.rm = T) - mean(Z_b, na.rm = T))/lambda_b
        
        # MUBM
        g_hat_m_b <- (phi_guess_b * sqrt(gamma_b) + gamma_b * (1 - phi_guess_b) * r_yz_b)/(phi_guess_b * r_yz_b + (1 - phi_guess_b) * sqrt(gamma_b))
        c_obs_m_b <- (sd(t_b$inc_adm_20192_std, na.rm = T) * sqrt(gamma_b)) / sd(Z_b, na.rm = T)
        rm_err_hat_b <- g_hat_m_b * c_obs_m_b * (mean(Z_b[S_b == 1], na.rm = T) - mean(Z_b, na.rm = T))/lambda_b
      }
      t_err_hat_b <- m_err_hat_b  + rm_err_hat_b
      boot_m[b] <- m_err_hat_b; boot_r[b] <- r_err_hat_b; boot_rm[b] <- rm_err_hat_b; boot_t[b] <-  t_err_hat_b
      boot_gamma[b] <- gamma_b
      
      # Absolute values
      abs_m[b] <- abs(m_err - m_err_hat_b)
      abs_r[b] <- abs(r_err - r_err_hat_b)
      abs_rm[b] <- abs(r_err - rm_err_hat_b)
      abs_t[b] <- abs(t_err - t_err_hat_b)
    }
    
    BootRes[[length(BootRes)+1]] <- data.frame(sample_pct = samps[i],
                                               abs_m_mean = mean(abs_m, na.rm = T), abs_m_lo = quantile(abs_m, 0.025, na.rm = T), abs_m_hi = quantile(abs_m, 0.975, na.rm = T), 
                                               abs_r_mean = mean(abs_r, na.rm = T), abs_r_lo = quantile(abs_r, 0.025, na.rm = T), abs_r_hi = quantile(abs_r, 0.975, na.rm = T),
                                               abs_rm_mean = mean(abs_rm, na.rm = T), abs_rm_lo = quantile(abs_rm, 0.025, na.rm = T), abs_rm_hi = quantile(abs_rm, 0.975, na.rm = T), 
                                               abs_t_mean = mean(abs_t, na.rm = T), abs_t_lo = quantile(abs_t, 0.025, na.rm = T), abs_t_hi = quantile(abs_t, 0.975, na.rm = T),
                                               gamma_mean = mean(boot_gamma, na.rm = T), gamma_sd = sd(boot_gamma, na.rm = T),
                                               gamma_lo = quantile(boot_gamma, 0.025, na.rm = T), gamma_hi = quantile(boot_gamma, 0.975, na.rm = T), stringsAsFactors = F)
    Res[i+3*(j-1),] <- c(m_err, r_err, t_err, m_err_hat, r_err_hat, rm_err_hat, t_err_hat, Fit)
  }
}

colnames(Res) <- c("m_err", "r_err", "t_err", "m_err_hat", "r_err_hat", "rm_err_hat", "t_err_hat", "chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic", "bic")
method <- c("5% - Monte Carlo", "10% - Monte Carlo", "20% - Monte Carlo",
                   "5% - Naïve", "10% - Naïve", "20% - Naïve",
                   "5% - Adjusted", "10% - Adjusted", "20% - Adjusted")
Res <- as.data.frame(cbind(method, Res))
Res[,-1] <- apply(Res[,-1], 2, as.numeric)

# Plot the results

Res$abe_m <- abs(Res$m_err - Res$m_err_hat)
Res$abe_r <- abs(Res$r_err-Res$r_err_hat)
Res$abe_rm <- abs(Res$r_err-Res$rm_err_hat)
Res$abe_t <- abs(Res$t_err-Res$t_err_hat)

Res_long <- reshape2::melt(Res, id.vars = "method", measure.vars = c("abe_m", "abe_r", "abe_rm", "abe_t"),
                 variable.name = "parameter", value.name = "abs_err")

Res_long <- Res_long %>% 
  separate_wider_delim(method, delim = " - ", names = c("Sample", "Step2")) %>% 
  mutate(Sample = factor(Sample, levels = c("5%", "10%", "20%"))) %>% 
  mutate(parameter = recode(as.character(Res_long$parameter), abe_m = "Measurement", abe_r = "Representation (MUB)", abe_rm = "Representation* (MUBM)", abe_t = "Total"))

Res_long$parameter <- factor(Res_long$parameter, levels = c("Measurement", "Representation (MUB)", "Representation* (MUBM)", "Total"))

Res_plot <- ggplot(Res_long, aes(x=fct_rev(Step2), y=abs_err, fill=fct_rev(parameter))) +
  geom_bar(stat="identity", position="dodge", alpha = 0.7) +
  scale_fill_brewer(palette = "Set1") + 
  coord_flip() +
  facet_grid(. ~ Sample, labeller = label_both) +
  guides(fill = guide_legend(reverse = T)) + 
  labs(y="Absolute Deviation", x="", title="", fill = "Error Type") + 
  theme(text = element_text(size = 20))
Res_plot
ggsave("Res_IT.png", plot = Res_plot, path = fold_graphs, width = 40, height = 15, units = "cm")

setwd(fold_graphs)
Res_tab <- Res %>% 
  select("m_err","r_err", "t_err", "method", "chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr") %>% 
  separate_wider_delim(method, delim = " - ", names = c("Sample", "Step 2"))
Res_tab <- Res_tab[, c(5,4,1,2,3,6,7,8,11,12,9,10)]
tab_fit <- xtable(Res_tab) 
print(tab_fit, file = "Table_fit_IT.tex", include.rownames = F)  

# Build long data with means and CI already aligned

BootRes_df <- do.call(rbind, BootRes)    # each entry was a data.frame; this row-binds
method2 <- c("5% - Monte Carlo","5% - Naïve","5% - Adjusted",
            "10% - Monte Carlo", "10% - Naïve","10% - Adjusted",
            "20% - Monte Carlo", "20% - Naïve","20% - Adjusted")
BootRes_df <- as.data.frame(BootRes_df, stringsAsFactors = FALSE)
BootRes_df <- as.data.frame(cbind(method2, BootRes_df))
BootRes_long <- reshape2::melt(BootRes_df, id.vars = "method2", measure.vars = c("abs_m_mean","abs_m_lo", "abs_m_hi", "abs_r_mean", "abs_r_lo", "abs_r_hi","abs_rm_mean",
                                                                      "abs_rm_lo", "abs_rm_hi", "abs_t_mean","abs_t_lo", "abs_t_hi"),
                           variable.name = "parameter", value.name = "abs_err") %>% 
  separate_wider_delim(method2, delim = " - ", names = c("Sample", "Step2")) %>% 
  mutate(Sample = factor(Sample, levels = c("5%", "10%", "20%"))) %>% 
  mutate(error_type = str_extract(parameter, "(?<=abs_)[a-z]+"), stat = str_extract(parameter, "(mean|lo|hi)")) %>%
  select(Step2, Sample, abs_err, error_type, stat) %>%
  pivot_wider(names_from = stat, values_from = abs_err) 

BootRes_long <- BootRes_long %>% 
  mutate(error_type = recode(as.character(BootRes_long$error_type), m = "Measurement", r = "Representation (MUB)", rm = "Representation* (MUBM)", t = "Total"))

boot_plot <- ggplot(BootRes_long, aes(x = fct_rev(Step2), y = mean, fill = fct_rev(error_type))) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.75) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.8), size = 0.3, colour = "gray40") +
  scale_fill_brewer(palette = "Set1") +
  coord_flip() +
  facet_grid(. ~ Sample, labeller = label_both) +
  guides(fill = guide_legend(reverse = T)) + 
  labs(x = "", y = "Absolute Error (mean ± 95% CI)",
       fill = "Error type") +
  theme_minimal(base_size = 16) +
  theme(text = element_text(size = 20),
    panel.border = element_rect(color = "gray70", fill = NA, size = 0.5),
    strip.background = element_rect(fill = "gray90", color = NA),
    panel.spacing = unit(0.5, "lines"), legend.position = "right"
  )
boot_plot

ggsave("Res_IT_bootstrap.png", plot = boot_plot,
       path = fold_graphs, width = 40, height = 15, units = "cm")

# Check the gamma conditions
gamma_tab <- BootRes_df[, c("method2", "sample_pct", 
                    "gamma_mean", "gamma_sd", "gamma_lo", "gamma_hi")]

# compute ratios
gamma_tab$ratio_mean <- (1 - gamma_tab$gamma_mean) / gamma_tab$gamma_sd
gamma_tab$ratio_lo   <- (1 - gamma_tab$gamma_lo)/gamma_tab$gamma_sd 
gamma_tab$ratio_hi   <- (1 - gamma_tab$gamma_hi)/gamma_tab$gamma_sd 

# round for readability
gamma_tab <- transform(gamma_tab,
                       gamma_mean = round(gamma_mean, 3),
                       gamma_sd   = round(gamma_sd, 3),
                       gamma_lo   = round(gamma_lo, 3),
                       gamma_hi   = round(gamma_hi, 3),
                       ratio_mean = round(ratio_mean, 3),
                       ratio_lo   = round(ratio_lo, 3),
                       ratio_hi   = round(ratio_hi, 3))


