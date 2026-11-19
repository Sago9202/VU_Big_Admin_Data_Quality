## --------------------------------------------------------------- ##
## Departure from Ignorable Sample Selection with Ordinal Outcomes ##
## Script 2 - Application                                          ##
## Written by: Santiago Gómez-Echeverry                            ##
## Last modified: 11/11/2026                                       ##
## --------------------------------------------------------------- ##

#### - (I) Working space and packages - ####

rm(list = ls()) # Remove everything from the working environment
packs <- c('data.table', 'ggplot2', 'depmixS4', 'tidyr', 'dplyr', 'forcats', 'reshape2', 'scales', 'MASS', 'polycor', 
           'stringr', 'ordinalNet', 'progress', 'knitr', 'MetBrewer')
ipacks <- packs  %in% rownames(installed.packages()) 
if(any(ipacks == F)){install.packages(packs[!ipacks])}
invisible(lapply(packs, library, character.only = T))
fold_data <- "C:/Users/sagom/Dropbox/PhD VU/Italian Data/Total Error"
fold_graphs <- "C:/Users/sagom/Dropbox/PhD VU/Papers/5 - Ordinal MUB/4 - Graphs and Tables"
setwd(fold_data)

plt_cont <- met.brewer("Benedictus", n=100)

#### - (II) Data arrangement - ####

full_dat <- fread("TotErr_Rel20.csv") # Entire data set
N <- nrow(full_dat)

# Arranging the auxiliary variables - Keep in mind that they should come from adm
arr_dat <- full_dat[full_dat$age31_12_rbiadm >= 18,]  # Only adults
arr_dat$gender_rbiadm <- factor(arr_dat$gender_rbiadm, levels = 1:2, labels = c('Male', 'Female'))
arr_dat$region_adm <- factor(arr_dat$region_adm)
arr_dat$nationality_rbiadm <- factor(arr_dat$nationality_rbiadm)
arr_dat$S <- ifelse(is.na(arr_dat$net_wage_class), NA, 1)
arr_dat$S[arr_dat$eligibility==1 & is.na(arr_dat$S)] <- 0
table(arr_dat$S, useNA = "ifany")
arr_dat$net_wage_class[arr_dat$net_wage_class==997] <- NA
arr_dat$net_wage_class[arr_dat$net_wage_class==998] <- NA

table(arr_dat$survey_year, arr_dat$net_wage_class)
# arr_dat$net_wage_class[(arr_dat$net_wage_class >=8 & arr_dat$survey_year>2017)] <- 8
# table(arr_dat$survey_year, arr_dat$net_wage_class)
K <- length(table(arr_dat$net_wage_class))
K
arr_dat$net_wage_class <- factor(arr_dat$net_wage_class)
is.factor(arr_dat$net_wage_class)
arr_dat$nm_Fcont <- str_count(arr_dat$ftpt12m_adm, "F")
arr_dat$nm_Pcont <- str_count(arr_dat$ftpt12m_adm, "P")
table(arr_dat$nm_Fcont)

edu_labels <- c("No higher than primary", "Lower-secondary", "Upper-secondary", "Bachelor's or short-cycle tertiary education",
                "Master's or equivalent level", "Doctoral or equivalent level")
# The education variables are: hatlevel_lfs for the survey and edu_isced_rbiadm in the admin data
arr_dat <- arr_dat %>% 
  mutate(edu_lev_lfs = case_match(hatlevel_lfs, 100 ~ 1L, c(303,200,343,353,393)  ~ 2L, c(304,400,344,354,440,450) ~ 3L,
                                  c(500, 540, 550, 600) ~ 4L, 700 ~ 5L, 800 ~ 6L)) %>% 
  mutate(edu_lev_lfs = factor(edu_lev_lfs, labels = edu_labels)) %>% 
  mutate(edu_lev_adm = case_match(edu_isced_rbiadm, c(1,2,3) ~ 1L, 4 ~ 2L, 5 ~ 3L, 6 ~ 4L, 7 ~ 5L, 8 ~ 6L)) %>% 
  mutate(edu_lev_adm = factor(edu_lev_adm, labels = edu_labels)) %>% 
  as.data.frame()

is.factor(arr_dat$edu_lev_lfs)         # Should be T
prop_true <- prop.table(table(arr_dat$edu_lev_adm)) # These are the true proportions we will check against with

padd <- arr_dat %>%
  filter(!is.na(survey_year), !is.na(net_wage_class)) %>%
  count(survey_year, net_wage_class) %>%
  group_by(survey_year) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = factor(survey_year), y = prop, fill = net_wage_class)) +
  geom_col(position = "fill", width = 0.7) + scale_fill_manual(values = met.brewer("Benedictus", n = 10)) +
  scale_y_continuous(labels = scales::percent) + labs(x = "Year", y = "Proportion",fill = "Net wage class") +
  theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank(), text = element_text(size = 20))

ggsave("padd.png", plot = padd, path = fold_graphs, width = 25, height = 20, units = "cm")

#### - (III) MUBM function - ####

mubp_k <- function(data, summary = c("arithmetic", "closed_geometric")){
  # 1. Split the data into sample and out of sample
  data    <- data[!is.na(data$X), ]
  s_data  <- data[which(data$S == 1), ]
  ns_data <- data[which(data$S == 0), ]
  s_data$Y <- droplevels(s_data$Y)
  
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
  valid <- !is.na(sigma0_uu) & sigma0_uu > 0
  m_valid <- sum(valid)
  m_dropped <- M - m_valid
  
  if (m_valid == 0) {
    out <- setNames(rep(NA_real_, K), paste0("C", 1:K))
    attr(out, "m_dropped") <- M
    attr(out, "mubp_phi0") <- attr(out, "mubp_phi1") <- rep(NA_real_, K)
    attr(out, "mubp_min")  <- attr(out, "mubp_max")  <- rep(NA_real_, K)
    warning("All draws of phi were invalid; check the proxy moments.")
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
  
  # 8. MUBP at the extremes
  mubp_at <- function(phi_vals) {
    g  <- (phi_vals + (1 - phi_vals) * rho_ux) / (phi_vals * rho_ux + (1 - phi_vals))
    mu <- g * ((mu0_x - mu1_x) / sqrt(sigma1_xx))
    s2 <- 1 + g^2 * ((sigma0_xx - sigma1_xx) / sigma1_xx)
    out_phi <- matrix(NA_real_, length(phi_vals), K)
    ok <- s2 > 0
    if (any(ok)) {
      z0 <- outer(1 / sqrt(s2[ok]), all_taus) - outer(mu[ok] / sqrt(s2[ok]), rep(1, K + 1))
      p0 <- pnorm(z0)
      v0 <- p0[, -1, drop = FALSE] - p0[, -(K + 1), drop = FALSE]
      vpop <- sweep((1 - pi_sample) * v0, 2, pi_sample * v_k_1, "+")
      out_phi[ok, ] <- sweep(-vpop, 2, v_k_1, "+")
    }
    out_phi
  }
  ext <- mubp_at(c(0, 1))
  mubp_draws <- sweep(-v_k_pop_mat, 2, v_k_1, "+")
  
  # 8. MUBP
  out <- v_k_1 - v_k_pop
  names(out) <- paste0("C", levels(s_data$Y))
  attr(out, "m_dropped") <- m_dropped
  attr(out, "summary") <- summary
  attr(out, "mubp_phi0") <- ext[1, ]
  attr(out, "mubp_phi1") <- ext[2, ]
  attr(out, "mubp_min")  <- apply(mubp_draws, 2, min)
  attr(out, "mubp_max")  <- apply(mubp_draws, 2, max)
  return(out)
}

#### - (IV) Estimation - ####

## - Obtain the proxy - ##
# We will consider two specifications: (1) standard ordinal regression and (2) lasso penalized ordinal regression

# # (1) Standard oprobit
# reg_proxy <- polr(net_wage_class ~ edu_lev_adm + age31_12_rbiadm + nationality_rbiadm + hhsize_adm  + gender_rbiadm + 
#                     region_adm + Sqroot_age_ES_adm + nm_Pcont + nm_Fcont, data = arr_dat, method = "probit")
# temp_proxy <- reg_proxy$lp                # Linear predictor
# y_miss <- reg_proxy$model$net_wage_class  # Y with missings
# cat("Correlation standard proxy:", polyserial(temp_proxy, y_miss))
# 
# aux_vars <- colnames(reg_proxy$model)[-1] # Auxiliary variables
# arr_dat <- arr_dat[complete.cases(arr_dat[, aux_vars]), ]
# aux_mat <- model.matrix(formula(reg_proxy)[-2], data = arr_dat)
# nrow(aux_mat)
# betas <- coef(reg_proxy)
# arr_dat$proxy1 <- aux_mat[, names(betas)] %*% betas 
# 
# # (2) Lasso oprobit
# preds <- c("edu_lev_adm", "age31_12_rbiadm", "nationality_rbiadm",  "hhsize_adm", "gender_rbiadm", "region_adm","Sqroot_age_ES_adm", "nm_Pcont", "nm_Fcont")
# response <- "net_wage_class"
# complete_idx <- complete.cases(arr_dat[, c(preds, response)])
# dat_clean <- arr_dat[complete_idx, ]
# form <- as.formula(paste("~", paste(preds, collapse = " + ")))
# X <- model.matrix(form, data = dat_clean)[, -1]  
# X <- as.matrix(X) 
# y <- dat_clean[[response]]
# 
# fit_proxy <- ordinalNet(x = X, y = y, family = "cumulative", link = "probit", alpha = 0, reverse = T)
# all_coef <- coef(fit_proxy)
# beta <- all_coef[-(1:n_thresh)]
# proxy <- as.vector(X %*% beta)
# cat("Correlation LASSO proxy:", polyserial(proxy, y))
# 
# arr_dat$proxy2 <- NA
# arr_dat$proxy2[complete_idx] <- proxy

yrs <- unique(arr_dat$survey_year)
nyrs <- length(yrs)
nyrs
Results <- vector(mode = "list", length = nyrs)

for (t in 1:nyrs){
  t_yr <- yrs[t]
  t_dat <- arr_dat[arr_dat$survey_year == t_yr & !is.na(arr_dat$S),]
  
  preds <- c("edu_lev_adm", "age31_12_rbiadm", "nationality_rbiadm", "hhsize_adm", "gender_rbiadm", "region_adm", "Sqroot_age_ES_adm", "nm_Pcont", "nm_Fcont")
  response <- "net_wage_class"
  complete_idx <- complete.cases(t_dat[, preds])
  fit_idx <- complete.cases(t_dat[, c(preds, response)])
  X <- model.matrix(as.formula(paste("~", paste(preds, collapse = " + "))), data = t_dat[complete_idx, ])[, -1]   
  X <- as.matrix(X)
  X_fit <- model.matrix(as.formula(paste("~", paste(preds, collapse = " + "))), data = t_dat[fit_idx, ])[, -1]
  X_fit <- as.matrix(X_fit)
  y <- t_dat$net_wage_class[fit_idx]
  fit_proxy <- ordinalNet(x = X_fit, y = y, family = "cumulative", link = "probit", alpha = 0, reverse = T)
  n_thresh <- length(levels(y)) - 1         
  beta <- coef(fit_proxy)[-(1:n_thresh)]
  common_cols <- intersect(colnames(X), colnames(X_fit))
  X <- X[, common_cols, drop = F]
  proxy <- as.vector(X %*% beta)
  proxy_fit <- as.vector(X_fit %*% beta) 
  y_obs_k <- droplevels(t_dat$net_wage_class[fit_idx])
  rho_ux <- polyserial(proxy_fit, y_obs_k)
  cat("Correlation LASSO proxy:", rho_ux, "\n")
  t_dat$X <- NA
  t_dat$X[complete_idx] <- proxy
  
  # 1. Ensure that the target variable is a factor
  t_dat$Y <- factor(t_dat$net_wage_class, levels = 1:K)
  
  # 2. Estimate the MUBP
  mubp_est <- mubp_k(t_dat, summary = "closed_geometric")
  m_dropped <- attr(mubp_est, "m_dropped")
  
  # 4. Store all the information
  Results[[t]] <- data.frame(t = t_yr, Cat = names(mubp_est), mubp_est = as.numeric(mubp_est), 
                             phi0 = attr(mubp_est, "mubp_phi0"), phi1 = attr(mubp_est, "mubp_phi1"), M_drop = m_dropped)
  cat("Year ", t_yr, " done \n")
}

app_results <- bind_rows(Results)

#### - (V) Graphs and Tables - ####

cor_df <- data.frame(t = 2017:2021, cor_val = c(0.7510348, 0.7307699, 0.7336984, 0.716343, 0.7262942)) %>% 
  mutate(label = sprintf("%.3f", cor_val))

# Ensure Cat is a factor with correct order (C1 to C8)
app_results$Cat <- factor(app_results$Cat, levels = paste0("C", 1:K))
app_results$mubp_est[app_results$Cat=="C9" & app_results$t==2017] <- NA
app_results$mubp_est[app_results$Cat=="C10" & app_results$t==2017] <- NA

# Define threshold (e.g., 0.05)
threshold <- 0.05

p_tile <- ggplot(app_results, aes(x = t, y = Cat, fill = mubp_est)) +
  geom_tile(aes(width = 1, height = 1), color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(mubp_est), "", percent(mubp_est, accuracy = 0.1))), color = "black", size = 6) +
  scale_fill_gradientn(colours = c(plt_cont[1], plt_cont[50], plt_cont[100]), limits  = c(-0.12, 0.12),
                       oob     = squish, na.value = "gray70", name    = "MUBP") +
  labs(x = "Year", y = "Income category") +
  theme_minimal() +
  scale_x_continuous(breaks = cor_df$t, sec.axis = sec_axis(transform = ~., breaks = cor_df$t, labels = cor_df$label, name = expression(rho[ux]))) +
  theme(panel.grid = element_blank(), text = element_text(size = 20))

sig_cells <- app_results %>% filter(!is.na(mubp_est), abs(mubp_est) > threshold)
if (nrow(sig_cells) > 0) {
  p_tile <- p_tile +
    geom_tile(data = sig_cells, aes(x = t, y = Cat, width = 1, height = 1), fill = NA, color = "black", linewidth = 1.5)
}
p_tile
# Manually save as 800x1200
ggsave("p9.png", plot = p_tile, path = fold_graphs, width = 20, height = 30, units = "cm")

## Aggregate results
w1_by_year <- app_results |> 
  mutate(cat_num = as.integer(sub("C", "", Cat))) |>  
  arrange(t, cat_num) |> 
  group_by(t) |> 
  summarise(W1 = sum(abs(cumsum(mubp_est)[-n()])), .groups = "drop")

tab <- app_results |> 
  group_by(t) |> 
  summarise(mean_abs = mean(abs(mubp_est)), max_abs = max(abs(mubp_est)), rms = sqrt(mean(mubp_est^2))) |> 
  left_join(w1_by_year, by = "t") |> 
  mutate(Year = t, `Mean Abs` = round(mean_abs, 4), `Max Abs` = round(max_abs, 4), RMS = round(rms, 4), `W1` = round(W1, 4)) |> 
  select(Year, `Mean Abs`, `Max Abs`, RMS, `W1`)

latex_table <- kable(tab, format = "latex", booktabs = T, caption = "Aggregated MUBP per year", label = "1", digits = 4)
latex_table
writeLines(latex_table, paste0(fold_graphs, "/Tab_1.tex"))

## Mean and ranges plot

int_df <- app_results |> 
  mutate(Cat = factor(Cat, levels = paste0("C", 1:K)), lo  = pmin(phi0, phi1), hi = pmax(phi0, phi1)) |> 
  filter(!(t == 2017 & Cat %in% c("C9", "C10"))) |> 
  left_join(cor_df, by = "t") |> 
  mutate(t_lab = paste0("atop('", t, "', rho[ux] == '", label, "')"))

cat_cols <- setNames(as.character(met.brewer("Benedictus", n = 10)), paste0("C", 1:10))

p_int <- ggplot(int_df, aes(y = Cat, colour = Cat)) +
  annotate("rect", xmin = -0.05, xmax = 0.05, ymin = -Inf, ymax = Inf, fill = "grey85") +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_segment(aes(x = lo, xend = hi, yend = Cat), linewidth = 1.2) +
  geom_point(aes(x = mubp_est, fill = Cat), shape = 21, colour = "grey20", size = 2.5, stroke = 0.4) +
  geom_text(aes(x = mubp_est, label = scales::percent(mubp_est, accuracy = 0.1)), colour = "grey20", size = 3, vjust = -0.9) +
  facet_wrap(~ t_lab, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = cat_cols, guide = "none") +
  scale_fill_manual(values = cat_cols, guide = "none") +
  scale_x_continuous(breaks = seq(-0.10, 0.10, by = 0.05), labels = scales::percent_format(accuracy = 1), guide = guide_axis(angle = 45)) +
  coord_cartesian(xlim = c(-0.15, 0.15)) +
  labs(x = "MUBP", y = "Income category") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA), strip.text = element_text(margin = margin(4, 0, 4, 0)))
p_int
ggsave("p9_intervals.png", plot = p_int, path = fold_graphs, width = 30, height = 15, units = "cm")
