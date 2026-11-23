## ------------------------------------------------- ##
## Bayesian Inference for proportions with ME and SE ##
## Script 2: Plots Simulation Study                  ##
## Written by: Santiago Gómez-Echeverry              ##
## Last update: 20/10/2026                           ##
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

res <- readRDS(file.path(fold_data, "results_all_models.rds"))

### - (II) RMSE and Coverage plots across models - ###

models     <- c("selection", "measurement", "joint")
model_labs <- c(selection = "Selection-only", measurement = "Measurement-only", joint = "Composite")
model_cols <- c(selection = "#2C7BB6", measurement = "#D7301F", joint = "#7B3294")
factors    <- c("J", "Delta_RY", "pi", "N", "delta_Z")
factor_lab <- c(J = "J[q]", Delta_RY = "Delta[Y[q]]", pi = "pi", N = "N", delta_Z = "delta[Z]")

base <- res %>%
  filter(prior == "base")

est <- bind_rows(base %>% filter(model == "selection",   variable == "ES") %>% mutate(quantity = "EC", true = true_EC),
                 base %>% filter(model == "measurement", variable == "EM") %>% mutate(quantity = "EC", true = true_EC),
                 base %>% filter(model == "joint",       variable == "EC") %>% mutate(quantity = "EC", true = true_EC),
                 # Comparsons with the separate errors
                 base %>% filter(variable == "ES", model %in% c("selection", "joint"))   %>% mutate(quantity = "ES", true = true_value),
                 base %>% filter(variable == "EM", model %in% c("measurement", "joint")) %>% mutate(quantity = "EM", true = true_value)) %>%
  mutate(model = factor(model, levels = models))

perf <- est %>%
  mutate(sq_err = (mean - true)^2, cover  = q5 <= true & true <= q95) %>%
  group_by(model, quantity, condition, across(all_of(factors))) %>%
  summarise(n_sim = n(), mse   = mean(sq_err), coverage = mean(cover), .groups   = "drop") %>%
  mutate(rmse      = sqrt(mse))

theme_res <- theme_classic(base_size = 12) +
  theme(strip.background = element_blank(), strip.placement = "outside",
        axis.text = element_text(size = 20), legend.text = element_text(size = 20),                       
        legend.key.size = unit(1.2, "lines"), 
        strip.text = element_text(face = "bold", size = 20), legend.position = "top", legend.title = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.spacing.x = unit(1, "lines"), plot.caption = element_text(colour = "grey40", hjust = 0))

to_long <- function(d) {
  d %>%
    pivot_longer(all_of(factors), names_to = "factor", values_to = "lev") %>%
    mutate(lev  = factor(as.character(lev), levels = as.character(sort(unique(lev)))),
           factor = factor(factor, levels = factors, labels = factor_lab[factors]))
}

R_rep    <- median(perf$n_sim)
cov_band <- 0.90 + c(-2, 2) * sqrt(0.9 * 0.1 / R_rep)

plot_perf<- function(d, cols = model_cols, labs_leg = model_labs) {
  dodge <- position_dodge(width = 0.8)
  dl    <- to_long(d)
  p_rmse <- ggplot(dl, aes(lev, rmse, colour = model, fill = model)) +
    geom_boxplot(position = dodge, width = 0.7, alpha = 0.25, outlier.size = 0.8, linewidth = 0.4) +
    facet_wrap(~ factor, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = cols, labels = labs_leg) +
    scale_fill_manual(values = cols, labels = labs_leg) +
    labs(x = NULL, y = "RMSE") + theme_res + theme(strip.text = element_blank())
  
  avg <- dl %>%
    group_by(model, factor, lev) %>%
    summarise(coverage = mean(coverage), .groups = "drop")
  
  p_cov <- ggplot(dl, aes(lev, coverage, colour = model, fill = model)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = cov_band[1], ymax = cov_band[2], fill = "grey88", colour = NA) +
    geom_hline(yintercept = 0.90, colour = "grey40", linewidth = 0.4) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.8, seed = 1),
               size = 0.9, alpha = 0.45, shape = 16) +
    geom_point(data = avg, position = dodge, shape = 21, size = 3, colour = "black", stroke = 0.6) +
    facet_wrap(~ factor, nrow = 1, scales = "free_x", strip.position = "bottom", labeller = label_parsed) +
    scale_colour_manual(values = cols, labels = labs_leg) + scale_fill_manual(values = cols, labels = labs_leg) +
    labs(x = NULL, y = "Coverage of the 90% interval") + theme_res + theme(legend.position = "none")
  
  p_rmse / p_cov
}

p1_comp <- perf |>  
  filter(quantity == "EC") |> 
  plot_perf()
ggsave("p1_comp.png", plot = p1_comp, path = fold_graphs, width = 50, height = 20, units = "cm")


p1_sel <- perf |>  
  filter(quantity == "ES") |> 
  plot_perf()
ggsave("p1_sel.png", plot = p1_sel, path = fold_graphs, width = 50, height = 20, units = "cm")

p1_mes <- perf |>  
  filter(quantity == "EM") |> 
  plot_perf()
ggsave("p1_mes.png", plot = p1_mes, path = fold_graphs, width = 50, height = 20, units = "cm")


### - (III) RMSE and Coverage plots across priors - ###

priors     <- c("weak", "base", "strong")
prior_labs <- c(weak = "Weak", base = "Base", strong = "Strong")
prior_cols <- c(weak = "#C2A5CF", base = "#9970AB", strong = "#40004B")

est <- res %>%
  filter(model == "joint", variable %in% c("EC", "ES", "EM")) %>%
  mutate(quantity = variable,
         true     = ifelse(variable == "EC", true_EC, true_value),
         model    = factor(prior, levels = priors))

perf <- est %>%
  mutate(sq_err = (mean - true)^2, cover  = q5 <= true & true <= q95) %>%
  group_by(model, quantity, condition, across(all_of(factors))) %>%
  summarise(n_sim = n(), mse   = mean(sq_err), coverage = mean(cover), .groups   = "drop") %>%
  mutate(rmse      = sqrt(mse))

p2_comp <- perf |>  
  filter(quantity == "EC") |> 
  plot_perf(cols = prior_cols, labs_leg = prior_labs)
ggsave("p2_comp.png", plot = p2_comp, path = fold_graphs, width = 50, height = 20, units = "cm")

p2_sel <- perf |>  
  filter(quantity == "ES") |> 
  plot_perf(cols = prior_cols, labs_leg = prior_labs)
ggsave("p2_sel.png", plot = p2_sel, path = fold_graphs, width = 50, height = 20, units = "cm")

p2_mes <- perf |>  
  filter(quantity == "EM") |> 
  plot_perf(cols = prior_cols, labs_leg = prior_labs)
ggsave("p2_mes.png", plot = p2_mes, path = fold_graphs, width = 50, height = 20, units = "cm")

### - (III) Allocaiton and Composite Error Bias - ###

f_q        <- 0.2       # sampling fraction in the design
priors     <- c("weak", "base", "strong")
prior_labs <- c(weak = "Weak", base = "Base", strong = "Strong")
prior_cols <- c(weak = "#C2A5CF", base = "#9970AB", strong = "#40004B")

## ---- 1. Bias per condition and prior: allocation and composite -------------
cond_err <- res %>%
  filter(model == "joint", variable %in% c("EM", "ES", "EC")) %>%
  select(prior, condition, J, Delta_RY, pi, N, delta_Z, draw, variable, mean,
         true_EM, true_ES, true_EC) %>%
  pivot_wider(names_from = variable, values_from = mean) %>%
  mutate(A  = ((ES - true_ES) - (EM - true_EM)) / 2,
         dC = EC - true_EC) %>%
  group_by(prior, condition, J, Delta_RY, pi, N, delta_Z) %>%
  summarise(true_ES = first(true_ES),
            Allocation = mean(A),
            Composite  = mean(dC),
            .groups = "drop") %>%
  mutate(share = true_ES / pmin(1 - pi, pi * (1 - f_q) / f_q)) %>%
  pivot_longer(c(Allocation, Composite), names_to = "metric", values_to = "bias") %>%
  mutate(prior  = factor(prior, levels = priors),
         metric = factor(metric, levels = c("Allocation", "Composite"),
                         labels = c("Bias of the allocation", "Bias of the composite error")),
         N_lab  = factor(paste0("N = ", N), levels = paste0("N = ", sort(unique(N)))))

## ---- 2. Plot ----------------------------------------------------------------
# Faint points: single conditions. Lines: smoothed average per prior (loess).
p_prior_dynamics <- ggplot(cond_err, aes(share, bias, colour = prior)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_point(alpha = 0.25, size = 1) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 1.1) +
  facet_grid(metric ~ N_lab, switch = "y") +              # fixed scales: same vertical axis in both rows
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_colour_manual(values = prior_cols, labels = prior_labs) +
  labs(x = expression("Intensity of selection:  " * E[S] / max ~ E[S]), y = NULL) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(), strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        legend.position = "top", legend.title = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.spacing = unit(1, "lines"), plot.caption = element_text(colour = "grey40", hjust = 0))

print(p_prior_dynamics)
ggsave("p3_priors.png", plot = p_prior_dynamics, path = fold_graphs, width = 30, height = 20, units = "cm")

## ---------------------------------------------------------------------- ##
## Convergence diagnostics                                                ##
##   1. Table: share of fits flagged by each criterion, per model / prior ##
##   2. Tile plot: share of fits flagged, criterion x model, per N        ##
##   3. Tile plots, one per model: J x N panels, Delta_RY x pi tiles      ##
##                                                                        ##
## A fit is flagged if any of:                                            ##
##   divergent transitions > 0, max R-hat > 1.01, min E-BFMI < 0.2,       ##
##   min bulk ESS < 400   (the criteria used in the simulation)           ##
## Tail ESS < 400 and max-treedepth hits are reported but not used for    ##
## the flag.                                                              ##
## ---------------------------------------------------------------------- ##
library(dplyr); library(tidyr); library(ggplot2)

# fold_data: folder with the simulation results (as in the simulation script)
diag <- readRDS(file.path(fold_data, "diagnostics_all_models.rds"))

factors    <- c("J", "Delta_RY", "pi", "N", "delta_Z")
factor_lab <- c(J = "J[q]", Delta_RY = "Delta[RY]", pi = "pi", N = "N", delta_Z = "delta[Z]")

# One series per model and prior
series_lev <- c("selection base", "measurement base", "joint weak", "joint base", "joint strong")
series_lab <- c("selection base"   = "Selection-only",
                "measurement base" = "Measurement-only",
                "joint weak"       = "Composite, weak prior",
                "joint base"       = "Composite, base prior",
                "joint strong"     = "Composite, strong prior")
series_col <- c("selection base" = "#2C7BB6", "measurement base" = "#D7301F",
                "joint weak" = "#C2A5CF", "joint base" = "#9970AB", "joint strong" = "#40004B")

d <- diag %>%
  mutate(series   = factor(paste(model, prior), levels = series_lev),
         failed   = coalesce(as.logical(failed), FALSE),
         f_div    = n_divergent > 0,
         f_rhat   = max_rhat > 1.01,
         f_ebfmi  = min_ebfmi < 0.2,
         f_bulk   = min_ess_bulk < 400,
         f_tail   = min_ess_tail < 400,
         f_tree   = n_max_treedepth > 0,
         flagged  = f_div | f_rhat | f_ebfmi | f_bulk)

## ---- 1. Table -------------------------------------------------------------------
# Percentages are over the fits that ran (failed fits are counted separately).
pct <- function(x) 100 * mean(x, na.rm = TRUE)

diag_table <- d %>%
  group_by(series) %>%
  summarise(Fits              = n(),
            Failed            = sum(failed, na.rm = TRUE),
            `Divergences (%)` = pct(f_div[!failed]),
            `R-hat > 1.01 (%)`= pct(f_rhat[!failed]),
            `E-BFMI < 0.2 (%)`= pct(f_ebfmi[!failed]),
            `Bulk ESS < 400 (%)` = pct(f_bulk[!failed]),
            `Tail ESS < 400 (%)` = pct(f_tail[!failed]),
            `Max treedepth (%)`  = pct(f_tree[!failed]),
            `Any flag (%)`    = pct(flagged[!failed]),
            `Max R-hat`       = max(max_rhat[!failed], na.rm = TRUE),
            `Median bulk ESS` = median(min_ess_bulk[!failed], na.rm = TRUE),
            `Median time (s)` = median(time_total_s[!failed], na.rm = TRUE),
            .groups = "drop") %>%
  mutate(series = series_lab[as.character(series)]) %>%
  dplyr::rename(Model = series)

print(as_tibble(diag_table), width = Inf)

# For the paper (LaTeX), if kableExtra is loaded:
# kableExtra::kbl(diag_table, format = "latex", booktabs = TRUE, digits = 2)

## ---- 2. Tile plot: criterion x model, by sample size --------------------------------
# One panel per N; within each panel, the share of fits meeting each criterion per model.
crit_lab <- c(f_div = "Divergences > 0", f_rhat = "R-hat > 1.01", f_ebfmi = "E-BFMI < 0.2",
              f_bulk = "Bulk ESS < 400", f_tail = "Tail ESS < 400", f_tree = "Max treedepth",
              flagged = "Any flag")

tile_crit <- d %>%
  filter(!failed) %>%
  pivot_longer(names(crit_lab), names_to = "crit", values_to = "flag") %>%
  group_by(series, N, crit) %>%
  summarise(pct = pct(flag), .groups = "drop") %>%
  mutate(crit   = factor(crit, levels = names(crit_lab), labels = crit_lab),
         series = factor(series, levels = rev(series_lev), labels = series_lab[rev(series_lev)]),
         N_lab  = factor(paste0("N = ", N), levels = paste0("N = ", sort(unique(N)))))

p_tile_crit <- ggplot(tile_crit, aes(crit, series, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f", pct), colour = pct > 0.5 * max(pct)), size = 3.2,
            show.legend = FALSE) +
  facet_wrap(~ N_lab, nrow = 1) +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white")) +
  scale_fill_gradient(low = "#fff5f0", high = "#a50f15", name = "Fits flagged (%)") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0),   # avoids overlapping labels
        strip.text = element_text(face = "bold"), strip.placement = "outside",
        panel.spacing.x = unit(1.2, "lines"))

print(p_tile_crit)

## ---- 3. Tile plots: where the flags occur, one figure per model -------------------
# Panels: J (rows) x N (columns). Within each panel: Delta_RY (horizontal) x pi (vertical).
# Each tile: share of fits with any flag over the 3 conditions with those five-factor
# levels except delta_Z (pooled), i.e. 3 x replications fits per tile.
# The colour scale is the same in all figures, so shades can be compared across models.
tile_grid <- d %>%
  filter(!failed) %>%
  group_by(series, J, N, Delta_RY, pi) %>%
  summarise(pct = pct(flagged), .groups = "drop") %>%
  mutate(Delta_RY = factor(Delta_RY), pi = factor(pi))
pct_max <- max(tile_grid$pct)

plot_flag_grid <- function(s) {
  ggplot(filter(tile_grid, series == s), aes(Delta_RY, pi, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.0f", pct), colour = pct > 0.5 * pct_max), size = 3,
              show.legend = FALSE) +
    facet_grid(J ~ N, switch = "y",
               labeller = labeller(J = as_labeller(function(x) sprintf("J[q] == %s", x), label_parsed),
                                   N = as_labeller(function(x) sprintf("N == %s", x), label_parsed))) +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white")) +
    scale_fill_gradient(low = "#fff5f0", high = "#a50f15", limits = c(0, pct_max),
                        name = "Fits flagged (%)") +
    labs(x = expression(Delta[Y[q]]), y = expression(pi), title = series_lab[s]) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"),
          strip.placement = "outside", panel.spacing = unit(0.8, "lines"),
          legend.position = "right")
}

for (s in series_lev) print(plot_flag_grid(s))

## ---- 4. Checks on the flagged selection-only fits ----------------------------------
# (a) How severe are the flags? Mostly just past the thresholds, or far beyond them?
d %>%
  filter(model == "selection", !failed, flagged) %>%
  summarise(n = n(),
            rhat_q50 = median(max_rhat),     rhat_q90 = quantile(max_rhat, 0.9),
            rhat_max = max(max_rhat),
            ess_q50  = median(min_ess_bulk), ess_q10  = quantile(min_ess_bulk, 0.1),
            ess_min  = min(min_ess_bulk)) %>%
  print()

# (b) Do the flagged fits give different estimates of E_S from the unflagged ones?
#     Needs the results file; same conditions and replications, flagged vs not.
res <- readRDS(file.path(fold_data, "results_all_models.rds"))
res %>%
  filter(model == "selection", variable == "ES") %>%
  inner_join(d %>% filter(model == "selection") %>% select(condition, draw, flagged),
             by = c("condition", "draw")) %>%
  group_by(N, J, flagged) %>%
  summarise(fits = n(), bias = mean(mean - true_value),
            rmse = sqrt(mean((mean - true_value)^2)), .groups = "drop") %>%
  print(n = Inf)