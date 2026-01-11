# #####################################
# Setup                               #
# Author: Abel Kjaersgaard            #
# Date: Wed Dec 17 2025               #
# #####################################

# 1. global variables ---------------------------------------
sims <- 10000 # number of simulated draws used to calculate threshold values

# 2. dependencies ---------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
})

# 3. helper functions ---------------------------------------
## function: calculate posterior from beta-binomial model ---------------------
rbetabinom <- function(samples, x, n, alpha, beta) {
  return(rbeta(n = samples, shape1 = x + alpha, shape2 = n - x + beta))
}

## function: comment -----------------------------
comment <- function(...) invisible(NULL)

## function: run threshold model on dataframe -----------------------------
calc_threshold <- function(
  n_sims = sims,
  pa_factor = 1,
  ps_factor = 1,
  prev_factor = 1,
  prev_true_col = "prev_true",
  prev_false_col = "prev_false",
  N_col = "mdrrr_tb",
  Tr_col = "lag_tr"
) {
  df <- dplyr::pick(everything())

  pmap(
    list(
      prev_true = df[[prev_true_col]],
      prev_false = df[[prev_false_col]],
      N = df[[N_col]],
      Tr = df[[Tr_col]],
      pa_factor = pa_factor,
      ps_factor = ps_factor,
      prev_factor = prev_factor,
      n_sims = n_sims
    ),
    threshold_model
  )
}

## function: plot output from threshold model -----------------------------
sa_plot <- function(df) {
  df |>
    ggplot(aes(
      x = year,
      y = thres_val,
      ymin = thres_val_lo,
      ymax = thres_val_hi,
      shape = sa,
      color = sa_change
    )) +
    geom_point() +
    geom_errorbar() +
    geom_hline(aes(yintercept = 1), color = "red", linetype = "dashed") +
    facet_wrap(~country, scales = "free_y", ncol = 3) +
    scale_y_log10() +
    scale_x_continuous(breaks = seq(2015, 2023, 2)) +
    scale_color_manual(values = c("black", "red")) +
    scale_shape_manual(name = "", values = c(16, 17)) +
    guides(color = "none") +
    labs(x = "Year", y = "Transmission threshold values") +
    theme_bw() +
    theme(legend.position = "top")
}

# 4. transmission threshold function ---------------------------------------
#fmt: skip
threshold_model <- function(
  # data inputs
  prev_true = NULL,    # number of BDQR for calculating BDQR prevalence 
  prev_false = NULL,   # number of BDQS for calculating BDQR prevalence
  N,                   # number of MDR/RR-TB 
  Tr,                  # number of MDR/RR-TB treated with BDQ

  # model parameters
  n_sims = sims,         # number of simulations to run per observation

  # sensitivity analysis parameters
  pa_factor = 1,       # multiplier for probability of acquired resistance
  ps_factor = 1,       # multiplier for probability of spontaneous resistance
  prev_factor = 1      # multiplier for prevalence of BDQR among MDR/RR-TB
) {
  
  # check inputs
  if (is.na(N) | is.na(Tr) | is.na(prev_true) | is.na(prev_false)) {
    stop("all inputs must not be NA")
  }

  # compute spontaneous resistance cases
  # beta-binomial model assuming flat prior (alpha = beta = 1). Successes and
  # trials (ps_x and ps_n) are calculated in 1_parameter_setup.R
  N_spontaneous <- rbetabinom(n_sims, ps_x, ps_n, 1, 1) * N * ps_factor

  # compute acquired resistance cases
  # abr_samples is a pre-computed vector of acquired resistance probabilities
  # (see 1_parameter_setup.R)
  abr_samples <- sample(abr_samples) * pa_factor
  N_acquired = abr_samples * (Tr - (N_spontaneous * (Tr / N)))

  # compute observed resistance
  x = round(prev_true)
  n = round(prev_true + prev_false)
  prev <- rbetabinom(n_sims, x, n, 1, 25)
  N_obs <- N * prev * prev_factor
  
  # compute threshold value
  thres_val = N_obs / (N_spontaneous + N_acquired)

  return(list(
    "thres_val" = median(thres_val),
    "thres_val_lo" = quantile(thres_val, 0.025),
    "thres_val_hi" = quantile(thres_val, 0.975),
    "N_spontaneous" = mean(N_spontaneous),
    "Ns_lo" = quantile(N_spontaneous, 0.025),
    "Ns_hi" = quantile(N_spontaneous, 0.975),
    "N_acquired" = mean(N_acquired),
    "Na_lo" = quantile(N_acquired, 0.025),
    "Na_hi" = quantile(N_acquired, 0.975),
    "N_obs" = mean(N_obs),
    "No_lo" = quantile(N_obs, 0.025),
    "No_hi" = quantile(N_obs, 0.975)
  ))
}
