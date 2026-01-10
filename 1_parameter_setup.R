# #####################################
# Parameter setup                     #
# Author: Abel Kjaersgaard            #
# Date: Wed Dec 17 2025               #
# #####################################

# 1. dependencies ---------------------------------------
suppressPackageStartupMessages({
  library(cmdstanr)
  library(tidybayes)
  library(tidyverse)
})


# 2. acquired resistance ---------------------------------------
## extract data from studies ---------------------------------------
mallick <- tribble(
  ~author        , ~year , ~country       , ~abr , ~follow_up , ~baseline ,
  "Conradie"     ,  2020 , "South Africa" ,    1 , NA         , 57 - 2    ,
  "Ghodousi"     ,  2022 , "Moldova"      ,    8 ,         30 ,        30 ,
  "Liu"          ,  2024 , "Georgia"      ,    8 ,         94 ,       277 ,
  "Nimmo"        ,  2022 , "South Africa" ,    6 , NA         ,        92 ,
  "Nimmo"        ,  2022 , "South Africa" ,    8 , NA         ,       207 ,
  "Diacon"       ,  2022 , "Moldova"      ,    1 ,         10 ,        79 ,
  "Pym"          ,  2024 , "Georgia"      ,   12 ,         24 ,       233 ,
  "Guglielmetti" ,  2022 , "South Africa" ,    1 , NA         ,        22 ,
  "Guglielmetti" ,  2022 , "South Africa" ,    0 ,         10 ,        10 ,
  "Diacon"       ,  2022 , "South Africa" ,    0 , NA         ,        20 ,
  "Kempker"      ,  2022 , "South Africa" ,    1 , NA         ,        62
)

hu <- tribble(
  ~author       , ~year , ~country       , ~abr , ~follow_up , ~baseline ,
  "Brown"       ,  2023 , "South Africa" ,    2 , NA         , 147 - 12  ,
  "O'Donnell"   ,  2022 , "South Africa" ,    3 , NA         ,        58 ,
  "Chesov"      ,  2022 , "Moldova"      ,    4 ,         26 ,        62 ,
  "Mikiashvili" ,  2024 , "Georgia"      ,    5 ,         21 ,       106 ,
  "Ismail"      ,  2022 , "South Africa" ,   16 ,        695 ,       695
)

abr <- mallick %>% # number of participants who acquired resistance in study
  bind_rows(hu) %>%
  pull(abr)

baseline <- mallick %>% # number of participants in study
  bind_rows(hu) %>%
  pull(baseline)

N <- length(abr) # number of studies

stan_data <- list(
  N = N,
  x = abr,
  n = baseline
)

model <- cmdstan_model("stan/hierarchical_binomial.stan")

fit <- model$sample(
  data = stan_data,
  seed = 1,
  chains = 1,
  iter_warmup = 1000,
  iter_sampling = sims
)

# abr_samples is a pre-computed vector of acquired resistance probabilities
abr_samples <- gather_draws(fit, p_mean_weighted)$.value

plot(density(abr_samples), main = "Probability of acquired BDQR")

print(paste0("mean pa: ", round(mean(abr_samples), 3)))
cat("95% CrI:", round(quantile(abr_samples, c(0.025, 0.975)), 3))

# 3. spontaneous resistance ---------------------------------------
## read data from stochastic simulation ---------------------------------------
bac_growth <- read_csv(
  "c_simulation/out.csv",
  skip = 1,
  col_names = c("sim", "timestep", "S", "R")
)
target <- 1e10 # max bacteria population

bac_growth <- bac_growth %>%
  group_by(sim) %>%
  summarise(
    mutants = {
      # interpolate mutants when population = target

      x0 <- max(which(S < target)) # find last time step before target
      x1 <- min(which(S >= target)) # find first time step after target

      y0 <- S[x0] # get susceptible population sizes
      y1 <- S[x1]

      # get size of partial time step to target (for linear interpolation)
      xt <- (target - y0) / (y1 - y0)

      r0 <- R[x0] # get mutant population sizes
      r1 <- R[x1]

      # increase mutant population using partial time step (linear
      # interpolation)
      r_estimate <- r0 + xt * (r1 - r0)

      r_estimate
    },
    .groups = "drop"
  )

bac_growth %>%
  mutate(per = mutants / target * 100) %>%
  ggplot(aes(per)) +
  geom_density() +
  geom_vline(aes(xintercept = 1), lty = "dashed", color = "red") +
  scale_x_log10()

prop <- mean(bac_growth$mutants / target > 0.01)
ps_n <- nrow(bac_growth)
ps_x <- sum(bac_growth$mutants / target > 0.01)
a <- 0.05

# qbeta(a / 2, x, n - x + 1) * 100
ps_lb <- qbeta(a / 2, ps_x + 1, ps_n - ps_x + 1)
# qbeta(1 - a / 2, x + 1, n - x) * 100
ps_ub <- qbeta(1 - a / 2, ps_x + 1, ps_n - ps_x + 1)

plot(
  density(rbeta(sims, ps_x + 1, ps_n - ps_x + 1)),
  main = "Distribution of spontaneous BDQR"
)
print(paste0("mean ps: ", prop))
cat("95% CrI:", round(c(ps_lb, ps_ub), 5))
