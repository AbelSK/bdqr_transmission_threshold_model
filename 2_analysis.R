# #####################################
# Analysis                            #
# Author: Abel Kjaersgaard            #
# Date: Wed Dec 17 2025               #
# #####################################

# 1. dependencies ---------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
})

source('0_global_setup.R')
source('1_clean_data_sheets.R')
source('1_data_setup.R')
source('1_parameter_setup.R')

# 2. run transmission threshold model ---------------------------------------
model_data <- notif_data_clean %>%
  left_join(
    rbind(
      south_africa_clean,
      mozambique_clean,
      china_clean,
      uzbekistan_clean,
      brazil_clean,
      kazakhstan_clean,
      india_clean,
      sierra_leone_clean,
      ethiopia_clean,
      japan_clean,
      georgia_clean,
      peru_clean
    )
  )

model_data <- model_data |>
  mutate(
    bdq_tr = case_when(
      # impute based on policy
      country == "South Africa" & year == 2022 ~ mdrrr_tb,

      country %in% c("Mozambique", "Sierra Leone") & year < 2017 ~ 0,

      # impute using lagging and leading year
      is.na(bdq_tr) &
        !is.na(lag(bdq_tr)) &
        !is.na(lead(bdq_tr)) ~
        (lag(bdq_tr) / lag(mdrrr_tb) + lead(bdq_tr) / lead(mdrrr_tb)) /
        2 *
        mdrrr_tb,

      # impute using lagging year only
      is.na(bdq_tr) & !is.na(lag(bdq_tr)) ~
        lag(bdq_tr) / lag(mdrrr_tb) * mdrrr_tb,

      # impute using leading year only
      is.na(bdq_tr) & !is.na(lead(bdq_tr)) ~
        lead(bdq_tr) / lead(mdrrr_tb) * mdrrr_tb,

      TRUE ~ bdq_tr
    )
  ) |>
  mutate(lag_tr = lag(bdq_tr, 1)) %>%
  mutate(lag_tr_rate = lag(bdq_tr / mdrrr_tb, 1)) %>%
  select(-bdq_tr) %>%
  drop_na()

model_data <- model_data |>
  bind_rows(timm2023_model_data) |> # add Timm et al multi-country data
  group_by(country, year) |> # sum observations by country and year
  summarise(
    mdrrr_tb = mean(mdrrr_tb),
    lag_tr = mean(lag_tr),
    lag_tr_rate = mean(lag_tr_rate),
    prev_true = sum(prev_true),
    prev_false = sum(prev_false)
  ) |>
  ungroup()

out <- model_data |>
  mutate(result = calc_threshold()) |>
  unnest_wider(result)

# 3. manuscript figures ---------------------------------------
## 1: model structure ---------------------------------------
comment('Figure 1 for the main text is generated outside of R')
## 2: BDQR prevalence ---------------------------------------
range(model_data$prev_false + model_data$prev_true)
range(model_data$prev_true / (model_data$prev_true + model_data$prev_false))

figure2 <- model_data %>%
  mutate(
    shape = case_when(
      country %in% c("Mozambique", "Peru", "Sierra Leone") ~ "g",
      country == "South Africa" & year == 2020 ~ "g",
      country == "Brazil" & year == 2022 ~ "g",
      TRUE ~ "p"
    )
  ) |>
  mutate(prev = prev_true / (prev_true + prev_false)) %>%
  ggplot(aes(year, prev)) +
  geom_line() +
  geom_point(aes(size = prev_false + prev_true, shape = shape), alpha = .7) +
  scale_x_continuous(
    limits = c(2015, 2025),
    breaks = seq(2015, 2025, 2),
    labels = \(x) {
      sprintf("'%02d", x %% 100)
    }
  ) +
  scale_shape_manual(values = c(15, 16), ) +
  labs(x = "Year", y = "BDQR prevalence") +
  theme_bw() +
  scale_size_continuous(
    name = "Sample size",
    breaks = c(10, 100, 500, 1000, 2000)
  ) +
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  facet_wrap(~country, ncol = 3) +
  guides(
    shape = "none"
  ) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

ggsave("figures/figure2.png", figure2, width = 6, height = 6, create.dir = TRUE)

## 3: Scenario analysis ---------------------------------------
N = 2000
M = 1000
scenario_data <- expand_grid(
  prev_true = seq(0.01, 0.05, 0.001) * M,
  prev_false = M,
  ps = 0.02 / 100,
  mdrrr_tb = N,
  lag_tr = seq(0.5, 0.9, 0.01) * N
) |>
  mutate(prev_false = M - prev_true) |>
  mutate(result = calc_threshold()) |>
  unnest_wider(result) |>
  mutate(lag_tr = lag_tr / N, prev = prev_true / (prev_true + prev_false))

eta = 0.05
# extract poitns need threhsold = 1
scenario_line <- scenario_data %>%
  filter(thres_val > 1 - eta & thres_val < 1 + eta)

figure3 <- scenario_data %>%
  mutate(color = thres_val > 1) %>%
  ggplot(aes(lag_tr, prev)) +
  geom_raster(aes(fill = thres_val)) +
  geom_smooth(
    data = scenario_line,
    method = "lm",
    se = FALSE,
    color = "black",
    aes(linetype = "dashed"),
    size = 0.5
  ) +
  labs(
    x = "MDR/RR-TB patients treated with BDQ last year",
    y = "BDQR prevalence among MDR/RR-TB"
  ) +
  scale_fill_gradient2(
    name = "Threshold value",
    low = "blue1",
    mid = "white",
    high = "red1",
    midpoint = 1
  ) +
  guides(fill = guide_colorbar(title.position = "top")) +
  guides(linetype = guide_legend(title.position = "top", title.hjust = 0.5)) +
  scale_linetype_manual(
    name = NULL,
    values = c("dashed" = "dashed"),
    labels = c("dashed" = "Threshold\nvalue = 1")
  ) +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw() +
  theme(
    plot.caption = element_text(hjust = 0, size = 7),
    legend.position = "top"
  )

ggsave("figures/figure3.png", figure3, width = 5, height = 5)

scenario_data |>
  filter(thres_val >= 1) |>
  filter(lag_tr == 0.5) |>
  pull(prev) |>
  min()

## 4: country analysis ---------------------------------------
figure4 <- out %>%
  ggplot(aes(x = year)) +
  geom_hline(
    aes(yintercept = 1, color = "values >1\nindicate transmission"),
    lty = "dashed"
  ) +
  geom_point(aes(y = thres_val, size = (prev_true + prev_false))) +
  geom_errorbar(aes(ymin = thres_val_lo, ymax = thres_val_hi)) +
  scale_x_continuous(breaks = seq(2015, 2023, 2), labels = \(x) {
    sprintf("'%02d", x %% 100)
  }) +
  scale_y_log10(label = function(x) format(x, scientific = FALSE)) +
  scale_size_continuous(
    range = c(0.2, 3.5),
    breaks = c(10, 100, 500, 1000, 2000),
    name = "BDQR prevalence\nsample size"
  ) +
  labs(x = "Year", y = "Transmissio threshold values") +
  guides(size = guide_legend(nrow = 2)) +
  scale_color_manual(name = "", values = "red") +
  facet_wrap(~country, scales = "free_y", ncol = 3) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.key.spacing.y = unit(0.0, "cm"),
    panel.grid.minor = element_blank()
  )

ggsave("figures/figure4.png", figure4, width = 6, height = 6)

# 4. appendix figures ---------------------------------------
## A1: imputation ---------------------------------------

figureA1 <- notif_data_clean |>
  filter(country %in% model_data$country) |>
  mutate(who_lag_tr_rate = lag(bdq_tr) / lag(mdrrr_tb)) |>
  filter(year >= 2015 & year <= 2024) |>
  left_join(
    model_data |>
      select(country, year, lag_tr_rate, prev_true)
  ) |>
  mutate(
    bdqr_prev_flag = ifelse(
      is.na(prev_true),
      "Missing BDQR prevalence",
      "Point has BDQR prevalence"
    ),
    bdq_tr_flag = ifelse(
      is.na(who_lag_tr_rate),
      "BDQ treatment is imputed",
      "WHO data"
    )
  ) |>
  mutate(
    # if WHO lag_tr_rate is missing, use imputed tr_rate
    lag_tr_rate = ifelse(is.na(who_lag_tr_rate), lag_tr_rate, who_lag_tr_rate)
  ) |>
  ggplot(aes(year, lag_tr_rate)) +
  geom_point(aes(shape = bdqr_prev_flag, color = bdq_tr_flag)) +
  geom_line() +
  facet_wrap(~country, scales = "free_y", ncol = 3) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(labels = \(x) {
    sprintf("'%02d", x %% 100)
  }) +
  theme_bw() +
  scale_color_manual(name = "", values = c("red", "black")) +
  scale_shape_manual(name = "", values = c(1, 19)) +
  labs(x = "Year", y = "BDQ treatment rate") +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    legend.spacing.y = unit(0.0, "cm"),
    legend.box = "vertical"
  )


ggsave("figures/figureA1.png", figureA1, height = 6, width = 6)

## A2: treatment vs prevalence ---------------------------------------
figureA2 <- out %>%
  mutate(prop = prev_true / (prev_true + prev_false) * 100) %>%
  ggplot(aes(x = year)) +

  geom_hline(aes(yintercept = 1), color = "red", lty = "dashed") +

  geom_point(aes(y = thres_val, color = "threshold value")) +
  geom_errorbar(aes(
    ymin = thres_val_lo,
    ymax = thres_val_hi,
    color = "threshold value"
  )) +

  geom_line(aes(y = prop, color = "BDQR prevalence")) +
  geom_line(aes(y = lag_tr / mdrrr_tb * 100, color = "treatment rate")) +
  geom_point(aes(y = prop, color = "BDQR prevalence")) +
  geom_point(aes(y = lag_tr / mdrrr_tb * 100, color = "treatment rate")) +

  scale_x_continuous(breaks = seq(2015, 2023, 2), labels = function(x) {
    sprintf("'%02d", x %% 100)
  }) +
  scale_y_log10(
    label = function(x) format(x, scientific = FALSE),
    sec.axis = sec_axis(
      ~.,
      name = "BDQR prevalence/treatment rate",
      label = function(x) format(paste0(round(x), "%"), scientific = FALSE)
    )
  ) +
  scale_color_manual(values = c("red", "black", "blue"), name = "") +
  labs(x = "Year", y = "Threshold of transmission value") +
  facet_wrap(~country, ncol = 3) +
  theme_bw() +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave("figures/figureA2.png", figureA2, width = 6, height = 6)

## A3-8: prev, pa, ps sensitivity analysis ---------------------------------------
sa_scenarios <- tibble(
  sa = c("pa x 2", "pa x 0.5", "ps x 2", "ps x 0.5", "prev x 2", "prev x 0.5"),
  ps_factor = c(1, 1, 2, 0.5, 1, 1),
  pa_factor = c(2, 0.5, 1, 1, 1, 1),
  prev_factor = c(1, 1, 1, 1, 2, 0.5)
)

for (i in 1:nrow(sa_scenarios)) {
  sa_plot_data <- model_data |>
    mutate(sa = sa_scenarios[[i, 1]]) |>
    mutate(
      result = calc_threshold(
        ps_factor = sa_scenarios[[i, 2]],
        pa_factor = sa_scenarios[[i, 3]],
        prev_factor = sa_scenarios[[i, 4]]
      )
    ) |>
    unnest_wider(result) |>
    bind_rows(out |> mutate(sa = "baseline")) |>
    left_join(
      out |>
        select(country, year, base_lo = thres_val_lo),
      by = c("country", "year")
    ) |>
    mutate(
      sa_change = (base_lo > 1 & thres_val_lo < 1) |
        (base_lo < 1 & thres_val_lo > 1)
    )

  figure_sa <- sa_plot_data |>
    mutate(year = ifelse(sa == "baseline", year - 0.2, year + 0.2)) |>
    filter(sa %in% c("baseline", sa_scenarios[[i, 1]])) |>
    sa_plot()

  ggsave(
    glue("figures/figureA{i+2}_", gsub(" ", "", sa_scenarios[[i, 1]]), ".png"),
    width = 6,
    height = 6
  )
}

## A9-10: WHO estimated vs notified  ---------------------------------------
model_data %>%
  left_join(est_data_clean) %>%
  mutate(difference = mdrrr_tb - e_inc_rr_num) %>%
  select(country, year, difference) %>%
  arrange(desc(difference))

figuresA9 <- model_data |>
  left_join(est_data_clean) |>
  mutate(sa = "Estimated MDR/RR-TB") |>
  mutate(result = calc_threshold(N_col = "e_inc_rr_num")) |>
  unnest_wider(result) |>
  bind_rows(out |> mutate(sa = "baseline")) |>
  left_join(
    out |>
      select(country, year, base_lo = thres_val_lo),
    by = c("country", "year")
  ) |>
  mutate(
    sa_change = (base_lo > 1 & thres_val_lo < 1) |
      (base_lo < 1 & thres_val_lo > 1)
  ) |>
  mutate(year = ifelse(sa == "baseline", year - 0.2, year + 0.2)) |>
  sa_plot()

ggsave("figures/figuresA9_est_mdrrr.png", figuresA9, width = 6, height = 6)

figuresA10 <- model_data |>
  left_join(est_data_clean) |>
  mutate(sa = "adjusted BDQ treatment rate") |>
  mutate(lag_tr = lag_tr * e_inc_rr_num / mdrrr_tb) |>
  mutate(result = calc_threshold()) |>
  unnest_wider(result) |>
  bind_rows(out |> mutate(sa = "baseline")) |>
  left_join(
    out |>
      select(country, year, base_lo = thres_val_lo),
    by = c("country", "year")
  ) |>
  mutate(
    sa_change = (base_lo > 1 & thres_val_lo < 1) |
      (base_lo < 1 & thres_val_lo > 1)
  ) |>
  mutate(year = ifelse(sa == "baseline", year - 0.2, year + 0.2)) |>
  mutate(
    sa = factor(sa, levels = c("baseline", "adjusted BDQ treatment rate"))
  ) |>
  sa_plot() +
  scale_shape_manual(
    name = "",
    values = c("baseline" = 16, "adjusted BDQ treatment rate" = 17)
  )

ggsave("figures/figuresA10_est_bdqtr.png", figuresA10, width = 6, height = 6)

## Table A3  ---------------------------------------
table_a3 <- out |>
  select(country, year, lag_tr, mdrrr_tb, prev_true, prev_false) |>
  left_join(
    notif_data_clean |>
      mutate(who_lag_tr_rate = lag(bdq_tr) / lag(mdrrr_tb)) |>
      mutate(
        bdq_tr_flag = ifelse(
          is.na(who_lag_tr_rate),
          "BDQ treatment is imputed",
          "WHO data"
        )
      ) |>
      select(country, year, bdq_tr_flag)
  ) |>
  mutate(lag_tr = round(lag_tr, 1)) |>
  mutate(
    lag_tr = ifelse(bdq_tr_flag == "WHO data", lag_tr, paste0("[", lag_tr, "]"))
  ) |>
  select(-bdq_tr_flag) |>
  mutate(
    country = ifelse(
      row_number() == 1,
      "Brazil",
      ifelse(country == lag(country), NA, country)
    )
  ) |>
  mutate(
    prev_true_fmt = paste0(
      round(prev_true, 1),
      " (",
      round(prev_true / (prev_true + prev_false) * 100, 1),
      "%)"
    ),
    prev_false_fmt = paste0(
      round(prev_false, 1),
      " (",
      round(prev_false / (prev_true + prev_false) * 100, 1),
      "%)"
    )
  ) |>
  select(-prev_true, -prev_false) |>
  rename_with(
    ~ c(
      "Country",
      "Year",
      "BDQ treated preceding year (imputed)",
      "MDR/RR-TB incidence",
      "BDQR (%)",
      "BDQS (%)"
    )
  )

write_csv(x = table_a3, file = "tableA3.csv")
