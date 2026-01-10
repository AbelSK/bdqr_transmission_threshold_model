# #####################################
# Clean data sheets                   #
# Author: Abel Kjaersgaard            #
# Date: Wed Dec 17 2025               #
# #####################################

# 1. dependencies ---------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
})

# 2. Clean WHO data ---------------------------------------
## read WHO data ---------------------------------------
# read notification data
notif_data <- read_csv(file = "who_data/tb_notifications.csv")
# read estimate data
est_data <- read_csv(
  file = "who_data/mdr_rr_tb_burden_estimates.csv"
)
dr_data <- read_csv(file = "who_data/tb_dr_surveillance.csv")

## clean WHO data ---------------------------------------
est_data_clean <- est_data %>%
  select(country, year, contains("num"))
#fmt: skip
dr_data_clean <- dr_data %>%
  select(
    country, year, rr_new, rr_rel,
    rr_ret, rr_fqr_bdqr_lzdr,
    rr_fqr_bdqr_lzds, rr_fqr_bdqr_lzdu
  ) %>%
  mutate(across(starts_with("rr_"), ~ ifelse(is.na(.x), 0, .x))) %>%
  transmute(
    country = country,
    year = year,
    rr = rr_new + rr_rel + rr_ret,
    all_bdqr = rr_fqr_bdqr_lzdr + rr_fqr_bdqr_lzds + rr_fqr_bdqr_lzdu
  )

notif_data_clean <- notif_data %>%
  select(country, year, mdrxdr_bdq_tx, conf_rrmdr) %>%
  rename_with(., ~ c("country", "year", "bdq_tr", "mdrrr")) %>%
  left_join(dr_data_clean, by = join_by(country, year)) %>%
  # some data was missing from one of the sheets
  mutate(mdrrr_tb = ifelse(is.na(mdrrr), rr, mdrrr)) %>%
  select(country, year, bdq_tr, mdrrr_tb)

notif_data_clean <- notif_data %>%
  select(country, year, mdrxdr_bdq_tx, conf_rrmdr) %>%
  rename_with(., ~ c("country", "year", "bdq_tr", "mdrrr")) %>%
  left_join(dr_data_clean, by = join_by(country, year)) %>%
  # some data was missing from one of the sheets
  mutate(mdrrr_tb = ifelse(is.na(mdrrr), rr, mdrrr)) %>%
  select(country, year, bdq_tr, mdrrr_tb)

# 3. clean data country data ---------------------------------------
## Mozambique ---------------------------------------
barilar2024 <- readxl::read_xlsx("bdqr_prev_data/barilar2024.xlsx")

barilar2024_clean <- barilar2024 %>%
  janitor::clean_names() %>%
  mutate(year = as.numeric(year_of_sample_collection)) %>%
  mutate(
    prediction_bdq = ifelse(prediction_bdq == "preXDR", "pre", prediction_bdq)
  ) %>%
  mutate(
    bdq = ifelse(
      str_detect(prediction_bdq, "BDQ") | str_detect(prediction_bdq, "XDR"),
      1,
      0
    )
  ) %>%
  reframe(
    country = "Mozambique",
    prev_true = sum(bdq),
    prev_false = length(bdq) - prev_true,
    .by = year
  ) %>%
  drop_na()

write_csv(barilar2024_clean, "bdqr_prev_data/barilar2024_clean.csv")

## South Africa ---------------------------------------
roberts2024 <- readxl::read_xlsx("bdqr_prev_data/roberts2024.xlsx", sheet = 4)
roberts2024_clean <- roberts2024 %>%
  janitor::clean_names() %>%
  filter(bdq_susceptibility != "NA") %>%
  mutate(
    bdq_susceptibility = case_when(
      bdq_susceptibility == "R" ~ "prev_true",
      TRUE ~ "prev_false"
    )
  ) %>%
  count(year, bdq_susceptibility) %>%
  pivot_wider(names_from = bdq_susceptibility, values_from = n) %>%
  mutate(country = "South Africa", .before = year)

write_csv(roberts2024_clean, "bdqr_prev_data/roberts2024_clean.csv")

## Sierra Leone ---------------------------------------
blankson2024 <- read_csv("bdqr_prev_data/blankson2024.csv")
blankson2024_clean <- blankson2024 %>%
  mutate(gBDQ = ifelse(str_detect(gBDQ, "Rv0678"), 1, 0)) %>%
  reframe(
    prev_true = sum(gBDQ, na.rm = TRUE),
    prev_false = length(gBDQ) - prev_true,
    .by = `year of isolation`
  ) %>%
  drop_na() %>% # 3 strains with missing years (all BDQS)
  transmute(
    country = "Sierra Leone",
    year = `year of isolation`,
    prev_true = prev_true,
    prev_false = prev_false
  )

write_csv(blankson2024_clean, "bdqr_prev_data/blankson2024_clean.csv")


## Timm et al. ---------------------------------------
timm2023 <- readxl::read_xlsx("bdqr_prev_data/timm2023.xlsx", sheet = 1)

timm2023_clean <- timm2023 %>%
  janitor::clean_names() %>%
  filter(original_tb_type != "DS") |>
  select(trial, country, year_of_randomisation, bedaquiline_mgit) %>%
  filter(trial != "STAND") %>%
  filter(bedaquiline_mgit != "NA") %>%
  mutate(
    bedaquiline_mgit = ifelse(
      bedaquiline_mgit == "<=0.125",
      0.125,
      bedaquiline_mgit
    )
  ) %>%
  mutate(bedaquiline_mgit = as.numeric(bedaquiline_mgit)) %>%
  mutate(res = ifelse(bedaquiline_mgit <= 1, "prev_false", "prev_true")) %>%
  count(country, year_of_randomisation, res) |>
  pivot_wider(names_from = res, values_from = n, values_fill = 0) |>
  mutate(
    country = case_when(
      country == "Tanzania" ~ "United Republic of Tanzania",
      country == "Moldova" ~ "Republic of Moldova",
      country == "Russia" ~ "Russian Federation",
      TRUE ~ country
    )
  ) |>
  rename(year = year_of_randomisation)

timm2023_model_data <- notif_data_clean |>
  left_join(timm2023_clean) |>
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

write_csv(timm2023_model_data, "bdqr_prev_data/timm2023_model_data.csv")
