# #####################################
# Data setup                          #
# Author: Abel Kjaersgaard            #
# Date: Wed Dec 17 2025               #
# #####################################

# 1. dependencies ---------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
})

# 2. country data ---------------------------------------
## Mozambique ---------------------------------------
mozambique_clean <- read_csv("bdqr_prev_data/barilar2024_clean.csv")

## Brazil ---------------------------------------
brazil_clean <- tibble(
  country = "Brazil", # Marcon et al 2025
  year = 2022,
  prev_false = 38 - 3,
  prev_true = 3
)

## South Africa ---------------------------------------
roberts2024_clean <- read_csv("bdqr_prev_data/roberts2024_clean.csv") |> tail(1)

moultrie2024_clean <- tibble(
  country = "South Africa",
  year = 2023,
  prev_false = 2308 - 149,
  prev_true = 149
)

# Ismail 2022 reports BDQR prevalence by guideline period, rather than anually.
# disaggregate numbers from each policy period, to produce annual numbers based
# on the number of months of sampling in each period and year

ismail2022 <- tribble(
  ~period     , ~start_date  , ~end_date    , ~susceptible , ~resistant , ~total ,
  "<2017"     , "2015-01-01" , "2016-01-01" ,         1074 ,         36 ,   1110 ,
  "<2017"     , "2016-01-01" , "2017-01-01" ,         1074 ,         36 ,   1110 ,
  "<2017"     , "2017-01-01" , "2017-06-01" ,         1074 ,         36 ,   1110 ,

  "2017-2018" , "2017-06-01" , "2018-01-01" ,          523 ,         25 ,    548 ,
  "2017-2018" , "2018-01-01" , "2018-06-01" ,          523 ,         25 ,    548 ,

  ">2018"     , "2018-06-01" , "2019-01-01" ,          350 ,         15 ,    365 ,
  ">2018"     , "2019-01-01" , "2019-07-31" ,          350 ,         15 ,    365
)

ismail2022_clean <- ismail2022 %>%
  mutate(end_date = as.Date(end_date), start_date = as.Date(start_date)) %>%
  mutate(months = interval(start_date, end_date) %/% months(1)) %>%
  mutate(year = c(2015, 2016, 2017, 2017, 2018, 2018, 2019)) %>%
  group_by(period) %>%
  mutate(total_months = sum(months)) %>%
  ungroup() %>%
  reframe(
    prev_false = sum(susceptible * months / total_months),
    prev_true = sum(resistant * months / total_months),
    .by = year
  ) %>%
  mutate(country = "South Africa", .before = 1)

south_africa_clean <- ismail2022_clean %>%
  bind_rows(roberts2024_clean) %>%
  bind_rows(moultrie2024_clean)

## Kazakhstan ---------------------------------------
# we disaggregate the total BDQR prevalence by the number of months of sampling
# in each year of data collection
kaz_int <- interval(as.Date("2020-09-01"), as.Date("2022-10-01")) %/%
  months(1)

kazakhstan_clean <- tribble(
  ~country     , ~year , ~prev_false              , ~prev_true       ,
  "Kazakhstan" ,  2020 , (137 - 3) * 3 / kaz_int  , 3 * 3 / kaz_int  ,
  "Kazakhstan" ,  2021 , (137 - 3) * 12 / kaz_int , 3 * 12 / kaz_int ,
  "Kazakhstan" ,  2022 , (137 - 3) * 10 / kaz_int , 3 * 10 / kaz_int
)

## China ---------------------------------------
chn_int <- interval(as.Date("2019-03-01"), as.Date("2020-07-01")) %/%
  months(1)

china_clean <- tribble(
  ~country , ~year , ~prev_false              , ~prev_true       ,
  "China"  ,  2019 , (205 - 9) * 10 / chn_int , 9 * 10 / chn_int , # Hu 2023
  "China"  ,  2020 , (205 - 9) * 6 / chn_int  , 9 * 6 / chn_int  ,
  "China"  ,  2021 , 245 - 5                  ,                5 , # Tong 2023
  "China"  ,  2022 , 263 - 4                  ,                4 # Li 2024
)

## Uzbekistan ---------------------------------------
uzbekistan_clean <- tibble(
  country = "Uzbekistan",
  year = 2019:2023,
  prev_true = c(6, 18, 14, 9, 14),
  prev_false = c(561 - 6 - 429, 683 - 18 - 1, 442 - 14, 446 - 9, 273 - 14)
)

## Sierra Leone ---------------------------------------
sierra_leone_clean <- read_csv("bdqr_prev_data/blankson2024_clean.csv")

## Ethiopia ---------------------------------------
eth_int <- interval(as.Date("2022-02-01"), as.Date("2024-08-01")) /
  months(1)

ethiopia_clean <- tibble(
  country = "Ethiopia",
  year = 2022:2024,
  prev_true = c(7 * 11, 7 * 12, 7 * 7) / eth_int,
  prev_false = c(468 * 11, 468 * 12, 468 * 7) / eth_int - prev_true
)

## India  ---------------------------------------
ind_int <- interval(as.Date("2019-04-01"), as.Date("2021-01-01")) /
  months(1)

india_clean <- tribble(
  ~country , ~year , ~prev_false              , ~prev_true       ,
  "India"  ,  2019 , (165 - 2) * 9 / ind_int  , 2 * 9 / ind_int  , # Padmapriyadarsini (2022)
  "India"  ,  2020 , (165 - 2) * 12 / ind_int , 2 * 12 / ind_int ,
  "India"  ,  2021 , (403 - 7) / 3            , 7 / 3            , # Padmapriyadarsini (2024)
  "India"  ,  2022 , (403 - 7) / 3            , 7 / 3            ,
  "India"  ,  2023 , (403 - 7) / 3            , 7 / 3            ,
)

## Georgia ---------------------------------------
geo_int <- interval(as.Date("2017-11-01"), as.Date("2021-01-01")) %/% months(1)

georgia_clean <- tribble(
  ~country  , ~year , ~prev_false       , ~prev_true       ,
  "Georgia" ,  2017 , 83 * 2 / geo_int  , 6 * 2 / geo_int  ,
  "Georgia" ,  2018 , 83 * 12 / geo_int , 6 * 12 / geo_int ,
  "Georgia" ,  2019 , 83 * 12 / geo_int , 6 * 12 / geo_int ,
  "Georgia" ,  2020 , 83 * 12 / geo_int , 6 * 12 / geo_int ,
)

## Japan ---------------------------------------
jpn_int <- interval(as.Date("2018-01-01"), as.Date("2022-10-01")) %/% months(1)

japan_clean <- tribble(
  ~country , ~year , ~prev_false       , ~prev_true       ,
  "Japan"  ,  2018 , 22 * 12 / jpn_int , 0 * 12 / jpn_int ,
  "Japan"  ,  2019 , 22 * 12 / jpn_int , 0 * 12 / jpn_int ,
  "Japan"  ,  2020 , 22 * 12 / jpn_int , 0 * 12 / jpn_int ,
  "Japan"  ,  2021 , 22 * 12 / jpn_int , 0 * 12 / jpn_int ,
  "Japan"  ,  2022 , 22 * 10 / jpn_int , 0 * 10 / jpn_int ,
)

## Peru ---------------------------------------
peru_int <- interval(as.Date("2017-05-01"), as.Date("2019-03-01")) %/% months(1)
peru_clean <- tribble(
  ~country , ~year , ~prev_false         , ~prev_true        ,
  "Peru"   ,  2017 , 165 * 8 / peru_int  , 6 * 8 / peru_int  ,
  "Peru"   ,  2018 , 165 * 12 / peru_int , 6 * 12 / peru_int ,
  "Peru"   ,  2019 , 165 * 2 / peru_int  , 6 * 2 / peru_int  ,
)

## Timm et al. ---------------------------------------
timm2023_model_data <- read_csv("bdqr_prev_data/timm2023_model_data.csv")
