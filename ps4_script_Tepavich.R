library(tidyverse)
library(fpp3)
library(BVAR)

liquor_raw <- read_csv(
  "data/MRTSSM4453USN.csv",
  na = ".",
  col_types = cols(
    observation_date = col_date(format = "%Y-%m-%d"),
    MRTSSM4453USN = col_double()
  )
)

liquor <- liquor_raw |>
  transmute(Month = yearmonth(observation_date), sales = MRTSSM4453USN) |>
  as_tsibble(index = Month) |>
  filter(Month <= yearmonth("2019 Dec"))

has_gaps(liquor)

liquor |>
  autoplot(sales) +
  labs(
    title = "Monthly Beer, Wine, and Liquor Store Sales",
    x = "Year",
    y = "Sales (Millions of Dollars)"
  )

liquor |>
  gg_season(log(sales)) +
  labs(
    title = "Seasonal Pattern of Log Liquor Sales",
    x = "Month",
    y = "Log(Sales)"
  )

#The data is upward trending with a monthly seasonal pattern that spikes every 
#January and then drops the next month.The seasonal swings grow with the level 
#of the series meaning taking logs would be a good idea.
