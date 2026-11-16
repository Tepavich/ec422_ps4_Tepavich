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

liquor_fit <- liquor |>
  model(arima = ARIMA(log(sales)),
        ets   = ETS(log(sales)))

liquor_fit |> select(arima) |> report()
liquor_fit |> select(ets) |> report()

#The ARIMA model selected was ARIMA(1,0,2)(2,1,1)[12] with drift. 
#This model includes one nonseasonal autoregressive term, two nonseasonal 
#moving average terms, two seasonal autoregressive terms, and one seasonal 
#moving average term. The model does not use regular differencing (d = 0), 
#but it does use one seasonal difference (D = 1) to account for the annual 
#seasonal pattern in monthly sales data. The ETS model selected was ETS(A,A,A),
#indicating additive error, additive trend, and additive seasonality.

liquor_fit |> select(arima) |> gg_tsresiduals()

liquor_fit |>
  select(arima) |>
  augment() |>
  features(.innov, ljung_box, lag = 24, dof = 6)

#The ljung-box test resulted in a p-value of approximately zero, meaning
#we reject the H0 that there is no autocorrelation which suggests the residuals 
#are not white noise.