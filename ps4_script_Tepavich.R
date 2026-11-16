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

liquor_cv <- liquor |>
  stretch_tsibble(.init = 120, .step = 1) |>
  filter(.id < max(.id))

liquor_cv_fit <- liquor_cv |>
  model(snaive = SNAIVE(sales),
        ets    = ETS(log(sales)),
        arima  = ARIMA(log(sales)))

saveRDS(liquor_cv_fit, "results/liquor_cv_fit.rds")

#216 forecast origins

liquor_cv_fc <- liquor_cv_fit |>
  forecast(h = 1)

liquor_cv_fc |>
  accuracy(liquor, measures = list(RMSE = RMSE, MAE = MAE)) |>
  knitr::kable(digits = 2)

#The first forecast: January 2002 
#The last forecast: December 2019
#The ARIMA model forecasts most accurately; both metrics agree 

liquor_errs <- liquor_cv_fc |>
  as_tibble() |>
  left_join(as_tibble(liquor) |> select(Month, actual = sales),
            by = "Month") |>
  mutate(e = actual - .mean) |>
  select(Month, .model, e) |>
  pivot_wider(names_from = .model, values_from = e)

forecast::dm.test(liquor_errs$arima, liquor_errs$ets, h = 1, power = 2)

forecast::dm.test(liquor_errs$arima, liquor_errs$snaive, h = 1, power = 2)

#DM stat (ARIMA & ETS): -2.9062
#p-value (ARIMA & ETS): 0.004042
#DM stat (ARIMA & snaive): -10.755
#p-value (ARIMA & snaive): < 2.2e-16
#Both tests reject H0 meaning neither the ARIMA and ETS nor the 
#ARIMA and snaive models have equal forecast accuracy
#If I were the liquor store chain, I would use the ARIMA model