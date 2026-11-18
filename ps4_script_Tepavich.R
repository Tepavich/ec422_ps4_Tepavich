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

##PART 3

liquor_combo_fc <- liquor_cv_fit |>
  mutate(combo = (arima + ets) / 2) |>
  forecast(h = 1)

liquor_combo_fc |>
  accuracy(liquor, measures = list(RMSE = RMSE, MAE = MAE)) |>
  knitr::kable(digits = 2)

#Combo has an RMSE and an MAE between the two models that go into it

cor(liquor_errs$arima, liquor_errs$ets)

#Combination didn't beat the ARIMA model because the ARIMA and ETS 
#model do not have equal forecast accuracy, and the correlation between
#the ARIMA error terms and the ETS error terms is quite high

##PART 4

qd_raw <- read_csv("data/fredqd_2026-09.csv", show_col_types = FALSE) |>
  slice(-(1:2))

growth <- function(x) 100 * difference(log(x))

qd <- qd_raw |>
  mutate(date = as.Date(sasdate, "%m/%d/%Y")) |>
  transmute(date,
            gdp    = growth(GDPC1),     # real GDP
            cons   = growth(PCECC96),   # real consumption
            inv    = growth(GPDIC1),    # real investment
            emp    = growth(PAYEMS),    # payroll employment
            infl   = growth(GDPCTPI),   # GDP price index
            unrate = difference(UNRATE),     # unemployment rate
            ffr    = difference(FEDFUNDS),   # federal funds rate
            gs10   = difference(GS10)) |>    # 10-year Treasury yield
  filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))

y <- qd |> select(-date) |> as.matrix()

nrow(y)
head(y)

N <- nrow(y)
origins <- 120:(N - 1)
actual <- y[121:N, ]

#We will produce 120 forecasts for each variable 
#First forecast: 1990 Q1
#Last forecast: 2019 Q4

var_ols_fc <- function(z, p) {
  e <- embed(z, p + 1)            # columns: y_t, then its lags 1 to p
  M <- ncol(z)                    # number of variables
  B <- coef(lm(e[, 1:M] ~ e[, -(1:M)]))
  x <- c(1, t(z[nrow(z):(nrow(z) - p + 1), ]))   # constant, latest p lags
  x %*% B
}

fc_var1 <- NULL
fc_var4 <- NULL
for (T in origins) {
  fc_var1 <- rbind(fc_var1, var_ols_fc(y[1:T, ], 1))
  fc_var4 <- rbind(fc_var4, var_ols_fc(y[1:T, ], 4))
}

rmse_by_var <- function(f) sqrt(colMeans((actual - f)^2))

rmse_by_var(fc_var1)
rmse_by_var(fc_var4)

#var(1) forecasts better for all the variables but inflation.
#The smaller model forecasts better because it has fewer parameters to 
#estimate and is less likely to overfit the data.

bvar_cv <- function(lambda) {
  pr <- bv_priors(hyper = "lambda",
                  mn = bv_mn(lambda = bv_lambda(mode = lambda, sd = 0.4,
                                                min = 0.999 * lambda,
                                                max = 1.001 * lambda),
                             b = 0))
  fc <- NULL
  for (T in origins) {
    fit_T <- bvar(y[1:T, ], lags = 4, priors = pr,
                  n_draw = 2000, n_burn = 1000, verbose = FALSE)
    x <- c(1, t(y[T:(T - 3), ]))         # constant and 4 lags of each variable
    fc <- rbind(fc, x %*% coef(fit_T))   # forecast of y_{T+1}
  }
  fc
}

set.seed(422)
fc_bvar_03 <- bvar_cv(0.3)
saveRDS(fc_bvar_03, "results/fc_bvar_03.rds")

set.seed(422)
fc_bvar_005 <- bvar_cv(0.05)

var_results <- tibble(
  Variable = colnames(y),
  `VAR(1)` = rmse_by_var(fc_var1),
  `VAR(4)` = rmse_by_var(fc_var4),
  `BVAR(4), lambda = 0.3` = rmse_by_var(fc_bvar_03),
  `BVAR(4), lambda = 0.05` = rmse_by_var(fc_bvar_005)
)

saveRDS(var_results, "results/var_results.rds")

var_results |> knitr::kable(digits = 3)


#The BVAR(4) with lambda = 0.3 generally performs better than the OLS VAR(4). 
#It has a lower RMSE for seven of the eight variables except GDP.
#Although both models use the same variables and four lags, they differ in how 
#they estimate their coefficients. The OLS VAR estimates coefficients without 
#shrinkage, while the BVAR uses a Minnesota prior to shrink coefficients. 
#This shrinkage helps reduce overfitting, particularly when estimating a 
#model with many parameters.

#The BVAR(4) with lambda = 0.3 performs better than VAR(1) for investment, 
#inflation, and the 10-year Treasury yield. However, VAR(1) performs better 
#for GDP, consumption, employment, unemployment, and the federal funds rate.
#This suggests that including additional lags and applying shrinkage can 
#improve forecasting for certain variables, but a simpler VAR(1) model can 
#still perform better for others.

#When lambda decreases from 0.3 to 0.05, the model applies less shrinkage, 
#allowing the estimated coefficients to move farther from their prior values. 
#The BVAR(4) with lambda = 0.05 performs worse than the lambda = 0.3 
#model for seven of the eight variables. The federal funds rate is the 
#exception, where lambda = 0.05 produces a lower RMSE.
#This suggests that stronger shrinkage generally improves 
#forecasting accuracy in this dataset, likely by reducing overfitting.
#Choosing the model with the lowest RMSE from this table and reporting that 
#same RMSE as an estimate of its future forecasting accuracy would likely not 
#be accurate. This is because the model was selected based on its 
#performance on the same cross-validation sample used to evaluate the models. 
#The winning model may have benefited from random variation in that particular 
#sample. A separate test sample or additional out-of-sample evaluation would 
#provide a more reliable estimate of future forecasting performance.

e_var4 <- actual[, "gdp"] - fc_var4[, 1]
e_bvar <- actual[, "gdp"] - fc_bvar_03[, 1]

forecast::dm.test(e_var4, e_bvar, h = 1, power = 2)

e_var4 <- actual[, "infl"] - fc_var4[, 5]
e_bvar <- actual[, "infl"] - fc_bvar_03[, 5]

forecast::dm.test(e_var4, e_bvar, h = 1, power = 2)

#For GDP growth, the p value is 0.03468, meaning we reject the null hypothesis
#that the two models have equal forecast accuracy; however, for inflation,
#the p value is 0.1489 meaning we fail to reject the null hypothesis that the 
#two models have equal forecast accuracy. This is inconsistent with my 
#findings in part 4.3b where the bvar forecasted more accurately for 
#inflation.

##PART 5

md_raw <- read_csv("data/fredmd_2026-09.csv", show_col_types = FALSE)
tcode <- unlist(md_raw[1, -1])
md <- md_raw[-1, ] |> mutate(date = as.Date(sasdate, "%m/%d/%Y"))

transform_series <- function(x, code) {
  d <- function(v) c(NA, diff(v))
  switch(as.character(code),
         "1" = x,                      # level
         "2" = d(x),                   # first difference
         "3" = d(d(x)),                # second difference
         "4" = log(x),                 # log
         "5" = 100 * d(log(x)),        # growth rate
         "6" = 100 * d(d(log(x))),     # change in the growth rate
         "7" = d(x / dplyr::lag(x) - 1))
}

fm <- map2_dfc(md[names(tcode)], tcode, transform_series) |>
  mutate(date = md$date, .before = 1) |>
  filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))

fm <- fm[, colSums(is.na(fm)) == 0]

dim(fm)

tcode[c("INDPRO", "PAYEMS", "CPIAUCSL", "FEDFUNDS")]

#INDPRO: 5 (growth rate, first difference) 
#PAYEMS: 5 (growth rate, first difference)
#CPIAUCSL: 5 (growth rate, second difference)
#FEDFUNDS: 2 (first difference)
#the FEDFUNDS series is the only one that is not transformed to a growth 
#rate because it is already a rate. The other three series are transformed 
#to growth rates to make them stationary and to focus on the percentage changes 
#rather than the absolute levels.

Z <- fm |> select(-date) |> as.matrix()

pc <- prcomp(scale(Z))

summary(pc)$importance[, 1:5]

#the first three factors explain 29.45% of the total variation in all series

F1 <- pc$x[, 1] * sign(cor(pc$x[, 1], fm$INDPRO))

tibble(date = fm$date, F1 = F1) |>
  ggplot(aes(date, F1)) +
  geom_line() +
  labs(title = "First Factor from FRED-MD", x = NULL, y = NULL) +
  theme_minimal()

sort(abs(cor(pc$x[, 1], Z))[1, ], decreasing = TRUE)[1:8]
sort(abs(cor(pc$x[, 2], Z))[1, ], decreasing = TRUE)[1:8]

#First Factor
#1. PAYEMS — 0.842
#2. IPMANSICS — 0.840
#3. USGOOD — 0.839
#4. INDPRO — 0.822
#5. MANEMP — 0.803
#6. IPFPNSS — 0.790
#7. CUMFNS — 0.778
#8. DMANEMP — 0.777
#Second factor
#1. CUSR0000SAC — 0.747
#2. DNDGRG3M086SBEA — 0.739
#3. CUSR0000SA0L2 — 0.728
#4. CPIAUCSL — 0.720
#5. CUSR0000SA0L5 — 0.703
#6. CPITRNSL — 0.681
#7. PCEPI — 0.680
#8. CPIULFSL — 0.657
#The first factor appears to represent real economic activity, particularly 
#employment and industrial production, while the second factor appears to 
#represent inflation and changes in consumer prices.

ip <- fm$INDPRO
n <- nrow(Z)
init <- which(fm$date == as.Date("1989-12-01"))

e_ar <- NULL
e_fa <- NULL
for (T in init:(n - 1)) {
  F_T <- prcomp(scale(Z[1:T, ]))$x[, 1:3]         # factors using data through T
  d <- data.frame(y_next = ip[2:T], y = ip[1:(T - 1)], F_T[1:(T - 1), ])
  new <- data.frame(y = ip[T], t(F_T[T, ]))
  e_ar <- c(e_ar, ip[T + 1] - predict(lm(y_next ~ y, d), new))
  e_fa <- c(e_fa, ip[T + 1] - predict(lm(y_next ~ ., d), new))
}

c(RMSE_AR = sqrt(mean(e_ar^2)), RMSE_FA = sqrt(mean(e_fa^2)))

#It would be cheating to estimate the factors once, using the full sample,
#and use those in the horse race because using the full sample to estimate the 
#factors introduces data leakage.

direct_race <- function(series, h = 12, r = 3) {
  lv <- md |>
    mutate(log_x = log(.data[[series]]),
           g = 100 * (log_x - lag(log_x))) |>
    filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))
  yh <- 100 * (lead(lv$log_x, h) - lv$log_x)      # growth over the next h months
  g  <- lv$g                                      # growth over the latest month
  
  e_ar <- NULL
  e_fa <- NULL
  for (T in init:(n - h)) {
    F_T <- prcomp(scale(Z[1:T, ]))$x[, 1:r, drop = FALSE]
    i <- 1:(T - h)
    d <- data.frame(yh = yh[i], g = g[i], F_T[i, , drop = FALSE])
    new <- data.frame(g = g[T], F_T[T, , drop = FALSE])
    e_ar <- c(e_ar, yh[T] - predict(lm(yh ~ g, d), new))
    e_fa <- c(e_fa, yh[T] - predict(lm(yh ~ ., d), new))
  }
  tibble(series, forecasts = length(e_ar),
         RMSE_AR = sqrt(mean(e_ar^2)), RMSE_FA = sqrt(mean(e_fa^2)),
         MAE_AR = mean(abs(e_ar)), MAE_FA = mean(abs(e_fa)),
         DM_pvalue = forecast::dm.test(e_ar, e_fa, h = h, power = 2,
                                       varestimator = "bartlett")$p.value)
}

direct_race("CPIAUCSL")

#The factors improve the forecast of 12-month CPI inflation. 
#The AR model produces an RMSE of 2.11, while the factor-augmented model 
#produces an RMSE of 1.89. This represents an approximately 10.4% 
#reduction in RMSE when the three factors are included. The MAE also decreases 
#from 1.58 to 1.35, further indicating improved forecasting accuracy. 
#The Diebold–Mariano test produces a p-value of 0.000443, allowing us to 
#reject the null hypothesis of equal predictive accuracy at the 5% significance
#level. Therefore, the factor-augmented model provides a statistically 
#significant improvement over the AR model in forecasting 12-month 
#CPI inflation.
