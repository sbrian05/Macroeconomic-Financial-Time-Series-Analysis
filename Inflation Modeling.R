
library(dplyr)

cpi_data <- read.csv("C:/Users/bryns/downloads/CPIAUCSL.csv")
cpi_data$DATE <- as.Date(cpi_data$DATE)

cpi_data <- cpi_data %>%
  arrange(DATE) %>%
  mutate(Inflation = 1200 * (log(CPIAUCSL) - lag(log(CPIAUCSL))))

filtered_data <- cpi_data %>%
  filter(DATE >= as.Date("1970-01-01") & DATE <= as.Date("2012-12-01"))

plot(filtered_data$DATE, filtered_data$Inflation, type = "l",
     xlab = "Date", ylab = "Annualized Inflation Rate",
     main = "Inflation Rate: 1970-2012")

filtered_data <- filtered_data %>%
  mutate(Diff_Inflation = Inflation - lag(Inflation)) %>%
  filter(!is.na(Diff_Inflation))

acf(filtered_data$Inflation, lag.max = 12, main = "ACF: Inflation")
acf(filtered_data$Diff_Inflation, lag.max = 12, main = "ACF: First Difference of Inflation")

plot(filtered_data$DATE, filtered_data$Diff_Inflation, type = "l",
     xlab = "Date", ylab = "First Difference of Inflation",
     main = "Δ Inflation: 1970-2012")

filtered_data <- filtered_data %>% mutate(Lag_Inflation = lag(Inflation)) %>% filter(!is.na(Lag_Inflation))
ols_model <- lm(Inflation ~ Lag_Inflation, data = filtered_data)
summary(ols_model)

filtered_data <- filtered_data %>% mutate(Lag_Inflation_2 = lag(Inflation, 2)) %>% filter(!is.na(Lag_Inflation_2))
ar1 <- lm(Inflation ~ Lag_Inflation, data = filtered_data)
ar2 <- lm(Inflation ~ Lag_Inflation + Lag_Inflation_2, data = filtered_data)
summary(ar1)
summary(ar2)

cat("Residual SE AR(1):", summary(ar1)$sigma, "\n")
cat("Residual SE AR(2):", summary(ar2)$sigma, "\n")

aic_bic <- data.frame(Lag = 0:8, AIC = NA, BIC = NA)

for (p in 0:8) {
  temp_data <- filtered_data
  for (i in 1:p) {
    temp_data <- temp_data %>% mutate(!!paste0("Lag", i) := lag(Inflation, i))
  }
  temp_data <- temp_data %>% filter(row_number() > p)

  if (p == 0) {
    model <- lm(Inflation ~ 1, data = temp_data)
  } else {
    formula <- as.formula(paste("Inflation ~", paste(paste0("Lag", 1:p), collapse = "+")))
    model <- lm(formula, data = temp_data)
  }

  n <- nrow(temp_data)
  k <- length(coef(model))
  rss <- sum(residuals(model)^2)
  aic_bic$AIC[p + 1] <- n * log(rss / n) + 2 * k
  aic_bic$BIC[p + 1] <- n * log(rss / n) + log(n) * k
}

print(aic_bic)
cat("Best lag by AIC:", which.min(aic_bic$AIC) - 1, "\n")
cat("Best lag by BIC:", which.min(aic_bic$BIC) - 1, "\n")

last_two <- tail(filtered_data$Inflation, 2)
predicted_inflation <- coef(ar2)[1] + coef(ar2)[2] * last_two[2] + coef(ar2)[3] * last_two[1]
cat("Predicted Inflation 2013:M01:", predicted_inflation, "\n")

adf_data <- filtered_data %>% mutate(Lag_Diff1 = lag(Diff_Inflation, 1), Lag_Diff2 = lag(Diff_Inflation, 2)) %>%
  filter(!is.na(Lag_Diff1), !is.na(Lag_Diff2))
adf_model <- lm(Diff_Inflation ~ Lag_Inflation + Lag_Diff1 + Lag_Diff2, data = adf_data)
summary(adf_model)

n <- nrow(filtered_data)
trim <- floor(0.15 * n)
start <- trim + 1
end <- n - trim
f_stats <- numeric(end - start + 1)

for (bp in start:end) {
  model_before <- lm(Inflation ~ Lag_Inflation + Lag_Inflation_2, data = filtered_data[1:bp, ])
  model_after <- lm(Inflation ~ Lag_Inflation + Lag_Inflation_2, data = filtered_data[(bp + 1):n, ])
  rss1 <- sum(residuals(model_before)^2)
  rss2 <- sum(residuals(model_after)^2)
  rss_combined <- rss1 + rss2

  full_model <- lm(Inflation ~ Lag_Inflation + Lag_Inflation_2, data = filtered_data)
  rss_full <- sum(residuals(full_model)^2)

  k <- length(coef(full_model))
  f_stats[bp - start + 1] <- ((rss_full - rss_combined) / k) / (rss_combined / (n - 2 * k))
}

QLR_stat <- max(f_stats)
cat("QLR Test Statistic:", QLR_stat, "\n")

forecast_indices <- which(filtered_data$DATE >= as.Date("2005-12-01"))
forecasts <- numeric(length(forecast_indices))

for (i in seq_along(forecast_indices)) {
  train_data <- filtered_data[1:(forecast_indices[i] - 1), ]
  model <- lm(Inflation ~ Lag_Inflation + Lag_Inflation_2, data = train_data)
  forecasts[i] <- predict(model, newdata = filtered_data[forecast_indices[i], ])
}

results <- data.frame(
  Date = filtered_data$DATE[forecast_indices],
  Actual = filtered_data$Inflation[forecast_indices],
  Forecast = forecasts
)

results$Error <- results$Actual - results$Forecast
mean_error <- mean(results$Error)
t_test <- t.test(results$Error, mu = 0)
cat("Mean Forecast Error:", mean_error, "\n")
print(t_test)

rmsfe <- sqrt(mean(results$Error^2))
cat("RMSFE:", rmsfe, "\n")

train_data <- filtered_data %>% filter(DATE <= as.Date("2005-12-01"))
final_model <- lm(Inflation ~ Lag_Inflation + Lag_Inflation_2, data = train_data)
cat("Residual SE of AR(2) model:", summary(final_model)$sigma, "\n")
