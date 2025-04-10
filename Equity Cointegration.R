library(tseries)  

cmcsa_data <- read.csv("C:/Users/bryns/downloads/CMCSA.csv")
chtr_data <- read.csv("C:/Users/bryns/downloads/CHTR.csv")


cmcsa_data$Date <- as.Date(cmcsa_data$Date)
chtr_data$Date <- as.Date(chtr_data$Date)


cmcsa_data$Return <- c(NA, diff(log(cmcsa_data$Adj.Close)))
chtr_data$Return <- c(NA, diff(log(chtr_data$Adj.Close)))


merged_data <- merge(cmcsa_data[, c("Date", "Return")], chtr_data[, c("Date", "Return")],
                     by = "Date", suffixes = c("_CMCSA", "_CHTR"))
merged_data <- na.omit(merged_data)


cointegration_model <- lm(Return_CMCSA ~ Return_CHTR, data = merged_data)
summary(cointegration_model)

merged_data$residuals <- resid(cointegration_model)

adf_result <- adf.test(merged_data$residuals)
cat("ADF Test on residuals:\n")
print(adf_result)

merged_data$ECT <- lag(merged_data$residuals, 1)

merged_data$Delta_CMCSA <- c(NA, diff(merged_data$Return_CMCSA))
merged_data$Delta_CHTR <- c(NA, diff(merged_data$Return_CHTR))
merged_data$Delta_CMCSA_Lag1 <- lag(merged_data$Delta_CMCSA, 1)
merged_data$Delta_CHTR_Lag1 <- lag(merged_data$Delta_CHTR, 1)


vecm_data <- na.omit(merged_data)

vecm_cm_model <- lm(Delta_CMCSA ~ ECT + Delta_CMCSA_Lag1 + Delta_CHTR_Lag1, data = vecm_data)
vecm_ch_model <- lm(Delta_CHTR ~ ECT + Delta_CMCSA_Lag1 + Delta_CHTR_Lag1, data = vecm_data)

cat("\nVECM Model for CMCSA:\n")
print(summary(vecm_cm_model))
cat("\nVECM Model for CHTR:\n")
print(summary(vecm_ch_model))

last_obs <- tail(vecm_data, 1)

next_return_cm <- coef(vecm_cm_model)[1] +
  coef(vecm_cm_model)["ECT"] * last_obs$ECT +
  coef(vecm_cm_model)["Delta_CMCSA_Lag1"] * last_obs$Delta_CMCSA_Lag1 +
  coef(vecm_cm_model)["Delta_CHTR_Lag1"] * last_obs$Delta_CHTR_Lag1

next_return_ch <- coef(vecm_ch_model)[1] +
  coef(vecm_ch_model)["ECT"] * last_obs$ECT +
  coef(vecm_ch_model)["Delta_CMCSA_Lag1"] * last_obs$Delta_CMCSA_Lag1 +
  coef(vecm_ch_model)["Delta_CHTR_Lag1"] * last_obs$Delta_CHTR_Lag1

cat("\nPredicted next period return for CMCSA:", next_return_cm, "\n")
cat("Predicted next period return for CHTR:", next_return_ch, "\n")

vecm_data$Forecast_CMCSA <- fitted(vecm_cm_model)
vecm_data$Forecast_CHTR <- fitted(vecm_ch_model)

errors_cm <- vecm_data$Delta_CMCSA - vecm_data$Forecast_CMCSA
errors_ch <- vecm_data$Delta_CHTR - vecm_data$Forecast_CHTR

rmsfe_cm <- sqrt(mean(errors_cm^2))
rmsfe_ch <- sqrt(mean(errors_ch^2))

cat("\nRMSFE for CMCSA:", rmsfe_cm, "\n")
cat("RMSFE for CHTR:", rmsfe_ch, "\n")
