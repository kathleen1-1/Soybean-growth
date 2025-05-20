# Load necessary libraries
# Data manipulation and base
library(forecast)       # Used for time series models (though not directly in your shared code)
library(lmtest)         # For gqtest(), dwtest(), coeftest()
library(sandwich)       # For robust standard errors (e.g., NeweyWest)
library(MASS)           # For stepAIC() function
library(car)            # For vif() - Variance Inflation Factor

# Visualization and regression effects
library(coefplot)       # For plotting regression coefficients
library(effects)        # For visualizing effects of predictors


# Load data
data <- read.csv("Data B.csv")  # Replace with your actual path

# Convert factors
data$Season <- as.factor(data$Season)
data$Cultivar <- as.factor(data$Cultivar)
data$Repetition <- as.factor(data$Repetition)

# Log transform target
data$log_GY <- log(data$GY)

#mlr
mlr_model <- lm(log_GY ~ Season + Cultivar + Repetition + PH + IFP + NLP + NGP + NGL + NS + MHG, data = data)
summary(mlr_model)
print(vif(mlr_model))
gqtest(mlr_model)
dwtest(mlr_model)

# --- Stepwise Regression ---
stepwise_model <- stepAIC(mlr_model, direction = "both")
summary(stepwise_model)

# Check for multicollinearity, heteroscedasticity, and autocorrelation for stepwise model
print(vif(stepwise_model))
gqtest(stepwise_model)
dwtest(stepwise_model)

# Check for multicollinearity, heteroscedasticity, and autocorrelation for stepwise model
print(vif(stepwise_model))
gqtest(stepwise_model)
dwtest(stepwise_model)

coeftest(stepwise_model, vcov = NeweyWest)

# Residual diagnostics for stepwise model
par(mfrow = c(2, 2))  # 2x2 layout
plot(stepwise_model)
par(mfrow = c(1, 1))  # reset


# Sorted by size, remove intercept for clarity
coefplot(stepwise_model, sort = "magnitude", intercept = FALSE, title = "Influence of Agronomic Variables")


library(car)
library(ggplot2)
library(dplyr)

# --- Function to prepare data ---
prepare_vif_data <- function(model) {
  vif_vals <- vif(model)
  if (is.matrix(vif_vals)) {
    adj_vif <- vif_vals[, "GVIF"]^(1 / (2 * vif_vals[, "Df"]))
    var_names <- rownames(vif_vals)
  } else {
    adj_vif <- vif_vals
    var_names <- names(vif_vals)
  }
  data.frame(Variable = var_names, VIF = adj_vif)
}

# --- MLR Model ---
vif_data_mlr <- prepare_vif_data(mlr_model)

ggplot(vif_data_mlr, aes(x = reorder(Variable, VIF), y = VIF)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
  geom_text(aes(label = round(VIF, 3)), hjust = -0.1, size = 4) +
  coord_flip() +
  labs(title = "Adjusted VIF - MLR Model", x = "Variable", y = "Adjusted VIF") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# --- Stepwise Model ---
vif_data_step <- prepare_vif_data(stepwise_model)

ggplot(vif_data_step, aes(x = reorder(Variable, VIF), y = VIF)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
  geom_text(aes(label = round(VIF, 3)), hjust = -0.1, size = 4) +
  coord_flip() +
  labs(title = "Adjusted VIF - Stepwise", x = "Variable", y = "Adjusted VIF") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Residuals vs Fitted
plot(stepwise_model, which = 1)

# Normal Q-Q plot
plot(stepwise_model, which = 2)

# Scale-Location plot (also for homoscedasticity)
plot(stepwise_model, which = 3)

# Residuals vs Leverage (for influential points)
plot(stepwise_model, which = 5)


