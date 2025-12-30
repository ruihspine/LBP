library(dplyr)

# Define the dataset with prevalence data
data <- data.frame(
  Age = c("95+ years", "65-69 years", "70-74 years", "75-79 years", "80-84 years", "85-89 years", "90-94 years"),
  Number = c(31895.8189, 1218160.7205, 941368.7506, 667776.4237, 464268.7814, 256695.5320, 108309.1373),
  Rate = c(668.2260, 471.0896, 503.1713, 525.5858, 549.9361, 590.3639, 642.4950)
)

# Define the population weights based on the percentage of each age group
population_data <- data.frame(
  Age = c("65-69 years", "70-74 years", "75-79 years", "80-84 years", "85-89 years", "90-94 years", "95+ years"),
  PopulationPercent = c(2.99, 2.27, 1.61, 1.11, 0.62, 0.26, 0.08)  # Values should match the age groups and be in the same order
)

# Merge datasets by Age
data <- merge(data, population_data, by = "Age")

# Calculate the age-standardized rate using the population percentages as weights
data <- mutate(data, WeightedRate = Rate * PopulationPercent)
total_population_percent <- sum(data$PopulationPercent)
age_standardized_rate <- sum(data$WeightedRate) / total_population_percent

# Calculate standard error, lower and upper 95% CI using the Number
data <- mutate(data,
               SE = Rate / sqrt(Number),
               LowerCI = Rate - 1.96 * SE,
               UpperCI = Rate + 1.96 * SE)

# Aggregate the 95% CI using the weighted sum of SEs
data <- mutate(data, WeightedSE = SE * PopulationPercent)
overall_SE <- sqrt(sum(data$WeightedSE^2))  # Weighted sum of SEs
lower_95CI <- age_standardized_rate - 1.96 * overall_SE
upper_95CI <- age_standardized_rate + 1.96 * overall_SE

# Output results
results <- data.frame(
  Age_Standardized_Rate = age_standardized_rate,
  Lower_95_CI = lower_95CI,
  Upper_95_CI = upper_95CI
)

print(results)
