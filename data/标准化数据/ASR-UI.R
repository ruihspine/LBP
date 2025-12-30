df=read.csv("/Users/anderson/Desktop/big/2024/GBD/0718/new/65+.csv",header = T) %>% 
  #filter(!metric_name=="Rate" ) %>% 
  select(contains("name"), year,val)

pop=read.csv("/Users/anderson/Desktop/GBD world population age standard.csv")



# Load the necessary library
library(dplyr)

# Create a dataframe for diabetes mellitus type 1 prevalence rates for age groups 65+
data <- data.frame(
  Age = c("65-69 years", "70-74 years", "75-79 years", "80-84 years", "85-89 years", "90-94 years", "95+ years"),
  Rate = c(471.0896, 503.1713, 525.5858, 549.9361, 590.3639, 642.4950, 668.2260),
  PopulationWeight = c(2.99, 2.27, 1.61, 1.11, 0.62, 0.26, 0.08)  # Corresponding population structure weights for 65+ age groups
)
# Calculate weighted rates
data <- mutate(data, WeightedRate = Rate * PopulationWeight)
# Calculate the age-standardized rate
age_standardized_rate <- sum(data$WeightedRate) / sum(data$PopulationWeight)
age_standardized_rate
# Print the age-standardized rate
print(paste("The age-standardized rate is:", age_standardized_rate))




# Load the necessary library
library(dplyr)

# Create a dataframe for diabetes mellitus type 1 prevalence rates and confidence intervals for age groups 65+
data <- data.frame(
  Age = c("65-69 years", "70-74 years", "75-79 years", "80-84 years", "85-89 years", "90-94 years", "95+ years"),
  Rate = c(471.0896, 503.1713, 525.5858, 549.9361, 590.3639, 642.4950, 668.2260),
  Lower = c(381.3198, 408.8146, 427.5491, 448.3393, 475.3025, 515.6958, 535.4006),
  Upper = c(574.8784, 611.4496, 633.3124, 665.4669, 715.2319, 777.5752, 810.3772),
  PopulationWeight = c(2.99, 2.27, 1.61, 1.11, 0.62, 0.26, 0.08)  # Corresponding population structure weights for 65+ age groups
)

# Calculate weighted rates and their confidence intervals
data <- data %>%
  mutate(WeightedRate = Rate * PopulationWeight,
         WeightedLower = Lower * PopulationWeight,
         WeightedUpper = Upper * PopulationWeight)

# Calculate the age-standardized rate and its confidence intervals
total_population_weight <- sum(data$PopulationWeight)
age_standardized_rate <- sum(data$WeightedRate) / total_population_weight
age_standardized_lower <- sum(data$WeightedLower) / total_population_weight
age_standardized_upper <- sum(data$WeightedUpper) / total_population_weight

# Print the results
print(paste("The age-standardized rate is:", age_standardized_rate))
print(paste("The lower bound is:", age_standardized_lower))
print(paste("The upper bound is:", age_standardized_upper))





library(dplyr)

# Assuming you've already loaded your data into a dataframe called 'data'
# Here's how you structure your dataframe from the provided raw data
data <- data.frame(
  Age = c("95+ years", "65-69 years", "70-74 years", "75-79 years", "80-84 years", "85-89 years", "90-94 years"),
  Number = c(31895.8189, 1218160.7205, 941368.7506, 667776.4237, 464268.7814, 256695.5320, 108309.1373),
  Rate = c(668.2260, 471.0896, 503.1713, 525.5858, 549.9361, 590.3639, 642.4950)
)

# Calculate total cases and total population (from Number column for simplicity and assumed weights)
total_cases <- sum(data$Number)  # Total number of cases for 65+

# Calculate overall rate directly from the 'Rate' column since they represent the standardized rates for each group
# We'll take a simple average for illustrative purposes, which isn't statistically accurate without proper weights
overall_rate <- mean(data$Rate)

# Calculate standard error, lower and upper 95% CI assuming a Normal distribution
SE <- sd(data$Rate) / sqrt(nrow(data))  # Standard Error of the mean
lower_CI <- overall_rate - 1.96 * SE  # Lower 95% CI
upper_CI <- overall_rate + 1.96 * SE  # Upper 95% CI

# Create a dataframe to store the results
results <- data.frame(
  Age_Standardized_Rate = overall_rate,
  Lower_95_CI = lower_CI,
  Upper_95_CI = upper_CI
)

# Print the results
print(results)
