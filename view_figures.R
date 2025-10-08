# Script to display figures in R
library(ggplot2)
library(readr)

# Read the data
power_data <- read_csv("power_results.csv")

# Create and display the plot
p <- ggplot(power_data, aes(x = effect_size, y = power)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 3) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray") +
  labs(
    x = "Effect Size",
    y = "Statistical Power",
    title = "Power Analysis: Effect Size vs Statistical Power",
    subtitle = "Dashed line shows 80% power threshold"
  ) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent)

# Display the plot
print(p)

# Keep the plot window open
readline("Press Enter to close the plot window...")