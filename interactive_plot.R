# Interactive R script - run this in RStudio or R console
library(ggplot2)
library(readr)

# Read data
power_data <- read_csv("power_results.csv")

# Create plot
p <- ggplot(power_data, aes(x = effect_size, y = power)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 3) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray") +
  labs(
    x = "Effect Size",
    y = "Statistical Power",
    title = "Power Analysis: Effect Size vs Statistical Power"
  ) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent)

# Display plot (will show in RStudio Plots pane or R graphics device)
p