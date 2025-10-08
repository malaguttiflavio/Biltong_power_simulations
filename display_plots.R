# Script to ensure plots are displayed properly
library(ggplot2)
library(readr)

# Read the data
power_data <- read_csv("power_results.csv")
print("Data loaded:")
print(power_data)

# Create the plot
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

# Method 1: Try to display directly
cat("Attempting to display plot...\n")
print(p)

# Method 2: Force display with explicit device
if(interactive()) {
  # If running interactively, try different display methods
  if(capabilities("X11")) {
    X11()
    print(p)
  } else if(.Platform$OS.type == "unix" && Sys.info()["sysname"] == "Darwin") {
    # macOS
    quartz()
    print(p)
  } else {
    # Windows
    windows()
    print(p)
  }
}

# Method 3: Save again with explicit dimensions
ggsave("current_power_plot.png", p, width = 10, height = 6, dpi = 150)
ggsave("current_power_plot.pdf", p, width = 10, height = 6)

cat("Plots saved as:\n")
cat("- current_power_plot.png\n")
cat("- current_power_plot.pdf\n")
cat("Files are located in:", getwd(), "\n")

# Method 4: Try to open the file automatically
if(.Platform$OS.type == "unix" && Sys.info()["sysname"] == "Darwin") {
  # macOS
  system("open current_power_plot.png")
} else if(.Platform$OS.type == "windows") {
  # Windows
  shell.exec("current_power_plot.png")
} else {
  # Linux
  system("xdg-open current_power_plot.png")
}