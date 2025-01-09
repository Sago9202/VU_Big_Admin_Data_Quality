# Load required library
library(ggplot2)
library(cols4all)
rm(list = ls())
sys <- Sys.info()
fold_graphs <- paste0("C:/Users/", sys[7], "/Dropbox/PhD VU/Papers/1 - Selection Bias Simulations/4 - Graphs & Tables")

# Simulate data
set.seed(42)
days <- 0:360

# Get Spectral colors
spectral_colors <- c4a("carto.pastel", n = 4)

simulate_response <- function(day, midpoint, growth_rate) {
  1 / (1 + exp(-growth_rate * (day - midpoint)))
}

# Parameters for each category
quarter <- c("Q1", "Q2", "Q3", "Q4")
midpoints <- c(45, 135, 225, 315)  # Midpoint for logistic growth
growth_rates <- c(0.3, 0.3, 0.3, 0.3)  # Growth rates

# Create data frame
data <- data.frame(
  Day = rep(days, times = 4),
  Response = c(
    simulate_response(days, midpoints[1], growth_rates[1]),
    simulate_response(days, midpoints[2], growth_rates[2]),
    simulate_response(days, midpoints[3], growth_rates[3]),
    simulate_response(days, midpoints[4], growth_rates[4])
  ),
  Quarter = factor(rep(quarter, each = length(days)))
)

# Define shaded regions for the quarters
quarters <- data.frame(
  xmin = c(0, 90, 180, 270),
  xmax = c(90, 180, 270, 360),
  label = c("Q2", "Q3", "Q4", "Q1")
)

data$Response[data$Quarter=="Q4" & data$Day<270] <- NA
data$Response[data$Quarter=="Q3" & data$Day<180] <- NA
data$Response[data$Quarter=="Q2" & data$Day<90] <- NA

# Create the plot
ggplot(data, aes(x = Day, y = Response, color = Quarter)) +
  geom_line(size = 0.7) +
  geom_point(size = 0.7) + 
  
  # Add vertical black lines for quarter boundaries
  geom_vline(xintercept = c(90, 180, 270, 360), color = "black") +
  
  # Add vertical dashed red lines
  geom_vline(xintercept = midpoints, linetype = "dashed", color = "darkred") +
  # Add shaded regions for quarters
  geom_rect(
    data = quarters,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = label),
    alpha = 0.15,
    inherit.aes = F
  ) + 
  
  # Customize axis labels and theme
  labs(
    x = "t - Days after the end of Yr1-Q1",
    y = "Completion rate", 
    color = "",
  ) +
  scale_color_manual(values = spectral_colors) + # Apply Spectral colors
  scale_fill_manual(values = spectral_colors, guide = "none") + # Apply Spectral colors to fills but hide the legend
  scale_x_continuous(breaks = c(0, 45, 90, 135, 180, 225, 270, 315, 360),
                     sec.axis = dup_axis(breaks = c(45, 135, 225, 315), labels = c("Yr1-Q2", "Yr1-Q3", "Yr1-Q4", "Yr2-Q1"), name = "Year-Quarter")) +
  theme_minimal() +
  theme(text = element_text(size = 15),
    axis.text.x.top = element_text(color = "darkblue"),
    axis.title.x.top = element_text(color = "darkblue"))
ggsave(filename = "Grad_Plot.png", path = fold_graphs, width = 30, height = 15, units = "cm")    



