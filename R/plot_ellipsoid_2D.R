library(ggplot2)

# Function to plot the ellipsoid for a given observation in 2D
plot_ellipsoid_2D <- function(obs_number, Ellipsoids) {
  # Select the ellipsoid points for the given observation
  ellipsoid_points <- Ellipsoids[[obs_number]]
  
  # Create a data frame with the two factors (dimensions)
  ellipsoid_2D <- data.frame(x = ellipsoid_points[, 1], y = ellipsoid_points[, 2])
  
  # Plot the ellipsoid using ggplot2
  p <- ggplot(ellipsoid_2D, aes(x = x, y = y)) +
    geom_point(color = "blue") +   # Plot points for the ellipsoid
    ggtitle(paste("2D Ellipsoid for Observation", obs_number)) +
    xlab("Factor 1") + 
    ylab("Factor 2") +
    theme_minimal()
  
  # Display the plot
  print(p)
}

