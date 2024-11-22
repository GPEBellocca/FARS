library(plotly)

plot_distribution <- function(Density_result, dates = NULL) {
  
  density_matrix <- Density_result$density_matrix  # Matrix of densities (T x est_points)
  t <- nrow(density_matrix)  # N observations
  #est_points <- seq(-20, 20, length.out = ncol(density_matrix))  
  est_points <- seq(min(density_matrix), max(density_matrix), length.out = ncol(density_matrix))
  
  # Generate default dates if not provided
  if (is.null(dates)) {
    dates <- seq_len(t)
  }
  
  # Convert matrix into a format suitable for 3D plotting
  density_data <- as.data.frame(expand.grid(est_points, dates))
  density_data$density <- as.vector(t(density_matrix))
  colnames(density_data) <- c("X", "Y", "Z")
  
  # Create 3D surface plot
  plot <- plot_ly(
    data = density_data,
    x = ~X,
    y = ~Y,
    z = ~Z,
    type = "scatter3d",
    mode = "lines",
    line = list(width = 2),
    color = ~Z  # Use Z (density) for color mapping
  ) %>%
    layout(
      scene = list(
        xaxis = list(title = "Estimated Points"),
        yaxis = list(title = "Time Observations"),
        zaxis = list(title = "Density")
      )
    )
  
  return(plot)
}
