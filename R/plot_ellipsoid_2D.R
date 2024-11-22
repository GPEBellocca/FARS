library(ggplot2)


plot_ellipsoid_2D <- function(obs_number, Ellipsoids) {
 
  ellipsoid_points <- Ellipsoids[[obs_number]]
  
  
  ellipsoid_2D <- data.frame(x = ellipsoid_points[, 1], y = ellipsoid_points[, 2])
  
  
  p <- ggplot(ellipsoid_2D, aes(x = x, y = y)) +
    geom_point(color = "blue") +   
    ggtitle(paste("2D Ellipsoid for Observation", obs_number)) +
    xlab("Factor 1") + 
    ylab("Factor 2") +
    theme_minimal()
  
 
  print(p)
}

