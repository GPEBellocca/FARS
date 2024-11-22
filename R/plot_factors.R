library(dplyr)
library(tidyr)

generate_factor_plot <- function(factor_name,data_plot,min,max) {
  
  plot_res <- data_plot %>% filter(Factors==factor_name) %>%
    ggplot() + geom_line(aes(x = Date,y = value),col="blue",alpha=0.5) +
    geom_ribbon(aes(x = as.Date(Date), ymax = UB,ymin=LB),alpha=0.3) +
    xlab("")+ 
    ylab("")+
    coord_cartesian(ylim = c(min, max))+ 
    theme_bw()+
    facet_wrap(vars(Factors),nrow = 3,) +
    theme(legend.box.spacing = unit(-10, "pt"),
          legend.margin=margin(0,0,0,0),
          legend.position="none", 
          legend.title = element_text(size=1)) 
   #labs(colour = paste0(""),title = paste0(factor_name,' factor'))
  
  #ggsave(filename = paste0("./plots/factors/",factor_name,".png"), plot = plot_res, width = 6, height = 4, dpi = 300)
  return(plot_res)
  
}


plot_factors <- function(MLDFM_result,dates = NULL) {
 
  
  Factors <- MLDFM_result$Factors
  Lambda <- MLDFM_result$Lambda
  Residuals <- MLDFM_result$Residuals
  
  # MSE CI
  t<-nrow(Residuals)
  N<-ncol(Residuals)
  
  PP<- solve((t(Lambda) %*% Lambda)/N)
  
  # Assumption uncorrelated idiosincratic
  sigma_e<-sum(diag(t(Residuals) %*% Residuals)/(N*t))
  Gamma<-sigma_e * (t(Lambda) %*% Lambda)/N
  SD<-sqrt(diag(PP %*% (Gamma %*% PP))/N)
  
 
  # Extract Factors
  factors<-Factors
  factors_df <- as.data.frame(factors)
  
  # Extract Factors names
  keys <- names(MLDFM_result$Factors_list)
  values <- unlist(MLDFM_result$Factors_list) 
  transformed_keys <- c()
  
  for (i in seq_along(keys)) {
   
    clean_key <- paste0("F", gsub("-", "", keys[i]))
    if(values[i]>1){
      repeated_keys <- paste0(clean_key,"n",seq_len(values[i]))
    }else{
      repeated_keys <- clean_key
    }
    transformed_keys <- c(transformed_keys, repeated_keys)
  }
  colnames(factors_df)<-transformed_keys
 
  
  # Add dates
  if (is.null(dates)) {
    dates <- 1:nrow(factors_df)
  }
  data_plot <- data.frame(Date = dates, factors_df)
  
  
  # add SD
  data_plot <- data_plot %>%
    pivot_longer(cols = -Date, values_to = "value", names_to = "Factors") %>%
    mutate(
      LB = value - 2 * SD[as.numeric(factor(Factors, levels = colnames(factors_df)))],
      UB = value + 2 * SD[as.numeric(factor(Factors, levels = colnames(factors_df)))]
    )
  
  
  # Compute scale 
  min <- min(data_plot$value, na.rm = TRUE)
  max <- max(data_plot$value, na.rm = TRUE)
 
 
 
  # Create factor plot
  plot_list <- list()
  for (factor_name in colnames(factors_df)) {
    p <- generate_factor_plot(factor_name,data_plot,min,max)
    plot_list[[length(plot_list) + 1]] <- p
  }
  
  return(plot_list)
  
}



