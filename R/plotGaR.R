# plotGaR
library(sn)

plotGaR <- function(dep_variable,distribution)
 {
  
  N_time <- length(dep_variable)
  
  QGaR.05 <- GaR(distribution, QTAU = 0.05) 
  QGaR.95 <- GaR(distribution, QTAU = 0.95) 
  
  print(QGaR.05)
  print(QGaR.95)
  
  
  axisFig<-c("2005Q3","2005Q4",
             "2006Q1","2006Q2","2006Q3","2006Q4",
             "2007Q1","2007Q2","2007Q3","2007Q4",
             "2008Q1","2008Q2","2008Q3","2008Q4",
             "2009Q1","2009Q2","2009Q3","2009Q4",
             "2010Q1","2010Q2","2010Q3","2010Q4",
             "2011Q1","2011Q2","2011Q3","2011Q4",
             "2012Q1","2012Q2","2012Q3","2012Q4",
             "2013Q1","2013Q2","2013Q3","2013Q4",
             "2014Q1","2014Q2","2014Q3","2014Q4",
             "2015Q1","2015Q2","2015Q3","2015Q4",
             "2016Q1","2016Q2","2016Q3","2016Q4",
             "2017Q1","2017Q2","2017Q3","2017Q4",
             "2018Q1","2018Q2","2018Q3","2018Q4",
             "2019Q1","2019Q2","2019Q3","2019Q4",
             "2020Q1")
  
  
  time <- 1:N_time
  MLGaRGiS <- data.frame(
    year=1:N_time,
    growth1=as.vector(dep_variable),
    GaR.05=as.vector(QGaR.05),
    GaR.95=as.vector(QGaR.95)
  )
  
  p <- ggplot(MLGaRGiS,aes(x=time,y=growth1)) + 
    geom_line() + theme_bw()+
    geom_line(aes(x=time,GaR.05),linewidth=0.1, colour="black",linetype = "dashed")+
    geom_line(aes(x=time,GaR.95),linewidth=0.1, colour="black",linetype = "dashed")+
    scale_y_continuous("Growth")+
    scale_x_continuous("",breaks=1:59, labels=axisFig)+
    theme(axis.text.x = element_text(angle = 90))
  
  fig <- ggplotly(p)
  fig
}

