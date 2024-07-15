# plotDensity
library(sn)

plotDensity <- function(z) {
  
  axis2<-c(2005,2005,
           2006,2006,2006,2006,
           2007,2007,2007,2007,
           2008,2008,2008,2008,
           2009,2009,2009,2009,
           2010,2010,2010,2010,
           2011,2011,2011,2011,
           2012,2012,2012,2012,
           2013,2013,2013,2013,
           2014,2014,2014,2014,
           2015,2015,2015,2015,
           2016,2016,2016,2016,
           2017,2017,2017,2017,
           2018,2018,2018,2018,
           2019,2019,2019,2019,
           2020)
  
  
  xd<-rep(axis2, times=1, each=512)    
  yd=seq(-40,10,abs((10+40)/511))
  y2<-rep(yd,59)
  mydata2<-cbind(t(t(xd)),y2,z)
  df = as.data.frame(mydata2)
  fig<-plot_ly(df, x = df$y2, y = df$V1, z = df$z, split = df$V1, 
               type = "scatter3d",mode = "lines", 
               lines=list( size=1 , opacity=0.5), color = ~z)
  axx <- list(title = "")
  axy <- list(title = "")
  axz <- list(title = "", range=c(0,0.5))
  fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  
  fig
}

