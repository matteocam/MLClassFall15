library(ggplot2)
library(plyr)

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# dataSum <- summarySE(data,measurevar = "r", groupvars = c("K", "DistanceType", "ND"))

plotData <- function()
{
  # Read file
  data <- read.csv("~/Proiects/ML_exp1/good-analysis.dat", sep="")
  
  # Make Summary
  dataSum <- summarySE(data,measurevar = "r", groupvars = c("K", "DistanceType", "ND"))

  # Filter by K

    mkPlot <- function(KTarget)
    {
      dataSum1 <- dataSum[dataSum$K == KTarget,]
      
      # Plot
      return (ggplot(dataSum1, aes(x=ND, y=r, colour=DistanceType, group = DistanceType)) + 
                        geom_errorbar(aes(ymin=r-se, ymax=r+se), width=5.0) +
                        geom_line() +
                        geom_point() )
    }
    
    plot1000 <- mkPlot(1000)
    ggsave("exp1-plot1000.pdf", width=9, height=5, dpi=150)
    
    plot5000 <- mkPlot(5000)
    ggsave("exp1-plot5000.pdf", width=9, height=5, dpi=150)
    
    plot10000 <- mkPlot(10000)
    ggsave("exp1-plot10000.pdf", width=9, height=5, dpi=150)
    
    #library(gridExtra)
    #grid.arrange(plot1000, plot5000, plot10000, nrow=3)
    #ggsave("analysis.pdf", width=9, height=5, dpi=150)
  
}
