#####
# Function for plotting density or cumulative density in variant and non-variant conditions
# The script is distributed in association with the ADPenetrance (https://github.com/ThomasPSpargo/adpenetrance/) approach to calculate genetic penetrance
# This allows testing for equal onset variability across two groups 
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk) 
# Please get in touch if you have any issues.
#####

#Purpose:
# This can be used to give an indication about whether penetrance estimates might be affected substantially by age of sampling,
# when the difference in variability between groups is large
# the main inputs are provided to: dataset OR groups and age onset should be provided

## INPUT:
#dataset      = a data.frame that must contain columns columns named "groups" and "ageonset" (additional columns will be ignodred).See groups and ageonset arguments for further details
#groups       = factor or character vector flagging to which group the corresponding 'ageonset' value corresponds. If this has factor class, factor order will be retained.
#ageonset     = numeric vector indicating the age at which disease onsets, becoming penetrant, for each person.
#center_on    = should be either "mean", "median", or a numeric between 0 and 1. These determine how age of onset values will be centered for each group. Numeric values are parsed by the quantile() function
#plotType     = character vector of either "density" or "cumulative density" and indicates whether to plot the density curve or cumulative density function for each group.
#...          = arguments passed to plot(); cannot currently accept x,y,col,ylim,xlim,xlab,type,xlab,ylab,sub since these are defined explicitly within the main function

## OUTPUT:
# Base R plot comparing variability in onset between two groups
# the plot 'sub' message indicates the relative difference in span of time between the first and third quartiles of onset.
#   Values closer to 1 represent a that the variability in age of onset is similar between groups.
#   Values <1 indicate that the group plotted in blue has a longer interquartile interval
#   Values >1 indicate that the group plotted in red has a a longer interquartile interval
#   A longer-form version of the variability in age of onset message is also printed to the console,
#   giving slightly more explicit detail about how the number relates to the two supplied groups

checkOnsetVariability <- function(dataset=NULL,groups=NULL,ageonset=NULL,centre_on="mean",plotType="density",...){
  
  if(!is.null(dataset)){
    data <- dataset[,c("groups","ageonset")]
  } else if (!is.null(groups) && !is.null(ageonset)){
    data <- data.frame(groups=groups,ageonset=ageonset)
  }
  
  
  if(is.factor(data$groups)){
    #If groups are provided as a factor, remove any empty levels and extract levels in order
    xlevels<- levels(droplevels(data$groups))
  } else {
    xlevels<- unique(data$groups)
  }
  #Split the dataset by group allocations
  data1 <- data[which(data$groups==xlevels[1]),]
  data2 <- data[which(data$groups==xlevels[2]),]
  
  #Function to prepare for density plotting, and calculate time interval between onset range (e.g. IQR)
  getDist <- function(data){
    if(centre_on=="mean"){
      #Identify mean as center point
      centrepoint <- mean(data$ageonset,na.rm = TRUE)
    } else if(centre_on=="median"){
      #Identify median as center point
      centrepoint <- median(data$ageonset,na.rm = TRUE)
    } else if(is.numeric(centre_on)){
      #Center on the numeric value indicated by quantile
      centrepoint <- quantile(data$ageonset,centre_on,na.rm=TRUE)
    }
    #Center age of onset about about the center point
    data$centred <- data$ageonset-centrepoint
    
    #Extract cumulative age of onset
    x <- sort(data$centred)
    y <- 1:length(x)/length(x)
    
    #Test difference in time across IQR for onset
    time1 <- 0.25
    startage<- x[which(abs(y - time1) == min(abs(y - time1)))]
    
    time2 <- 0.75
    endage<- x[which(abs(y - time2) == min(abs(y - time2)))]
    
    agediff <- endage-startage
    if(plotType=="density"){
      #Calculate density function around center point
      dense<-density(data$centred)
    } else if (plotType=="cumulative"){
      
      dense<- cbind(x,y)
    }
    
    return(list(dense=dense,agediff_time=agediff))
  }
  Dist1 <- getDist(data1)
  Dist2 <- getDist(data2)
  
  #Calculate density function about 
  d1 <- Dist1$dense
  d2 <- Dist2$dense
  
  #obtain time difference between density range
  time_g1 <- Dist1$agediff_time
  time_g2 <- Dist2$agediff_time
  
  t_diff <- time_g1/time_g2
  
  t_diff_message <- paste0("[",round(t_diff,3), " relative difference in onset variability]")
  
  message(paste0("The time taken for disease occurences across the interquartile range (IQR) to onset in group ",xlevels[2]," is ", round(t_diff,3), " times shorter than in group ",xlevels[1],".\n"))
  
  if(is.numeric(centre_on)){
    #Rename center on for labelling
    centre_on <- paste0(centre_on*100,"% quantile")
  }
  
  
  if(plotType=="density"){
  #Determine plotting limits according to limits from each density function
  maxy<- max(d1$y,d2$y)
  rangex<- range(d1$x,d2$x)
  
  #Plot with base R plotting functions
  
  
  #Main plot
  plot(d1,
       ...,
       col="red",
       ylim=c(0,maxy),
       xlim=rangex,
       xlab=paste0("Age of onset, centered by group ", centre_on),
       sub=t_diff_message
       )
  #Overlay 2nd density curve
  lines(d2,
        col="blue",
        lty="dashed")
  
  #Include legend
  legend("topright", legend = xlevels,
        col = c("red", "blue"), 
        lty = c("solid","dashed"), cex = 0.6)
  
  } else if(plotType=="cumulative"){
    
    rangex<- range(d1[,"x"],d2[,"x"])
    
  
    #Main plot
    plot(x=d1[,"x"],
         y=d1[,"y"],
         ...,
         col="red",
         ylim=c(0,1),
         xlim=rangex,
         type="l",
         xlab=paste0("Age of onset, centered by group ", centre_on),
         ylab="Cumulative density",
         sub=t_diff_message
    )
    #Overlay 2nd density curve
    lines(x=d2[,"x"],
          y=d2[,"y"],
          col="blue",
          lty="dashed")
    
    #Include legend
    legend("bottomright", legend = xlevels,
           col = c("red", "blue"), 
           lty = c("solid","dashed"), cex = 0.6)
  }
  #Add vline at y intercept
  abline(v=0,lty="dotted")
  
}