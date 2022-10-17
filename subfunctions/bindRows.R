#####
# Specify bindRows function, which generates a matrix by combining a list of vectors or matrices with varying length/ncol;
# duplicate the final value/column of shorter list elements until max length/ncol is reached for all elements and rows can be bound
# Written to prepare simulated data for use in validating the ADPenetrance (https://github.com/ThomasPSpargo/adpenetrance/) approach to calculate genetic penetrance
# Author: Thomas Spargo (thomas.spargo@kcl.ac.uk) 
# Please get in touch if you have any issues.
####

## INPUT:
#x              = A list of vectors or matrices
#isMatrix       = Logical defaults to FALSE. Specify TRUE if x is a list of matrices. If FALSE, assumes a list of vectors as input
#silent         = Logical, defaults to FALSE. Specify TRUE to silence warning about dropping of NULL list elements

## OUTPUT:
# A single matrix combining all elements of an input list rowwise. The final element (if vector) or column (if matrix) of shorter list elements are duplicated until all elements have equal length/ncol and rows can be bound.

bindRows  <- function(x,isMatrix=FALSE,silent=FALSE){
  
  #Drop any null elements of the list
  checknull<- sapply(x,is.null)
  if(any(checknull)){
    x <- x[!checknull]
    if(silent==FALSE){
    warning("Certain elements of the list passed to bindRows were null. Elements ",paste0(which(checknull==TRUE),collapse=", ")," have been dropped")
    }
  }
  #Default method, assuming list elements are vectors
  if(isMatrix==FALSE){
    all_lengths<- unique(sapply(x,length)) #Check all elements have the same length
    if(length(all_lengths)>1){ #If not all datasets are the right length, duplicate the final element of the short ones
      maximum<- max(all_lengths)
      for(j in 1:length(x)){
        len<- length(x[[j]])
        if(len<maximum){
          times <-maximum-len
          x[[j]] <- c(x[[j]],rep(x[[j]][len],times))
        }
      }
    }
    x <- matrix(unlist(x), nrow=length(x), byrow=TRUE) #Unlist into matrix
  } else if(isMatrix==TRUE){
    all_lengths<- unique(sapply(x,ncol)) #Check all elements have the same length
    if(length(all_lengths)>1){ #If not all datasets are the right length, duplicate the final element of the short ones
      maximum<- max(all_lengths)
      for(j in 1:length(x)){
        col<- ncol(x[[j]])
        if(col<maximum){
          times <-maximum-col
          x[[j]] <- cbind(x[[j]],x[[j]][,rep(col,each=times)])                
        }
      }
    }
    #Recursively bind together list elements
    x_mat <- NULL
    for(j in 1:length(x)){x_mat <- rbind(x_mat,x[[j]])}
    x <- x_mat
  }
  
  return(x)
}

############ END FUNCTION SPECIFICATION ############

# #Quick functionality test
# 
# #Generate var.Char
# f=0.75
# g=0.00
# numsteps=10
# onsetRateDiff=1
# var.Char<- varChars(f,g,numsteps,onsetRateDiff) #Demo varchars settings
# 
# 
# #Short lapply genFamily for two example populations
# final_time=FALSE
# stepHazard=FALSE
# eldestAt0=FALSE
# return_indivs=FALSE
# 
# aTime <- lapply(0:5,genFamily,group="var",var.Char=var.Char,eldestAt0=eldestAt0,final_time=final_time,stepHazard=stepHazard,return_indivs=return_indivs)
# 
# bTime <- lapply(0:5,genFamily,group="novar",var.Char=var.Char,eldestAt0=eldestAt0,final_time=final_time,stepHazard=stepHazard,return_indivs=return_indivs)
# 
# #Combine list elements for each of aTime and bTime with bindRows (vector combining)
# #Use default cbind function to append a first column flagging the different populations
# aBound<- cbind("var",bindRows(aTime))
# bBound<- cbind("novar",bindRows(bTime))
# 
# #Combine aBound and bBound with bindRows (matrix combining)
# allBound<- bindRows(list(aBound,bBound),isMatrix = TRUE)


