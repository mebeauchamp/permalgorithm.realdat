##################################################################################################
##################################################################################################
###
### Code for the function 'permalgorithm.realDat'
###
### Code by: Marie-Eve Beauchamp                                            
### Last update: May 18, 2023                                              
###
##################################################################################################
##################################################################################################

###################################################################################################
## CODE FOR THE MAIN FUNCTION: permalgorithm.realDat
###################################################################################################

permalgorithm.realDat <- function(data, id, start, stop, covariates, eventTimes, betas){
  
  ### Verification that inputs are appropriate ###
  
  if (length(eventTimes) > length(unique(data[, id]))){
    stop("The number of events to assign must be smaller or equal to the number of subjects.")
  }
  if (max(eventTimes) > max(data[, stop])){
    stop("Some values in eventTimes are higher than the maximum follow-up time of subjects.") 
  }
  if (min(eventTimes) <= max(by(data[, start], data[, id], min))){
    stop("Some values in eventTimes are earlier than the beginning of follow-up of some subjects.") 
  }

  eventTimes <- sort(eventTimes)
  ids <- unique(data[, id])
  
  ### Calculate the numerator of sampling probability for subjects in risk set of each event time ###
    # numRS: 1 line per event time (increasing time order), 1 column per subject in 'data'.
    # NA assigned to subjects not in the risk set of an event time.

  if (max(data[, stop] - data[, start]) == 1){
    numRS <- do.call('cbind', by(data, data[, id], 
                                 numRS.1linePerUnit.fct, stop, covariates, eventTimes, betas, 
                                 simplify = FALSE))
  } else {
    numRS <- do.call('cbind', by(data, data[, id], 
                                 numRS.intervals.fct, start, stop, covariates, eventTimes, betas, 
                                 simplify = FALSE))
  }
  colnames(numRS) <- paste0('id.', ids)

  ### Assign event times ###
  
  ev.assigned.ids <- rep(NA, times = length(eventTimes)) # Store id being assigned an event
  
  # For each event time, in increasing time order
  for (i in 1:length(eventTimes)){

    if (sum(!is.na(numRS[i, ])) == 0){
      stop("Events assigned (in increasing time order) up to the ", i-1, "th out of ", 
      length(eventTimes), " events. 
  Cannot assign the remaining events because no subjects are available in their risk sets.")
    }
    
    # Sampling probabilities of being assigned the current event time (see equation (2) in 
    # Sylvestre and Abrahamowicz (2008))
    probPA <- numRS[i, ] / sum(numRS[i, ], na.rm = TRUE) 
    
    # Randomly sample 1 subject available in current risk set 
    if (sum(!is.na(numRS[i, ])) == 1){    
      ev.assigned.ids[i] <- ids[!is.na(probPA)]
    } else {
      ev.assigned.ids[i] <- sample(ids[!is.na(probPA)], size = 1, prob = probPA[!is.na(probPA)])
    }
    
    # Remove the subject selected from risk sets for event times to be subsequently assigned
    if (i < length(eventTimes)){
      numRS[(i+1):length(eventTimes), which(ids == ev.assigned.ids[i])] <- NA
    }
  }

  ### Construct the generated dataset ###
  
  # Notes: Subjects assigned an event finish their follow-up at the assigned event times.
  #        For those not assigned an event, their follow-up ends at same time as in 'data'.
  
  data.GEN <- do.call('rbind', by(data, data[, id], cutFUPev.fct, id, start, stop, eventTimes, 
                                  ev.assigned.ids))
  return(data.GEN)
}


###################################################################################################
## CODE FOR INTERNAL FUNCTIONS 
###################################################################################################

#-----------------------------------------------------------------------------------------------
# Function evaluating if each item in the vector 'p' is included in the interval 'int' delimited  
# by int[1] and int[2] 
#-----------------------------------------------------------------------------------------------

in.interval.fct <- function(int, p){
  p > int[1] & p <= int[2]
}

#------------------------------------------------------------------------------------------------
# Functions to calculate the numerators of sampling probabilities for subjects in the risk set of   
# each event time (NA assigned to subjects not in a risk set)  
#------------------------------------------------------------------------------------------------

### Version of the function for 'data' including time intervals longer than 1 time unit 
  # (function generalized to follow-up starting at other value than 0)

numRS.intervals.fct <- function(curSubj, start, stop, covariates, eventTimes, betas){
  # Arguments:
    # curSubj: current subject on whom the function is applied
  
  # Indicates if (TRUE/FALSE) each event time is included in follow-up of 'curSubj' 
  InRS <- eventTimes <= max(curSubj[, stop])
    # Note: the line above will need to be changed if delayed entry or gaps time intervals are allowed
  
  # Identifies the line in 'curSubj' for each event time during follow-up of 'curSubj' 
  tmp <- t(apply(curSubj[, c(start, stop)],  1, in.interval.fct, p = eventTimes[InRS]))
  tmp2 <- tmp * 1 * 1:nrow(curSubj)
  tmp2[tmp2 == 0] <- NA
  evLineNo.curSubj <- tmp2[!is.na(tmp2)]
  
  # Numerators of sampling probabilities for each event time in 'eventTimes'; NA indicated if 
  # 'curSubj' is not in the risk set for that event time  
  num <- rep(NA, times = length(InRS))
  num[InRS] <- exp(as.matrix(curSubj[evLineNo.curSubj, covariates]) %*%
                     as.matrix(betas, nrow = length(betas)))
    # Note: ties in event times are repeated above
  return(num)
}

### Version of the function for 'data' with the format of 1 line per time unit 
  # (function generalized to follow-up starting at other value than 0)

numRS.1linePerUnit.fct <- function(curSubj, stop, covariates, eventTimes, betas){
  # Arguments:
    # curSubj: current subject on whom the function is applied
  
  # Indicates if (TRUE/FALSE) each event time is included in follow-up of 'curSubj'
  InRS <- eventTimes <= max(curSubj[, stop])
    # Note: the line above will need to be changed if delayed entry or gaps time intervals are allowed
  
  # Numerators of sampling probabilities for each event time in 'eventTimes'; NA indicated if 
  # 'curSubj' is not in the risk set for that event time 
  num <- rep(NA, times = length(InRS))
  num[InRS] <- exp(as.matrix(curSubj[eventTimes[InRS] - (curSubj[1, stop] - 1), covariates]) %*%
                     as.matrix(betas, nrow = length(betas)))
    # Bote: "+ curSubj[1, stop] - 1" accounts for the fact that 1st stop value can be different from 0.
    # Note: ties in event times are repeated above.
  return(num)
}

#----------------------------------------------------------
# Function to cut the follow-up at the assigned event times
#----------------------------------------------------------

cutFUPev.fct <- function(curSubj, id, start, stop, eventTimes, ev.assigned.ids){
  # Arguments:
    # curSubj: current subject on whom the function is applied
    # ev.assigned.ids: id of all subjects being assigned an event
  
  # If 'curSubj' is assigned an event
  if (curSubj[1, id] %in% ev.assigned.ids){
    
    # Identifies the line of 'curSujb' in which falls the assigned event time
    tmp <- t(apply(curSubj[, c(start, stop)],  1, 
                   in.interval.fct, p = eventTimes[curSubj[1, id] == ev.assigned.ids]))
    evLine <- max(tmp * 1 * 1:nrow(curSubj))
    
    curSubj <- curSubj[1:evLine, ]
    curSubj[nrow(curSubj), stop] <- eventTimes[curSubj[1, id] == ev.assigned.ids]
    curSubj[, 'Event.NEW'] <- c(rep(0, times=(evLine - 1)), 1)
    curSubj[, 'Fup.NEW'] <- curSubj[nrow(curSubj), stop] 
    return(curSubj)
    
  # If 'curSubj' is not assigned an event
  } else {
    curSubj[, 'Event.NEW'] <- 0
    curSubj[, 'Fup.NEW'] <- curSubj[nrow(curSubj), stop]
    return(curSubj) 
  }
}
