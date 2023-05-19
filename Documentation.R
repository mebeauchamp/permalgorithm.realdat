##################################################################################################
##################################################################################################
###
### Documentation for the function 'permalgorithm.realDat'
###
### Code by: Marie-Eve Beauchamp                                            
### Last update: May 19, 2023                                              
###
##################################################################################################
##################################################################################################

#--------------------------------------------------------------------------------------------------
# permalgorithm.realDat  Generates a Dataset in Which Event Times Are Conditional On Time-Dependent 
#                        Covariates For Which Values Are Known Up To Different Follow-Up Times 
#                        Across Subjects, As In Most Real-World Datasets
#--------------------------------------------------------------------------------------------------

### Description ###

  # The 'permalgorithm.realDat' function generates a dataset in which event times are conditional       
  # on a user-specified list of covariates, some or all of which can be time-dependent. This
  # function was created to permit the use of covariate values taken from a real-world dataset 
  # with time-dependent covariate values known up to different time points for different subjects. 
  # Note that time-dependent covariate values from a real-world data are often unknown after the 
  # event or censoring time of a subject. In contrast, the 'permalgorithm' function in the 
  # 'PermAlgo' package requires that the matrix of covariate values (argument 'Xmat') passed to   
  # the function includes covariate values for each subject up to the maximum length of follow-up  
  # (argument 'maxTime'). 


### Usage ###

# permalgorithm.realDat(data, id, start, stop, covariates, eventTimes, betas)


### Arguments ###

  # data:   is a data frame in the counting process format where every line represents a time
  #         interval which can either correspond to 1 time unit or several time units during which 
  #         all covariate values for a given subject remain constant. The data frame must include 
  #         columns indicating: the subject identifier, the beginning and end of time intervals,     
  #         and the covariates on which the hazard of events depends. The length of follow-up can  
  #         vary across subjects. However, all subjects must start their follow-up at the same  
  #         time (i.e. delayed entry is not possible).
  # id:     is a character string representing the name of the column in 'data' identifying the 
  #         subjects, e.g. "id".
  # start:  is a character string representing the name of the column in 'data' identifying the   
  #         beginning of time intervals for the counting process format, e.g. "start".
  # stop:   is a character string representing the name of the column in 'data' identifying the end   
  #         of the time interval for the counting process format, e.g. "stop".
  # covariates: is a vector of character strings representing the names of covariates in 'data' on   
  #         which the hazard of events depends, e.g. c("dose", "sex", "age"). Covariates cannot be  
  #         factors; instead, users need to code them with binary indicators. 
  # eventTimes: is the vector of individual event times to be assigned. The length of 'eventTimes' 
  #         must be smaller or equal to the number of subjects in 'data'. Event times must be 
  #         smaller or equal to the maximum follow-up time of subjects.
  # betas:  is a vector of regression coefficients (log hazard) that represent the magnitude of the
  #         relationship between each covariate and the hazard of an event. The length of 'betas' 
  #         must correspond to the length of the vector 'covariates'.


### Details ###

  # The algorithm implemented in the function 'permalgorithm.realDat' is an adaptation of the 
  # permutational algorithm described in Sylvestre and Abrahamowicz (2008) and implemented in the
  # 'permalgorithm' function of the 'PermAlgo' package. This adaption allows using covariate  
  # values from a real-world dataset with time-dependent covariate values known up to different   
  # time points for different subjects. 
  
  # The current adaptation performs the matching of s event times with s vectors of covariates  
  # values. The number of events and individual event times are fixed by the user. They can either
  # be i) extracted from a real-word dataset (e.g., the same dataset as for the argument 'data'),     
  # or ii) generated from a user-defined distribution. The number of events (s) must be smaller or  
  # equal to the number of subjects (n) in 'data'.
  
  # Event times are assigned in increasing time order. Each event time is randomly matched, with 
  # weighted sampling, to one vector of covariate values (corresponding to one subject) among  
  # subjects available in the corresponding risk set. The risk set for an event time is defined as
  # the set of subjects: a) for whom the covariate values are known up to this event time, and 
  # b) who have not been selected yet for any earlier event. 
  
  # If there are no subjects available in the risk set of an event time to be assigned, then the 
  # function stops and returns an error message. This situation may occur if the number of events  
  # is too large with respect to the known follow-up of subjects (person-time in the cohort), or   
  # it could also depend on the distribution of event times to be assigned. 
  
  # After the s event times from the argument 'eventTimes' are assigned to s subjects, the other     
  # n-s subjects in 'data' are censored at the end of their known follow-up in 'data'. Note that     
  # this differs from the 'permalgorithm' function, in which n observed times (event or censoring   
  # times) are randomly assigned to the n subjects in the dataset generated. 
  
  # Factor variables are not allowed in 'data'. Instead, users need to code them with binary 
  # indicators.


### Value ###

  # The function 'permalgorithm.realDat' returns a data frame including the same columns as in   
  # 'data' and with the following two additional columns:
  # Event.NEW: indicator of generated event. 'Event.NEW'=1 when an event occurs and 0 otherwise.
  # Fup.NEW:   individual follow-up time for the corresponding subject in the generated data.
  
  # The columns indicating the beginning and end of time intervals in the returned data frame have  
  # the same names as in 'data' (i.e. arguments 'start' and 'stop'). However, the end of time   
  # intervals are adjusted to event times assigned. All other columns in the 'data' argument 
  # remain unchanged in the returned data frame. Note that if the column with event status from 
  # the original dataset is included in 'data', then this column will still be included in the 
  # data output. Therefore, to avoid confusion with 'Event.NEW', it may be advisable to exclude 
  # the original event status from 'data'.
  
  # The number of subjects in the generated dataset is identical to 'data'.


### References ###

  # This algorithm is a variation of the permutational algorithm described and validated in  
  # Sylvestre and Abrahamowicz (2008), and originally presented by Abrahamowicz, MacKenzie and 
  # Esdaile (1996) and MacKenzie and Abrahamowicz (2002).
  
  # Please reference the manuscript by Sylvestre and Abrahamowicz (2008), cited below, if this 
  # function is used in any published material.
  
  # Sylvestre M.-P., Abrahamowicz M. (2008) Comparison of algorithms to generate event times 
  # conditional on time-dependent covariates. Statistics in Medicine 27(14):2618–34.
  # 
  # Abrahamowicz M., MacKenzie T., Esdaile J.M. (1996) Time-dependent hazard ratio: modelling and
  # hypothesis testing with application in lupus nephritis. JASA 91:1432–9.
  # 
  # MacKenzie T., Abrahamowicz M. (2002) Marginal and hazard ratio specific random data generation:
  # Applications to semi-parametric bootstrapping. Statistics and Computing 12(3):245–252.


### Examples ###

# Set the path to the current folder
setwd("C://path//to//your//directory")

# Source the code for the function 'permalgorithm.realDat'
source("permalgorithm.RealDat.R")

# Load the dataset (called 'data_original') used for the examples below
load("data.RData")
head(data_original)

  # Description of variables in 'data_original':
    # id: patient identifier
    # start: start time of the interval
    # stop: stop time of the interval
    # event: event status 
    # fup: follow-up time for the current patient
    # age: patient's age
    # sex: patient's sex
    # dose: time-dependent value of the current dose of a drug

# Maximum follow-up time varies across subjects, from 26 to 366 days
summary(by(data_original$stop, data_original$id, max))

## EXAMPLE 1: Generate a new dataset using the same event times as in the original dataset

# Extract event times from the original datasets
evt_original <- data_original$stop[data_original$event == 1]

data_new <- permalgorithm.realDat(data = data_original, id = "id", start = "start", stop = "stop", 
                                  covariates = c("age", "sex", "dose"), eventTimes = evt_original, 
                                  betas = log(c(1.04, 1.21, 1.10)))

head(data_new)
  # Note: the columns 'event' and 'fup' correspond to the original dataset and are irrelevant for 
  # the new dataset generated. Event status and follow-up for the generated dataset are given by 
  # the columns 'Event.NEW' and 'Fup.NEW'.

## EXAMPLE 2: Generated a dataset using newly generated event times

# Identify the cohort entry time ('start'=1 for all subjects) and maximum follow-up time ('stop'=366)
# in the original dataset. Therefore, the generated event times must be after 1 and up to 366. 
summary(by(data_original$start, data_original$id, min))
summary(by(data_original$stop, data_original$id, max))

# Generate 100 event times during the study follow-up
evt_new <- ceiling(runif(n = 100, min = 1, max = 366))

data_new2 <- permalgorithm.realDat(data = data_original, id = "id", start = "start", stop = "stop",
                                   covariates = c("age", "sex", "dose"), eventTimes = evt_new,
                                   betas = log(c(1.04, 1.21, 1.10)))
