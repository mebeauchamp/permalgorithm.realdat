The R function 'permalgorithm.realDat' generates a dataset in which event times are conditional on a user-specified list of covariates, some or all of which can be time-dependent. This function was created to permit the use of covariate values taken from a real-world dataset with time-dependent covariate values known up to different time points for different subjects. Note that time-dependent covariate values from a real-world data are often unknown after the event or censoring time of a subject. In contrast, the 'permalgorithm' function in the R package 'PermAlgo' requires that the matrix of covariate values (argument 'Xmat') passed to the function includes covariate values for each subject up to the maximum length of follow-up (argument 'maxTime'). 

For questions or comments about the code contact Marie-Eve Beauchamp (marie-eve.beauchamp at rimuhc.ca).

## Content

#### `permalgorithm.RealDat.R`
Code of the function 'permalgorithm.realDat' and of the internal functions.

#### `Documentation.R`
Documentation for the function 'permalgorithm.realDat', including examples of its use.

#### `data.RData`
File including the two datasets used in the examples presented in `Documentation.R`. 
