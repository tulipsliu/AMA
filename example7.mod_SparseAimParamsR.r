setsigmaGlobally <- function(newValue){
sigma<<- newValue;
} 

setsigmaLocally <- function(paramVector, newValue){
#Requires paramVector to have names
paramVector['sigma'] <<- newValue 
 }

setdeltaGlobally <- function(newValue){
delta<<- newValue;
} 

setdeltaLocally <- function(paramVector, newValue){
#Requires paramVector to have names
paramVector['delta'] <<- newValue 
 }

setlambdaGlobally <- function(newValue){
lambda<<- newValue;
} 

setlambdaLocally <- function(paramVector, newValue){
#Requires paramVector to have names
paramVector['lambda'] <<- newValue 
 }

setrhoGlobally <- function(newValue){
rho<<- newValue;
} 

setrhoLocally <- function(paramVector, newValue){
#Requires paramVector to have names
paramVector['rho'] <<- newValue 
 }

setgampiGlobally <- function(newValue){
gampi<<- newValue;
} 

setgampiLocally <- function(paramVector, newValue){
#Requires paramVector to have names
paramVector['gampi'] <<- newValue 
 }

setAllParamsGlobally <- function(newParamVector){
	 sigma<<- paramVector1
	 delta<<- paramVector2
	 lambda<<- paramVector3
	 rho<<- paramVector4
	 gampi<<- paramVector5
} 

updateParamVector <- function(oldParamVector, newParamVector){
	 oldParamVector <<- newParamVector 
 }
