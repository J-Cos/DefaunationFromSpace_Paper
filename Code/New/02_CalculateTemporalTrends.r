library(terra)
library(tidyverse)
library(Kendall)

# functions
fun_kendall <- function(x){ return(unlist(Kendall::MannKendall(x)))
}

# calculate
FRIPchange<-rast("Outputs/FRIPchange.tif") 

mk_results<-terra::app(FRIPchange, fun_kendall)
trends<-mask(x=mk_results$tau, mask=mk_results$tau!=1, maskvalues=0)

writeRaster(trends, "Outputs/FRIPtrend.tif", overwrite=TRUE)



