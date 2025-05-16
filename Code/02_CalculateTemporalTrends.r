library(terra)
library(tidyverse)
library(Kendall)

# functions
fun_kendall <- function(x){ return(unlist(Kendall::MannKendall(x)))
}

# calculate
FRIPchange_l<-lapply(list.files("Outputs/FRIPchange", full.names=TRUE), function(x){rast(x)})

trends_l<-lapply(
    FRIPchange_l, 
    function(FRIPchange){
        mk_results<-terra::app(FRIPchange, fun_kendall)
        trends<-mask(x=mk_results$tau, mask=mk_results$tau!=1, maskvalues=0)
        varnames(trends)<-varnames(FRIPchange)
        return(trends)
    }
)

lapply(trends_l, function(x){ writeRaster(x, paste0("Outputs/FRIPtrend/", varnames(x)[1], ".tif"), overwrite=TRUE)  })



