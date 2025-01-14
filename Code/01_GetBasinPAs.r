library(terra)
library(tidyverse)


v<-vect(list(
    vect("Data/WDPA_Nov2024_Public_shp_0"),
    vect("Data/WDPA_Nov2024_Public_shp_1"),
    vect("Data/WDPA_Nov2024_Public_shp_2"))
)

af<-vect("Data/hybas_af_lev02_v1c")
congo<-af[(af$HYBAS_ID==1020018110)]
sa<-vect("Data/hybas_sa_lev02_v1c")
amazon<-sa[(sa$SUB_AREA==max(sa$SUB_AREA))]
basins<-vect(list(congo, amazon))

x <- crop(v, basins)
writeVector(basins, "Outputs/Basins", overwrite=TRUE)
writeVector(x, "Outputs/BasinPAs", overwrite=TRUE)