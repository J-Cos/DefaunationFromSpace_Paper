library(terra)
library(tidyverse)

# 1) Load FRIP Data
rawFRIP<-rast("Data/NppCorr_5km_includingNonsignif.tif")$correlation
rawFRIPchange<-rast("Data/NppCorr_5km_Annual.tif") 

# 2) Get Two Basins in FRIP projection
af<-vect("Data/hybas_af_lev02_v1c")
congo<-af[(af$HYBAS_ID==1020018110)]
sa<-vect("Data/hybas_sa_lev02_v1c")
amazon<-sa[(sa$SUB_AREA==max(sa$SUB_AREA))]
basins<-vect(list(congo, amazon)) %>%
    project(., crs(rawFRIP))

# 3) Crop FRIP data to the basins and make a list
FRIP_l <-list(
    crop(rawFRIP, basins[1], mask=TRUE) %>% trim,
    crop(rawFRIP, basins[2], mask=TRUE) %>% trim)

# 4) Crop countries to the available FRIP data and reproject
countries_congo<-vect("Data/world-administrative-boundaries") %>%
    project(., crs(FRIP)) %>%
    crop(FRIP_l[[1]])
countries_amazon<-vect("Data/world-administrative-boundaries") %>%
    project(., crs(FRIP)) %>%
    crop(FRIP_l[[2]])
countries<-vect(c(countries_congo, countries_amazon))

# 5) Crop PAs to the available FRIP data and reproject
rawPAs<-vect(list(
    vect("Data/WDPA_Nov2024_Public_shp_0"),
    vect("Data/WDPA_Nov2024_Public_shp_1"),
    vect("Data/WDPA_Nov2024_Public_shp_2"))
)
PAs_congo<-rawPAs %>%
    project(., crs(FRIP)) %>%
    crop(FRIP_l[[1]])
PAs_amazon<-rawPAs%>%
    project(., crs(FRIP)) %>%
    crop(FRIP_l[[2]])
PAs<-vect(c(PAs_congo, PAs_amazon))

# 6) Convert FRIP list into raster and crop FRIPchange
FRIPchange <-merge(
    crop(rawFRIPchange, basins[1], mask=TRUE) %>% trim,
    crop(rawFRIPchange, basins[2], mask=TRUE) %>% trim)
FRIP <- merge(FRIP_l[[1]], FRIP_l[[2]])

# 7) Save all outputs
writeVector(basins, "Outputs/Basins", overwrite=TRUE)
writeVector(countries, "Outputs/BasinCountries", overwrite=TRUE)
writeVector(PAs, "Outputs/BasinPAs", overwrite=TRUE)

writeRaster(FRIP, "Outputs/FRIP.tif", overwrite=TRUE)
writeRaster(FRIPchange, "Outputs/FRIPchange.tif", overwrite=TRUE)

