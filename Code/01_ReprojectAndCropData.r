library(terra)
library(tidyverse)

# 1) Load FRIP Data
rawFRIP_l<-lapply(list.files("Data/FRIP", full.names=TRUE), function(x){rast(x)$correlation})
rawFRIPchange_l<-lapply(list.files("Data/FRIP_Annual", full.names=TRUE), function(x){rast(x)})
preds_l<-lapply(list.files("Data/PredictorStack", full.names=TRUE), function(x){rast(x)})


# 2) Get Two Basins in FRIP projection
af<-vect("Data/hybas_af_lev02_v1c")
congo<-af[(af$HYBAS_ID==1020018110)]
sa<-vect("Data/hybas_sa_lev02_v1c")
amazon<-sa[(sa$SUB_AREA==max(sa$SUB_AREA))]
basins<-vect(list(congo, amazon)) %>%
    project(., crs(rawFRIP_l[[1]]))

# 3) Crop FRIP data to the basins and make a list
cropFRIP<-function(rawFRIP){
    list(
        crop(rawFRIP, basins[1], mask=TRUE) %>% trim,
        crop(rawFRIP, basins[2], mask=TRUE) %>% trim)
}
seperateFRIP_l <- lapply(rawFRIP_l, cropFRIP)
seperateFRIPchange_l <- lapply(rawFRIPchange_l, cropFRIP)


# 4) Crop countries to the available FRIP data and reproject
countries_congo<-vect("Data/world-administrative-boundaries") %>%
    project(., crs(rawFRIP_l[[1]])) %>%
    crop(seperateFRIP_l[[1]][[1]])
countries_amazon<-vect("Data/world-administrative-boundaries") %>%
    project(., crs(rawFRIP_l[[1]])) %>%
    crop(seperateFRIP_l[[1]][[2]])
countries<-vect(c(countries_congo, countries_amazon))

# 5) Crop PAs to the available FRIP data and reproject
rawPAs<-vect(list(
    vect("Data/WDPA_Nov2024_Public_shp_0"),
    vect("Data/WDPA_Nov2024_Public_shp_1"),
    vect("Data/WDPA_Nov2024_Public_shp_2"))
)
PAs_congo<-rawPAs %>%
    project(., crs(rawFRIP_l[[1]])) %>%
    crop(seperateFRIP_l[[1]][[1]])
PAs_amazon<-rawPAs%>%
    project(., crs(rawFRIP_l[[1]])) %>%
    crop(seperateFRIP_l[[1]][[2]])
PAs<-vect(c(PAs_congo, PAs_amazon))


# 6) Convert FRIP list into raster and crop FRIPchange
FRIP_l <- lapply(seperateFRIP_l, function(x){merge(x[[1]], x[[2]])})
FRIPchange_l <- lapply(seperateFRIPchange_l, function(x){merge(x[[1]], x[[2]])})

# 6) collate existing metrics and align with each FRIP resolution
addMetricsToFRIP<-function(FRIP){
    access<-rast("Data/accessibility.tif")  %>% crop(ext(FRIP))
    access_reprojected <- project(access, FRIP)
    ext(access_reprojected)<-ext(FRIP)
    names(access_reprojected)<-"access"


    DI<-rast("Data/DefInd.tif") %>% crop(ext(FRIP))
    crs(DI)<-"+proj=moll"
    DI_reprojected<-project(DI, FRIP)
    ext(DI_reprojected)<-ext(FRIP)
    names(DI_reprojected)<-"DI"

    DIlarge<-rast("Data/largeDefInd.tif") %>% crop(ext(FRIP))
    crs(DIlarge)<-"+proj=moll"
    DIlarge_reprojected<-project(DIlarge, FRIP)
    ext(DIlarge_reprojected)<-ext(FRIP)
    names(DIlarge_reprojected)<-"DIlarge"

    #get bogoni
    Bogoni<-rast("Data/Bogoni_DI.tif")
    Bogoni_reprojected<-project(Bogoni, FRIP)
    ext(Bogoni_reprojected)<-ext(FRIP)
    Bogoni_reprojected<-Bogoni_reprojected %>% crop(FRIP)

    names(Bogoni_reprojected)<-"DI_Bogoni"

    combinedRast<-c(FRIP, access_reprojected, DI_reprojected, DIlarge_reprojected, Bogoni_reprojected)
    combinedRast_masked<-mask(combinedRast, combinedRast$correlation)
    varnames(combinedRast_masked)<-varnames(FRIP)
    return(combinedRast_masked)
}
FRIPandMetrics_l<-lapply(FRIP_l, addMetricsToFRIP)


# 7) align predictorstack to FRIP
alignedPreds_l<-list()
for (i in 1:length(FRIP_l)){  alignedPreds_l[[i]]<-mask(crop(preds_l[[i]],  FRIP_l[[i]]),  FRIP_l[[i]])  }

# 8) Save all outputs
writeVector(basins, "Outputs/Basins", overwrite=TRUE)
writeVector(countries, "Outputs/BasinCountries", overwrite=TRUE)
writeVector(PAs, "Outputs/BasinPAs", overwrite=TRUE)

lapply(FRIPandMetrics_l, function(x){ writeRaster(x, paste0("Outputs/FRIP/", varnames(x)[1], ".tif"), overwrite=TRUE)  })
lapply(FRIPchange_l, function(x){ writeRaster(x, paste0("Outputs/FRIPchange/", varnames(x)[1], ".tif"), overwrite=TRUE)  })
lapply(alignedPreds_l, function(x){ writeRaster(x, paste0("Outputs/Predictors/", varnames(x)[1], ".tif"), overwrite=TRUE)  })

