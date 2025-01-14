library(terra)
library(tidyverse)
library(tidyterra)
library(ggsignif)
library(multcompView)


cor<-rast("Data/NppCorr_5km_includingNonsignif.tif") 
basins<-vect("Outputs/Basins")
countries<-vect("Data/world-administrative-boundaries")
PAs<-vect("Outputs/BasinPAs")

#basin raster
basin_rast<-as.factor(terra::rasterize(x=project(basins, crs(cor)), y=cor, "PFAF_ID"))

#country raster
countries_rast<-terra::rasterize(x=project(countries, crs(cor)), y=cor, "name")

#PA raster by protection binary
PAs$IUCN_number<-as.numeric(factor(PAs$IUCN_CAT, levels = c("Ia", "Ib", "II", "III", "IV", "V", "VI", "Not Applicable", "Not Assigned", "Not Reported")))
PA_status_rast<-terra::rasterize(x=project(PAs, crs(cor)), y=cor, field="IUCN_number", fun="min")
PA_status_rast[is.na(PA_status_rast)]<-11 # remove those with uncertain status

PA_binary_rast<-terra::rasterize(x=project(PAs, crs(cor)), y=cor, field="IUCN_number", fun="min")
PA_binary_rast$protection_binary<-as.factor(!is.na(PA_binary_rast$IUCN_number)) # protected = TRUE

#get different levels of protection
Protected<-(PA_status_rast$IUCN_number<4)*1
Unprotected<-(PA_status_rast$IUCN_number==11)*2
Uncertain<-(PA_status_rast$IUCN_number>3 & PA_status_rast$IUCN_number<11)*3

categories<-Protected+Unprotected+Uncertain
names(categories)<-"protectionCat"


PA_binary_rast$IUCN_number_factor<-as.factor(PA_status_rast$IUCN_number) 

PAs_merged<-aggregate(PAs) %>% disagg()
PAs_merged$size<-expanse(PAs_merged)
PAsize<-terra::rasterize(x=project(PAs_merged, crs(cor)), y=cor, field="size", fun="max")
# get pa by size



model_rast<-c(cor, basin_rast, countries_rast, PA_binary_rast$protection_binary, PA_binary_rast$IUCN_number_factor, categories, PAsize)
model_rast_agg <- aggregate(model_rast, fact=15, fun="modal", na.rm=TRUE)
model_rast_agg$correlation <- aggregate(model_rast$correlation, fact=15, fun="mean", na.rm=TRUE)

model_df<-as.data.frame(model_rast_agg, xy=TRUE, na.rm=TRUE) %>% as_tibble()

#model_df<-model_df[model_df$name %in% names(which(table(model_df$name)>90)),]

library(nlme)
subset<-1171

#filter out uncertain protection category

data.spatialCor.gls <- gls(correlation ~ size*PFAF_ID, data = model_df,
    method = "REML", verbose=TRUE, subset=1:subset)
data.spatialCor.glsExp <- gls(correlation ~ size*PFAF_ID, data = model_df,
    correlation = corExp(form = ~x + y, nugget = TRUE),
    method = "REML", verbose=TRUE, subset=1:subset)
#11.27

plot(residuals(data.spatialCor.gls, type = "normalized") ~
    fitted(data.spatialCor.gls))
plot(data.spatialCor.gls)
plot(nlme:::Variogram(data.spatialCor.gls, form = ~x +
    y, resType = "normalized"))


plot(residuals(data.spatialCor.glsExp, type = "normalized") ~
    fitted(data.spatialCor.glsExp))

AIC(data.spatialCor.gls, data.spatialCor.glsExp)

summary(data.spatialCor.gls)
summary(data.spatialCor.glsExp)