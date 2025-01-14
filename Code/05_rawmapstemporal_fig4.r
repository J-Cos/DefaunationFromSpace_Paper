library(terra)
library(tidyverse)
library(tidyterra)
library(Kendall)

# functions
fun_kendall <- function(x){ return(unlist(Kendall::MannKendall(x)))
}

cor<-rast("Data/NppCorr_5km_Annual.tif") 
basins<-vect("Outputs/Basins")
countries<-vect("Data/world-administrative-boundaries")

mk_cor<-terra::app(cor, fun_kendall)
mk_cor<- mask(x=mk_cor, mask=mk_cor$tau==1, maskvalues=1)
writeRaster(mk_cor, "Outputs/KendallResults.tif")


trends<-mask(x=mk_cor$tau, mask=mk_cor$sl<0.05, maskvalues=0)


a<-crop(trends, project(basins[2], crs(cor)), mask=TRUE) %>% trim
c<-crop(trends, project(basins[1], crs(cor)), mask=TRUE) %>% trim


fig4a<-ggplot() +
    geom_spatraster(data = a, aes(fill=tau)) +
    geom_spatvector(data=crop(project(countries, crs(a)), ext(a)), fill=NA, color="black") +
    scale_fill_gradient2(limits = c(-1,1), na.value = "transparent", high=scales::muted("red"), low=scales::muted("blue"), mid="light grey", name="Change in FRIP\n(2020-2023)")+
    theme_classic()

fig4b<-ggplot() +
    geom_spatraster(data = c, aes(fill=tau)) +
    geom_spatvector(data=crop(project(countries, crs(c)), ext(c)), fill=NA, color="black") +
    scale_fill_gradient2(limits = c(-1,1), na.value = "transparent", high=scales::muted("red"), low=scales::muted("blue"), mid="light grey", name="Change in FRIP\n(2020-2023)")+
    theme_classic()



#make combined figure
cowplot::plot_grid(fig4a, fig4b, labels = c('A', 'B'), label_size = 12, ncol=1, rel_heights = c(1.025, 1))
ggsave("Figures/Figure4.png", height=13.4, width=10)



