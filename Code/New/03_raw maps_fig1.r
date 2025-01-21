library(terra)
library(tidyverse)
library(tidyterra)

#functions
makeDataFigure<-function(plotDat, fillColumn){
    p<-ggplot() +
        geom_spatraster(data = plotDat, aes(fill=.data[[fillColumn]])) +
        geom_spatvector(data=crop(countries, ext(plotDat)), fill=NA, color="black") +
        theme_classic()
    if (fillColumn=="correlation") { p<- p + scale_fill_gradient2(na.value = "transparent", high=scales::muted("red"), low=scales::muted("blue"), mid="light grey", name="Flooding\nrole in\nproductivity", limits = c(-1,1)) }
    else if (fillColumn=="tau") {   p<-  p + scale_fill_gradient2(na.value = "transparent", high=scales::muted("orange"), low=scales::muted("green"), mid="light grey", name="Trend in\nflooding\nrole in\nproductivity\n(2000-2023)", limits = c(-1,1)) }
    return(p)
}

FRIP<-rast("Outputs/FRIP.tif") 
FRIPtrend<-rast("Outputs/FRIPtrend.tif")
basins<-vect("Outputs/Basins")
countries<-vect("Outputs/BasinCountries")


fig1a<-makeDataFigure(trim(crop(FRIP, basins[1])), "correlation")
fig1b<-makeDataFigure(crop(FRIP, basins[2]), "correlation")
fig1c<-makeDataFigure(aggregate(trim(crop(FRIP, basins[1])), 10, na.rm=TRUE), "correlation")
fig1d<-makeDataFigure(aggregate(crop(FRIP, basins[2]), 10, na.rm=TRUE), "correlation")

fig1e<-makeDataFigure(trim(crop(FRIPtrend, basins[1])), "tau")
fig1f<-makeDataFigure(crop(FRIPtrend, basins[2]), "tau")
fig1g<-makeDataFigure(aggregate(trim(crop(FRIPtrend, basins[1])), 10, na.rm=TRUE), "tau")
fig1h<-makeDataFigure(aggregate(crop(FRIPtrend, basins[2]), 10, na.rm=TRUE), "tau")


#make combined figure
cowplot::plot_grid(fig1a, fig1b, fig1c, fig1d, fig1e, fig1f, fig1g, fig1h, labels = c('A', 'B', "C", "D", "E", "F", "G", "H"), label_size = 12, ncol=2)#, rel_heights = c(1.1, 1))
ggsave("Figures/Figure1.png", height=15, width=13)