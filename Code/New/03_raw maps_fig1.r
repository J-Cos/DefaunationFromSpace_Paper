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
    else if (fillColumn=="tau") {   p<-  p + scale_fill_gradient2(na.value = "transparent", high=scales::muted("orange"), low=scales::muted("green"), mid="light grey", name="Trend in\nflooding\nrole in\nproductivity\n(2000-2023)", limits = c(-0.8,0.8)) }
    else if (fillColumn=="DIunderestimate") {   p<-  p + scale_fill_gradient2(na.value = "transparent", high=scales::muted("yellow"), low=scales::muted("blue"), mid="light grey", name="Potential\ndefaunation\nmisestimation\n(obs.-pred.)", limits = c(-1.2,1.2)) }
    #else if (fillColumn=="DIunderestimate") {   p<-  p + viridis::  scale_fill_viridis(option="cividis", na.value = "transparent", name="Potential\nDefaunation\nMisestimation\n(Observation\n-\nPrediction)", limits = c(-1.2,1.2))}
    return(p)
}

FRIP<-rast("Outputs/FRIP.tif") 
FRIPtrend<-rast("Outputs/FRIPtrend.tif")
DIunderestimate<-rast("Outputs/DIunderestimate.tif")
basins<-vect("Outputs/Basins")
countries<-vect("Outputs/BasinCountries")



fig1a<-makeDataFigure(trim(crop(FRIP, basins[1])), "correlation")+guides(fill="none")
fig1b<-makeDataFigure(crop(FRIP, basins[2]), "correlation")
fig1c<-makeDataFigure(aggregate(trim(crop(FRIP, basins[1])), 10, na.rm=TRUE), "correlation") +guides(fill="none")
fig1d<-makeDataFigure(aggregate(crop(FRIP, basins[2]), 10, na.rm=TRUE), "correlation")

fig1e<-makeDataFigure(trim(crop(FRIPtrend, basins[1])), "tau")+guides(fill="none")
fig1f<-makeDataFigure(crop(FRIPtrend, basins[2]), "tau")
fig1g<-makeDataFigure(aggregate(trim(crop(FRIPtrend, basins[1])), 10, na.rm=TRUE), "tau")+guides(fill="none")
fig1h<-makeDataFigure(aggregate(crop(FRIPtrend, basins[2]), 10, na.rm=TRUE), "tau")


fig_DIc<-makeDataFigure(trim(crop(DIunderestimate, basins[1])), "DIunderestimate")+guides(fill="none")
fig_DIa<-makeDataFigure(trim(crop(DIunderestimate, basins[2])), "DIunderestimate")



#make combined figure
#cowplot::plot_grid(fig1a, fig1b, fig1c, fig1d, fig1e, fig1f, fig1g, fig1h, labels = c('A', 'B', "C", "D", "E", "F", "G", "H"), label_size = 12, ncol=2)#, rel_heights = c(1.1, 1))

cowplot::plot_grid(fig1c, fig1d, fig_DIc, fig_DIa, fig1g, fig1h, labels = c('A', 'B', "C", "D", "E", "F"), label_size = 12, ncol=2, rel_widths = c(1, 1.3))

ggsave("Figures/Figure1.png", height=15, width=15)

cowplot::plot_grid(fig1a, fig1b, fig1e, fig1f, labels = c('A', 'B', "C", "D"), label_size = 12, ncol=2, rel_widths = c(1, 1.25))
ggsave("Figures/FigureS1.png", height=10, width=15)
