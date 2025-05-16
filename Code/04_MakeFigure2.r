library(terra)
library(tidyverse)
library(tidyterra)

#functions
makeDataFigure<-function(plotDat, fillColumn){
    p<-ggplot() +
        geom_spatraster(data = plotDat, aes(fill=.data[[fillColumn]])) +
        geom_spatvector(data=crop(countries, ext(plotDat)), fill=NA, color="black") +
        theme_classic()
    if (fillColumn=="correlation") { p<- p + scale_fill_gradient2(na.value = "transparent", high=scales::muted("red"), low=scales::muted("blue"), mid="light grey", name="Flooding role\nin productivity", limits = c(-0.8,0.8)) }
    else if (fillColumn=="tau") {   p<-  p + scale_fill_gradient2(na.value = "transparent", high=scales::muted("orange"), low=scales::muted("green"), mid="light grey", name="Trend in\nflooding role\nin productivity\n(2000-2023)", limits = c(-0.8,0.8)) }
    else if (fillColumn=="DIunderestimate") {   p<-  p + scale_fill_gradient2(na.value = "transparent", high=scales::muted("yellow"), low=scales::muted("blue"), mid="light grey", name="Potential\ndefaunation\nmisestimation\n(obs.-pred.)", limits = c(-1.2,1.2)) }
    #else if (fillColumn=="DIunderestimate") {   p<-  p + viridis::  scale_fill_viridis(option="cividis", na.value = "transparent", name="Potential\nDefaunation\nMisestimation\n(Observation\n-\nPrediction)", limits = c(-1.2,1.2))}
    return(p)
}

FRIP_l<-lapply(list.files("Outputs/FRIP", full.names=TRUE), function(x){rast(x)$correlation})
FRIPtrend_l<-lapply(list.files("Outputs/FRIPtrend", full.names=TRUE), function(x){rast(x)})
basins<-vect("Outputs/Basins")
countries<-vect("Outputs/BasinCountries")
panel_Bogoni<-readRDS("Outputs/panel_Bogoni.RDS")


fig1a<-makeDataFigure(trim(crop(FRIP_l[[5]], basins[1])), "correlation")+guides(fill="none")
fig1b<-makeDataFigure((crop(FRIP_l[[5]], basins[2])), "correlation")
fig1c<-makeDataFigure(trim(crop(FRIPtrend_l[[5]], basins[1])), "tau")+guides(fill="none")
fig1d<-makeDataFigure((crop(FRIPtrend_l[[5]], basins[2])), "tau")


fig1e<-makeDataFigure(trim(crop(FRIP_l[[2]], basins[1])), "correlation")+guides(fill="none")
fig1f<-makeDataFigure((crop(FRIP_l[[2]], basins[2])), "correlation")
fig1g<-makeDataFigure(trim(crop(FRIPtrend_l[[2]], basins[1])), "tau")+guides(fill="none")
fig1h<-makeDataFigure((crop(FRIPtrend_l[[2]], basins[2])), "tau")



#make combined figures for state
toprow<-cowplot::plot_grid(fig1a, fig1b, labels = c('A', 'B'), label_size = 12, ncol=2, rel_widths = c(1, 1.22))
all<-cowplot::plot_grid(toprow, panel_Bogoni+geom_vline(xintercept=25, color="red", linetype=2), labels = c('', 'C'), label_size = 12, ncol=1)
ggsave("Figures/Figure2.png", height=8, width=12, bg="white")


toprow<-cowplot::plot_grid(fig1e, fig1f, labels = c('A', 'B'), label_size = 12, ncol=2, rel_widths = c(1, 1.07))
all<-cowplot::plot_grid(toprow, panel_Bogoni+geom_vline(xintercept=100, color="red", linetype=2), labels = c('', 'C'), label_size = 12, ncol=1)
ggsave("Figures/FigureS1.png", height=8, width=12, bg="white")


#make combined figures for trends
toprow<-cowplot::plot_grid(fig1c, fig1d, labels = c('A', 'B'), label_size = 12, ncol=2, rel_widths = c(1, 1.15))
ggsave("Figures/Figure5.png", height=5, width=16, bg="white")


toprow<-cowplot::plot_grid(fig1g, fig1h, labels = c('A', 'B'), label_size = 12, ncol=2, rel_widths = c(1, 1))
ggsave("Figures/FigureS7.png", height=5, width=16, bg="white")
