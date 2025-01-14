library(terra)
library(tidyverse)
library(tidyterra)


cor<-rast("Data/NppCorr_5km_includingNonsignif.tif") 
basins<-vect("Outputs/Basins")
countries<-vect("Data/world-administrative-boundaries")

a<-crop(cor, project(basins[2], crs(cor)), mask=TRUE) %>% trim
c<-crop(cor, project(basins[1], crs(cor)), mask=TRUE) %>% trim


fig1a<-ggplot() +
    geom_spatraster(data = a, aes(fill=correlation)) +
    geom_spatvector(data=crop(project(countries, crs(a)), ext(a)), fill=NA, color="black") +
    scale_fill_gradient2(na.value = "transparent", high=scales::muted("red"), low=scales::muted("blue"), mid="light grey", name="Flooding\nrole in\nproductivity")+
    theme_classic()

fig1b<-ggplot() +
    geom_spatraster(data = c, aes(fill=correlation)) +
    geom_spatvector(data=crop(project(countries, crs(c)), ext(c)), fill=NA, color="black") +
    scale_fill_gradient2(na.value = "transparent", high=scales::muted("red"), low=scales::muted("blue"), mid="light grey", name="Flooding\nrole in\nproductivity")+
    theme_classic()



#make combined figure
cowplot::plot_grid(fig1a, fig1b, labels = c('A', 'B'), label_size = 12, ncol=1, rel_heights = c(1.025, 1))
ggsave("Figures/Figure1.png", height=13.4, width=10)