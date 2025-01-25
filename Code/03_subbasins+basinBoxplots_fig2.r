library(terra)
library(tidyverse)
library(tidyterra)
library(ggsignif)


cor<-rast("Data/NppCorr_5km_includingNonsignif.tif") 
basins<-vect("Outputs/Basins")
countries<-vect("Data/world-administrative-boundaries")
subbasins_congo<-crop(vect("Data/hybas_af_lev08_v1c"), basins[1])
subbasins_amazon<-crop(vect("Data/hybas_sa_lev08_v1c"), basins[2])


a<-crop(cor, project(basins[2], crs(cor)), mask=TRUE) %>% trim
c<-crop(cor, project(basins[1], crs(cor)), mask=TRUE) %>% trim



vals<-terra::extract( cor, subbasins_congo, weights=TRUE, ID=TRUE) %>%
    filter(!is.na(correlation)) %>%
    filter(weight>0.99) %>%
    group_by(ID) %>%
    summarise(mean=mean(correlation), n=n()) %>%
    filter(n>5) #%>%
    #mutate(mean=case_when(mean<0~0, .default =mean))
values(subbasins_congo)[ vals$ID, "mean"]<-vals$mean
subbasins_congo<-crop(subbasins_congo, ext(project(c, crs(subbasins_congo))))

fig2a<-ggplot() +
    geom_spatvector(data=subbasins_congo, aes(fill=mean), color=NA)+
    geom_spatvector(data=crop(countries, ext(subbasins_congo)), color="black", linewidth=1, fill=NA)+
    scale_fill_gradient2(limits = c(-1,1), na.value = "transparent", high=scales::muted("red"), low=scales::muted("blue"), mid="light grey", name="Flooding\nrole in\nproductivity")+
    theme_classic()


vals<-terra::extract( cor, subbasins_amazon, weights=TRUE, ID=TRUE) %>%
    filter(!is.na(correlation)) %>%
    filter(weight>0.99) %>%
    group_by(ID) %>%
    summarise(mean=mean(correlation), n=n()) %>%
    filter(n>5) #%>%
    #mutate(mean=case_when(mean<0~0, .default =mean))
values(subbasins_amazon)[ vals$ID, "mean"]<-vals$mean
subbasins_amazon<-crop(subbasins_amazon, ext(project(a, crs(subbasins_amazon))))

fig2b<-ggplot() +
    geom_spatvector(data=subbasins_amazon, aes(fill=mean), color=NA)+
    geom_spatvector(data=crop(countries, ext(subbasins_amazon)), color="black", linewidth=1, fill=NA)+
    scale_fill_gradient2(limits = c(-1,1), na.value = "transparent", high=scales::muted("red"), low=scales::muted("blue"), mid="light grey", name="Flooding\nrole in\nproductivity")+
    theme_classic()


df<- rbind(
    data.frame("Mean"=subbasins_amazon$mean, "Basin"="Amazon"),
    data.frame("Mean"=subbasins_congo$mean, "Basin"="Congo")
    ) %>%
    as_tibble %>%
    dplyr::filter(!is.na(Mean))



lm(Mean~Basin, df) %>% summary

fig2c<-ggplot(data=df, aes(y=Basin, x=Mean, shape=Basin))+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.6), alpha=0.25, size=1)+
    geom_vline(xintercept=0, linetype=2, alpha=0.5)+
    geom_boxplot(outliers=FALSE, color="black", fill=NA)+
    geom_signif(
        y_position = c(0.8), xmin = c(1), xmax = c(2),
        annotation = c("***"), tip_length = 0, color="black"
    ) +
    theme_classic()+
    guides(shape="none")+
    xlab("Flooding role in productivity")+ylab("")



#make combined figure
cowplot::plot_grid(fig2a, fig2b, fig2c, labels = c('A', 'B', "C"), label_size = 12, ncol=1, rel_heights = c(0.98, 1, 0.2))
ggsave("Figures/Figure2.png", height=14.3, width=10)
