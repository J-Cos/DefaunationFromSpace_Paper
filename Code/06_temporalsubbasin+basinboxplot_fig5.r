library(terra)
library(tidyverse)
library(tidyterra)
library(ggsignif)


mk_cor<-rast("Outputs/KendallResults.tif")
basins<-vect("Outputs/Basins")
countries<-vect("Data/world-administrative-boundaries")
subbasins_congo<-crop(vect("Data/hybas_af_lev08_v1c"), basins[1])
subbasins_amazon<-crop(vect("Data/hybas_sa_lev08_v1c"), basins[2])


a<-crop(mk_cor, project(basins[2], crs(mk_cor)), mask=TRUE) %>% trim
c<-crop(mk_cor, project(basins[1], crs(mk_cor)), mask=TRUE) %>% trim



vals<-terra::extract( mk_cor, subbasins_congo, weights=TRUE, ID=TRUE) %>%
    filter(!is.na(tau)) %>%
    filter(weight>0.99) %>%
    group_by(ID) %>%
    summarise(mean=mean(tau), n=n()) %>%
    filter(n>10) #%>% 
    #mutate(mean=case_when(mean<0~0, .default =mean))
values(subbasins_congo)[ vals$ID, "mean"]<-vals$mean
subbasins_congo<-crop(subbasins_congo, ext(project(c, crs(subbasins_congo))))

fig5a<-ggplot() +
    geom_spatvector(data=subbasins_congo, aes(fill=mean), color=NA)+
    geom_spatvector(data=crop(countries, ext(subbasins_congo)), color="black", linewidth=1, fill=NA)+
    scale_fill_gradient2(limits = c(-1,1), na.value = "transparent", high=scales::muted("red"), low=scales::muted("blue"), mid="light grey", name="Change in FRIP\n(2020-2023)")+
    theme_classic()


vals<-terra::extract( mk_cor, subbasins_amazon, weights=TRUE, ID=TRUE) %>%
    filter(!is.na(tau)) %>%
    filter(weight>0.99) %>%
    group_by(ID) %>%
    summarise(mean=mean(tau), n=n()) %>%
    filter(n>10) #%>%
    #mutate(mean=case_when(mean<0~0, .default =mean))
values(subbasins_amazon)[ vals$ID, "mean"]<-vals$mean
subbasins_amazon<-crop(subbasins_amazon, ext(project(a, crs(subbasins_amazon))))

fig5b<-ggplot() +
    geom_spatvector(data=subbasins_amazon, aes(fill=mean), color=NA)+
    geom_spatvector(data=crop(countries, ext(subbasins_amazon)), color="black", linewidth=1, fill=NA)+
    scale_fill_gradient2(limits = c(-1,1), na.value = "transparent", high=scales::muted("red"), low=scales::muted("blue"), mid="light grey", name="Change in FRIP\n(2020-2023)")+
    theme_classic()


df<- rbind(
    data.frame("Mean"=subbasins_amazon$mean, "Basin"="Amazon"),
    data.frame("Mean"=subbasins_congo$mean, "Basin"="Congo")
    ) %>%
    as_tibble %>%
    dplyr::filter(!is.na(Mean))



lm(Mean~Basin, df) %>% summary

fig5c<-ggplot(data=df, aes(y=Basin, x=Mean, shape=Basin))+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.6), alpha=0.25, size=1)+
    geom_vline(xintercept=0, linetype=2, alpha=0.5)+
    geom_boxplot(outliers=FALSE, color="black", fill=NA)+
    geom_signif(
        y_position = c(0.32), xmin = c(1), xmax = c(2),
        annotation = c("***"), tip_length = 0, color="black"
    ) +
    theme_classic()+
    guides(shape="none")+
    xlab("Change in FRIP\n(2020-2023)")+ylab("")



#make combined figure
cowplot::plot_grid(fig5a, fig5b, fig5c, labels = c('A', 'B', "C"), label_size = 12, ncol=1, rel_heights = c(0.98, 1, 0.2))
ggsave("Figures/Figure5.png", height=14.3, width=10)
