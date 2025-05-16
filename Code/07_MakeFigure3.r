library(terra)
library(tidyverse)
library(tidyterra)

# 1) load data --------------------------------------------------
FRIP<-lapply(list.files("Outputs/FRIP", full.names=TRUE), function(x){rast(x)$correlation})[[1]]
basins<-vect("Outputs/Basins")
countries<-vect("Outputs/BasinCountries")
PAs<-vect("Outputs/BasinPAs")


PAs_df<-readRDS("Outputs/PlottingData_PAs.RDS")
countries_df<-readRDS("Outputs/PlottingData_bcp.RDS")



# 2) prep data for plotting -----------------------------------------------
PAs_amazon<- crop(PAs, basins[2])
PAs_amazon_Focal<-PAs_amazon[expanse(PAs_amazon, unit="km")>40000]# & !PAs_reprojected$IUCN_CAT %in% c("Not Reported", "Not Applicable") ]
PAs_amazon_Focal<-PAs_amazon_Focal[PAs_amazon_Focal$NAME!="Central Amazon Conservation Complex"]

PAs_congo<- crop(PAs, basins[1])
PAs_congo_Focal<-PAs_congo[expanse(PAs_congo, unit="km")>20000]# & !PAs_reprojected$IUCN_CAT %in% c("Not Reported", "Not Applicable") ]
PAs_congo_Focal<-PAs_congo_Focal[PAs_congo_Focal$NAME!="Salonga"]

PAs_Focal<-vect(c(PAs_amazon_Focal, PAs_congo_Focal))




# 3) plot countries -----------------------------------------------------

df<-countries_df %>% filter(name %in% names(which(table(countries_df$name)>50)))

unique(df$Letters) %>% sort
colorRampPalette(c("grey", scales::muted("red")))(5)
group_by(df, Letters) %>% summarise(mean=mean(correlation)) %>% arrange(mean)
myPalette = c(
    'b' =   "#BEBEBE",
    'bc' = "#B4A4A4", 
    'abc' = "#AA8A8A", 
    'a' = "#A07171", 
    'x' = "#965757",
    'de' = "#8C3D3D", 
    'e' = "#832424")

fig_c1<-ggplot(data=df, aes(y=reorder(name, correlation, "median"), x=correlation, color=Letters))+
    geom_vline(xintercept=0, linetype=2, alpha=0.5)+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.6), aes(x=correlation), alpha=1, size=1)+
    geom_boxplot(aes(x=correlation), outliers=FALSE, fill=NA, color="black")+
    geom_text(data=group_by(df, name) %>% summarise(Letters=first(Letters), correlation=median(correlation)), aes(label=Letters), x=1.1)+
    xlim(-1.1, 1.2)+
    theme_classic()+
    guides(color="none")+
    #scale_fill_manual(values=myPalette)+
    scale_color_manual(values=myPalette)+
    xlab("Flooding role in productivity")+ylab("")

countriesForPlotting<-countries
values(countriesForPlotting)<-values(countries) %>%
    as_tibble %>%
    left_join(., summarise(group_by(df, treatment), Letters=first(Letters)), join_by(name==treatment))


plotDat<-crop(FRIP, basins[1])
fig_c2<-ggplot() +
    geom_spatvector(data=crop(countriesForPlotting, ext(plotDat)), aes(fill=Letters), linewidth=1, color="black")+
    scale_fill_manual(values=myPalette, na.value = "transparent")+
    guides(fill="none")+
    theme_classic()


plotDat<-crop(FRIP, basins[2])
fig_c3<-ggplot() +
    geom_spatvector(data=crop(countriesForPlotting, ext(plotDat)), aes(fill=Letters), linewidth=1, color="black")+
    scale_fill_manual(values=myPalette, na.value = "transparent")+
    guides(fill="none")+
    theme_classic()

# 4) plot PAS FRIP -------------------------------------------

df<-PAs_df %>% filter(NAME %in% names(which(table(PAs_df$NAME)>2)))  %>% 
  mutate(across('NAME', str_replace,  "Área De Proteção Ambiental Do Arquipélago Do Marajó", "Área De Proteção Ambiental\nDo Arquipélago Do Marajó"))


unique(df$Letters) %>% sort
colorRampPalette(c("grey", scales::muted("red")))(7)
group_by(df, Letters) %>% summarise(mean=mean(correlation)) %>% arrange(mean)
myPalette = c(
    'd' =   "#BEBEBE",
    'ad' = "#B4A4A4", 
    'abd' = "#AA8A8A", 
    'abcd' = "#A07171", 
    'bc' = "#965757",
    'c' = "#8C3D3D", 
    'abc' = "#832424")


fig_pa1<-ggplot(data=df, aes(y=reorder(NAME, correlation, "median"), x=correlation, color=Letters))+
    geom_vline(xintercept=0, linetype=2, alpha=0.5)+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.6), aes(x=correlation), alpha=1, size=1)+
    geom_boxplot(aes(x=correlation), outliers=FALSE, fill=NA, color="black")+
    geom_text(data=group_by(df, NAME) %>% summarise(Letters=first(Letters), correlation=median(correlation)), aes(label=Letters), x=1.1)+
    xlim(-1.1, 1.2)+
    theme_classic()+
    guides(color="none")+
    #scale_fill_manual(values=myPalette)+
    scale_color_manual(values=myPalette)+
    xlab("Flooding role in productivity")+ylab("")

PAs_wLetters<-PAs_Focal
values(PAs_wLetters)<-values(PAs_Focal) %>%
    as_tibble %>%
    left_join(., summarise(group_by(df, treatment), Letters=first(Letters)), join_by(NAME==treatment))

plotDat<-crop(FRIP, basins[1])
fig_pa2<-ggplot() +
    geom_spatvector(data=crop(PAs_wLetters[!is.na(PAs_wLetters$Letters)], ext(plotDat)), aes(fill=Letters), linewidth=1, color="black")+
    scale_fill_manual(values=myPalette, na.value = "transparent")+
    geom_spatvector(data=crop(countries, ext(plotDat)), fill=NA, linewidth=1, color="black")+
    guides(fill="none")+
    theme_classic()

plotDat<-crop(FRIP, basins[2])
fig_pa3<-ggplot() +
    geom_spatvector(data=crop(PAs_wLetters[!is.na(PAs_wLetters$Letters)], ext(plotDat)), aes(fill=Letters), linewidth=1, color="black")+
    scale_fill_manual(values=myPalette, na.value = "transparent")+
    geom_spatvector(data=crop(countries, ext(plotDat)), fill=NA, linewidth=1, color="black")+
    guides(fill="none")+
    theme_classic()

# 4) plot PAS FRIP trends-------------------------------------------

df<-PAs_df %>% filter(NAME %in% names(which(table(PAs_df$NAME)>2))) %>% 
  mutate(across('NAME', str_replace,  "Área De Proteção Ambiental Do Arquipélago Do Marajó", "Área De Proteção Ambiental\nDo Arquipélago Do Marajó"))


 
myPalette = c(
    'Decreasing' =   scales::muted("green"),
    'Increasing' = scales::muted("orange"), 
    'None' = "grey")


fig_trend1<-ggplot(data=df, aes(y=reorder(NAME, tau, "median"), x=tau, color=change))+
    geom_vline(xintercept=0, linetype=2, alpha=0.5)+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.6), aes(x=tau), alpha=1, size=1)+
    geom_boxplot(aes(x=tau), outliers=FALSE, fill=NA, color="black")+
    theme_classic()+
    guides(color="none")+
    #scale_fill_manual(values=myPalette)+
    scale_color_manual(values=myPalette)+
    xlab("Trend in flooding role in\nproductivity (2000-2023)")+ylab("")

PAs_wLetters<-PAs_Focal
values(PAs_wLetters)<-values(PAs_Focal) %>%
    as_tibble %>%
    left_join(., summarise(group_by(df, NAME), Letters=first(change)))

plotDat<-crop(FRIP, basins[1])
fig_trend2<-ggplot() +
    geom_spatvector(data=crop(PAs_wLetters[!is.na(PAs_wLetters$Letters)], ext(plotDat)), aes(fill=Letters), linewidth=1, color="black")+
    scale_fill_manual(values=myPalette, na.value = "transparent")+
    geom_spatvector(data=crop(countries, ext(plotDat)), fill=NA, linewidth=1, color="black")+
    guides(fill="none")+
    theme_classic()

plotDat<-crop(FRIP, basins[2])
fig_trend3<-ggplot() +
    geom_spatvector(data=crop(PAs_wLetters[!is.na(PAs_wLetters$Letters)], ext(plotDat)), aes(fill=Letters), linewidth=1, color="black")+
    scale_fill_manual(values=myPalette, na.value = "transparent")+
    geom_spatvector(data=crop(countries, ext(plotDat)), fill=NA, linewidth=1, color="black")+
    guides(fill="none")+
    theme_classic()

#
cowplot::plot_grid(
    fig_c2, fig_c3, fig_c1, fig_pa2, fig_pa3, fig_pa1, fig_trend2, fig_trend3, fig_trend1, 
    labels = c("A", "B",'C', "D", "E", "F", "G", "H", "I"), 
    label_size = 12, ncol=3, rel_widths = c(1, 1, 1.5))
ggsave("Figures/Figure4.png", height=10, width=14, bg="white")


