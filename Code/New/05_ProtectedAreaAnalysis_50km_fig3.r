# explore negret et al 2020 method

library(terra)
library(tidyverse)
library(multcompView)
library(tidyterra)

generate_label_df <- function(TUKEY, variable){
 
     # Extract labels and factor levels from Tukey post-hoc 
     Tukey.levels <- TUKEY[[variable]][,4] 
     Tukey.levels<-Tukey.levels[!is.na(Tukey.levels)]
     Tukey.labels <- data.frame(multcompView::multcompLetters(x=Tukey.levels)['Letters'])
     
     #I need to put the labels in the same order as in the boxplot :
     Tukey.labels$treatment=rownames(Tukey.labels)
     Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
     return(Tukey.labels)
     }
 
# 0) set parameters
aggregationFactor<-10

# 1) load data 
FRIP<-rast("Outputs/FRIP.tif") 
basins<-vect("Outputs/Basins")
countries<-vect("Outputs/BasinCountries")
PAs<-vect("Outputs/BasinPAs")

# 2) make single raster
# aggregate FRIP
FRIPagg<-aggregate(FRIP, aggregationFactor, na.rm=TRUE)

#basin raster
PAs_amazon<- crop(PAs, basins[2])
PAs_amazon_Focal<-PAs_amazon[expanse(PAs_amazon, unit="km")>40000]# & !PAs_reprojected$IUCN_CAT %in% c("Not Reported", "Not Applicable") ]
PAs_amazon_Focal<-PAs_amazon_Focal[PAs_amazon_Focal$NAME!="Central Amazon Conservation Complex"]

PAs_congo<- crop(PAs, basins[1])
PAs_congo_Focal<-PAs_congo[expanse(PAs_congo, unit="km")>20000]# & !PAs_reprojected$IUCN_CAT %in% c("Not Reported", "Not Applicable") ]
PAs_congo_Focal<-PAs_congo_Focal[PAs_congo_Focal$NAME!="Salonga"]

PAs_Focal<-vect(c(PAs_amazon_Focal, PAs_congo_Focal))

PAs_rast<-(terra::rasterize(PAs_Focal, FRIPagg, "NAME"))


r<-c(FRIPagg, PAs_rast)

# 3) get spatial autocorrelation layers
r$I3<-focal(r$correlation, w=3, fun=mean, na.rm=TRUE)
r$I5<-focal(r$correlation, w=5, fun=mean, na.rm=TRUE)
r$I7<-focal(r$correlation, w=7, fun=mean, na.rm=TRUE)
r$I9<-focal(r$correlation, w=9, fun=mean, na.rm=TRUE)
r$I11<-focal(r$correlation, w=11, fun=mean, na.rm=TRUE)



# 4) get modelling Dataframe
df<-as.data.frame(r, cells=TRUE, na.rm=TRUE)

df<-df %>%
    as_tibble %>%
    filter (NAME%in% names(which(table(df$NAME)>5))) %>% 
    mutate(NAME = str_replace_all(NAME, "-", ":"))


# 5) modelling basins
# n.b. this is competed seperately due to nesting of country in basin making group compariosns non-estimable if modelled together

# basin + protection
m1<-aov(correlation~NAME, df)
m2<-aov(correlation~NAME+I3, df) #
m3<-aov(correlation~NAME+I3+I5, df)
m4<-aov(correlation~NAME+I3+I5+I7, df)
m5<-aov(correlation~NAME+I3+I5+I7+I9, df)
m6<-aov(correlation~NAME+I3+I7, df)
m7<-aov(correlation~NAME+I3+I9, df)
m8<-aov(correlation~NAME+I9, df) #
m9<-aov(correlation~NAME+I3+I11, df)
m10<-aov(correlation~NAME+I11, df) #
m11<-aov(correlation~NAME+I7, df) #
m12<-aov(correlation~NAME+I5, df) #
m13<-aov(correlation~NAME+I3+I7+I9, df)
m14<-aov(correlation~NAME+I5+I7, df)


# basin* protection

MuMIn::model.sel(
    m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14)

#m3 is best
summary(m3)
TukeyHSD(m3, "NAME")


# 7) plot model results for PAs
m<-m3
# country wise plot
LABELS <- generate_label_df(TukeyHSD(m, "NAME") , "NAME")

df_wLetters<-df %>% mutate(treatment=NAME) %>% left_join(., LABELS) %>% mutate(NAME = str_replace_all(NAME, "-", ":"))

myPalette = c(
    'e' =   "#AD002AFF",
    'af' = "#925E9FFF", 
    'a' = "#ADB6B6FF", 
    'b' = "#42B540FF", 
    'bd' = "#00468BFF",
    'ef' = "#ED0000FF",
    'cd' = "#FDAF91FF",
    'bd' = "#ED0000FF",
    'cd' = "#FDAF91FF",
    'ac' = "orange",
    'bcd' ="cyan")

fig3c<-ggplot(data=df_wLetters, aes(y=reorder(NAME, correlation), x=correlation, fill=Letters))+
    geom_vline(xintercept=0, linetype=2, alpha=0.5)+
    geom_boxplot(outliers=FALSE, color="black")+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.6), alpha=1, size=0.5)+
    theme_classic()+
    #guides(color="none")+
    scale_fill_manual(values=myPalette)+
    xlab("Flooding role in productivity")+ylab("")

PAsForPlotting<-PAs_Focal
values(PAsForPlotting)<-values(PAs_Focal) %>%
    as_tibble %>%
    mutate(NAME = str_replace_all(NAME, "-", ":")) %>%
    left_join(., LABELS, join_by(NAME==treatment))
PAsForPlotting<-PAsForPlotting[!is.na(values(PAsForPlotting)$Letters)]

plotDat<-crop(FRIP, basins[1])
fig3a<-ggplot() +
    geom_spatvector(data=crop(PAsForPlotting, ext(plotDat)), aes(fill=Letters), linewidth=1, color="black")+
    geom_spatvector(data=crop(countries, ext(plotDat)), fill=NA, linewidth=1, color="black")+
    scale_fill_manual(values=myPalette, na.value = "transparent")+
    guides(fill="none")+
    theme_classic()


plotDat<-crop(FRIP, basins[2])
fig3b<-ggplot() +
    geom_spatvector(data=crop(PAsForPlotting, ext(plotDat)), aes(fill=Letters), linewidth=1, color="black")+
    geom_spatvector(data=crop(countries, ext(plotDat)), fill=NA, linewidth=1, color="black")+
    scale_fill_manual(values=myPalette, na.value = "transparent")+
    guides(fill="none")+
    theme_classic()


fig3column<-cowplot::plot_grid(fig3a, fig3b, labels = c('A', 'B'), label_size = 12, ncol=1, rel_heights = c(0.825, 1))

cowplot::plot_grid(fig3column, fig3c, labels = c("", 'C'), label_size = 12, ncol=2, rel_widths = c(1, 1))
ggsave("Figures/Figure3.png", height=15, width=18)

