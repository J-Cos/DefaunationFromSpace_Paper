# explore negret et al 2020 method

library(terra)
library(tidyverse)
library(multcompView)
library(tidyterra)

#checking no autocorrelation in residuals
testResidualAutocorrelation<-function(model, type="pvalue"){
    df$resid<-m$resid
    n <- r$correlation
    n[df$cell] <-df$resid

    ac<-autocor( n, global=TRUE, method="moran")

    m <- sapply(1:99, function(i) {
        tempRast<-n
        tempRast[!is.na(tempRast)]<- sample(n[!is.na(n)])
        autocor( tempRast, global=TRUE, method="moran")
    })
    pval <- sum(m >= ac) / 100
    if (type=="pvalue") { return(pval) } 
    else if (type=="plot") { return(plot(n)) } 
}

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
aggregationFactor<-20

# 1) load data 
FRIP<-rast("Outputs/FRIP.tif") 
basins<-vect("Outputs/Basins")
countries<-vect("Outputs/BasinCountries")
PAs<-vect("Outputs/BasinPAs")

# 2) make single raster
# aggregate FRIP
FRIPagg<-aggregate(FRIP, aggregationFactor, na.rm=TRUE)

#basin raster
basin_rast<-as.factor(terra::rasterize(basins, FRIPagg, "PFAF_ID"))

#country raster
countries_rast<-(terra::rasterize(countries, FRIPagg, "name"))

# protected raster
protected<-rasterize(PAs, FRIPagg, "IUCN_CAT")
protected$protected<-!is.na(protected)

r<-c(FRIPagg, basin_rast, countries_rast, protected$protected)

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
    filter (name %in% names(which(table(df$name)>5))) %>%
    mutate(protected=as.factor(protected))

m<-aov(correlation~PFAF_ID+name+protected, df)
TukeyHSD(m)
summary(lm(correlation~0+PFAF_ID+name+protected, df))

# 5) modelling basins
# n.b. this is competed seperately due to nesting of country in basin making group compariosns non-estimable if modelled together

# basin + protection
m1<-aov(correlation~PFAF_ID+protected, df)
m2<-aov(correlation~PFAF_ID+protected+I3, df) #
m3<-aov(correlation~PFAF_ID+protected+I3+I5, df)
m4<-aov(correlation~PFAF_ID+protected+I3+I5+I7, df)
m5<-aov(correlation~PFAF_ID+protected+I3+I5+I7+I9, df)
m6<-aov(correlation~PFAF_ID+protected+I3+I7, df)
m7<-aov(correlation~PFAF_ID+protected+I3+I9, df)
m8<-aov(correlation~PFAF_ID+protected+I9, df) #
m9<-aov(correlation~PFAF_ID+protected+I3+I11, df)
m10<-aov(correlation~PFAF_ID+protected+I11, df) #
m11<-aov(correlation~PFAF_ID+protected+I7, df) #
m12<-aov(correlation~PFAF_ID+protected+I5, df) #
m13<-aov(correlation~PFAF_ID+protected+I3+I7+I9, df)
m14<-aov(correlation~PFAF_ID+protected+I5+I7, df)


# basin* protection
m1i<-aov(correlation~PFAF_ID*protected, df)
m2i<-aov(correlation~PFAF_ID*protected+I3, df) #
m3i<-aov(correlation~PFAF_ID*protected+I3+I5, df)
m4i<-aov(correlation~PFAF_ID*protected+I3+I5+I7, df)
m5i<-aov(correlation~PFAF_ID*protected+I3+I5+I7+I9, df)
m6i<-aov(correlation~PFAF_ID*protected+I3+I7, df)
m7i<-aov(correlation~PFAF_ID*protected+I3+I9, df)
m8i<-aov(correlation~PFAF_ID*protected+I9, df) #
m9i<-aov(correlation~PFAF_ID*protected+I3+I11, df)
m10i<-aov(correlation~PFAF_ID*protected+I11, df) #
m11i<-aov(correlation~PFAF_ID*protected+I7, df) #
m12i<-aov(correlation~PFAF_ID*protected+I5, df) #
m13i<-aov(correlation~PFAF_ID*protected+I3+I7+I9, df)
m14i<-aov(correlation~PFAF_ID*protected+I5+I7, df)

MuMIn::model.sel(
    m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14,
    m1i, m2i, m3i, m4i, m5i, m6i, m7i, m8i, m9i, m10i, m11i, m12i, m13i, m14i)

#m4 is best
summary(m4)
TukeyHSD(m4, "PFAF_ID")
TukeyHSD(m4, "protected")

testResidualAutocorrelation(m4, "plot")
testResidualAutocorrelation(m4, "pvalue")



# 6) modelling countries

#country+protected
m1<-aov(correlation~name+protected, df)
m2<-aov(correlation~name+protected+I3, df) #
m3<-aov(correlation~name+protected+I3+I5, df)
m4<-aov(correlation~name+protected+I3+I5+I7, df)
m5<-aov(correlation~name+protected+I3+I5+I7+I9, df)
m6<-aov(correlation~name+protected+I3+I7, df)
m7<-aov(correlation~name+protected+I3+I9, df)
m8<-aov(correlation~name+protected+I9, df) #
m9<-aov(correlation~name+protected+I3+I11, df)
m10<-aov(correlation~name+protected+I11, df) #
m11<-aov(correlation~name+protected+I7, df) #
m12<-aov(correlation~name+protected+I5, df) #
m13<-aov(correlation~name+protected+I3+I7+I9, df)
m14<-aov(correlation~name+protected+I5+I7, df)

# country* protection
m1i<-aov(correlation~name*protected, df)
m2i<-aov(correlation~name*protected+I3, df) #
m3i<-aov(correlation~name*protected+I3+I5, df)
m4i<-aov(correlation~name*protected+I3+I5+I7, df)
m5i<-aov(correlation~name*protected+I3+I5+I7+I9, df)
m6i<-aov(correlation~name*protected+I3+I7, df)
m7i<-aov(correlation~name*protected+I3+I9, df)
m8i<-aov(correlation~name*protected+I9, df) #
m9i<-aov(correlation~name*protected+I3+I11, df)
m10i<-aov(correlation~name*protected+I11, df) #
m11i<-aov(correlation~name*protected+I7, df) #
m12i<-aov(correlation~name*protected+I5, df) #
m13i<-aov(correlation~name*protected+I3+I7+I9, df)
m14i<-aov(correlation~name*protected+I5+I7, df)

MuMIn::model.sel(
    m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14,
    m1i, m2i, m3i, m4i, m5i, m6i, m7i, m8i, m9i, m10i, m11i, m12i, m13i, m14i)

#m3 is best
summary(m3)
TukeyHSD(m3, "name")
TukeyHSD(m3, "protected")

testResidualAutocorrelation(m3, "plot")
testResidualAutocorrelation(m3, "pvalue")

# 7) plot model results for countries
m<-m3
# country wise plot
LABELS <- generate_label_df(TukeyHSD(m, "name") , "name")

df<-df %>% mutate(treatment=name) %>% left_join(., LABELS)

myPalette = c(
    'abc' =   "#AD002AFF",
    'ab' = "#925E9FFF", 
    'a' = "#ADB6B6FF", 
    'c' = "#42B540FF", 
    'bc' = "#00468BFF")#
    #'abd' = "#ED0000FF",
    #'c' = "#FDAF91FF")

fig3c<-ggplot(data=df, aes(y=reorder(name, correlation), x=correlation, fill=Letters))+
    geom_vline(xintercept=0, linetype=2, alpha=0.5)+
    geom_boxplot(outliers=FALSE, color="black")+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.6), alpha=1, size=0.5)+
    theme_classic()+
    #guides(color="none")+
    #scale_fill_manual(values=myPalette)+
    xlab("Flooding role in productivity")+ylab("")

countriesForPlotting<-countries
values(countriesForPlotting)<-values(countries) %>%
    as_tibble %>%
    left_join(., LABELS, join_by(name==treatment))


plotDat<-crop(FRIP, basins[1])
fig3a<-ggplot() +
    geom_spatvector(data=crop(countriesForPlotting, ext(plotDat)), aes(fill=Letters), linewidth=1, color="black")+
    scale_fill_manual(values=myPalette, na.value = "transparent")+
    guides(fill="none")+
    theme_classic()


plotDat<-crop(FRIP, basins[2])
fig3b<-ggplot() +
    geom_spatvector(data=crop(countriesForPlotting, ext(plotDat)), aes(fill=Letters), linewidth=1, color="black")+
    scale_fill_manual(values=myPalette, na.value = "transparent")+
    guides(fill="none")+
    theme_classic()


fig3column<-cowplot::plot_grid(fig3a, fig3b, labels = c('A', 'B'), label_size = 12, ncol=1, rel_heights = c(0.825, 1))

cowplot::plot_grid(fig3column, fig3c, labels = c("", 'C'), label_size = 12, ncol=2, rel_widths = c(1, 1))
ggsave("Figures/Figure2.png", height=15, width=18)

