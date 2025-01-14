library(terra)
library(tidyverse)
library(tidyterra)
library(ggsignif)
library(multcompView)


cor<-rast("Data/NppCorr_5km_includingNonsignif.tif") 
basins<-vect("Outputs/Basins")
countries<-vect("Data/world-administrative-boundaries")
subbasins_congo<-crop(vect("Data/hybas_af_lev08_v1c"), basins[1])
PAs<-vect("Outputs/BasinPAs")


c<-crop(cor, project(basins[1], crs(cor)), mask=TRUE) %>% trim



## congo
PAs_reprojected<- crop(project(PAs, crs(c)), ext(c))
PAsFocal<-PAs_reprojected[expanse(PAs_reprojected, unit="km")>5000 ]#& !PAs_reprojected$IUCN_CAT %in% c("Not Reported", "Not Applicable") ]
PAsFocal<-PAsFocal[!PAsFocal$NAME %in% 
    c(  "Salonga",  
        "Dja",
        "Forest Massif of Odzala-Kokoua",
        "Garamba",
        "Kahuzi-Biega",
        "Odzala Kokoua",
        "Okapis",
        "Réserve du triangle de la Ngiri",
        "Tumba-Lediima"

    )
    ]
plot(PAsFocal, col="black")


vals<-terra::extract( cor, PAsFocal, weights=TRUE, ID=TRUE) %>%
    filter(!is.na(correlation)) %>%
    filter(weight>0.99) %>%
    group_by(ID) %>%
    summarise(mean=mean(correlation), n=n()) %>%
    filter(n>10) #%>%
    #mutate(mean=case_when(mean<0~0, .default =mean))
values(PAsFocal)[ vals$ID, "mean"]<-vals$mean
PAsFocal<-crop(PAsFocal, ext(project(c, crs(PAsFocal))))

figpanel_unused<-ggplot() +
    geom_spatvector(data=PAsFocal, aes(fill=mean), color="black")+
    geom_spatvector(data=crop(countries, ext(project(PAsFocal, crs(countries)))), color="black", linewidth=2, fill=NA)+
    scale_fill_gradient2(limits = c(-1,1), na.value = "transparent", high=scales::muted("red"), low=scales::muted("blue"), mid="light grey", name="Flooding\nrole in\nproductivity")+
    theme_classic()

#########################################
name_df<-data.frame("ID" = seq(1:length(PAsFocal$NAME)), "NAME"=PAsFocal$NAME)
df<-terra::extract( cor, PAsFocal, weights=TRUE, ID=TRUE, xy=TRUE) %>%
    left_join(name_df) %>%
    filter(!is.na(correlation)) %>%
    filter(weight>0.99)

summary_df<-df %>% group_by(NAME) %>%
    summarise(mean=mean(correlation, na.rm=TRUE), n=n())

df<-df %>%
    filter(!duplicated(paste0(df$x, df$y))) %>%
    #filter(!NAME %in% summary_df$NAME[summary_df$n>100]) %>% 
    mutate(NAME = str_replace_all(NAME, "-", ":"))


m<-aov(correlation~NAME, df) 
summary(m)
tHSD<-TukeyHSD(m, conf.level=0.05)
df_network<-tHSD$NAME %>% as.data.frame %>%rownames_to_column %>% as_tibble %>% filter(`p adj`<0.001)

# I need to group the treatments that are not different each other together.
generate_label_df <- function(TUKEY, variable){
 
     # Extract labels and factor levels from Tukey post-hoc 
     Tukey.levels <- TUKEY[[variable]][,4]
     Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
     
     #I need to put the labels in the same order as in the boxplot :
     Tukey.labels$treatment=rownames(Tukey.labels)
     Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
     return(Tukey.labels)
     }
 
# Apply the function on my dataset
LABELS <- generate_label_df(tHSD , "NAME")

df<-df %>% left_join(., LABELS, by=join_by(NAME==treatment))


fig9a<-ggplot(data=df, aes(y=reorder(NAME, correlation), x=correlation, fill=Letters))+
    geom_boxplot(outliers=FALSE)+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.6), alpha=1, size=0.5)+
    geom_vline(xintercept=0, linetype=2, alpha=0.5)+
    theme_classic()+
    viridis::scale_fill_viridis(discrete = TRUE) +    
    #guides(color="none")+
    xlab("Flooding role in productivity")+ylab("")



PAs_wLetters<-PAs
values(PAs_wLetters)<-values(PAs) %>%
    as_tibble %>%
    mutate(NAME = str_replace_all(NAME, "-", ":")) %>%
    left_join(., LABELS, join_by(NAME==treatment))


fig9b<-ggplot() +
    geom_spatvector(data=crop(PAs_wLetters[!is.na(PAs_wLetters$Letters)], ext(subbasins_congo)), aes(fill=Letters, color=), linewidth=1, color="black")+
    viridis::scale_fill_viridis(discrete = TRUE) +    
    geom_spatvector(data=crop(countries, ext(subbasins_congo)), fill=NA, linewidth=1, color="black")+
    #guides(fill="none")+
    theme_classic()

cowplot::plot_grid(fig9b, fig9a, labels = c("A", 'B'), label_size = 12, ncol=1, rel_widths = c(2, 1))
ggsave("Figures/Figure9.png", height=20, width=13)






########### - alternativ emodelling approach accounting for spatial autocorrelation

library(nlme)

#filter out uncertain protection category

data.spatialCor.gls <- gls(correlation ~ NAME, data = df,
    method = "REML", verbose=TRUE)
data.spatialCor.glsExp <- gls(correlation ~ NAME, data = df,
    correlation = corExp(form = ~x + y, nugget = TRUE),
    method = "REML", verbose=TRUE)
data.spatialCor.glsGaus <- gls(correlation ~ NAME, data = df,
    correlation = corGaus(form = ~x + y, nugget = TRUE),
    method = "REML")
data.spatialCor.glsLin <- gls(correlation ~ NAME, data = df,
    correlation = corLin(form = ~x + y, nugget = TRUE),
    method = "REML")
data.spatialCor.glsRatio <- gls(correlation ~ NAME, data = df,
    correlation = corRatio(form = ~x + y, nugget = TRUE),
    method = "REML")
data.spatialCor.glsSpher <- gls(correlation ~ NAME, data = df,
    correlation = corSpher(form = ~x + y, nugget = TRUE),
    method = "REML")


AIC(data.spatialCor.gls, data.spatialCor.glsExp, data.spatialCor.glsGaus,
    data.spatialCor.glsLin, data.spatialCor.glsRatio,
    data.spatialCor.glsSpher)


plot(residuals(data.spatialCor.gls, type = "normalized") ~
    fitted(data.spatialCor.gls))
plot(data.spatialCor.gls)
plot(nlme:::Variogram(data.spatialCor.gls, form = ~x +
    y, resType = "normalized"))


plot(residuals(data.spatialCor.glsExp, type = "normalized") ~
    fitted(data.spatialCor.glsExp))


summary(data.spatialCor.gls)
summary(data.spatialCor.glsExp)
summary(data.spatialCor.glsRatio)