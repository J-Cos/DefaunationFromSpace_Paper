library(terra)
library(tidyverse)
library(tidyterra)
library(ggsignif)
library(multcompView)


cor<-rast("Data/NppCorr_5km_includingNonsignif.tif") 
basins<-vect("Outputs/Basins")
countries<-vect("Data/world-administrative-boundaries")
subbasins_congo<-crop(vect("Data/hybas_af_lev08_v1c"), basins[1])
subbasins_amazon<-crop(vect("Data/hybas_sa_lev08_v1c"), basins[2])


a<-crop(cor, project(basins[2], crs(cor)), mask=TRUE) %>% trim
c<-crop(cor, project(basins[1], crs(cor)), mask=TRUE) %>% trim



## amazon
countries_reprojected<- crop(project(countries, crs(subbasins_amazon)), ext(subbasins_amazon))

subbasins_by_country_l<-list()
vals_by_country_l<-list()
for (i in 1:length(countries_reprojected) ) {
    subbasins_by_country_l[[countries_reprojected$name[i]]]<-crop(subbasins_amazon, countries_reprojected[i])


    print(paste0("Starting extraction for ", countries_reprojected$name[i]))

    vals_by_country_l[[countries_reprojected$name[i]]]<-terra::extract( cor, subbasins_by_country_l[[countries_reprojected$name[i]]], weights=TRUE, ID=TRUE) %>%
        filter(!is.na(correlation)) %>%
        filter(weight>0.99) %>%
        group_by(ID) %>%
        summarise(mean=mean(correlation), n=n()) %>%
        filter(n>10) #%>%
        #mutate(mean=case_when(mean<0~0, .default =mean))
    print(paste0("finised extraction for ", countries_reprojected$name[i]))

}

df_amazon<-bind_rows(vals_by_country_l, .id = "id")


## congo
countries_reprojected<- crop(project(countries, crs(subbasins_congo)), ext(subbasins_congo))

subbasins_by_country_l<-list()
vals_by_country_l<-list()
for (i in 1:length(countries_reprojected) ) {
    subbasins_by_country_l[[countries_reprojected$name[i]]]<-crop(subbasins_congo, countries_reprojected[i])


    print(paste0("Starting extraction for ", countries_reprojected$name[i]))

    vals_by_country_l[[countries_reprojected$name[i]]]<-terra::extract( cor, subbasins_by_country_l[[countries_reprojected$name[i]]], weights=TRUE, ID=TRUE) %>%
        filter(!is.na(correlation)) %>%
        filter(weight>0.99) %>%
        group_by(ID) %>%
        summarise(mean=mean(correlation), n=n()) %>%
        filter(n>10) #%>%
        #mutate(mean=case_when(mean<0~0, .default =mean))
    print(paste0("finised extraction for ", countries_reprojected$name[i]))

}

df_congo<-bind_rows(vals_by_country_l, .id = "id")




df<-rbind(df_amazon, df_congo) %>%
    filter()
 
summary_df<-df %>% group_by(id) %>%
    summarise(mean=mean(mean), n=n())

df<-filter(df, !id %in% summary_df$id[summary_df$n<5])


basin_df<-data.frame(
    id=c( "Ecuador", "Venezuela", "Brazil" ,   "Peru" , "Colombia",   "Bolivia"    ,                     
        "Congo"   ,     "Cameroon"  , "Central African Republic"    ,     "Democratic Republic of the Congo" , "Gabon" ),
    basin=c(rep("Amazon", 6), rep("Congo", 5)))

df<-left_join(df, basin_df)


m<-aov(mean~basin/id, df) 
summary(m)
TukeyHSD(m) %>% plot(las=1)
df_network<-TukeyHSD(m)['basin:id'] %>% as.data.frame %>%rownames_to_column %>% as_tibble %>% filter(`basin.id.p.adj`<0.05)

# I need to group the treatments that are not different each other together.
generate_label_df <- function(TUKEY, variable){
 
     # Extract labels and factor levels from Tukey post-hoc 
     Tukey.levels <- TUKEY[[variable]][,4] 
     Tukey.levels<-Tukey.levels[!is.na(Tukey.levels)]
     Tukey.labels <- data.frame(multcompLetters(x=Tukey.levels)['Letters'])
     
     #I need to put the labels in the same order as in the boxplot :
     Tukey.labels$treatment=rownames(Tukey.labels)
     Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
     return(Tukey.labels)
     }
 
# Apply the function on my dataset
LABELS <- generate_label_df(TukeyHSD(m) , "basin:id")

df<-df %>% mutate(treatment=paste0(basin, ":", id)) %>% left_join(., LABELS)

myPalette = c(
    'abcd' =  "#ADB6B6FF", 
    'abc' = "#925E9FFF", 
    'd' = "#00468BFF",
    'bc' = "#42B540FF", 
    'ad' = "#AD002AFF", 
    'abd' = "#ED0000FF",
    'c' = "#FDAF91FF")


fig3c<-ggplot(data=df, aes(y=reorder(id, mean), x=mean, fill=Letters))+
    geom_vline(xintercept=0, linetype=2, alpha=0.5)+
    geom_boxplot(outliers=FALSE, color="black")+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.6), alpha=1, size=1)+
    theme_classic()+
    #guides(color="none")+
    scale_fill_manual(values=myPalette)+
    xlab("Flooding role in productivity")+ylab("")




values(countries)<-values(countries) %>%
    as_tibble %>%
    left_join(., mutate(LABELS, name=str_split_i(treatment, ":", i=2)))

fig3a<-ggplot() +
    geom_spatvector(data=crop(countries, ext(subbasins_amazon)), aes(fill=Letters), linewidth=1, color="black")+
    scale_fill_manual(values=myPalette, na.value = "transparent")+
    guides(fill="none")+
    theme_classic()


fig3b<-ggplot() +
    geom_spatvector(data=crop(countries, ext(subbasins_congo)), aes(fill=Letters), linewidth=1, color="black")+
    scale_fill_manual(values=myPalette, na.value = "transparent")+
    guides(fill="none")+
    theme_classic()

fig3column<-cowplot::plot_grid(fig3a, fig3b, labels = c('A', 'B'), label_size = 12, ncol=1, rel_heights = c(0.825, 1))

cowplot::plot_grid(fig3column, fig3c, labels = c("", 'C'), label_size = 12, ncol=2, rel_widths = c(1, 1))
ggsave("Figures/Figure3.png", height=15, width=18)

