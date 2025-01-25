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
PAs<-vect("Outputs/BasinPAs")


a<-crop(cor, project(basins[2], crs(cor)), mask=TRUE) %>% trim
c<-crop(cor, project(basins[1], crs(cor)), mask=TRUE) %>% trim



## amazon
PAs_reprojected<- crop(project(PAs, crs(a)), ext(a))
PAsFocal<-PAs_reprojected[expanse(PAs_reprojected, unit="km")>40000]# & !PAs_reprojected$IUCN_CAT %in% c("Not Reported", "Not Applicable") ]
PAsFocal<-PAsFocal[PAsFocal$NAME!="Central Amazon Conservation Complex"]
plot(PAsFocal, col="black")

subbasins_by_PA_l<-list()
vals_by_PA_l<-list()
for (i in 1:length(PAsFocal) ) {
    subbasins_by_PA_l[[PAsFocal$NAME[i]]]<-crop(project(subbasins_amazon, crs(PAsFocal)), PAsFocal[i])


    print(paste0("Starting extraction for ", PAsFocal$NAME[i]))

    vals_by_PA_l[[PAsFocal$NAME[i]]]<-terra::extract( cor, subbasins_by_PA_l[[PAsFocal$NAME[i]]], weights=TRUE, ID=TRUE) %>%
        filter(!is.na(correlation)) %>%
        filter(weight>0.5) %>%
        group_by(ID) %>%
        summarise(mean=mean(correlation), n=n()) %>%
        filter(n>10) #%>%
        #mutate(mean=case_when(mean<0~0, .default =mean))
    print(paste0("finised extraction for ", PAsFocal$NAME[i]))

}

df_amazon<-bind_rows(vals_by_PA_l, .id = "id")


## congo
PAs_reprojected<- crop(project(PAs, crs(c)), ext(c))
PAsFocal<-PAs_reprojected[expanse(PAs_reprojected, unit="km")>20000 ]#& !PAs_reprojected$IUCN_CAT %in% c("Not Reported", "Not Applicable") ]
PAsFocal<-PAsFocal[PAsFocal$NAME!="Salonga"]
plot(PAsFocal, col="black")

subbasins_by_PA_l<-list()
vals_by_PA_l<-list()
for (i in 1:length(PAsFocal) ) {
    subbasins_by_PA_l[[PAsFocal$NAME[i]]]<-crop(project(subbasins_congo, crs(PAsFocal)), PAsFocal[i])


    print(paste0("Starting extraction for ", PAsFocal$NAME[i]))

    vals_by_PA_l[[PAsFocal$NAME[i]]]<-terra::extract( cor, subbasins_by_PA_l[[PAsFocal$NAME[i]]], weights=TRUE, ID=TRUE) %>%
        filter(!is.na(correlation)) %>%
        filter(weight>0.5) %>%
        group_by(ID) %>%
        summarise(mean=mean(correlation), n=n()) %>%
        filter(n>10) #%>%
        #mutate(mean=case_when(mean<0~0, .default =mean))
    print(paste0("finised extraction for ", PAsFocal$NAME[i]))

}

df_congo<-bind_rows(vals_by_PA_l, .id = "id")




df<-rbind(df_amazon, df_congo) %>%
    filter()
 
summary_df<-df %>% group_by(id) %>%
    summarise(mean=mean(mean), n=n())

df<-filter(df, !id %in% summary_df$id[summary_df$n<5]) %>% mutate(id = str_replace_all(id, "-", ":"))


m<-aov(mean~id, df) 
summary(m)
TukeyHSD(m) %>% plot(las=1)
df_network<-TukeyHSD(m)$id %>% as.data.frame %>%rownames_to_column %>% as_tibble %>% filter(`p adj`<0.05)

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
LABELS <- generate_label_df(TukeyHSD(m) , "id")

df<-df %>% left_join(., LABELS, by=join_by(id==treatment))


myPalette = c(
    'bc' = "#42B540FF", 
    'ad' = "#925E9FFF", 
    'bd' = "#00468BFF",
    'a' = "#ADB6B6FF", 
    'b' = "#AD002AFF", 
    'abcd' = "#ED0000FF",
    'c' = "#FDAF91FF")


fig7c<-ggplot(data=df, aes(y=reorder(id, mean), x=mean, fill=Letters))+
    geom_boxplot(outliers=FALSE)+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.6), alpha=1, size=1)+
    geom_vline(xintercept=0, linetype=2, alpha=0.5)+
    theme_classic()+
    scale_fill_manual(values=myPalette)+
    #guides(color="none")+
    xlab("Flooding role in productivity")+ylab("")



PAs_wLetters<-PAs
values(PAs_wLetters)<-values(PAs) %>%
    as_tibble %>%
    left_join(., LABELS, join_by(NAME==treatment))

fig7a<-ggplot() +
    geom_spatvector(data=crop(PAs_wLetters[!is.na(PAs_wLetters$Letters)], ext(subbasins_amazon)), aes(fill=Letters), linewidth=1, color="black")+
    scale_fill_manual(values=myPalette, na.value = "transparent")+
    geom_spatvector(data=crop(countries, ext(subbasins_amazon)), fill=NA, linewidth=1, color="black")+
    #guides(fill="none")+
    theme_classic()


fig7b<-ggplot() +
    geom_spatvector(data=crop(PAs_wLetters[!is.na(PAs_wLetters$Letters)], ext(subbasins_congo)), aes(fill=Letters), linewidth=1, color="black")+
    scale_fill_manual(values=myPalette, na.value = "transparent")+
    geom_spatvector(data=crop(countries, ext(subbasins_congo)), fill=NA, linewidth=1, color="black")+
    #guides(fill="none")+
    theme_classic()

fig7column<-cowplot::plot_grid(fig7a, fig7b, labels = c('A', 'B'), label_size = 12, ncol=1, rel_heights = c(0.825, 1))

cowplot::plot_grid(fig7column, fig7c, labels = c("", 'C'), label_size = 12, ncol=2, rel_widths = c(1, 1))
ggsave("Figures/Figure7.png", height=15, width=19)

