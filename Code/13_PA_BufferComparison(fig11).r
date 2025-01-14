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
PAs_agg<-aggregate(PAs_reprojected)# fore erasing purposes

a<-crop(cor, project(basins[2], crs(cor)), mask=TRUE) %>% trim
c<-crop(cor, project(basins[1], crs(cor)), mask=TRUE) %>% trim



## amazon #############################################
#first pas
PAs_reprojected<- crop(project(PAs, crs(a)), ext(a))
PAsFocal<-PAs_reprojected[expanse(PAs_reprojected, unit="km")>40000]# & !PAs_reprojected$IUCN_CAT %in% c("Not Reported", "Not Applicable") ]
PAsFocal<-PAsFocal[PAsFocal$NAME!="Central Amazon Conservation Complex"]
plot(PAsFocal, col="black")

subbasins_by_PA_l<-list()
vals_by_PA_l<-list()
for (i in 1:length(PAsFocal) ) {
    subbasins_by_PA_l[[PAsFocal$NAME[i]]]<-crop(project(subbasins_amazon, crs(PAsFocal)), PAsFocal[i])
    BufferForThisPA <-erase(erase(buffer(PAsFocal[i], 10^5), PAsFocal[i]), PAs_agg)
    subbasins_by_PA_l[[ paste0( PAsFocal$NAME[i], " unprotected buffer") ]]<-crop(project(subbasins_amazon, crs(PAsFocal)), BufferForThisPA)

    print(paste0("Starting extraction for ", PAsFocal$NAME[i]))

    vals_by_PA_l[[PAsFocal$NAME[i]]]<-terra::extract( cor, subbasins_by_PA_l[[PAsFocal$NAME[i]]], weights=TRUE, ID=TRUE) %>%
        filter(!is.na(correlation)) %>%
        filter(weight>0.5) %>%
        group_by(ID) %>%
        summarise(mean=mean(correlation), n=n()) %>%
        filter(n>10) #%>%
        #mutate(mean=case_when(mean<0~0, .default =mean))
    print(paste0("finised extraction for ", PAsFocal$NAME[i]))

    print(paste0("Starting extraction for buffer of ", PAsFocal$NAME[i]))

    vals_by_PA_l[[  paste0( PAsFocal$NAME[i], " unprotected buffer")    ]]<-terra::extract( cor, subbasins_by_PA_l[[ paste0( PAsFocal$NAME[i], " unprotected buffer") ]], weights=TRUE, ID=TRUE) %>%
        filter(!is.na(correlation)) %>%
        filter(weight>0.5) %>%
        group_by(ID) %>%
        summarise(mean=mean(correlation), n=n()) %>%
        filter(n>10) #%>%
        #mutate(mean=case_when(mean<0~0, .default =mean))
    print(paste0("finised extraction for buffer of ", PAsFocal$NAME[i]))
}

df_amazon<-bind_rows(vals_by_PA_l, .id = "id")


## congo####################################
#first pas
PAs_reprojected<- crop(project(PAs, crs(c)), ext(c))
PAsFocal<-PAs_reprojected[expanse(PAs_reprojected, unit="km")>20000 ]#& !PAs_reprojected$IUCN_CAT %in% c("Not Reported", "Not Applicable") ]
PAsFocal<-PAsFocal[PAsFocal$NAME!="Salonga"]
plot(PAsFocal, col="black")

subbasins_by_PA_l<-list()
vals_by_PA_l<-list()
for (i in 1:length(PAsFocal) ) {
    subbasins_by_PA_l[[PAsFocal$NAME[i]]]<-crop(project(subbasins_congo, crs(PAsFocal)), PAsFocal[i])
    BufferForThisPA <-erase(erase(buffer(PAsFocal[i], 10^5), PAsFocal[i]), PAs_agg)
    subbasins_by_PA_l[[ paste0( PAsFocal$NAME[i], " unprotected buffer") ]]<-crop(project(subbasins_congo, crs(PAsFocal)), BufferForThisPA)

    print(paste0("Starting extraction for ", PAsFocal$NAME[i]))

    vals_by_PA_l[[PAsFocal$NAME[i]]]<-terra::extract( cor, subbasins_by_PA_l[[PAsFocal$NAME[i]]], weights=TRUE, ID=TRUE) %>%
        filter(!is.na(correlation)) %>%
        filter(weight>0.5) %>%
        group_by(ID) %>%
        summarise(mean=mean(correlation), n=n()) %>%
        filter(n>10) #%>%
        #mutate(mean=case_when(mean<0~0, .default =mean))
    print(paste0("finised extraction for ", PAsFocal$NAME[i]))

    print(paste0("Starting extraction for buffer of ", PAsFocal$NAME[i]))

    vals_by_PA_l[[  paste0( PAsFocal$NAME[i], " unprotected buffer")    ]]<-terra::extract( cor, subbasins_by_PA_l[[ paste0( PAsFocal$NAME[i], " unprotected buffer") ]], weights=TRUE, ID=TRUE) %>%
        filter(!is.na(correlation)) %>%
        filter(weight>0.5) %>%
        group_by(ID) %>%
        summarise(mean=mean(correlation), n=n()) %>%
        filter(n>10) #%>%
        #mutate(mean=case_when(mean<0~0, .default =mean))
    print(paste0("finised extraction for buffer of ", PAsFocal$NAME[i]))
}

df_congo<-bind_rows(vals_by_PA_l, .id = "id")


df<-rbind(df_amazon, df_congo)
 
summary_df<-df %>% group_by(id) %>%
    summarise(mean=mean(mean), n=n())

df<-filter(df, !id %in% summary_df$id[summary_df$n<10]) %>% mutate(id = str_replace_all(id, "-", ":"))


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




ggplot(data=df, aes(y=reorder(id, mean), x=mean, fill=Letters))+
    geom_boxplot(outliers=FALSE)+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.6), alpha=1, size=1)+
    geom_vline(xintercept=0, linetype=2, alpha=0.5)+
    theme_classic()+
    ggsci::scale_fill_npg()+
    #guides(color="none")+
    xlab("Flooding role in productivity")+ylab("")


ggsave("Figures/Figure11.png", height=15, width=19)

