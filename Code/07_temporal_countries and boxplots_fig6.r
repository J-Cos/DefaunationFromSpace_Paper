library(terra)
library(tidyverse)
library(tidyterra)
library(ggsignif)
library(multcompView)


mk_cor<-rast("Outputs/KendallResults.tif")
basins<-vect("Outputs/Basins")
countries<-vect("Data/world-administrative-boundaries")
subbasins_congo<-crop(vect("Data/hybas_af_lev08_v1c"), basins[1])
subbasins_amazon<-crop(vect("Data/hybas_sa_lev08_v1c"), basins[2])


a<-crop(mk_cor, project(basins[2], crs(mk_cor)), mask=TRUE) %>% trim
c<-crop(mk_cor, project(basins[1], crs(mk_cor)), mask=TRUE) %>% trim



## amazon
countries_reprojected<- crop(project(countries, crs(subbasins_amazon)), ext(subbasins_amazon))

subbasins_by_country_l<-list()
vals_by_country_l<-list()
for (i in 1:length(countries_reprojected) ) {
    subbasins_by_country_l[[countries_reprojected$name[i]]]<-crop(subbasins_amazon, countries_reprojected[i])


    print(paste0("Starting extraction for ", countries_reprojected$name[i]))

    vals_by_country_l[[countries_reprojected$name[i]]]<-terra::extract( mk_cor, subbasins_by_country_l[[countries_reprojected$name[i]]], weights=TRUE, ID=TRUE) %>%
        filter(!is.na(tau)) %>%
        filter(weight>0.99) %>%
        group_by(ID) %>%
        summarise(mean=mean(tau), n=n()) %>%
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

    vals_by_country_l[[countries_reprojected$name[i]]]<-terra::extract( mk_cor, subbasins_by_country_l[[countries_reprojected$name[i]]], weights=TRUE, ID=TRUE) %>%
        filter(!is.na(tau)) %>%
        filter(weight>0.99) %>%
        group_by(ID) %>%
        summarise(mean=mean(tau), n=n()) %>%
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


m<-aov(mean~id, df) 
m_reduced<-aov(mean~0, df) 
anova(m_reduced, m)

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
    'a' = "#42B540FF", 
    'abc' = "#925E9FFF", 
    'abcd' = "#00468BFF",
    'ac' = "#ADB6B6FF", 
    'ad' = "#AD002AFF", 
    'b' = "#ED0000FF",
    'bc' = "#FDAF91FF", 
    'd' = "#0099B4FF")

fig6c<-ggplot(data=df, aes(y=reorder(id, mean), x=mean, color=Letters))+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.6), alpha=1, size=2)+
    geom_vline(xintercept=0, linetype=2, alpha=0.5)+
    geom_boxplot(outliers=FALSE, color="black", fill=NA)+
    theme_classic()+
    scale_color_manual(values=myPalette)+
    #guides(color="none")+
    xlab("Change in FRIP(2020-2023)")+ylab("")


values(countries)<-values(countries) %>%
    as_tibble %>%
    left_join(., LABELS, join_by(name==treatment))

fig6a<-ggplot() +
    geom_spatvector(data=crop(countries, ext(subbasins_amazon)), aes(fill=Letters), linewidth=1, color="black")+
    guides(fill="none")+
    scale_fill_manual(values=myPalette, na.value = "transparent")+
    theme_classic()


fig6b<-ggplot() +
    geom_spatvector(data=crop(countries, ext(subbasins_congo)), aes(fill=Letters), linewidth=1, color="black")+
    guides(fill="none")+
    scale_fill_manual(values=myPalette, na.value = "transparent")+
    theme_classic()

fig6column<-cowplot::plot_grid(fig6a, fig6b, labels = c('A', 'B'), label_size = 12, ncol=1, rel_heights = c(0.825, 1))

cowplot::plot_grid(fig6column, fig6c, labels = c("", 'C'), label_size = 12, ncol=1, rel_heights = c(1, 0.4))
ggsave("Figures/Figure6.png", height=18, width=13)

