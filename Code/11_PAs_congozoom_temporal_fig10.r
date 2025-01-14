library(terra)
library(tidyverse)
library(tidyterra)
library(ggsignif)
library(multcompView)


mk_cor<-rast("Outputs/KendallResults.tif")
basins<-vect("Outputs/Basins")
countries<-vect("Data/world-administrative-boundaries")
subbasins_congo<-crop(vect("Data/hybas_af_lev08_v1c"), basins[1])
PAs<-vect("Outputs/BasinPAs")


c<-crop(mk_cor, project(basins[1], crs(mk_cor)), mask=TRUE) %>% trim



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


vals<-terra::extract( mk_cor, PAsFocal, weights=TRUE, ID=TRUE) %>%
    filter(!is.na(tau)) %>%
    filter(weight>0.99) %>%
    group_by(ID) %>%
    summarise(mean=mean(tau), n=n()) %>%
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
df<-terra::extract( mk_cor, PAsFocal, weights=TRUE, ID=TRUE) %>%
    left_join(name_df) %>%
    filter(!is.na(tau)) %>%
    filter(weight>0.99)

summary_df<-df %>% group_by(NAME) %>%
    summarise(mean=mean(tau, na.rm=TRUE), n=n())

df<-filter(df, !NAME %in% summary_df$NAME[summary_df$n<5]) %>% mutate(NAME = str_replace_all(NAME, "-", ":"))


m<-aov(tau~NAME, df) 
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


fig10a<-ggplot(data=df, aes(y=reorder(NAME, tau), x=tau, fill=Letters))+
    geom_boxplot(outliers=FALSE)+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.6), alpha=1, size=1)+
    geom_vline(xintercept=0, linetype=2, alpha=0.5)+
    theme_classic()+
    viridis::scale_fill_viridis(discrete = TRUE) +    
    #guides(color="none")+
    xlab("Change in FRIP(2020-2023)")+ylab("")



PAs_wLetters<-PAs
values(PAs_wLetters)<-values(PAs) %>%
    as_tibble %>%
    mutate(NAME = str_replace_all(NAME, "-", ":")) %>%
    left_join(., LABELS, join_by(NAME==treatment))


fig10b<-ggplot() +
    geom_spatvector(data=crop(PAs_wLetters[!is.na(PAs_wLetters$Letters)], ext(subbasins_congo)), aes(fill=Letters, color=), linewidth=1, color="black")+
    viridis::scale_fill_viridis(discrete = TRUE) +    
    geom_spatvector(data=crop(countries, ext(subbasins_congo)), fill=NA, linewidth=1, color="black")+
    #guides(fill="none")+
    theme_classic()

cowplot::plot_grid(fig10b, fig10a, labels = c("A", 'B'), label_size = 12, ncol=1, rel_widths = c(2, 1))
ggsave("Figures/Figure10.png", height=20, width=13)

