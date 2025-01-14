library(terra)
library(tidyverse)
library(tidyterra)


cor100<-rast("Data/NppCorr_5km.tif") # 5km loaded as cor100 for spped as was working late
#plot(cor5)

#cor100<-rast("Data/NppCorr_100km.tif")
plot(cor100)


PAs<-vect("Outputs/BasinPAs")
basins<-vect("Outputs/Basins")


PA_congo<-crop(PAs, basins[1]) 
PA_amazon<-crop(PAs, basins[2]) 

protectedStatus<-c("Ia", "Ib", "II")

unprotected_c<-erase(basins[1], PAs)
semiprotected_c<-PA_congo[(!PA_congo$IUCN_CAT %in% protectedStatus)]
protected_c<-PA_congo[(PA_congo$IUCN_CAT %in% protectedStatus)]

unprotected_a<-erase(basins[2], PAs)
semiprotected_a<-PA_amazon[(!PA_amazon$IUCN_CAT %in% protectedStatus)]
protected_a<-PA_amazon[(PA_amazon$IUCN_CAT %in% protectedStatus)]

ggplot() +
    geom_spatvector(data=protected_a, fill="green", color=NA ) +
    geom_spatvector(data=unprotected_a, fill="red", color=NA ) +
    geom_spatvector(data=protected_c, fill="green", color=NA ) +
    geom_spatvector(data=unprotected_c, fill="red", color=NA ) +
    geom_spatraster(data = cor100, aes(fill=correlation))+
    viridis::scale_fill_viridis(na.value = "transparent")
ggsave("Figures/FullMap.png")


d<-rbind(
    cbind(
        rbind(
            cbind(terra::extract( cor100, unprotected_a, weights=TRUE), class="unprotected"),
            cbind(terra::extract( cor100, semiprotected_a, weights=TRUE), class="semiprotected"),
            cbind(terra::extract( cor100, protected_a, weights=TRUE), class="protected")
        ), 
        basin="amazon"
        ),
        cbind(
        rbind(
            cbind(terra::extract( cor100, unprotected_c, weights=TRUE), class="unprotected"),
            cbind(terra::extract( cor100, semiprotected_c, weights=TRUE), class="semiprotected"),
            cbind(terra::extract( cor100, protected_c, weights=TRUE), class="protected")
        ), 
        basin="congo"
        )
)


#variance tests - neds different input
congo<-d %>%
    as_tibble() %>%
    filter(weight>0.99) %>% 
    filter(!is.na(correlation)) %>%
    filter(class!="semiprotected") %>%
    filter(basin=="congo") 


    shapiro.test(congo$correlation[1:5000])
    var.test(correlation~class, congo ,  ratio=1, alternative = "less")


ggplot(data=congo, aes(x=class, y=correlation, color=basin))+
    #geom_hline(data=av, aes(yintercept=median, color=basin), linetype=2, alpha=0.5)+
    geom_boxplot(data=congo, outliers=FALSE)+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha=0.1)+
    theme_classic()

amazon<-d %>%
    as_tibble() %>%
    filter(weight>0.99) %>% 
    filter(!is.na(correlation)) %>%
    filter(class!="semiprotected") %>%
    filter(basin=="amazon") 

    shapiro.test(amazon$correlation[1:5000])
    #fligner.test(correlation~class, .)
    var.test(correlation~class, . ,  alternative = "less")


fligner.test(correlation~class, datAll)



#uordered factor
dat<- d %>% filter(weight>0.99) %>% filter(!is.na(correlation)) %>% filter(correlation<0)
ggplot()+
    geom_boxplot(data=dat, aes(x=class, y=correlation, color=basin))

lm(correlation~class*basin, dat) %>% summary

#ordered factor
classes <- c("unprotected","semiprotected","protected")
dat2<- dat %>%    mutate(class= factor(class, levels = classes, ordered = T) )

lm(correlation~class*basin, dat2) %>% summary



###make composite boxplot
pos<- d %>% filter(weight>0.99) %>% filter(!is.na(correlation)) %>% filter(correlation>0)
neg<- d %>% filter(weight>0.99) %>% filter(!is.na(correlation)) %>% filter(correlation<0)
composite_data<-rbind(pos, neg)%>%
    mutate(dir=correlation>0)

av<-composite_data %>%
    group_by(dir, basin) %>%
    summarise(median=median(correlation))


ggplot(data=composite_data, aes(x=class, y=correlation, color=basin, shape=dir))+
    geom_hline(data=av, aes(yintercept=median, color=basin), linetype=2, alpha=0.5)+
    geom_boxplot(data=composite_data, outliers=FALSE)+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha=0.1)+
    theme_classic()
ggsave("Figures/DemonstrationPlot.png")
    





#############
d_basins<-terra::extract( cor100, basins, weights=TRUE) %>% filter(weight>0.99) %>% filter(!is.na(correlation))%>% filter(correlation<0)

head(d_basins)

ggplot()+
    geom_boxplot(data=d_basins, aes(x=as.factor(ID), y=correlation))


