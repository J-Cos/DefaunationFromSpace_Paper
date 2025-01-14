library(terra)
library(tidyverse)
library(tidyterra)
library(Kendall)

# functions
fun_kendall <- function(x){ return(unlist(Kendall::MannKendall(x)))
}

getDataframe<-function(input){
    d<-rbind(
        cbind(
            rbind(
                cbind(terra::extract( input, unprotected_a, weights=TRUE), class="unprotected"),
                cbind(terra::extract( input, semiprotected_a, weights=TRUE), class="semiprotected"),
                cbind(terra::extract( input, protected_a, weights=TRUE), class="protected")
            ), 
            basin="amazon"
            ),
            cbind(
            rbind(
                cbind(terra::extract( input, unprotected_c, weights=TRUE), class="unprotected"),
                cbind(terra::extract( input, semiprotected_c, weights=TRUE), class="semiprotected"),
                cbind(terra::extract( input, protected_c, weights=TRUE), class="protected")
            ), 
            basin="congo"
            )
    )
    return(d)
}

# 1) get all polygon data
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



# 2) load and preprocess data
cor<-rast("Data/NppCorr_5km_Annual.tif") 
plot(cor)

mk_cor<-terra::app(cor, fun_kendall)
mk_cor<- mask(x=mk_cor, mask=mk_cor$tau==1, maskvalues=1)
plot(mk_cor$tau)

# get only those pixels that always have positive or negative correlations
# this is a very strict criteria a looser critera: pos<-(mean(cor)>0) - gives 
# qualtitatively similar results
pos<-(min(cor)>0)
neg<-(max(cor)<0) 



# 3) exploatory analysis of significant trends only

trends<-mask(x=mk_cor$tau, mask=mk_cor$sl<0.05, maskvalues=0)
plot(trends)

pos_trends<-mask(x=mk_cor, mask=pos, maskvalues=0)
plot(pos_trends)

neg_trends<-mask(x=trends, mask=neg, maskvalues=0)
plot(neg_trends)

plotDat<-crop(trends, terra::project(basins[2], terra::crs(trends)), mask=TRUE) %>% trim()
ggplot() +
    geom_spatraster(data = plotDat) +
   # geom_spatvector(data=crop(project(PAs, crs(trends)), ext(plotDat)), fill=NA ) +
    scale_fill_gradient2(na.value = "transparent", mid="light grey")+
    theme_classic()

subbasins_congo<-crop(vect("Data/hybas_af_lev08_v1c"), basins[1])
subbasins_amazon<-crop(vect("Data/hybas_sa_lev08_v1c"), basins[2])
    vals<-terra::extract( mk_cor$tau, subbasins_amazon, weights=TRUE, ID=TRUE) %>%
        filter(!is.na(tau)) %>%
        filter(weight>0.99) %>%
        group_by(ID) %>%
        summarise(mean=mean(tau), n=n())
        #mutate(mean=case_when(mean<0~0, .default =mean))
    values(subbasins_amazon)[ vals$ID, "mean"]<-vals$mean
    values(subbasins_amazon)[ vals$ID, "n"]<-vals$n

ggplot() +
    geom_spatvector(data=subbasins_amazon, aes(fill=mean))+
    scale_fill_gradient2(na.value = "light grey", mid="light grey")



######
d_pos<-getDataframe(pos_trends)
dat_pos<- d_pos %>% filter(weight>0.99) %>% filter(!is.na(tau)) #%>% filter(tau>0)
classes <- c("unprotected","semiprotected","protected")
dat_pos2<- dat_pos %>%    mutate(class= factor(class, levels = classes, ordered = T) )

ggplot()+
    geom_violin(data=dat_pos, aes(x=class, y=tau, color=basin))

lm(tau~class+basin, dat_pos) %>% summary

lm(tau~class+basin, dat_pos2) %>% summary

######
d_neg<-getDataframe(neg_trends)
dat_neg<- d_neg %>% filter(weight>0.99) %>% filter(!is.na(tau)) #%>% filter(tau>0)
classes <- c("unprotected","semiprotected","protected")
dat_neg2<- dat_neg %>%    mutate(class= factor(class, levels = classes, ordered = T) )

ggplot()+
    geom_violin(data=dat_neg, aes(x=class, y=tau, color=basin))

lm(tau~class+basin, dat_neg) %>% summary

lm(tau~class+basin, dat_neg2) %>% summary

###########################################


amazon<-crop(pos_trends, terra::project(basins[2], terra::crs(trends)), mask=TRUE)
ggplot() +
    geom_spatraster(data = amazon) +
    geom_spatvector(data=PA_amazon, fill=NA ) 


congo<-crop(pos_trends, terra::project(basins[1], terra::crs(trends)), mask=TRUE)

ggplot() +
    geom_spatraster(data = congo) +
    geom_spatvector(data=PA_congo, fill=NA ) 


# 4) main analysis of all data (including non-significant)
#############################
#all data
AllDat<-getDataframe(mk_cor) %>%
    filter(weight>0.99)
str(AllDat)

AllDat %>%
    as_tibble() %>%
    filter(!is.na(tau)) %>%
    #filter(sl<0.05) %>%
    filter(class!="semiprotected") %>%
    ggplot(data=. , aes(x=class, color=basin, y=tau))+
        geom_boxplot()+
        geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha=0.1)

m1<-AllDat %>%
    as_tibble() %>%
    filter(!is.na(tau)) %>%
    #filter(sl<0.05) %>%
    filter(class!="semiprotected") %>%
    lm (tau~class*basin, .) 
m2<-AllDat %>%
    as_tibble() %>%
    filter(!is.na(tau)) %>%
    #filter(sl<0.05) %>%
    filter(class!="semiprotected") %>%
    lm (tau~class+basin, .) 

summary(m1)
summary(m2)
AIC(m1, m2)

##############################
# subset to negative and positive correlations
# this makes more sense as the FNP effect of a positive or negative trend
# depends on whether the correlation was positive or negative in the first place


pos_trends<-mask(x=mk_cor, mask=pos, maskvalues=0)
ggplot() +
    geom_spatvector(data=protected_a, fill="green", color=NA ) +
    geom_spatvector(data=unprotected_a, fill="red", color=NA ) +
    geom_spatvector(data=protected_c, fill="green", color=NA ) +
    geom_spatvector(data=unprotected_c, fill="red", color=NA ) +
    geom_spatraster(data = pos_trends, aes(fill=tau))+
    viridis::scale_fill_viridis(na.value = "transparent")
ggsave("Figures/TemporalTrendsPlot.png")

neg_trends<-mask(x=mk_cor, mask=neg, maskvalues=0)
plot(neg_trends)


AllPos<-getDataframe(pos_trends) %>%
    filter(weight>0.99)
str(AllPos)

AllPos %>%
    as_tibble() %>%
    filter(!is.na(tau)) %>%
    #filter(sl<0.05) %>%
    filter(class!="semiprotected") %>%
    ggplot(data=., aes(x=class, y=tau, color=basin))+
        geom_hline(yintercept=0)+
        geom_boxplot(outliers=FALSE)+
        geom_jitter(position = position_jitterdodge(jitter.width = 0.3), alpha=0.1)+
        theme_classic()+
        ggtitle("Mann-kendall trends in FF-NPP correlation for pixels \n that show a positive correlation in each year from 2020-2023")
ggsave("Figures/TemporalTrendsPlot.png")

#subsample?
#test<-AllPos[sample(1:nrow(AllPos),nrow(AllPos)*0.2 ),]

m1<-AllPos %>%
    as_tibble() %>%
    filter(!is.na(tau)) %>%
    #filter(sl<0.05) %>%
    filter(class!="semiprotected") %>%
    lm (tau~class*basin, .) 

m2<-AllPos %>%
    as_tibble() %>%
    filter(!is.na(tau)) %>%
    #filter(sl<0.05) %>%
    filter(class!="semiprotected") %>%
    lm (tau~class+basin, .) 

m3<-AllPos %>%
    as_tibble() %>%
    filter(!is.na(tau)) %>%
    #filter(sl<0.05) %>%
    filter(class!="semiprotected") %>%
    lm (tau~class:basin, .) 

summary(m1)
AIC(m1, m2, m3)
summary(m2)
plot(m1)

# the negative models show nothing the negative correlations should
#potential be considered largely meaningless
AllNeg<-getDataframe(neg_trends) %>%
    filter(weight>0.99)
str(AllNeg)

AllNeg %>%
    as_tibble() %>%
    filter(!is.na(tau)) %>%
    #filter(sl<0.05) %>%
    filter(class!="semiprotected") %>%
    ggplot(data=. , aes(x=class, color=basin, y=tau))+
        geom_boxplot()+
        geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha=0.1)

AllNeg %>%
    as_tibble() %>%
    filter(!is.na(tau)) %>%
    #filter(sl<0.05) %>%
    filter(class!="semiprotected") %>%
    lm (tau~class*basin, .) %>% 
    summary