
library(terra)
library(tidyverse)
library(tidyterra)


cor100<-rast("Data/NppCorr_5km.tif") # 5km loaded as cor100 for spped as was working late
#plot(cor5)

#cor100<-rast("Data/NppCorr_100km.tif")
plot(cor100)

countries<-vect("Data/world-administrative-boundaries")
basins<-vect("Outputs/Basins")
subbasins_congo<-crop(vect("Data/hybas_af_lev06_v1c"), basins[1])
subbasins_amazon<-crop(vect("Data/hybas_sa_lev06_v1c"), basins[2])

PAs<-vect("Outputs/BasinPAs")
PA_congo<-crop(PAs, basins[1]) 
PA_amazon<-crop(PAs, basins[2]) 

lowprotectedStatus<-c("V", "VI")

unprotected_c<-erase(subbasins_congo, PAs)
unprotected_c<-unprotected_c[expanse(unprotected_c, unit="km")>100]
#semiprotected_c<-PA_congo[(!PA_congo$IUCN_CAT %in% protectedStatus)]
#protected_c<-PA_congo[(PA_congo$IUCN_CAT %in% protectedStatus)]
protected_c<-PA_congo[!PA_congo$IUCN_CAT %in% lowprotectedStatus]


unprotected_a<-erase(subbasins_amazon, PAs)
unprotected_a<-unprotected_a[expanse(unprotected_a, unit="km")>100]
#semiprotected_a<-PA_amazon[(!PA_amazon$IUCN_CAT %in% protectedStatus)]
#protected_a<-PA_amazon[(PA_amazon$IUCN_CAT %in% protectedStatus)]
protected_a<-PA_amazon[!PA_amazon$IUCN_CAT %in% lowprotectedStatus]

d<-crop(cor100, project(basins[2], crs(cor100))) %>% trim

ggplot() +
    geom_spatvector(data=crop(project(PAs, crs(d)), ext(d)), linewidth=1.2, alpha=0.1, fill="green" ) +
    geom_spatraster(data = d) +
    geom_spatvector(data=crop(project(countries, crs(d)), ext(d)), linewidth=2, fill=NA, color="black") +
    scale_fill_gradient2(na.value = "transparent", mid="light grey")+
    theme_classic()


getAveragePerArea<-function(cat){
    terra::extract( cor100, cat, weights=TRUE, ID=TRUE) %>%
        filter(!is.na(correlation)) %>%
        filter(weight>0.99) %>%
        group_by(ID) %>%
        summarise(mean=mean(correlation))
}

d<-rbind(
    cbind(
        rbind(
            cbind(getAveragePerArea(unprotected_a), class="unprotected"),
            #cbind(getAveragePerArea(semiprotected_a), class="semiprotected"),
            cbind(getAveragePerArea(protected_a), class="protected")
        ), 
        basin="amazon"
        ),
        cbind(
        rbind(
            cbind(getAveragePerArea(unprotected_c), class="unprotected"),
            #cbind(getAveragePerArea(semiprotected_c), class="semiprotected"),
            cbind(getAveragePerArea(protected_c), class="protected")
        ), 
        basin="congo"
        )
)

d_med<-rbind(
    cbind(
        rbind(
            cbind(getAveragePerArea(unprotected_a), class="unprotected"),
            cbind(getAveragePerArea(semiprotected_a), class="semiprotected"),
            cbind(getAveragePerArea(protected_a), class="protected"),
            cbind(getAveragePerArea(subbasins_amazon), class="subbasins")

        ), 
        basin="amazon"
        ),
        cbind(
        rbind(
            cbind(getAveragePerArea(unprotected_c), class="unprotected"),
            cbind(getAveragePerArea(semiprotected_c), class="semiprotected"),
            cbind(getAveragePerArea(protected_c), class="protected"),
            cbind(getAveragePerArea(subbasins_congo), class="subbasins")
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
    fligner.test(correlation~class, congo)

    var.test(correlation~class, congo ,  ratio=1, alternative = "less")


ggplot(data=amazon, aes(x=class, y=correlation, color=basin))+
    #geom_hline(data=av, aes(yintercept=median, color=basin), linetype=2, alpha=0.5)+
    geom_boxplot(data=amazon, outliers=FALSE)+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2), alpha=0.1)+
    theme_classic()

amazon<-d %>%
    as_tibble() %>%
    filter(weight>0.99) %>% 
    filter(!is.na(correlation)) %>%
    filter(class!="semiprotected") %>%
    filter(basin=="amazon") 

    shapiro.test(amazon$correlation[1:5000])
    fligner.test(correlation~class, amazon)
    var.test(correlation~class, . ,  alternative = "less")


fligner.test(correlation~class, datAll)


d %>%
    as_tibble() %>%
    filter(weight>0.99) %>% 
    filter(!is.na(correlation)) %>%
    filter(class!="semiprotected") %>%
    mutate(poscor=sqrt(correlation^2)) %>%
    ggplot(data=., aes(x=class, y=poscor, color=basin))+
        #geom_hline(data=av, aes(yintercept=median, color=basin), linetype=2, alpha=0.5)+
        geom_boxplot(outliers=FALSE)+
        #geom_violin()+
        geom_jitter(position = position_jitterdodge(jitter.width = 0.3), alpha=0.1)+
        #geom_violin()+
        theme_classic()
ggsave("Figures/NegBinomPlot.png")

    d %>%
    as_tibble() %>%
    filter(weight>0.99) %>% 
    filter(!is.na(correlation)) %>%
    filter(class=="unprotected") %>%
    fligner.test(correlation~basin, .)

    d %>%
    as_tibble() %>%
    filter(weight>0.99) %>% 
    filter(!is.na(correlation)) %>%
    filter(class=="protected") %>%
    fligner.test(correlation~basin, .)




library(MASS)
#subsample
#test<-d[sample(1:nrow(d),nrow(d)*0.2 ),]

negbinom<-d %>%
    as_tibble() %>%
    filter(weight>0.99) %>% 
    filter(!is.na(correlation)) %>%
    filter(class!="semiprotected") %>%
    mutate(poscor=sqrt(correlation^2))
m1 <- glm.nb(poscor ~ class+basin, data = negbinom)
m2 <- glm.nb(poscor ~ class*basin, data = negbinom)
m3 <- glm.nb(poscor ~ class, data = negbinom)
m4 <- glm.nb(poscor ~ basin, data = negbinom)

AIC(m1, m2, m3, m4)
summary(m1)
summary(m2)

#test if negbinom better than poisson
m3 <- glm(poscor ~ class+basin, data = negbinom, family = "poisson")
pchisq(2 * (logLik(m1) - logLik(m3)), df = 1, lower.tail = FALSE)

(est <- cbind(Estimate = coef(m1), confint(m1)))
exp(est)


#########new
negbinom<-d_med %>%
    as_tibble() %>%
    filter(class=="subbasins") %>%
    mutate(poscor=sqrt(mean^2))

ggplot(negbinom, aes(y=mean, x=class, color=basin))+
    geom_boxplot(outliers=FALSE)+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2))+
    theme_classic()

m1 <- lm(mean ~ class+basin, data = negbinom)
m2 <- lm(mean ~ class*basin, data = negbinom)
m3 <- lm(mean ~ class, data = negbinom)
m4 <- lm(mean ~ basin, data = negbinom)

MuMIn::model.sel(m1, m2, m3, m4)
summary(m1)
summary(m2)
summary(m3)
summary(m4)








### testing in small subbasins

subbasins_congo<-crop(vect("Data/hybas_af_lev08_v1c"), basins[1])
subbasins_amazon<-crop(vect("Data/hybas_sa_lev08_v1c"), basins[2])

getAveragePerArea<-function(cat){
    terra::extract( cor100, cat, weights=TRUE, ID=TRUE) %>%
        filter(!is.na(correlation)) %>%
        filter(weight>0.99) %>%
        group_by(ID) %>%
        summarise(mean=mean(correlation), n=n())
}

d<-rbind(
    cbind(
        rbind(
            cbind(getAveragePerArea(subbasins_amazon), class="subbasins")

        ), 
        basin="amazon"
        ),
        cbind(
        rbind(
            cbind(getAveragePerArea(subbasins_congo), class="subbasins")
        ), 
        basin="congo"
        )
)

m <- lm(mean ~ basin, data = filter(d, n>20))

summary(m)

ggplot(filter(d, n>20), aes(y=mean, x=basin, color=basin))+
    geom_jitter(position = position_jitterdodge(jitter.width = 1), aes(size=log(n)), alpha=0.4)+
    geom_boxplot(outliers=FALSE, color="black", fill=NA)+
    theme_classic()

m<-data.frame()
for (i in 1:20){
    m[i, "est"] <-  lm(mean ~ basin, data = filter(d, n>i))$coef[2]
    m[i, "pval"] <-  summary(lm(mean ~ basin, data = filter(d, n>i)))$coef[2,4]
    m[i, "filter"]<-i
    m[i, "n"]<-sum(d$n>i)
}

ggplot(m, aes(x=filter, y=est))+
    geom_point(aes( color=pval>0.05, size=log(n)))+
    geom_smooth(method="lm", color="black")+
    theme_classic()




#########################
    vals<-terra::extract( cor100, subbasins_congo, weights=TRUE, ID=TRUE) %>%
        filter(!is.na(correlation)) %>%
        filter(weight>0.99) %>%
        group_by(ID) %>%
        summarise(mean=mean(correlation), n=n()) %>%
        filter(n>5) #%>%
        #mutate(mean=case_when(mean<0~0, .default =mean))
    values(subbasins_congo)[ vals$ID, "mean"]<-vals$mean

ggplot() +
    geom_spatvector(data=subbasins_congo, aes(fill=mean))+
    viridis::scale_fill_viridis()+
    geom_spatvector(data=crop(countries, ext(subbasins_congo)), color="black", linewidth=2, fill=NA)+
    geom_spatvector(data= allPA, color=NA, fill="black", linetwidth=1, alpha=0.5)


allPA<-aggregate(PA_congo)
allPA<-aggregate(PA_congo[expanse(PA_congo, unit="km")>10000])
allPA<-aggregate(PA_congo[expanse(PA_congo, unit="km")>10000 & !PA_congo$IUCN_CAT %in% c("V", "VI")])

ggplot() +
    geom_spatvector(data= PA_congo[expanse(PA_congo, unit="km")>10000 & !PA_congo$IUCN_CAT %in% c("V", "VI")], color=NA, aes(fill=IUCN_CAT), linetwidth=1, alpha=0.5)

#########################
    vals<-terra::extract( cor100, subbasins_amazon, weights=TRUE, ID=TRUE) %>%
        filter(!is.na(correlation)) %>%
        filter(weight>0.99) %>%
        group_by(ID) %>%
        summarise(mean=mean(correlation), n=n()) %>%
        filter(n>5) #%>%
        #mutate(mean=case_when(mean<0~0, .default =mean))
    values(subbasins_amazon)[ vals$ID, "mean"]<-vals$mean

ggplot() +
    geom_spatvector(data=subbasins_amazon, aes(fill=mean))+
    viridis::scale_fill_viridis()+
    geom_spatvector(data=crop(countries, ext(subbasins_amazon)), color="black", linewidth=2, fill=NA)+
    geom_spatvector(data= allPA, color=NA, fill="black", linetwidth=1, alpha=0.5)

allPA<-aggregate(PA_amazon)
allPA<-aggregate(PA_amazon[expanse(PA_amazon, unit="km")>10000])
allPA<-aggregate(PA_amazon[expanse(PA_amazon, unit="km")>10000 & !PA_amazon$IUCN_CAT %in% c("V", "VI")])

ggplot() +
    geom_spatvector(data= PA_amazon[expanse(PA_amazon, unit="km")>10000 & !PA_amazon$IUCN_CAT %in% c("V", "VI")], color=NA, aes(fill=IUCN_CAT), linetwidth=1, alpha=0.5)




######################################

protectedSubbasins<-crop(subbasins_amazon, allPA)
otherSubbasins<-subbasins_amazon[!subbasins_amazon$HYBAS_ID %in% protectedSubbasins$HYBAS_ID]
