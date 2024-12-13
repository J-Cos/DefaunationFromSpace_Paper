
library(terra)
library(tidyverse)
library(tidyterra)


cor100<-rast("Data/NppCorr_5km_includingNonsignif.tif") # 5km loaded as cor100 for spped as was working late
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
    geom_spatraster(data = cor100) #+
   # geom_spatvector(data=unprotected_a)+
    #geom_spatvector(data=protected_c)



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
AIC(m1, m2)
summary(m1)

#test if negbinom better than poisson
m3 <- glm(poscor ~ class+basin, data = negbinom, family = "poisson")
pchisq(2 * (logLik(m1) - logLik(m3)), df = 1, lower.tail = FALSE)

(est <- cbind(Estimate = coef(m1), confint(m1)))
exp(est)