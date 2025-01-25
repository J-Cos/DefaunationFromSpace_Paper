# explore accesibility

library(terra)
library(tidyverse)
library(multcompView)
library(tidyterra)

# 1) load data 
FRIP<-rast("Outputs/FRIP.tif") 
basins<-vect("Outputs/Basins")
countries<-vect("Outputs/BasinCountries")
PAs<-vect("Outputs/BasinPAs")

metrics<-rast("Outputs/existingMetrics.tif")

# 2) make single raster

#basin raster
basin_rast<-as.factor(terra::rasterize(basins, FRIP, "PFAF_ID"))

#country raster
countries_rast<-(terra::rasterize(countries, FRIP, "name"))

# protected raster
protected<-rasterize(PAs, FRIP, "IUCN_CAT")
protected$protected<-!is.na(protected)

r<-c(FRIP, basin_rast, countries_rast, protected$protected, metrics)

# 3) fit models across aggregation factors
weights_list<-mod1_list<-mod2_list<-list()
for (i in 1:50) {
    ragg<-aggregate(r, i, na.rm=TRUE)

    # 4) get modelling Dataframe
    df<-as.data.frame(ragg, cells=TRUE, na.rm=TRUE)

    df<-df %>%
        as_tibble %>%
        #filter (name %in% names(which(table(df$name)>50))) %>%
        mutate(protected=as.factor(protected))

    m1<-lm(correlation~accessibility, df)
    m2<-lm(correlation~DI, df)

    weights_list[[i]]<-MuMIn::model.sel(m1, m2) %>% as.data.frame %>% rownames_to_column("model") %>% select(model, weight)

    mod1_list[[i]]<-m1
    mod2_list[[i]]<-m2
}


# 4) get results dataframe
getDfFromModList<-function(mods_list){
    Res_df<-lapply(mods_list, function(x){rownames_to_column(as.data.frame( confint(x) ), "NAME")}) %>% 
        bind_rows( .id="aggFactor") %>% 
        as_tibble %>% 
        mutate(
            aggFactor=(5*as.numeric(aggFactor))^2,
            lwr=`2.5 %`,
            upr=`97.5 %`) %>%
        mutate(signif=lwr>0 | upr<0) 
    return(Res_df)
}

Res_df<-rbind(
        cbind(getDfFromModList(mod1_list), "model"=1),
        cbind(getDfFromModList(mod2_list), "model"=2)
    ) %>% 
    filter(NAME!="(Intercept)") %>% 
    mutate(NAME=as.factor(NAME)) %>%
    as_tibble


# 5) print out results for text
mod2_list[[10]] %>% summary
mod2_list[[10]] %>% confint
MaxAggFactorOfFinding<-which(!filter(Res_df, NAME=="DI")$signif)[1]-1
(MaxAggFactorOfFinding*5)^2

#6) make predicted FRIP raster
r10<-aggregate(r, 10, na.rm=TRUE)
FRIPpred<-predict(r10, mod2_list[[10]]) %>%
    crop(., r10$correlation, mask=TRUE)

DIunderestimate<-r10$correlation- FRIPpred
names(DIunderestimate)<-"DIunderestimate"

writeRaster(DIunderestimate, "Outputs/DIunderestimate.tif", overwrite=TRUE)

# 6) make supplementary figure
Strip_names <- c(
                    `accessibility` = "Inaccessibility (travel time to city in minutes)- Weiss et al. (2018)",
                    `DI` = "Defauantion Index - Benítez-López et al. (2019)")

access_df<-filter(Res_df, NAME=="accessibility")
panel_access<-ggplot(access_df, aes(x=aggFactor, alpha=signif))+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=(5*10)^2, color="red")+
    geom_hline(yintercept=  mean(c(access_df$upr, access_df$lwr)), linetype=2 )+
    geom_linerange(aes(ymin = lwr, ymax = upr))+
    facet_wrap(~NAME, labeller = as_labeller(Strip_names))+
    theme_classic()+
    guides(alpha=FALSE)+
    ylab("Coefficient estimate")+
    xlab("Pixel size (km2)")

DI_df<-filter(Res_df, NAME=="DI")
panel_DI<-ggplot(DI_df, aes(x=aggFactor, alpha=signif))+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=(5*10)^2, color="red")+
    geom_hline(yintercept=  mean(c(DI_df$upr, DI_df$lwr)), linetype=2 )+
    geom_linerange(aes(ymin = lwr, ymax = upr))+
    facet_wrap(~NAME, labeller = as_labeller(Strip_names))+
    theme_classic()+
    guides(alpha=FALSE)+
    ylab("Coefficient estimate")+
    xlab("")



panel_weight<-bind_rows(weights_list, .id="aggFactor")  %>% 
        as_tibble %>% 
        mutate(
            aggFactor=(5*as.numeric(aggFactor))^2) %>%
        filter(model=="m2") %>%
    ggplot(., aes(x=aggFactor, y=weight))+
        geom_line()+
        geom_hline(yintercept=0.5, linetype=3)+
        ylim(0,1)+
        theme_classic()+
        ylab("Akaike weight of evidence in favour of Defaunation Index model")+
        xlab("")

fig<-cowplot::plot_grid(panel_DI, panel_weight, panel_access, labels = c('A', 'B', 'C'), label_size = 12, ncol=1)
ggsave("Figures/FigureS2.png", height=18, width=10)

