
library(terra)
library(tidyverse)
library(multcompView)
library(tidyterra)

# functions
getDfFromModList<-function(mods_list){
    Res_df<-lapply(mods_list, function(x){rownames_to_column(as.data.frame( confint(x) ), "NAME")}) %>% 
        bind_rows( .id="aggFactor") %>% 
        as_tibble %>% 
        mutate(
            aggFactor=(as.numeric(aggFactor)),
            lwr=`2.5 %`,
            upr=`97.5 %`) %>%
        mutate(signif=lwr>0 | upr<0) 
    return(Res_df)
}

################
# set results resolution
################
res<-"25000"

# 1) load data 
FRIPandMetrics_l<-lapply(list.files("Outputs/FRIP", full.names=TRUE), function(x){rast(x)})
FRIPtrend_l<-lapply(list.files("Outputs/FRIPtrend", full.names=TRUE), function(x){rast(x)})
r_l<-list()
for (i in 1:length(FRIPandMetrics_l)){
    r_l[[i]]<-c(FRIPandMetrics_l[[i]], FRIPtrend_l[[i]])
}

basins<-vect("Outputs/Basins")
countries<-vect("Outputs/BasinCountries")
PAs<-vect("Outputs/BasinPAs")



# Amazon - Bogoni vs Benitez
# 2) fit models across aggregation factors ----------------------------------------------------------------
aggFactors<-unlist(lapply(FRIPandMetrics_l, function(x){ as.character(str_sub(varnames(x)[1], 6) )}))

weights_list<-mod1_list<-mod2_list<-df_list<-list()
for (i in 1:length(aggFactors)) {

    ragg<-r_l[[i]] 
    ragg <- mask(ragg, is.na(ragg$DI_Bogoni), maskvalue=TRUE)

    # 4) get modelling Dataframe
    df<-as.data.frame(ragg, cells=TRUE, na.rm=TRUE) %>%
        as_tibble# %>%
        #filter (name %in% names(which(table(df$name)>50))) 

    m1<-lm(correlation~DI, df)
    m2<-lm(correlation~DI_Bogoni, df)

    weights_list[[aggFactors[i]]]<-MuMIn::model.sel(m1, m2) %>% as.data.frame %>% rownames_to_column("model") %>% select(model, weight)

    mod1_list[[aggFactors[i]]]<-m1
    mod2_list[[aggFactors[i]]]<-m2
    df_list[[aggFactors[i]]]<-df
}

# 4) get results dataframe

Res_df<-rbind(
        cbind(getDfFromModList(mod1_list), "model"=1),
        cbind(getDfFromModList(mod2_list), "model"=2)
    ) %>% 
    filter(NAME!="(Intercept)") %>% 
    mutate(NAME=as.factor(NAME)) %>%
    mutate(aggFactor=as.numeric(aggFactor)/1000) %>%
    as_tibble %>%
    arrange((aggFactor))


# 5) print out results for text
mod2_list[[res]] %>% summary
mod2_list[[res]] %>% confint
MaxAggFactorOfFinding<-which(!filter(Res_df, NAME=="DI")$signif)[1]-1
lapply(mod1_list, function(x){summary(x)$r.squared})


# 6) make supplementary figure
Strip_names <- c(
                    `DI_Bogoni` = "Defaunation Index- Bogoni et al. (2020)",
                    `DI` = "Defauantion Index - Benítez-López et al. (2019)")

Bogoni_df<-filter(Res_df, NAME=="DI_Bogoni")
panel_Bogoni<-ggplot(Bogoni_df, aes(x=aggFactor, alpha=signif))+
    geom_hline(yintercept=0)+
    #geom_vline(xintercept=25000, color="red")+
    geom_hline(yintercept=  mean(c(Bogoni_df$upr, Bogoni_df$lwr)), linetype=2 )+
    geom_linerange(aes(ymin = lwr, ymax = upr))+
    facet_wrap(~NAME, labeller = as_labeller(Strip_names))+
    theme_classic()+
    guides(alpha=FALSE)+
    ylab("Coefficient estimate")+
    xlab("Pixel size (km)")

Benitez_df<-filter(Res_df, NAME=="DI")
panel_DI<-ggplot(Benitez_df, aes(x=aggFactor, alpha=signif))+
    geom_hline(yintercept=0)+
    #geom_vline(xintercept=25000, color="red")+
    geom_hline(yintercept=  mean(c(Benitez_df$upr, Benitez_df$lwr)), linetype=2 )+
    geom_linerange(aes(ymin = lwr, ymax = upr))+
    facet_wrap(~NAME, labeller = as_labeller(Strip_names))+
    theme_classic()+
    guides(alpha=FALSE)+
    ylab("Coefficient estimate")+
    xlab("Pixel size (km)")



panel_weight<-bind_rows(weights_list, .id="aggFactor")  %>% 
        as_tibble %>% 
        mutate(
            aggFactor=as.numeric(aggFactor)/1000) %>%
        filter(model=="m2") %>%
    ggplot(., aes(x=aggFactor, y=weight))+
        geom_line()+
        geom_hline(yintercept=0.5, linetype=3)+
        ylim(0,1)+
        theme_classic()+
        ylab("Akaike weight")+
        facet_wrap(~"Evidence in favour of Bogoni et al. Defaunation Index Model")+
        xlab("")

fig<-cowplot::plot_grid(panel_Bogoni+xlab(""), panel_weight, panel_DI, labels = c('A', 'B', 'C'), label_size = 12, ncol=1)
fig

ggsave("Figures/FigureS2.png", height=10, width=10)

saveRDS(panel_Bogoni, "Outputs/panel_Bogoni.RDS")
