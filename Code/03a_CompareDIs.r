
library(terra)
library(tidyverse)
library(multcompView)
library(tidyterra)

# functions
getDfFromModList<-function(mods_list){
    Res_df<-lapply(mods_list, function(x){
            cbind(
                rownames_to_column(as.data.frame( confint(x) ), "NAME"),
                "R2"=summary(x)$r.squared)
        }) %>% 
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

mod_list<-df_list<-list()
for (i in 1:length(aggFactors)) {

    ragg<-r_l[[i]] 
    ragg <- mask(ragg, is.na(ragg$DI_Bogoni), maskvalue=TRUE)

    # 4) get modelling Dataframe
    df<-as.data.frame(ragg, cells=TRUE, na.rm=TRUE) %>%
        as_tibble# %>%
        #filter (name %in% names(which(table(df$name)>50))) 

    m1<-lm(DI_Bogoni~DI, df)


    mod_list[[aggFactors[i]]]<-m1
    df_list[[aggFactors[i]]]<-df
}


# 4) get results dataframe

Res_df<-getDfFromModList(mod_list) %>% 
    filter(NAME!="(Intercept)") %>% 
    mutate(NAME=as.factor(NAME)) %>%
    mutate(aggFactor=as.numeric(aggFactor)/1000) %>%
    as_tibble %>%
    arrange((aggFactor))

#3) compare DIs
bind_rows(mod_list, .id="aggFactor") %>%
mutate(aggFactor=as.numeric(aggFactor)/1000) %>%
ggplot(., aes(x=DI, y=DI_Bogoni, group=as.factor(aggFactor), color=aggFactor))+
    geom_point(alpha=0.1)+
    geom_smooth(method="lm")+
    facet_wrap(~aggFactor)

lapply(df_list, function(item){
    summary(lm(DI~DI_Bogoni, item))$r.squared
    }
) %>% unlist%>% summary

################################
# 6) make supplementary figure
Strip_names <- c(
                    `DI_Bogoni` = "Defaunation Index- Bogoni et al. (2020)",
                    `DI` = "Defauantion Index - Benítez-López et al. (2019)")




p_coef<-ggplot(Res_df, aes(x=aggFactor))+
    #geom_hline(yintercept=0)+
    #geom_vline(xintercept=25000, color="red")+
    geom_hline(yintercept=  mean(c(Res_df$upr, Res_df$lwr)), linetype=2 )+
    geom_linerange(aes(ymin = lwr, ymax = upr))+
    #facet_wrap(~NAME, labeller = as_labeller(Strip_names))+
    theme_classic()+
    ylab("Coefficient\nestimate\n")+
    xlab("")

p_r2<-ggplot(Res_df, aes(x=aggFactor))+
    #geom_hline(yintercept=0)+
    #geom_vline(xintercept=25000, color="red")+
    geom_hline(yintercept=  mean(Res_df$R2), linetype=2 )+
    geom_line(aes(y=R2))+
    #facet_wrap(~NAME, labeller = as_labeller(Strip_names))+
    theme_classic()+
    ylab("Proportion\nvariance\nexplained")+
    xlab("Pixel size (km)")






#make maps
#functions
makeDataFigure<-function(plotDat, fillColumn){
    p<-ggplot() +
        geom_spatraster(data = plotDat, aes(fill=.data[[fillColumn]])) +
        geom_spatvector(data=crop(countries, ext(plotDat)), fill=NA, color="black") +
        theme_classic() + 
        viridis::scale_fill_viridis(na.value="transparent")
    return(p)
}

MapBenitez<-makeDataFigure(trim(crop(FRIPandMetrics_l[[5]], basins[2])), "DI")+ labs(fill='Benitez-Lopez\net al.\nDefaunation\nIndex') 

MapBogoni<-makeDataFigure(trim(crop(FRIPandMetrics_l[[5]], basins[2])), "DI_Bogoni")+ labs(fill='Bogoni et al.\nDefaunation\nIndex') 


#make combined figures for state
toprow<-cowplot::plot_grid(MapBenitez, MapBogoni, labels = c('A', 'B'), label_size = 12, ncol=2, rel_widths=c(1.05, 1))
all<-cowplot::plot_grid(
    toprow, 
    p_coef+geom_vline(xintercept=25, color="red", linetype=2), 
    p_r2+geom_vline(xintercept=25, color="red", linetype=2), 
    labels = c('', 'C', "D"), label_size = 12, ncol=1, rel_heights=c(1, 0.3, 0.3))

ggsave("Figures/FigureS9.png", height=10, width=18, bg="white")
