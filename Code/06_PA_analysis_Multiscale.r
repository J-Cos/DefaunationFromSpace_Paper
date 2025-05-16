# explore no spatial autocorrelation control PA-wise trends
library(terra)
library(tidyverse)
library(multcompView)
library(tidyterra)

source("Code/Functions.r")


# 1) load data 
FRIP_l<-lapply(list.files("Outputs/FRIP", full.names=TRUE), function(x){rast(x)$correlation})
FRIPtrend_l<-lapply(list.files("Outputs/FRIPtrend", full.names=TRUE), function(x){rast(x)})
basins<-vect("Outputs/Basins")
PAs<-vect("Outputs/BasinPAs")

# 2) fit models across aggregation factors ----------------------------------------------------------------
aggFactors<-unlist(lapply(FRIP_l, function(x){ as.character(str_sub(varnames(x)[1], 6) )}))

#basin raster
PAs_amazon<- crop(PAs, basins[2])
PAs_amazon_Focal<-PAs_amazon[expanse(PAs_amazon, unit="km")>40000]# & !PAs_reprojected$IUCN_CAT %in% c("Not Reported", "Not Applicable") ]
PAs_amazon_Focal<-PAs_amazon_Focal[PAs_amazon_Focal$NAME!="Central Amazon Conservation Complex"]

PAs_congo<- crop(PAs, basins[1])
PAs_congo_Focal<-PAs_congo[expanse(PAs_congo, unit="km")>20000]# & !PAs_reprojected$IUCN_CAT %in% c("Not Reported", "Not Applicable") ]
PAs_congo_Focal<-PAs_congo_Focal[PAs_congo_Focal$NAME!="Salonga"]

PAs_Focal<-vect(c(PAs_amazon_Focal, PAs_congo_Focal))


r_list<-list()
df_list<-list()
mods_list<-modstrend_list<-list()
for (i in 1:length(aggFactors)) {

    FRIPagg<-FRIP_l[[i]]
    FRIPtrendagg<-FRIPtrend_l[[i]]

    PAs_rast<-rasteriseAndMask(PAs_Focal, FRIPagg, "NAME")

    r_list[[aggFactors[i]]]<-c(FRIPagg, FRIPtrendagg, PAs_rast)

    print(paste0("Completed all rasterisation at aggregation factor ", aggFactors[i]))


    # 4) get modelling Dataframe
    df<-as.data.frame( r_list[[aggFactors[i]]], cells=TRUE, na.rm=TRUE)

    df_list[[aggFactors[i]]]<-df %>%
        as_tibble %>%
       # filter (NAME%in% names(which(table(df$NAME)>5))) %>% 
        mutate(NAME = as.factor(str_replace_all(NAME, "-", ":")))

    mods_list[[aggFactors[i]]]<-aov(correlation~NAME, df_list[[aggFactors[i]]])
    modstrend_list[[aggFactors[i]]]<-lm(tau~0+NAME, df_list[[aggFactors[i]]])

    print(paste0("Completed at aggregation factor ", aggFactors[i]))
}


# 3) get results dataframe #--------------------------------------------------------------------
# FRIPtrends
Res_df<-lapply(mods_list, function(x){ TukeyHSD(x) %>% 
                                            lapply(., as.data.frame) %>%
                                            bind_rows() %>%
                                            rownames_to_column("Comparison") %>%
                                            as_tibble }) %>%
    bind_rows( .id="aggFactor") %>% 
    as_tibble %>%    
    mutate( 
        signif=`p adj`<0.05,
        aggFactor=as.numeric(aggFactor))  %>%
    filter(Comparison %in% c(
        "Salonga National Park-Ngiri"  ,
        "Salonga National Park-Grands affluents"   ,
        "Yanomami-Menkragnotí"   
    )) %>%
    arrange((aggFactor))


#FRIP
Restrend_df<-lapply(modstrend_list, function(x){rownames_to_column(as.data.frame( confint(x) ), "NAME")}) %>% 
    bind_rows( .id="aggFactor") %>% 
    as_tibble %>% 
    mutate(
        aggFactor=as.numeric(aggFactor),
        lwr=`2.5 %`,
        upr=`97.5 %`) %>%
    mutate(signif=lwr>0 | upr<0) %>%
    mutate(NAME=str_replace_all(NAME, "NAME", ""))%>%
    arrange((aggFactor))


################
# set results resolution
################
res<-"25000"

# 4) print out results for text --------------------------------------------------------------------
TukeyHSD(mods_list[[res]])
TukeyHSD(mods_list[[res]])$NAME[c(47,63, 65),]


MaxAggFactorOfFinding<-which(!filter(Res_df, Comparison== "Salonga National Park-Ngiri" )$signif)[1]-1
(MaxAggFactorOfFinding*5)^2

MaxAggFactorOfFinding<-which(!filter(Res_df, Comparison=="Salonga National Park-Grands affluents" )$signif)[1]-1
(MaxAggFactorOfFinding*5)^2

MaxAggFactorOfFinding<-which(!filter(Res_df, Comparison==    "Yanomami-Menkragnotí"  )$signif)[1]-1
(MaxAggFactorOfFinding*5)^2

summary(modstrend_list[[res]])
confint(modstrend_list[[res]])


MaxAggFactorOfFinding<-which(!filter(Restrend_df, NAME== "Salonga National Park" )$signif)[1]-1
(sum(filter(Restrend_df, NAME== "Salonga National Park" )$signif)*5)^2 # have to sum as all are significant

MaxAggFactorOfFinding<-which(!filter(Restrend_df, NAME== "Menkragnotí" )$signif)[1]-1
(MaxAggFactorOfFinding*5)^2

MaxAggFactorOfFinding<-which(!filter(Restrend_df, NAME== "Alto Rio Negro" )$signif)[1]-1
(MaxAggFactorOfFinding*5)^2

MaxAggFactorOfFinding<-which(!filter(Restrend_df, NAME== "Rio Negro" )$signif)[1]-1
(MaxAggFactorOfFinding*5)^2



# 5) get dataframe for main figure ----------------------------------------------------------------------------

LABELS <- generate_label_df(TukeyHSD(mods_list[[res]]) , "NAME")
df<-df_list[[res]] %>% mutate(treatment=NAME) %>% left_join(., LABELS)

trendSignif_df<-as.data.frame(summary(modstrend_list[[res]])$coef[,c(1,4)]) %>%
    rownames_to_column("NAME") %>%
    mutate(
        NAME=str_replace_all(NAME, "NAME", ""),
        signif=`Pr(>|t|)`<0.05) %>%
    mutate(change = case_when(
        !signif ~ "None",
        signif & Estimate>0     ~ "Increasing",
        signif & Estimate<0     ~ "Decreasing"
  )) %>%
  select(NAME, change)

df_wTrendSignifs <-left_join(df, trendSignif_df)

saveRDS(df_wTrendSignifs, "Outputs/PlottingData_PAs.RDS")

# 6) make supplementary figure ------------------------------------------------------
meanLines<- Res_df %>%
    pivot_longer(c(lwr, upr)) %>%
    group_by(Comparison) %>%
    summarise(mean=mean(value, na.rm=TRUE))

Res_df_wMeans <-left_join(Res_df, meanLines)


figS4<-ggplot(Res_df_wMeans, aes(x=aggFactor/1000))+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=(25), color="red", linetype=2)+
    geom_hline(aes(yintercept=  mean), linetype=2 )+
    geom_linerange(aes(ymin = lwr, ymax = upr, alpha=signif))+
    facet_wrap(~Comparison, ncol=1)+#, labeller = as_labeller(Strip_names))+
    theme_classic()+
    guides(alpha=FALSE)+
    ylab("Coefficient estimate")+
    xlab("Pixel size (km)")

figS4
ggsave("Figures/FigureS6.png", height=10, width=10)

meanLinestrend<- Restrend_df %>%
    pivot_longer(c(lwr, upr)) %>%
    group_by(NAME) %>%
    summarise(mean=mean(value, na.rm=TRUE))

Restrend_df_wMeans <-left_join(Restrend_df, meanLinestrend)

figS5<-ggplot(Restrend_df_wMeans, aes(x=aggFactor/1000))+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=(25), color="red", linetype=2)+
    geom_hline(aes(yintercept=  mean), linetype=2 )+
    geom_linerange(aes(ymin = lwr, ymax = upr, alpha=signif))+
    facet_wrap(~NAME, ncol=3)+#, labeller = as_labeller(Strip_names))+
    theme_classic()+
    guides(alpha=FALSE)+
    ylab("Coefficient estimate")+
    xlab("Pixel size (km)")

figS5
ggsave("Figures/FigureS8.png", height=10, width=12)
