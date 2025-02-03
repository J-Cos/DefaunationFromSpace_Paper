# explore no spatial autocorrelation control PA-wise trends
library(terra)
library(tidyverse)
library(multcompView)
library(tidyterra)

source("Code/Functions.r")


# 1) load data 
FRIP<-rast("Outputs/FRIP.tif") 
FRIPtrend<-rast("Outputs/FRIPtrend.tif") 
basins<-vect("Outputs/Basins")
PAs<-vect("Outputs/BasinPAs")

# 2) fit models across aggregation factors ----------------------------------------------------------------


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
for (aggregationFactor in 1:50) {

    FRIPagg<-aggregate(FRIP, aggregationFactor, na.rm=TRUE)
    FRIPtrendagg<-aggregate(FRIPtrend, aggregationFactor, na.rm=TRUE)

    PAs_rast<-rasteriseAndMask(PAs_Focal, FRIPagg, "NAME")

    r_list[[aggregationFactor]]<-c(FRIPagg, FRIPtrendagg, PAs_rast)

    print(paste0("Completed all rasterisation at aggregation factor ", aggregationFactor))


    # 4) get modelling Dataframe
    df<-as.data.frame( r_list[[aggregationFactor]], cells=TRUE, na.rm=TRUE)

    df_list[[aggregationFactor]]<-df %>%
        as_tibble %>%
       # filter (NAME%in% names(which(table(df$NAME)>5))) %>% 
        mutate(NAME = as.factor(str_replace_all(NAME, "-", ":")))

    mods_list[[aggregationFactor]]<-aov(correlation~NAME, df_list[[aggregationFactor]])
    modstrend_list[[aggregationFactor]]<-lm(tau~0+NAME, df_list[[aggregationFactor]])

    print(paste0("Completed at aggregation factor ", aggregationFactor))
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
        aggFactor=(as.numeric(aggFactor)*5) ^2   )  %>%
    filter(Comparison %in% c(
        "Salonga National Park-Ngiri"  ,
        "Salonga National Park-Grands affluents"   ,
        "Yanomami-Menkragnotí"   
    ))


#FRIP
Restrend_df<-lapply(modstrend_list, function(x){rownames_to_column(as.data.frame( confint(x) ), "NAME")}) %>% 
    bind_rows( .id="aggFactor") %>% 
    as_tibble %>% 
    mutate(
        aggFactor=(5*as.numeric(aggFactor))^2,
        lwr=`2.5 %`,
        upr=`97.5 %`) %>%
    mutate(signif=lwr>0 | upr<0) %>%
    mutate(NAME=str_replace_all(NAME, "NAME", ""))



# 4) print out results for text --------------------------------------------------------------------
TukeyHSD(mods_list[[10]])

MaxAggFactorOfFinding<-which(!filter(Res_df, Comparison== "Salonga National Park-Ngiri" )$signif)[1]-1
(MaxAggFactorOfFinding*5)^2

MaxAggFactorOfFinding<-which(!filter(Res_df, Comparison=="Salonga National Park-Grands affluents" )$signif)[1]-1
(MaxAggFactorOfFinding*5)^2

MaxAggFactorOfFinding<-which(!filter(Res_df, Comparison==    "Yanomami-Menkragnotí"  )$signif)[1]-1
(MaxAggFactorOfFinding*5)^2

summary(modstrend_list[[10]])
confint(modstrend_list[[10]])


MaxAggFactorOfFinding<-which(!filter(Restrend_df, NAME== "Salonga National Park" )$signif)[1]-1
(sum(filter(Restrend_df, NAME== "Salonga National Park" )$signif)*5)^2 # have to sum as all are significant

MaxAggFactorOfFinding<-which(!filter(Restrend_df, NAME== "Menkragnotí" )$signif)[1]-1
(MaxAggFactorOfFinding*5)^2

MaxAggFactorOfFinding<-which(!filter(Restrend_df, NAME== "Alto Rio Negro" )$signif)[1]-1
(MaxAggFactorOfFinding*5)^2

MaxAggFactorOfFinding<-which(!filter(Restrend_df, NAME== "Rio Negro" )$signif)[1]-1
(MaxAggFactorOfFinding*5)^2



# 5) get dataframe for main figure ----------------------------------------------------------------------------

LABELS <- generate_label_df(TukeyHSD(mods_list[[10]]) , "NAME")
df<-df_list[[10]] %>% mutate(treatment=NAME) %>% left_join(., LABELS)

trendSignif_df<-as.data.frame(summary(modstrend_list[[10]])$coef[,c(1,4)]) %>%
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


figS4<-ggplot(Res_df_wMeans, aes(x=aggFactor))+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=(5*10)^2, color="red")+
    geom_hline(aes(yintercept=  mean), linetype=2 )+
    geom_linerange(aes(ymin = lwr, ymax = upr, alpha=signif))+
    facet_wrap(~Comparison, ncol=1)+#, labeller = as_labeller(Strip_names))+
    theme_classic()+
    guides(alpha=FALSE)+
    ylab("Coefficient estimate")+
    xlab("Pixel size (km2)")

figS4
ggsave("Figures/FigureS4.png")

meanLinestrend<- Restrend_df %>%
    pivot_longer(c(lwr, upr)) %>%
    group_by(NAME) %>%
    summarise(mean=mean(value, na.rm=TRUE))

Restrend_df_wMeans <-left_join(Restrend_df, meanLinestrend)

figS5<-ggplot(Restrend_df_wMeans, aes(x=aggFactor))+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=(5*10)^2, color="red")+
    geom_hline(aes(yintercept=  mean), linetype=2 )+
    geom_linerange(aes(ymin = lwr, ymax = upr, alpha=signif))+
    facet_wrap(~NAME)+#, labeller = as_labeller(Strip_names))+
    theme_classic()+
    guides(alpha=FALSE)+
    ylab("Coefficient estimate")+
    xlab("Pixel size (km2)")

figS5
ggsave("Figures/FigureS5.png", height=10, width=10)
