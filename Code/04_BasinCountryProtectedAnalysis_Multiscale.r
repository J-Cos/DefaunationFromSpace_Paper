# explore no spatial autocorrelation control basin/country/protected-wise trends
library(terra)
library(tidyverse)
library(multcompView)
library(tidyterra)

source("Code/Functions.r")

# 1) load data 
FRIP<-rast("Outputs/FRIP.tif") 
basins<-vect("Outputs/Basins")
countries<-vect("Outputs/BasinCountries")
PAs<-vect("Outputs/BasinPAs")

# 2) fit models across aggregation factors ----------------------------------------------------------------

#aggregate PAs first for substantially increaased efficiency
PAs_agg<-aggregate(PAs, dissolve=TRUE)
PAs_agg$IUCN_CAT<-"All"

r_list<-list()
df_list<-list()
mods_list<-list()
for (aggregationFactor in 1:50) {
    # make single raster
    # aggregate FRIP
    FRIPagg<-aggregate(FRIP, aggregationFactor, na.rm=TRUE)

    #basin raster
    basin_rast<-rasteriseAndMask(basins, FRIPagg, "PFAF_ID")
    print(paste0("Completed basin rasterisation at aggregation factor ", aggregationFactor))

    #country raster
    countries_rast<-rasteriseAndMask(countries, FRIPagg, "name")
    print(paste0("Completed country rasterisation at aggregation factor ", aggregationFactor))

    # protected raster
    protected_rast<-rasteriseAndMask(PAs_agg, FRIPagg, "IUCN_CAT")
    protected_rast$protected<-!is.na(protected_rast)

    r_list[[aggregationFactor]]<-c(FRIPagg, basin_rast, countries_rast, protected_rast$protected)

    print(paste0("Completed all rasterisation at aggregation factor ", aggregationFactor))


    # get modelling Dataframe
    df<-as.data.frame(r_list[[aggregationFactor]], cells=TRUE, na.rm=TRUE)

    df_list[[aggregationFactor]]<-df %>%
        as_tibble %>%
        #filter (name %in% names(which(table(df$name)>50))) %>%
        mutate(
            protected=as.factor(protected),
            PFAF_ID=as.factor(PFAF_ID))

    mods_list[[aggregationFactor]]<-aov(correlation~PFAF_ID+name+protected, df_list[[aggregationFactor]])
    print(paste0("Completed at aggregation factor ", aggregationFactor))
}


# 3) get results dataframe #--------------------------------------------------------------------
# FRIP
Res_df<-lapply(mods_list, function(x){ TukeyHSD(x) %>% 
                                            lapply(., as.data.frame) %>%
                                            bind_rows() %>%
                                            rownames_to_column("Comparison") %>%
                                            as_tibble }) %>%
    bind_rows( .id="aggFactor") %>% 
    as_tibble %>%    
    mutate( 
        signif=`p adj`<0.05,
        aggFactor=(as.numeric(aggFactor)*5) ^2   ) 

KeyComparison_df<-Res_df %>%
    filter(Comparison %in% 
        c(
            "62-13",
            "Peru-Brazil",
            "Gabon-Democratic Republic of the Congo",
            "TRUE-FALSE"
        )
        )


# 4) print out results for text --------------------------------------------------------------------
TukeyHSD(mods_list[[10]])

MaxAggFactorOfFinding<-which(!filter(Res_df, Comparison=="62-13")$signif)[1]-1
(MaxAggFactorOfFinding*5)^2

MaxAggFactorOfFinding<-which(!filter(Res_df, Comparison=="Peru-Brazil")$signif)[1]-1
(MaxAggFactorOfFinding*5)^2

MaxAggFactorOfFinding<-which(!filter(Res_df, Comparison=="Gabon-Democratic Republic of the Congo")$signif)[1]-1
(MaxAggFactorOfFinding*5)^2

MaxAggFactorOfFinding<-which(!filter(Res_df, Comparison=="TRUE-FALSE")$signif)[1]-1
(MaxAggFactorOfFinding*5)^2





# 5) get dataframe for main figure ----------------------------------------------------------------------------

LABELS <- generate_label_df(TukeyHSD(mods_list[[10]], "name") , "name")

df<-df_list[[10]] %>% mutate(treatment=name) %>% left_join(., LABELS)
saveRDS(df, "Outputs/PlottingData_bcp.RDS")

# 6) make supplementary figure
Strip_names <- c(
            "62-13" = "Amazon-Congo",
            "Peru-Brazil" =             "Peru-Brazil",
            "Gabon-Democratic Republic of the Congo" =             "Gabon-Democratic Republic of the Congo",
            "TRUE-FALSE" = "Protected-Unprotected"
                    )

meanLines<- KeyComparison_df %>%
    pivot_longer(c(lwr, upr)) %>%
    group_by(Comparison) %>%
    summarise(mean=mean(value))

KeyComparison_df_wMeans <-left_join(KeyComparison_df, meanLines)

figS<-ggplot(KeyComparison_df_wMeans, aes(x=aggFactor))+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=(5*10)^2, color="red")+
    geom_hline(aes(yintercept=  mean), linetype=2 )+
    geom_linerange(aes(ymin = lwr, ymax = upr, alpha=signif))+
    facet_wrap(~Comparison, labeller = as_labeller(Strip_names))+
    theme_classic()+
    guides(alpha=FALSE)+
    ylab("Coefficient estimate")+
    xlab("Pixel size (km2)")

figS
ggsave("Figures/FigureS3.png", height=10, width=10)
