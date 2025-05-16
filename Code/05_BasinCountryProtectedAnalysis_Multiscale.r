# explore no spatial autocorrelation control basin/country/protected-wise trends
library(terra)
library(tidyverse)
library(multcompView)
library(tidyterra)

source("Code/Functions.r")

# 1) load data 
FRIP_l<-lapply(list.files("Outputs/FRIP", full.names=TRUE), function(x){rast(x)$correlation})
basins<-vect("Outputs/Basins")
countries<-vect("Outputs/BasinCountries")
PAs<-vect("Outputs/BasinPAs")

# 2) fit models across aggregation factors ----------------------------------------------------------------
aggFactors<-unlist(lapply(FRIP_l, function(x){ as.character(str_sub(varnames(x)[1], 6) )}))

#aggregate PAs first for substantially increaased efficiency
PAs_agg<-aggregate(PAs, dissolve=TRUE)
PAs_agg$IUCN_CAT<-"All"

r_list<-list()
df_list<-list()
mods_list<-list()
for (i in 1:length(aggFactors)) {
    # make single raster
    # aggregate FRIP
    FRIPagg<-FRIP_l[[i]]

    #basin raster
    basin_rast<-rasteriseAndMask(basins, FRIPagg, "PFAF_ID")
    print(paste0("Completed basin rasterisation at aggregation factor ", aggFactors[i]))

    #country raster
    countries_rast<-rasteriseAndMask(countries, FRIPagg, "name")
    print(paste0("Completed country rasterisation at aggregation factor ", aggFactors[i]))

    # protected raster
    protected_rast<-rasteriseAndMask(PAs_agg, FRIPagg, "IUCN_CAT")
    protected_rast$protected<-!is.na(protected_rast)

    r_list[[aggFactors[i]]]<-c(FRIPagg, basin_rast, countries_rast, protected_rast$protected)

    print(paste0("Completed all rasterisation at aggregation factor ", aggFactors[i]))


    # get modelling Dataframe
    df<-as.data.frame(r_list[[aggFactors[i]]], cells=TRUE, na.rm=TRUE)

    df_list[[aggFactors[i]]]<-df %>%
        as_tibble %>%
        #filter (name %in% names(which(table(df$name)>50))) %>%
        mutate(
            protected=as.factor(protected),
            PFAF_ID=as.factor(PFAF_ID))

    mods_list[[aggFactors[i]]]<-aov(correlation~PFAF_ID+name+protected, df_list[[aggFactors[i]]])
    print(paste0("Completed at aggregation factor ", aggFactors[i]))
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
        aggFactor=as.numeric(aggFactor)) %>%
    arrange((aggFactor))

KeyComparison_df<-Res_df %>%
    filter(Comparison %in% 
        c(
            "62-13",
            "Peru-Brazil",
            "Gabon-Democratic Republic of the Congo",
            "TRUE-FALSE"
        )
        )


################
# set results resolution
################
res<-"50000"

# 4) print out results for text --------------------------------------------------------------------
TukeyHSD(mods_list[[res]])

KeyComparison_df %>% filter(aggFactor==res)

MaxAggFactorOfFinding<-which(!filter(Res_df, Comparison=="62-13")$signif)[1]-1

MaxAggFactorOfFinding<-which(!filter(Res_df, Comparison=="Peru-Brazil")$signif)[1]-1

MaxAggFactorOfFinding<-which(!filter(Res_df, Comparison=="Gabon-Democratic Republic of the Congo")$signif)[1]-1

MaxAggFactorOfFinding<-which(!filter(Res_df, Comparison=="TRUE-FALSE")$signif)[1]-1





# 5) get dataframe for main figure ----------------------------------------------------------------------------

LABELS <- generate_label_df(TukeyHSD(mods_list[[res]], "name") , "name")

df<-df_list[[res]] %>% mutate(treatment=name) %>% left_join(., LABELS)
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

figS<-ggplot(KeyComparison_df_wMeans, aes(x=aggFactor/1000))+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=(50), color="red", linetype=2)+
    geom_hline(aes(yintercept=  mean), linetype=2 )+
    geom_linerange(aes(ymin = lwr, ymax = upr, alpha=signif))+
    facet_wrap(~Comparison, labeller = as_labeller(Strip_names))+
    theme_classic()+
    guides(alpha=FALSE)+
    ylab("Coefficient estimate")+
    xlab("Pixel size (km)")

figS
ggsave("Figures/FigureS5.png", height=10, width=10)
