##################
# part 1
##################

library(terra)
library(tidyverse)
library(tidyterra)

# load data
# 1) load data 
FRIPandMetrics_l<-lapply(list.files("Outputs/FRIP", full.names=TRUE), function(x){rast(x)})
preds_l<-lapply(list.files("Outputs/Predictors", full.names=TRUE), function(x){rast(x)})
basins<-vect("Outputs/Basins")
countries<-vect("Outputs/BasinCountries")



# get list of combined rasters at each resolution
FullStackMain_l<-FullStackRCV_l<-FullStackSCV_l<-list()
for (i in 1:length(FRIPandMetrics_l)){ 
    FullStackMain_l[[i]]<-c(FRIPandMetrics_l[[i]], preds_l[[i]]) 
    }

# for each raster split into tiles
    # Make a 3x3 tile “raster” over each basin
    tile_r1 <- rast( nrow=3, ncol=3, ext=ext(basins[1]))
    tile_r2 <- rast( nrow=3, ncol=3, ext=ext(basins[2]))
    tiles1 <- as.polygons(tile_r1, dissolve = FALSE)
    tiles2 <- as.polygons(tile_r2, dissolve = FALSE)
    tiles<-vect(c(tiles1, tiles2))

tvdf_l<-list()
for (i in 1:length(FullStackMain_l)){ 
    # your SpatRaster
    r <- FullStackMain_l[[i]]
    tvdf_l[[i]] <- terra::extract(r, tiles, ID=TRUE) %>%  select(-DI_Bogoni) %>% filter(complete.cases(.)) %>% as_tibble
    tvdf_l[[i]] %>% group_by(ID) %>% summarise(n()) %>% print(n=99)
    }


# make dummy with no signal to confirm OSC is real evidence
#tvdf_l<-lapply(tvdf_l, function(item){ item$correlation<-rnorm(dim(item)[1], 0, 1) ; return(item)} )



# train residual predicting models on
mod_df<-data.frame(r2=NA, r2orig=NA, res=NA, basin=NA, r2diff=NA)
#param_df<-data.frame(id=NA, res=NA,  int=NA,  elevation_mean=NA, sur_refl_b02_mean=NA, #depth_mean=NA, hnd_mean=NA, hnd_variance=NA, depth_variance=NA)
param_df<-data.frame(id=NA, res=NA,  int=NA,ndvi_mean=NA, hnd_mean=NA)#, #depth_mean=NA, hnd_mean=NA, hnd_variance=NA, depth_variance=NA)

 fit_df<-data.frame(id=NA, res=NA,  int=NA, denoised=NA, p=NA, cor=NA)

row<-1
for (i in c(1, 3:20)){ 
    for (thisID in unique(tvdf_l[[i]]$ID)) {
        train_df<-tvdf_l[[i]] %>% filter(ID == thisID)
        validate_df<-tvdf_l[[i]] %>% filter(ID != thisID)


        train_df$resid<-lm(DIlarge~correlation, train_df) %>%resid
        denoisingModel<-lm(resid~ndvi_mean+hnd_mean, train_df)
        #denoisingModel<-lm(resid~elevation_mean+sur_refl_b02_mean+hnd_mean+ hnd_variance+ depth_variance, train_df)

        validate_df$denoised<-validate_df$correlation-predict( denoisingModel, validate_df)

        m1<-lm(DIlarge~denoised, validate_df)
        m2<-lm(DIlarge~correlation, validate_df)
   
        mod_df[row,]<-c(as.numeric(summary(m1)$r.squared), as.numeric(summary(m2)$r.squared), varnames(FRIPandMetrics_l[[i]]),thisID , as.numeric(summary(m1)$r.squared)- as.numeric(summary(m2)$r.squared))
        param_df[row,1:5]<-c(row, varnames(FRIPandMetrics_l[[i]]), denoisingModel$coef)
        fit_df[row,1:6]<-c(row, varnames(FRIPandMetrics_l[[i]]), m1$coef, summary(m1)$coef[2,4], lm(validate_df$DIlarge~ validate_df$denoised)$coef[2])

        row<-row+1
    }
}

p1<-mod_df %>%
    as_tibble  %>%
    mutate(r2=as.numeric(r2)) %>%
    mutate(r2diff=as.numeric(r2diff)) %>%
    mutate(r2orig=as.numeric(r2orig)) %>%
    mutate(res=as.numeric(substr(res, 6,20))/1000) %>%
    mutate(r2ratio=r2/r2orig) %>%
    ggplot(., aes(x=res, y=(r2)))+
        geom_boxplot(aes(group=res))+
        geom_hline(yintercept=0, linetype=2)+
        ylab("R2 of denoised FRIP") + 
        xlab("Pixel size (km)")+
        theme_classic()


p2<-mod_df %>%
    as_tibble  %>%
    mutate(r2=as.numeric(r2)) %>%
    mutate(r2diff=as.numeric(r2diff)) %>%
    mutate(r2orig=as.numeric(r2orig)) %>%
    mutate(res=as.numeric(substr(res, 6,20))/1000) %>%
    mutate(r2ratio=r2/r2orig) %>%
    ggplot(., aes(x=res, y=log10(r2ratio), group=res))+
        geom_boxplot( outliers=FALSE)+
        geom_jitter(width=0.3, alpha=0.5, size=2)+
        geom_hline(yintercept=0, linetype=2)+
        ylab("Log ratio denoised R2:original R2") + 
        xlab("Pixel size (km)")+
        theme_classic()

Strip_names <- c(
                    `ndvi_mean` = "NDVI",
                    `hnd_mean` = "Height above nearest drainage")

p3<-param_df %>%
    as_tibble  %>%
    select(-int) %>%
    mutate(res=as.numeric(substr(res, 6,20))/1000) %>%
    pivot_longer(-c(id, res), names_to="name", values_to="value") %>%
    mutate(value=as.numeric(value)) %>%
    ggplot(., aes(y=(value), x=res))+
        geom_boxplot(aes(group=res))+
        geom_hline(yintercept=0, linetype=2)+
        facet_wrap(~name, scales="free", ncol=1, labeller = as_labeller(Strip_names))+
        ylab("Variable coefficient") + 
        xlab("Pixel size (km)")+
        theme_classic()


p4<-fit_df %>%
    as_tibble  %>%
    select(-int) %>%
    mutate(res=as.numeric(substr(res, 6,20))/1000) %>%
    mutate(denoised=as.numeric(denoised)) %>%
    mutate(p=as.numeric(p)) %>%
    mutate(cor=as.numeric(cor)) %>%
    ggplot(., aes(y=(cor), x=res))+
        geom_boxplot(aes(group=res))+
        geom_hline(yintercept=0, linetype=2)+
        ylab("Denoised FRIP - Defaunation index correlation") + 
        xlab("Pixel size (km)")+
        theme_classic()


p2
ggsave("Figures/Figure3.png", height=5, width=10, bg="white")



cowplot::plot_grid(
    p1, p3,
    labels = c("A", "B"), 
    label_size = 12, ncol=1, rel_widths = c(0.7, 1))
ggsave("Figures/FigureS4.png", height=10, width=13, bg="white")








################################
# END OF SCRIPT
###############################


























#########################################
# basin test train split
#########################################


basintvdf_l<-list()
for (i in 1:length(FullStackMain_l)){ 
    # your SpatRaster
    r <- FullStackMain_l[[i]]
    basintvdf_l[[i]] <- terra::extract(r, basins, ID=TRUE) %>%  select(-DI_Bogoni) %>% filter(complete.cases(.)) %>% as_tibble
    basintvdf_l[[i]] %>% group_by(ID) %>% summarise(n()) %>% print(n=99)
    }


mod_df<-data.frame(r2=NA, r2orig=NA, res=NA, basin=NA, r2diff=NA)
#param_df<-data.frame(id=NA, res=NA,  int=NA,  elevation_mean=NA, sur_refl_b02_mean=NA, #depth_mean=NA,  hnd_mean=NA, hnd_variance=NA, depth_variance=NA)
param_df<-data.frame(id=NA, res=NA,  int=NA, hnd_mean=NA,ndvi_mean=NA,  hnd_variance=NA, ndvi_variance=NA)#, #depth_mean=NA, hnd_mean=NA, hnd_variance=NA, depth_variance=NA)
 fit_df<-data.frame(id=NA, res=NA,  int=NA, denoised=NA, p=NA, cor=NA)

row<-1
for (i in 1:length(FullStackMain_l)){ 
    for (thisID in unique(basintvdf_l[[i]]$ID)) {
        train_df<-basintvdf_l[[i]] %>% filter(ID == thisID)
        validate_df<-basintvdf_l[[i]] %>% filter(ID != thisID)


        train_df$resid<-lm(DIlarge~correlation, train_df) %>%resid
        denoisingModel<-lm(resid~ hnd_mean+ndvi_mean+hnd_variance+ndvi_variance, train_df)
       # denoisingModel<-lm(resid~elevation_mean+sur_refl_b02_mean+hnd_mean+ hnd_variance+ depth_variance, train_df)

        validate_df$denoised<-validate_df$correlation-predict( denoisingModel, validate_df)

        m1<-lm(DIlarge~denoised, validate_df)
        m2<-lm(DIlarge~correlation, validate_df)
   
        mod_df[row,]<-c(as.numeric(summary(m1)$r.squared), as.numeric(summary(m2)$r.squared), varnames(FRIPandMetrics_l[[i]]),thisID , as.numeric(summary(m1)$r.squared)- as.numeric(summary(m2)$r.squared))
        param_df[row,1:7]<-c(row, varnames(FRIPandMetrics_l[[i]]), denoisingModel$coef)
        fit_df[row,1:6]<-c(row, varnames(FRIPandMetrics_l[[i]]), m1$coef, summary(m1)$coef[2,4], cor(validate_df$DIlarge, validate_df$denoised))


        row<-row+1
    }
}

p1<-mod_df %>%
    as_tibble  %>%
    mutate(r2=as.numeric(r2)) %>%
    mutate(r2diff=as.numeric(r2diff)) %>%
    mutate(r2orig=as.numeric(r2orig)) %>%
    mutate(res=as.numeric(substr(res, 6,20))/1000) %>%
    mutate(r2ratio=r2/r2orig) %>%
    ggplot(., aes(x=res, y=(r2)))+
        geom_point(aes(group=res, color=basin))+
        geom_hline(yintercept=0, linetype=2)+
        ylab("R2 denoised FRIP") + 
        xlab("Pixel resolution (km)")+
        theme_classic()


p2<-mod_df %>%
    as_tibble  %>%
    mutate(r2=as.numeric(r2)) %>%
    mutate(r2diff=as.numeric(r2diff)) %>%
    mutate(r2orig=as.numeric(r2orig)) %>%
    mutate(res=as.numeric(substr(res, 6,20))/1000) %>%
    mutate(r2ratio=r2/r2orig) %>%
    ggplot(., aes(x=res, y=log10(r2ratio)))+
        geom_point(aes(group=res, color=basin))+
        geom_hline(yintercept=0, linetype=2)+
        ylab("log10 ratio of denoised to original R2") + 
        xlab("Pixel resolution (km)")+
        theme_classic()


p3<-param_df %>%
    as_tibble  %>%
    select(-int) %>%
    mutate(res=as.numeric(substr(res, 6,20))/1000) %>%
    pivot_longer(-c(id, res), names_to="name", values_to="value") %>%
    mutate(value=as.numeric(value)) %>%
    ggplot(., aes(y=(value), x=res))+
        geom_point(aes(group=res))+
        geom_hline(yintercept=0, linetype=2)+
        facet_wrap(~name, scales="free", ncol=1)+ 
        ylab("Variable coefficient") + 
        xlab("Pixel resolution (km)")+
        theme_classic()


p4<-fit_df %>%
    as_tibble  %>%
    select(-int) %>%
    mutate(res=as.numeric(substr(res, 6,20))/1000) %>%
    mutate(denoised=as.numeric(denoised)) %>%
    mutate(p=as.numeric(p)) %>%
    mutate(cor=as.numeric(cor)) %>%
    ggplot(., aes(y=(cor), x=res))+
        geom_point(aes(group=res))+
        geom_hline(yintercept=0, linetype=2)+
        ylab("Denoised FRIP - Defaunation index correlation") + 
        xlab("Pixel resolution (km)")+
        theme_classic()



cowplot::plot_grid(
    p1, p2, p3, p4,
    labels = c("A", "B",'C', "D"), 
    label_size = 12, ncol=1, rel_widths = c(1, 1, 1.5))


ggsave("Figures/DenoiseDemoFigure.png", height=20, width=13, bg="white")




fit_df %>%
    as_tibble  %>%
    select(-int) %>%
    mutate(res=as.numeric(substr(res, 6,20))) %>%
    mutate(denoised=as.numeric(denoised)) %>%
    mutate(p=as.numeric(p)) %>%
    mutate(cor=as.numeric(cor)) %>%
    ggplot(., aes(x=(cor), y=res))+
        geom_boxplot(aes(group=res)) +
        geom_vline(xintercept=0)


    #full raster (dor application at other resolutions)
    # RCV train dfs
    # SCV train dfs (both basins)

# apply denoising (predicted residual subtraction) to:
    #each raster using other resolutions
    #RCV validation data setes 
    #SCV validation data sets (both basins)

# save:
    # denosied FRIP values as raster layers
r<-        FullStackMain_l[[11]]

r$denoised<-r$correlation-terra::predict(r, denoisingModel)
plot(c(r$denoised, mask(1-r$DIlarge, r$denoised)))
    ggplot() +
        geom_spatraster(data=r, aes(fill=denoised)) +
        geom_spatvector(data=crop(countries, ext(r)), fill=NA, color="black") +
        theme_classic()+
        scale_fill_gradient2(na.value = "transparent", high=scales::muted("red"), low=scales::muted("blue"), mid="light grey", name="Flooding role\nin productivity")#, limits = c(-1,1))



##################
# part 2
##################

#calculate r2 and make matrix figure

# make supplementary signal strengthening figure and 





