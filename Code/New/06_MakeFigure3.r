#plotting for bcp plots

df<-df %>% filter(name %in% names(which(table(df$name)>50)))

myPalette = c(
    'b' =   "#AD002AFF",
    'ab' = "#925E9FFF", 
    'abc' = "#ADB6B6FF", 
    'acd' = "#42B540FF", 
    'cd' = "#00468BFF")#
    #'abd' = "#ED0000FF",
    #'c' = "#FDAF91FF")

fig_bcp<-list()
fig_bcp[[1]]<-ggplot(data=df, aes(y=reorder(name, correlation, "median"), x=correlation, color=Letters))+
    geom_vline(xintercept=0, linetype=2, alpha=0.5)+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.6), alpha=1, size=1)+
    geom_boxplot(outliers=FALSE, fill=NA, color="black")+
    theme_classic()+
    guides(color="none")+
    scale_fill_manual(values=myPalette)+
    scale_color_manual(values=myPalette)+
    xlab("Flooding role in productivity")+ylab("")

countriesForPlotting<-countries
values(countriesForPlotting)<-values(countries) %>%
    as_tibble %>%
    left_join(., LABELS, join_by(name==treatment))


plotDat<-crop(FRIP, basins[1])
fig_bcp[[2]]<-ggplot() +
    geom_spatvector(data=crop(countriesForPlotting, ext(plotDat)), aes(fill=Letters), linewidth=1, color="black")+
    scale_fill_manual(values=myPalette, na.value = "transparent")+
    guides(fill="none")+
    theme_classic()


plotDat<-crop(FRIP, basins[2])
fig_bcp[[3]]<-ggplot() +
    geom_spatvector(data=crop(countriesForPlotting, ext(plotDat)), aes(fill=Letters), linewidth=1, color="black")+
    scale_fill_manual(values=myPalette, na.value = "transparent")+
    guides(fill="none")+
    theme_classic()

saveRDS(fig_bcp, "Outputs/fig_bcp.RDS")
