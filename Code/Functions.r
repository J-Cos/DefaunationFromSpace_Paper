# functions

# turn spatvectors into rasters with any cells not 99% within a single polygon masked
# used for model tests
rasteriseAndMask<-function(v, rastToMatch, col){
    r_unmasked<-terra::rasterize(v, rastToMatch, field=col)
    r_coverMask<-(terra::rasterize(v, rastToMatch, cover=TRUE, by=col) >0.99) %>% any(na.rm=TRUE)
    r<-mask(r_unmasked, r_coverMask, maskvalues=FALSE)
    return(r)
}

# create significance indicating labels from TukeyHSD tests
# used to create dataframes for joint maps and box plots figure
generate_label_df <- function(TUKEY, variable){
 
     # Extract labels and factor levels from Tukey post-hoc 
     Tukey.levels <- TUKEY[[variable]][,4] 
     Tukey.levels<-Tukey.levels[!is.na(Tukey.levels)]
     Tukey.labels <- data.frame(multcompView::multcompLetters(x=Tukey.levels)['Letters'])
     
     #I need to put the labels in the same order as in the boxplot :
     Tukey.labels$treatment=rownames(Tukey.labels)
     Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
     return(Tukey.labels)
     }
