library(tidyverse)
continents<-data.frame(
    ISO3=c(  'BRA', "COL", "ECU", "GUF", "PER", "VEN", "BOL", "GUY", "SUR",
        "AGO", "CAF", "CMR", "COD", "COG", "GAB", "COG;CMR;CAF" ),
    continent= c(rep("SouthAmerica", 9), rep("Africa", 7)))


data<-read.csv("Downloads/PaAndMeans(3).csv") %>% left_join(continents)

d<- data %>% filter(GIS_AREA>1000) %>% filter(!is.na(mean))

ggplot(data=d, aes(y=mean, x=GIS_AREA, color=IUCN_CAT))+
    geom_point()+
    geom_smooth(method="lm")+
    facet_wrap(~IUCN_CAT)

library(lme4)
lmerTest::lmer(mean~GIS_AREA+(1|ISO3), data=d) %>% summary

table(d$IUCN_CAT)

#############


d<- data %>% filter(GIS_AREA>100) %>% filter(!is.na(mean)) %>% filter(IUCN_CAT %in% c("Ia", "Ib", "II"))

ggplot(data=d, aes(y=mean, x=GIS_AREA, color=continent))+
    geom_point()+
    geom_smooth(method="lm")

lmerTest::lmer(mean~GIS_AREA+(1|co), data=d) %>% summary

lm(mean~continent*GIS_AREA, data=d) %>% summary
ggplot(data=d, aes(y=mean, x=GIS_AREA, color=ISO3))+
    geom_point()+
    facet_wrap(~ISO3)


d %>% filter(GIS_AREA>30000)
d %>% filter(continent=="Africa") %>% filter(GIS_AREA>10000)




states<-read.csv("Downloads/StatesAndMeans(1).csv")



s<- states %>% filter(!is.na(mean)) 
largeCountries<-names(table(s$ADM0_NAME))[table(s$ADM0_NAME)>=3]
s<-s%>% filter(ADM0_NAME %in% largeCountries)


ggplot(data=s, aes(x=mean, y=ADM0_NAME))+
    geom_point()

s %>% filter(ADM0_NAME=="Brazil") %>% arrange(desc(mean))