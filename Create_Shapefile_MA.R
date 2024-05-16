
## Aggregates births and heat data at the ZCTA level and creates a shapefile for mapping

library(dplyr)
library(sf)
library(arrow)

shp <- read_sf("Data/shp/tl_2020_us_zcta510.shp")%>%
  mutate(zip=as.numeric(ZCTA5CE10))%>%
  filter(zip >= 27006 & zip <= 28909)

data <- read_parquet("Data/Deliv_Heatwaves.parquet")

datazip <- data %>%
  group_by(zip)%>%
  summarize(n=n(),
            hwct=sum(heatwave),
            smmct=sum(smm21),
            gdmct=sum(gdm),
            hdpct=sum(hdp),
            mdpct=sum(mdp),
            ptbct=sum(ptb),
            smict=sum(smi),
            pmadct=sum(pmad),
            subct=sum(sub),
            smmrt=(smmct/n),
            gdmrt=(gdmct/n),
            hdprt=(hdpct/n),
            mdprt=(mdpct/n),
            ptbrt=(ptbct/n),
            smirt=(smict/n),
            pmadrt=(pmadct/n),
            subrt=(subct)/n)

shp <- shp %>%
  as.data.frame()%>%
  select(-ZCTA5CE10, -GEOID10, -CLASSFP10, -MTFCC10, -FUNCSTAT10, -ALAND10, -AWATER10, -INTPTLAT10, -INTPTLON10)

sf <- shp %>%
  left_join(datazip, by=c('zip'))

write_sf(sf, "Data/Outcomes_HW_ZCTA.shp")

