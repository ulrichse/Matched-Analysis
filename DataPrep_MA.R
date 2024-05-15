

## Use this code to create a data frame for use in the matched analysis  
## This code will create a data frame consisting of daily observations at the ZCTA level for 
## the selected time period (May-Sept, 2016-2019), exposure metric (heatwaves) and health outcome(s) (births, PTB, hypertension etc.)
## This code also includes a yearly ZCTA-level population metric. This example uses maternal population
## (female ages 18 to 44) for 2016-2019, set as the 'denom' variable. 

# The should have the following columns after running this code: 

## "zip": county zip code
## "date": date of each day
## "year": year for each day
## "doy": day of year, used to find similar non-heatwave days in other years
## "dow": day of week
## "heatwave": an indicator variable for heatwave exposure, "heatwave" is 1 
##          if a day is exposed to storm and 0 for non-heatwave exposed days (or
##          other binary exposure of interest)
## "denom": relevant denominator for the health outcome of interest; this example uses the total female population
##          between age 18 and 44 for each ZCTA and year
## "births: count of deliveries in a ZCTA on a given day in the study period
## "SMM21, GDM, HDP, SMI, PTB, MDP, PMAD, SUB": outcomes included in delivery data 


# Read in packages and set working directory

library(dplyr)
library(arrow)
library(lubridate)
library(data.table)
library(haven)
library(sjPlot)
library(MASS)

setwd("~/EHF Matched Analysis")
getwd()

#### DATA PREP ####

#HEATWAVES
## Read in Heatwave_Metrics.parquet file and restrict to warm months (May to September) of 2016 to 2019
## This file contains daily, ZCTA-level heatwave data. We use the Heatwave column in this analysis, which 
## indicates whether a day within a ZCTA was a heatwave or not (1/0), based on EHF calculations. 
## Further documentation of EHF heatwave data is available here: https://github.com/wertisml/Heatwave
## The 'Zip' column in the Heatwave_Metrics.parquet file represent ZCTAs, so we do not need to crosswalk. 

data <- read_parquet("Data/Heatwave_Metrics.parquet")%>%
  setDT()%>%
  filter(Zip >= 27006 & Zip <= 28909)%>%  # Filter for NC zip codes
  rename(date=Date, # Rename date column to lowercase
         heatwave=Heatwave, # Select Heatwave column as our heatwave variable
         zip=Zip)%>% # Rename zip code column to lowercase
  mutate(year=year(date),
         doy=yday(date),
         dow=wday(date))%>%
  dplyr::select(date, heatwave, year, doy, dow, zip)%>% # Select variables of interest
  mutate(heatwave=as.integer(heatwave)) %>% # Set as an integer just in case
  filter(date >= "2016-01-01" & date <= "2019-12-31")%>% # Filter for years of interest
  filter(month(date) >= 5 & month(date) <= 9) # Filter for months of interest

#DELIVERIES
## deliveries.parquet is individual delivery data with date and zip code. 
## Includes 0/1 indicators for GDM, HDP, PTB, SMI, PMAD, MDP, SMM21, and SUB for each delivery. 
## We will restrict deliveries to May to September of 2016 to 2019. 
## ptzip is the patient's zip code of residence. We will need to crosswalk to ZCTA values. 
## We will aggregate the deliveries to create daily, ZCTA-level counts for births and outcomes.

deliv <- read_parquet("Data/Deliveries_Outcomes.parquet")%>% # Read in deliveries file
  setDT()%>%
  mutate(date=as.Date(admitdt), # Create date column
         zip=as.numeric(ptzip))%>% # This is not in ZCTA form so we will have to crosswalk it below
  filter(zip >= 27006 & zip <= 28909)%>% # Filter for NC zip codes
  dplyr::select(zip, date, SMM21, GDM, HDP, SMI, PTB, MDP, PMAD, SUB)%>% # Select outcomes of interest
  filter(month(date) >= 5 & month(date) <= 9)%>% # Filter for warm months
  filter(year(date) >= 2016 & year(date) <= 2019) # Filter for years of interest

# Read in crosswalk file to convert from zip code to ZCTA
cw <- read.csv("Data/Crosswalk_NC.csv")%>%
  setDT()%>%
  rename(zip=ZIP)

# Crosswalk the delivery data 
deliv <- deliv %>%
  left_join(cw, by=c('zip'))%>%
  dplyr::select(-zip, ZIP_TYPE)%>%
  rename(zip=ZCTA) #The 'zip' variable in the deliv dataframe will now represent ZCTAs, not zip codes. There should be ~804 unique ZCTAs. 

# Aggregate deliveries and select outcomes at the daily ZCTA level
delivzip <- deliv %>%
  group_by(zip, date)%>%
  summarize(births=n(), 
            smm21=sum(SMM21),
            gdm=sum(GDM), 
            hdp=sum(HDP), 
            smi=sum(SMI),
            ptb=sum(PTB),
            mdp=sum(MDP),
            pmad=sum(PMAD),
            sub=sum(SUB))

# Read in maternal population file to use as denom variable in models
pop <- read.csv("Data/TotalPop_F_18_44_2016_2021.csv")%>%
  mutate(denom=as.numeric(denom))%>%
  filter(year >= 2016 & year <= 2019)

#### MERGE DATASETS ####

data <- data %>%
  left_join(delivzip, by=c('zip', 'date')) %>% # Merge delvieries by ZCTA and data
  left_join(pop, by=c('zip', 'year')) # Merge population dataframe

data[is.na(data)] <- 0 # Replace health outcome NA values with 0

# Examine the count of heatwave events for the study period for each ZCTA and remove ZCTAs without heatwave events
datagrp <- data %>%
  group_by(zip)%>%
  summarize(hw=sum(heatwave))
View(datagrp) # Visually inspect for 0s 

datagrp <- datagrp %>% # Create list of ZCTAs without heatwaves
  filter(hw==0)
no_hw <- datagrp$zip 

data <- data %>% # Remove ZCTAs without heatwaves from the data
  dplyr::filter(!zip %in% no_hw)

#Filter out heatwaves with lag periods outside of May-September (heatwaves occurring during the last week of September)
data <- data %>% # Create sept_heatwave variable that =1 if heatwave occurs during September 23-30
  mutate(sept_heatwave=ifelse(heatwave==1 & date %in% c("2016-09-23","2016-09-24","2016-09-25", "2016-09-26", "2016-09-27", "2016-09-28", "2016-09-29", "2016-09-30", 
                                                        "2017-09-23","2017-09-24","2017-09-25", "2017-09-26", "2017-09-27", "2017-09-28", "2017-09-29", "2017-09-30", 
                                                        "2018-09-23","2018-09-24","2018-09-25", "2018-09-26", "2018-09-27", "2018-09-28", "2018-09-29", "2018-09-30", 
                                                        "2019-09-23","2019-09-24","2019-09-25", "2019-09-26", "2019-09-27", "2019-09-28", "2019-09-29", "2019-09-30"), 1, 0))

data <- data %>%
  mutate(heatwave=ifelse(heatwave==1 & sept_heatwave==0, 1, 0)) %>% # Remove late September heatwaves from the dataset
  dplyr::select(-sept_heatwave)

write_parquet(data, "Deliv_Heatwaves.parquet")
