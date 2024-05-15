
## Code to calculate cumulative lag (lag0 to lag3) relative risk for outcomes for all deliveries and demogaphic subgroups
## Uses generalized nonlinear model

##### FULL #####
#HEATWAVES
data <- read_parquet("Data/Heatwave_Metrics.parquet")%>%
  setDT()%>%
  filter(Zip >= 27006 & Zip <= 28909)%>%
  rename(date=Date,
         heatwave=Heatwave,
         zip=Zip)%>%
  mutate(X=row_number(),
         year=year(date),
         doy=yday(date),
         dow=wday(date))%>%
  dplyr::select(date, heatwave, year, doy, dow, zip)%>%
  mutate(heatwave=as.integer(heatwave)) %>%
  filter(date >= "2016-01-01" & date <= "2019-12-31")%>%
  filter(month(date) >= 5 & month(date) <= 9)

datagrp <- data %>%
  group_by(zip)%>%
  summarize(hw=sum(heatwave))

#DELIVERIES
deliv <- read_parquet("Data/Deliv_Heatwaves.parquet")%>%
  setDT()%>%
  mutate(date=as.Date(admitdt),
         zip=as.numeric(ptzip))%>%
  filter(zip >= 27006 & zip <= 28909)%>%
  dplyr::select(zip, date, SMM21, GDM, HDP, SMI, PTB, MDP, PMAD, SUB)%>%
  filter(month(date) >= 5 & month(date) <= 9)%>%
  filter(year(date) >= 2016 & year(date) <= 2019)

#POP
pop <- read.csv("Data/TotalPop_F_18_44_2016_2021.csv")%>%
  mutate(denom=as.numeric(denom))

#CROSSWALK
cw <- read.csv("Data/Crosswalk_NC.csv")%>%
  setDT()%>%
  rename(zip=ZIP)

### MERGE DATASETS ###

deliv <- deliv %>%
  left_join(cw, by=c('zip'))%>%
  dplyr::select(-zip, ZIP_TYPE)%>%
  rename(zip=ZCTA)

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

data <- data %>%
  left_join(delivzip, by=c('zip', 'date'))

data[is.na(data)] <- 0

data <- data %>%
  left_join(pop, by=c('zip', 'year'))

table(is.na(data))

data <- data %>% # Remove zip codes where denom = NA
  dplyr::filter(zip != 28583 & zip != 28647 & zip != 28018) 

datagrp <- data %>%
  group_by(zip)%>%
  summarize(hw=sum(heatwave))

datagrp <- datagrp %>%
  filter(hw==0)
no_hw <- datagrp$zip

data <- data %>%
  dplyr::filter(!zip %in% no_hw)

data <- data %>%
  mutate(sept_heatwave=ifelse(heatwave==1 & date %in% c("2016-09-23","2016-09-24","2016-09-25", "2016-09-26", "2016-09-27", "2016-09-28", "2016-09-29", "2016-09-30", 
                                                        "2017-09-23","2017-09-24","2017-09-25", "2017-09-26", "2017-09-27", "2017-09-28", "2017-09-29", "2017-09-30", 
                                                        "2018-09-23","2018-09-24","2018-09-25", "2018-09-26", "2018-09-27", "2018-09-28", "2018-09-29", "2018-09-30", 
                                                        "2019-09-23","2019-09-24","2019-09-25", "2019-09-26", "2019-09-27", "2019-09-28", "2019-09-29", "2019-09-30"), 1, 0))

data <- data %>%
  mutate(old_heatwave=heatwave)

data <- data %>%
  mutate(heatwave=ifelse(old_heatwave==1 & sept_heatwave==0, 1, 0))%>%
  dplyr::select(-old_heatwave, -sept_heatwave)

### MATCH HEATWAVE DAYS TO NON-EXPOSED DAYS ###

dat<-data
head(dat)

zip_list <- unique(dat$zip) 
setorder(dat, zip, date)

gc()

### MATCHING PART ###

for (i in 1:length(zip_list)) {
  df <- subset(dat, zip == zip_list[i])
  
  # exclude the 3 days within any other heatwave
  df$time <- 1:nrow(df)
  cand_control <- unique(c(which(df$heatwave == 1), which(df$heatwave == 1) + 1, which(df$heatwave == 1) - 1))
  df$cand_control <- TRUE
  df$cand_control[cand_control[cand_control <= nrow(df)]] <- FALSE # Ensure cand_control indices are within df's range
  
  case_dates <- subset(df, heatwave == 1)
  control_dates <- subset(df, heatwave == 0)
  
  for (j in 1:nrow(case_dates)) {
    # choose lags (lagged 0 to 3)
    lag_dates <- case_dates[j, ]$date + 0:3
    lag_case <- subset(df, date %in% lag_dates)
    
    # choose 10 comparable unexposed days for each heatwave-exposed day
    control_range <- case_dates[j, ]$doy + -3:3
    control_subset <- subset(control_dates,
                             control_dates$year != case_dates[j, ]$year &
                               doy %in% control_range &
                               cand_control)
    controls <- dplyr::sample_n(control_subset, 3) #Changed from 10 to 3
    
    # choose lagged days for selected unexposed days
    la_con <- c(1:7)
    for (p in 1:length(la_con)) {
      lag_control_dates <- controls$date + la_con[p]  # Fixed index from j to p
      lag_control_each <- subset(df, date %in% lag_control_dates)
      
      if (p == 1) {
        lag_control <- lag_control_each
      } else {
        lag_control <- rbind(lag_control, lag_control_each)
      }
    }
    j_stratum <- rbind(lag_case, controls, lag_control)
    stratum <- paste("stratum", j, sep = ".")
    j_stratum$stratum <- stratum
    status <- c(rep("case", nrow(lag_case)), rep("control", nrow(controls)), rep("control", nrow(lag_control)))
    j_stratum$status <- status
    # Adjusted the length calculation for the lag column
    lag <- c(rep(0:3, length.out = nrow(lag_case)), rep(0, length.out = nrow(controls)), rep(c(1:3), length.out = nrow(lag_control)))
    j_stratum$lag <- lag
    
    if (j == 1) {
      new_df <- j_stratum
    } else {
      new_df <- rbind(new_df, j_stratum)
    }
  }
  if (i == 1) {
    matched_df <- new_df
  } else {
    matched_df <- rbind(matched_df, new_df)
  }
}

gc()

# We use "dlnm" package to generate the distributed lag function for "heatwave" 
library(lme4); library(dlnm); library(dplyr); library(stats)

for(i in 1:length(zip_list)){
  orig_dat <- subset(dat, zip == zip_list[i])
  match_dat <- subset(matched_df, zip == zip_list[i])
  
  orig_cb <- dlnm::crossbasis(orig_dat$heatwave, lag = c(0,3),
                              argvar = list(fun = "lin"),
                              arglag = list(fun = "integer"))
  obs_n <- nrow(orig_dat)
  orig_cb_matr <- as.data.frame(subset(orig_cb, nrow = obs_n))
  orig_cb_matr$date <- orig_dat$date
  matched_date <- match_dat %>% dplyr::select(date)
  matched_cb_matr <- orig_cb_matr %>%
    dplyr::right_join(matched_date, by = "date") %>%
    dplyr::select(-date) %>% as.matrix()
  
  if(i == 1){
    matched_cb_matrix <- matched_cb_matr
  }else{
    matched_cb_matrix <- rbind(matched_cb_matrix, matched_cb_matr)
  }
  # add attributes
  matched_dim <- dim(matched_cb_matrix)
  attr <- attributes(orig_cb)
  attr$dim <- matched_dim
  matched_cb <- matched_cb_matrix
  attributes(matched_cb) <- attr
}

gc()

# Check for incomplete strata

table(matched_df$lag)

restricted_df <- matched_df %>%
  group_by(zip, stratum)%>%
  summarize(n=n())
setorder(restricted_df, zip, stratum)

### GET CUMULATIVE LAG ###

outcomes <- c("gdm", "ptb", "hdp", "smi", "pmad", "mdp", "sub", "smm21")
outcome_data <- list()

for (outcome in outcomes) {
  formula <- as.formula(paste0(outcome, " ~ matched_cb + factor(dow) + year * zip"))
  
  fit <- gnm::gnm(formula,
                  eliminate = factor(zip), family = quasipoisson(link = "log"),
                  data = matched_df, offset = log(denom))
  
  pred <- dlnm::crosspred(matched_cb, fit, at = 1)
  
  over_rr <- sum(pred$matRRfit) / 4
  
  library(msm)
  estvar <- pred$vcov
  estmean <- c(pred$coefficients)
  
  over_rr_se <- msm::deltamethod(~ (exp(x1) + exp(x2) + exp(x3) + exp(x4)) / 4, 
                                 estmean, estvar)
  over_rr_low <- over_rr / exp(1.96 * over_rr_se)
  over_rr_high <- over_rr * exp(1.96 * over_rr_se)
  
  # Create a data frame for the current outcome
  outcome_df <- data.frame(outcome = outcome,
                           over_rr = over_rr,
                           over_rr_se = over_rr_se,
                           over_rr_low = over_rr_low,
                           over_rr_high = over_rr_high,
                           estvar = estvar,
                           estmean = estmean)
  
  # Store the data frame in the list
  outcome_data[[paste0(outcome, "_cumlag")]] <- outcome_df
}

data_frames <- unname(outcome_data)
merged_df <- do.call(rbind, data_frames)
print(merged_df)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="All Deliveries")
write.csv(merged_df, "Results/cumulative_lag_full.csv")

##### BLACK #####
#HEATWAVES
data <- read_parquet("Data/Heatwave_Metrics.parquet")%>%
  setDT()%>%
  filter(Zip >= 27006 & Zip <= 28909)%>%
  rename(date=Date,
         heatwave=Heatwave,
         zip=Zip)%>%
  mutate(X=row_number(),
         year=year(date),
         doy=yday(date),
         dow=wday(date))%>%
  dplyr::select(date, heatwave, year, doy, dow, zip)%>%
  mutate(heatwave=as.integer(heatwave)) %>%
  filter(date >= "2016-01-01" & date <= "2019-12-31")%>%
  filter(month(date) >= 5 & month(date) <= 9)

datagrp <- data %>%
  group_by(zip)%>%
  summarize(hw=sum(heatwave))

#DELIVERIES
deliv <- read_parquet("Data/Deliv_Heatwaves.parquet")%>%
  setDT()%>%
  mutate(date=as.Date(admitdt),
         zip=as.numeric(ptzip))%>%
  filter(zip >= 27006 & zip <= 28909)%>%
  filter(Race=="Black") %>%
  dplyr::select(zip, date, SMM21, GDM, HDP, SMI, PTB, MDP, PMAD, SUB)%>%
  filter(month(date) >= 5 & month(date) <= 9)%>%
  filter(year(date) >= 2016 & year(date) <= 2019)

#POP
pop <- read.csv("Data/TotalPop_F_18_44_2016_2021.csv")%>%
  mutate(denom=as.numeric(denom))

#CROSSWALK
cw <- read.csv("Data/Crosswalk_NC.csv")%>%
  setDT()%>%
  rename(zip=ZIP)

# Merge datasets

deliv <- deliv %>%
  left_join(cw, by=c('zip'))%>%
  dplyr::select(-zip, ZIP_TYPE)%>%
  rename(zip=ZCTA)

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

data <- data %>%
  left_join(delivzip, by=c('zip', 'date'))

data[is.na(data)] <- 0

data <- data %>%
  left_join(pop, by=c('zip', 'year'))

table(is.na(data))

data <- data %>% # Remove zip codes where denom = NA
  dplyr::filter(zip != 28583 & zip != 28647 & zip != 28018) 

datagrp <- data %>%
  group_by(zip)%>%
  summarize(hw=sum(heatwave))

datagrp <- datagrp %>%
  filter(hw==0)
no_hw <- datagrp$zip

data <- data %>%
  dplyr::filter(!zip %in% no_hw)

data <- data %>%
  mutate(sept_heatwave=ifelse(heatwave==1 & date %in% c("2016-09-23","2016-09-24","2016-09-25", "2016-09-26", "2016-09-27", "2016-09-28", "2016-09-29", "2016-09-30", 
                                                        "2017-09-23","2017-09-24","2017-09-25", "2017-09-26", "2017-09-27", "2017-09-28", "2017-09-29", "2017-09-30", 
                                                        "2018-09-23","2018-09-24","2018-09-25", "2018-09-26", "2018-09-27", "2018-09-28", "2018-09-29", "2018-09-30", 
                                                        "2019-09-23","2019-09-24","2019-09-25", "2019-09-26", "2019-09-27", "2019-09-28", "2019-09-29", "2019-09-30"), 1, 0))

data <- data %>%
  mutate(old_heatwave=heatwave)

data <- data %>%
  mutate(heatwave=ifelse(old_heatwave==1 & sept_heatwave==0, 1, 0))%>%
  dplyr::select(-old_heatwave, -sept_heatwave)

# Match days

dat<-data
head(dat)

zip_list <- unique(dat$zip) 
setorder(dat, zip, date)

gc()

# Match function

for (i in 1:length(zip_list)) {
  df <- subset(dat, zip == zip_list[i])
  
  # exclude the 3 days within any other heatwave
  df$time <- 1:nrow(df)
  cand_control <- unique(c(which(df$heatwave == 1), which(df$heatwave == 1) + 1, which(df$heatwave == 1) - 1))
  df$cand_control <- TRUE
  df$cand_control[cand_control[cand_control <= nrow(df)]] <- FALSE # Ensure cand_control indices are within df's range
  
  case_dates <- subset(df, heatwave == 1)
  control_dates <- subset(df, heatwave == 0)
  
  for (j in 1:nrow(case_dates)) {
    # choose lags (lagged 0 to 3)
    lag_dates <- case_dates[j, ]$date + 0:3
    lag_case <- subset(df, date %in% lag_dates)
    
    # choose 10 comparable unexposed days for each heatwave-exposed day
    control_range <- case_dates[j, ]$doy + -3:3
    control_subset <- subset(control_dates,
                             control_dates$year != case_dates[j, ]$year &
                               doy %in% control_range &
                               cand_control)
    controls <- dplyr::sample_n(control_subset, 3) #Changed from 10 to 3
    
    # choose lagged days for selected unexposed days
    la_con <- c(1:7)
    for (p in 1:length(la_con)) {
      lag_control_dates <- controls$date + la_con[p]  # Fixed index from j to p
      lag_control_each <- subset(df, date %in% lag_control_dates)
      
      if (p == 1) {
        lag_control <- lag_control_each
      } else {
        lag_control <- rbind(lag_control, lag_control_each)
      }
    }
    j_stratum <- rbind(lag_case, controls, lag_control)
    stratum <- paste("stratum", j, sep = ".")
    j_stratum$stratum <- stratum
    status <- c(rep("case", nrow(lag_case)), rep("control", nrow(controls)), rep("control", nrow(lag_control)))
    j_stratum$status <- status
    # Adjusted the length calculation for the lag column
    lag <- c(rep(0:3, length.out = nrow(lag_case)), rep(0, length.out = nrow(controls)), rep(c(1:3), length.out = nrow(lag_control)))
    j_stratum$lag <- lag
    
    if (j == 1) {
      new_df <- j_stratum
    } else {
      new_df <- rbind(new_df, j_stratum)
    }
  }
  if (i == 1) {
    matched_df <- new_df
  } else {
    matched_df <- rbind(matched_df, new_df)
  }
}

gc()

# We use "dlnm" package to generate the distributed lag function for "heatwave" 
library(lme4); library(dlnm); library(dplyr); library(stats)

for(i in 1:length(zip_list)){
  orig_dat <- subset(dat, zip == zip_list[i])
  match_dat <- subset(matched_df, zip == zip_list[i])
  
  orig_cb <- dlnm::crossbasis(orig_dat$heatwave, lag = c(0,3),
                              argvar = list(fun = "lin"),
                              arglag = list(fun = "integer"))
  obs_n <- nrow(orig_dat)
  orig_cb_matr <- as.data.frame(subset(orig_cb, nrow = obs_n))
  orig_cb_matr$date <- orig_dat$date
  matched_date <- match_dat %>% dplyr::select(date)
  matched_cb_matr <- orig_cb_matr %>%
    dplyr::right_join(matched_date, by = "date") %>%
    dplyr::select(-date) %>% as.matrix()
  
  if(i == 1){
    matched_cb_matrix <- matched_cb_matr
  }else{
    matched_cb_matrix <- rbind(matched_cb_matrix, matched_cb_matr)
  }
  # add attributes
  matched_dim <- dim(matched_cb_matrix)
  attr <- attributes(orig_cb)
  attr$dim <- matched_dim
  matched_cb <- matched_cb_matrix
  attributes(matched_cb) <- attr
}

gc()

# Check for incomplete strata

table(matched_df$lag)

restricted_df <- matched_df %>%
  group_by(zip, stratum)%>%
  summarize(n=n())
setorder(restricted_df, zip, stratum)

### GET CUMULATIVE LAG ###

outcomes <- c("gdm", "ptb", "hdp", "smi", "pmad", "mdp", "sub", "smm21")
outcome_data <- list()

for (outcome in outcomes) {
  formula <- as.formula(paste0(outcome, " ~ matched_cb + factor(dow) + year * zip"))
  
  fit <- gnm::gnm(formula,
                  eliminate = factor(zip), family = quasipoisson(link = "log"),
                  data = matched_df, offset = log(denom))
  
  pred <- dlnm::crosspred(matched_cb, fit, at = 1)
  
  over_rr <- sum(pred$matRRfit) / 4
  
  library(msm)
  estvar <- pred$vcov
  estmean <- c(pred$coefficients)
  
  over_rr_se <- msm::deltamethod(~ (exp(x1) + exp(x2) + exp(x3) + exp(x4)) / 4, 
                                 estmean, estvar)
  over_rr_low <- over_rr / exp(1.96 * over_rr_se)
  over_rr_high <- over_rr * exp(1.96 * over_rr_se)
  
  # Create a data frame for the current outcome
  outcome_df <- data.frame(outcome = outcome,
                           over_rr = over_rr,
                           over_rr_se = over_rr_se,
                           over_rr_low = over_rr_low,
                           over_rr_high = over_rr_high,
                           estvar = estvar,
                           estmean = estmean)
  
  # Store the data frame in the list
  outcome_data[[paste0(outcome, "_cumlag")]] <- outcome_df
}

data_frames <- unname(outcome_data)
merged_df <- do.call(rbind, data_frames)
print(merged_df)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Black")
write.csv(merged_df, "Results/cumulative_lag_Black.csv")

#### HISPANIC ####

data <- read_parquet("Data/Heatwave_Metrics.parquet")%>%
  setDT()%>%
  filter(Zip >= 27006 & Zip <= 28909)%>%
  rename(date=Date,
         heatwave=Heatwave,
         zip=Zip)%>%
  mutate(X=row_number(),
         year=year(date),
         doy=yday(date),
         dow=wday(date))%>%
  dplyr::select(date, heatwave, year, doy, dow, zip)%>%
  mutate(heatwave=as.integer(heatwave)) %>%
  filter(date >= "2016-01-01" & date <= "2019-12-31")%>%
  filter(month(date) >= 5 & month(date) <= 9)

datagrp <- data %>%
  group_by(zip)%>%
  summarize(hw=sum(heatwave))

#DELIVERIES
deliv <- read_parquet("Data/Deliv_Heatwaves.parquet")%>%
  setDT()%>%
  mutate(date=as.Date(admitdt),
         zip=as.numeric(ptzip))%>%
  filter(zip >= 27006 & zip <= 28909)%>%
  filter(Ethnicity=="Hispanic") %>%
  dplyr::select(zip, date, SMM21, GDM, HDP, SMI, PTB, MDP, PMAD, SUB)%>%
  filter(month(date) >= 5 & month(date) <= 9)%>%
  filter(year(date) >= 2016 & year(date) <= 2019)

#POP
pop <- read.csv("Data/TotalPop_F_18_44_2016_2021.csv")%>%
  mutate(denom=as.numeric(denom))

#CROSSWALK
cw <- read.csv("Data/Crosswalk_NC.csv")%>%
  setDT()%>%
  rename(zip=ZIP)

# Merge data

deliv <- deliv %>%
  left_join(cw, by=c('zip'))%>%
  dplyr::select(-zip, ZIP_TYPE)%>%
  rename(zip=ZCTA)

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

data <- data %>%
  left_join(delivzip, by=c('zip', 'date'))

data[is.na(data)] <- 0

data <- data %>%
  left_join(pop, by=c('zip', 'year'))

table(is.na(data))

data <- data %>% # Remove zip codes where denom = NA
  dplyr::filter(zip != 28583 & zip != 28647 & zip != 28018) 

datagrp <- data %>%
  group_by(zip)%>%
  summarize(hw=sum(heatwave))

datagrp <- datagrp %>%
  filter(hw==0)
no_hw <- datagrp$zip

data <- data %>%
  dplyr::filter(!zip %in% no_hw)

data <- data %>%
  mutate(sept_heatwave=ifelse(heatwave==1 & date %in% c("2016-09-23","2016-09-24","2016-09-25", "2016-09-26", "2016-09-27", "2016-09-28", "2016-09-29", "2016-09-30", 
                                                        "2017-09-23","2017-09-24","2017-09-25", "2017-09-26", "2017-09-27", "2017-09-28", "2017-09-29", "2017-09-30", 
                                                        "2018-09-23","2018-09-24","2018-09-25", "2018-09-26", "2018-09-27", "2018-09-28", "2018-09-29", "2018-09-30", 
                                                        "2019-09-23","2019-09-24","2019-09-25", "2019-09-26", "2019-09-27", "2019-09-28", "2019-09-29", "2019-09-30"), 1, 0))

data <- data %>%
  mutate(old_heatwave=heatwave)

data <- data %>%
  mutate(heatwave=ifelse(old_heatwave==1 & sept_heatwave==0, 1, 0))%>%
  dplyr::select(-old_heatwave, -sept_heatwave)

# Match days

dat<-data
head(dat)

zip_list <- unique(dat$zip) 
setorder(dat, zip, date)

gc()

#Match days

for (i in 1:length(zip_list)) {
  df <- subset(dat, zip == zip_list[i])
  
  # exclude the 3 days within any other heatwave
  df$time <- 1:nrow(df)
  cand_control <- unique(c(which(df$heatwave == 1), which(df$heatwave == 1) + 1, which(df$heatwave == 1) - 1))
  df$cand_control <- TRUE
  df$cand_control[cand_control[cand_control <= nrow(df)]] <- FALSE # Ensure cand_control indices are within df's range
  
  case_dates <- subset(df, heatwave == 1)
  control_dates <- subset(df, heatwave == 0)
  
  for (j in 1:nrow(case_dates)) {
    # choose lags (lagged 0 to 3)
    lag_dates <- case_dates[j, ]$date + 0:3
    lag_case <- subset(df, date %in% lag_dates)
    
    # choose 10 comparable unexposed days for each heatwave-exposed day
    control_range <- case_dates[j, ]$doy + -3:3
    control_subset <- subset(control_dates,
                             control_dates$year != case_dates[j, ]$year &
                               doy %in% control_range &
                               cand_control)
    controls <- dplyr::sample_n(control_subset, 3) #Changed from 10 to 3
    
    # choose lagged days for selected unexposed days
    la_con <- c(1:7)
    for (p in 1:length(la_con)) {
      lag_control_dates <- controls$date + la_con[p]  # Fixed index from j to p
      lag_control_each <- subset(df, date %in% lag_control_dates)
      
      if (p == 1) {
        lag_control <- lag_control_each
      } else {
        lag_control <- rbind(lag_control, lag_control_each)
      }
    }
    j_stratum <- rbind(lag_case, controls, lag_control)
    stratum <- paste("stratum", j, sep = ".")
    j_stratum$stratum <- stratum
    status <- c(rep("case", nrow(lag_case)), rep("control", nrow(controls)), rep("control", nrow(lag_control)))
    j_stratum$status <- status
    # Adjusted the length calculation for the lag column
    lag <- c(rep(0:3, length.out = nrow(lag_case)), rep(0, length.out = nrow(controls)), rep(c(1:3), length.out = nrow(lag_control)))
    j_stratum$lag <- lag
    
    if (j == 1) {
      new_df <- j_stratum
    } else {
      new_df <- rbind(new_df, j_stratum)
    }
  }
  if (i == 1) {
    matched_df <- new_df
  } else {
    matched_df <- rbind(matched_df, new_df)
  }
}

gc()

# We use "dlnm" package to generate the distributed lag function for "heatwave" 
library(lme4); library(dlnm); library(dplyr); library(stats)

for(i in 1:length(zip_list)){
  orig_dat <- subset(dat, zip == zip_list[i])
  match_dat <- subset(matched_df, zip == zip_list[i])
  
  orig_cb <- dlnm::crossbasis(orig_dat$heatwave, lag = c(0,3),
                              argvar = list(fun = "lin"),
                              arglag = list(fun = "integer"))
  obs_n <- nrow(orig_dat)
  orig_cb_matr <- as.data.frame(subset(orig_cb, nrow = obs_n))
  orig_cb_matr$date <- orig_dat$date
  matched_date <- match_dat %>% dplyr::select(date)
  matched_cb_matr <- orig_cb_matr %>%
    dplyr::right_join(matched_date, by = "date") %>%
    dplyr::select(-date) %>% as.matrix()
  
  if(i == 1){
    matched_cb_matrix <- matched_cb_matr
  }else{
    matched_cb_matrix <- rbind(matched_cb_matrix, matched_cb_matr)
  }
  # add attributes
  matched_dim <- dim(matched_cb_matrix)
  attr <- attributes(orig_cb)
  attr$dim <- matched_dim
  matched_cb <- matched_cb_matrix
  attributes(matched_cb) <- attr
}

gc()

# Check for incomplete strata

table(matched_df$lag)

restricted_df <- matched_df %>%
  group_by(zip, stratum)%>%
  summarize(n=n())
setorder(restricted_df, zip, stratum)

### GET CUMULATIVE LAG ###

outcomes <- c("gdm", "ptb", "hdp", "smi", "pmad", "mdp", "sub", "smm21")
outcome_data <- list()

for (outcome in outcomes) {
  formula <- as.formula(paste0(outcome, " ~ matched_cb + factor(dow) + year * zip"))
  
  fit <- gnm::gnm(formula,
                  eliminate = factor(zip), family = quasipoisson(link = "log"),
                  data = matched_df, offset = log(denom))
  
  pred <- dlnm::crosspred(matched_cb, fit, at = 1)
  
  over_rr <- sum(pred$matRRfit) / 4
  
  library(msm)
  estvar <- pred$vcov
  estmean <- c(pred$coefficients)
  
  over_rr_se <- msm::deltamethod(~ (exp(x1) + exp(x2) + exp(x3) + exp(x4)) / 4, 
                                 estmean, estvar)
  over_rr_low <- over_rr / exp(1.96 * over_rr_se)
  over_rr_high <- over_rr * exp(1.96 * over_rr_se)
  
  # Create a data frame for the current outcome
  outcome_df <- data.frame(outcome = outcome,
                           over_rr = over_rr,
                           over_rr_se = over_rr_se,
                           over_rr_low = over_rr_low,
                           over_rr_high = over_rr_high,
                           estvar = estvar,
                           estmean = estmean)
  
  # Store the data frame in the list
  outcome_data[[paste0(outcome, "_cumlag")]] <- outcome_df
}

data_frames <- unname(outcome_data)
merged_df <- do.call(rbind, data_frames)
print(merged_df)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Hispanic")
write.csv(merged_df, "Results/cumulative_lag_Hispanic.csv")

#### WHITE ####

data <- read_parquet("Data/Heatwave_Metrics.parquet")%>%
  setDT()%>%
  filter(Zip >= 27006 & Zip <= 28909)%>%
  rename(date=Date,
         heatwave=Heatwave,
         zip=Zip)%>%
  mutate(X=row_number(),
         year=year(date),
         doy=yday(date),
         dow=wday(date))%>%
  dplyr::select(date, heatwave, year, doy, dow, zip)%>%
  mutate(heatwave=as.integer(heatwave)) %>%
  filter(date >= "2016-01-01" & date <= "2019-12-31")%>%
  filter(month(date) >= 5 & month(date) <= 9)

datagrp <- data %>%
  group_by(zip)%>%
  summarize(hw=sum(heatwave))

#DELIVERIES
deliv <- read_parquet("Data/Deliv_Heatwaves.parquet")%>%
  setDT()%>%
  mutate(date=as.Date(admitdt),
         zip=as.numeric(ptzip))%>%
  filter(zip >= 27006 & zip <= 28909)%>%
  filter(Race=="White") %>%
  dplyr::select(zip, date, SMM21, GDM, HDP, SMI, PTB, MDP, PMAD, SUB)%>%
  filter(month(date) >= 5 & month(date) <= 9)%>%
  filter(year(date) >= 2016 & year(date) <= 2019)

#POP
pop <- read.csv("Data/TotalPop_F_18_44_2016_2021.csv")%>%
  mutate(denom=as.numeric(denom))

#CROSSWALK
cw <- read.csv("Data/Crosswalk_NC.csv")%>%
  setDT()%>%
  rename(zip=ZIP)

# Merge data

deliv <- deliv %>%
  left_join(cw, by=c('zip'))%>%
  dplyr::select(-zip, ZIP_TYPE)%>%
  rename(zip=ZCTA)

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

data <- data %>%
  left_join(delivzip, by=c('zip', 'date'))

data[is.na(data)] <- 0

data <- data %>%
  left_join(pop, by=c('zip', 'year'))

table(is.na(data))

data <- data %>% # Remove zip codes where denom = NA
  dplyr::filter(zip != 28583 & zip != 28647 & zip != 28018) 

datagrp <- data %>%
  group_by(zip)%>%
  summarize(hw=sum(heatwave))

datagrp <- datagrp %>%
  filter(hw==0)
no_hw <- datagrp$zip

data <- data %>%
  dplyr::filter(!zip %in% no_hw)

data <- data %>%
  mutate(sept_heatwave=ifelse(heatwave==1 & date %in% c("2016-09-23","2016-09-24","2016-09-25", "2016-09-26", "2016-09-27", "2016-09-28", "2016-09-29", "2016-09-30", 
                                                        "2017-09-23","2017-09-24","2017-09-25", "2017-09-26", "2017-09-27", "2017-09-28", "2017-09-29", "2017-09-30", 
                                                        "2018-09-23","2018-09-24","2018-09-25", "2018-09-26", "2018-09-27", "2018-09-28", "2018-09-29", "2018-09-30", 
                                                        "2019-09-23","2019-09-24","2019-09-25", "2019-09-26", "2019-09-27", "2019-09-28", "2019-09-29", "2019-09-30"), 1, 0))

data <- data %>%
  mutate(old_heatwave=heatwave)

data <- data %>%
  mutate(heatwave=ifelse(old_heatwave==1 & sept_heatwave==0, 1, 0))%>%
  dplyr::select(-old_heatwave, -sept_heatwave)

# Match days

dat<-data
head(dat)

zip_list <- unique(dat$zip) 
setorder(dat, zip, date)

gc()

#Match days

for (i in 1:length(zip_list)) {
  df <- subset(dat, zip == zip_list[i])
  
  # exclude the 3 days within any other heatwave
  df$time <- 1:nrow(df)
  cand_control <- unique(c(which(df$heatwave == 1), which(df$heatwave == 1) + 1, which(df$heatwave == 1) - 1))
  df$cand_control <- TRUE
  df$cand_control[cand_control[cand_control <= nrow(df)]] <- FALSE # Ensure cand_control indices are within df's range
  
  case_dates <- subset(df, heatwave == 1)
  control_dates <- subset(df, heatwave == 0)
  
  for (j in 1:nrow(case_dates)) {
    # choose lags (lagged 0 to 3)
    lag_dates <- case_dates[j, ]$date + 0:3
    lag_case <- subset(df, date %in% lag_dates)
    
    # choose 10 comparable unexposed days for each heatwave-exposed day
    control_range <- case_dates[j, ]$doy + -3:3
    control_subset <- subset(control_dates,
                             control_dates$year != case_dates[j, ]$year &
                               doy %in% control_range &
                               cand_control)
    controls <- dplyr::sample_n(control_subset, 3) #Changed from 10 to 3
    
    # choose lagged days for selected unexposed days
    la_con <- c(1:7)
    for (p in 1:length(la_con)) {
      lag_control_dates <- controls$date + la_con[p]  # Fixed index from j to p
      lag_control_each <- subset(df, date %in% lag_control_dates)
      
      if (p == 1) {
        lag_control <- lag_control_each
      } else {
        lag_control <- rbind(lag_control, lag_control_each)
      }
    }
    j_stratum <- rbind(lag_case, controls, lag_control)
    stratum <- paste("stratum", j, sep = ".")
    j_stratum$stratum <- stratum
    status <- c(rep("case", nrow(lag_case)), rep("control", nrow(controls)), rep("control", nrow(lag_control)))
    j_stratum$status <- status
    # Adjusted the length calculation for the lag column
    lag <- c(rep(0:3, length.out = nrow(lag_case)), rep(0, length.out = nrow(controls)), rep(c(1:3), length.out = nrow(lag_control)))
    j_stratum$lag <- lag
    
    if (j == 1) {
      new_df <- j_stratum
    } else {
      new_df <- rbind(new_df, j_stratum)
    }
  }
  if (i == 1) {
    matched_df <- new_df
  } else {
    matched_df <- rbind(matched_df, new_df)
  }
}

gc()

# We use "dlnm" package to generate the distributed lag function for "heatwave" 
library(lme4); library(dlnm); library(dplyr); library(stats)

for(i in 1:length(zip_list)){
  orig_dat <- subset(dat, zip == zip_list[i])
  match_dat <- subset(matched_df, zip == zip_list[i])
  
  orig_cb <- dlnm::crossbasis(orig_dat$heatwave, lag = c(0,3),
                              argvar = list(fun = "lin"),
                              arglag = list(fun = "integer"))
  obs_n <- nrow(orig_dat)
  orig_cb_matr <- as.data.frame(subset(orig_cb, nrow = obs_n))
  orig_cb_matr$date <- orig_dat$date
  matched_date <- match_dat %>% dplyr::select(date)
  matched_cb_matr <- orig_cb_matr %>%
    dplyr::right_join(matched_date, by = "date") %>%
    dplyr::select(-date) %>% as.matrix()
  
  if(i == 1){
    matched_cb_matrix <- matched_cb_matr
  }else{
    matched_cb_matrix <- rbind(matched_cb_matrix, matched_cb_matr)
  }
  # add attributes
  matched_dim <- dim(matched_cb_matrix)
  attr <- attributes(orig_cb)
  attr$dim <- matched_dim
  matched_cb <- matched_cb_matrix
  attributes(matched_cb) <- attr
}

gc()

# Check for incomplete strata

table(matched_df$lag)

restricted_df <- matched_df %>%
  group_by(zip, stratum)%>%
  summarize(n=n())
setorder(restricted_df, zip, stratum)

### GET CUMULATIVE LAG ###

outcomes <- c("gdm", "ptb", "hdp", "smi", "pmad", "mdp", "sub", "smm21")
outcome_data <- list()

for (outcome in outcomes) {
  formula <- as.formula(paste0(outcome, " ~ matched_cb + factor(dow) + year * zip"))
  
  fit <- gnm::gnm(formula,
                  eliminate = factor(zip), family = quasipoisson(link = "log"),
                  data = matched_df, offset = log(denom))
  
  pred <- dlnm::crosspred(matched_cb, fit, at = 1)
  
  over_rr <- sum(pred$matRRfit) / 4
  
  library(msm)
  estvar <- pred$vcov
  estmean <- c(pred$coefficients)
  
  over_rr_se <- msm::deltamethod(~ (exp(x1) + exp(x2) + exp(x3) + exp(x4)) / 4, 
                                 estmean, estvar)
  over_rr_low <- over_rr / exp(1.96 * over_rr_se)
  over_rr_high <- over_rr * exp(1.96 * over_rr_se)
  
  # Create a data frame for the current outcome
  outcome_df <- data.frame(outcome = outcome,
                           over_rr = over_rr,
                           over_rr_se = over_rr_se,
                           over_rr_low = over_rr_low,
                           over_rr_high = over_rr_high,
                           estvar = estvar,
                           estmean = estmean)
  
  # Store the data frame in the list
  outcome_data[[paste0(outcome, "_cumlag")]] <- outcome_df
}

data_frames <- unname(outcome_data)
merged_df <- do.call(rbind, data_frames)
print(merged_df)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="White")
write.csv(merged_df, "Results/cumulative_lag_White.csv")

#### AGE ####

data <- read_parquet("Data/Heatwave_Metrics.parquet")%>%
  setDT()%>%
  filter(Zip >= 27006 & Zip <= 28909)%>%
  rename(date=Date,
         heatwave=Heatwave,
         zip=Zip)%>%
  mutate(X=row_number(),
         year=year(date),
         doy=yday(date),
         dow=wday(date))%>%
  dplyr::select(date, heatwave, year, doy, dow, zip)%>%
  mutate(heatwave=as.integer(heatwave)) %>%
  filter(date >= "2016-01-01" & date <= "2019-12-31")%>%
  filter(month(date) >= 5 & month(date) <= 9)

datagrp <- data %>%
  group_by(zip)%>%
  summarize(hw=sum(heatwave))

#DELIVERIES
deliv <- read_parquet("Data/Deliv_Heatwaves.parquet")%>%
  setDT()%>%
  mutate(date=as.Date(admitdt),
         zip=as.numeric(ptzip))%>%
  filter(zip >= 27006 & zip <= 28909)%>%
  filter(Age=="35-39" | Age=="40+") %>%
  dplyr::select(zip, date, SMM21, GDM, HDP, SMI, PTB, MDP, PMAD, SUB)%>%
  filter(month(date) >= 5 & month(date) <= 9)%>%
  filter(year(date) >= 2016 & year(date) <= 2019)

#POP
pop <- read.csv("Data/TotalPop_F_18_44_2016_2021.csv")%>%
  mutate(denom=as.numeric(denom))

#CROSSWALK
cw <- read.csv("Data/Crosswalk_NC.csv")%>%
  setDT()%>%
  rename(zip=ZIP)

# Merge data

deliv <- deliv %>%
  left_join(cw, by=c('zip'))%>%
  dplyr::select(-zip, ZIP_TYPE)%>%
  rename(zip=ZCTA)

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

data <- data %>%
  left_join(delivzip, by=c('zip', 'date'))

data[is.na(data)] <- 0

data <- data %>%
  left_join(pop, by=c('zip', 'year'))

table(is.na(data))

data <- data %>% # Remove zip codes where denom = NA
  dplyr::filter(zip != 28583 & zip != 28647 & zip != 28018) 

datagrp <- data %>%
  group_by(zip)%>%
  summarize(hw=sum(heatwave))

datagrp <- datagrp %>%
  filter(hw==0)
no_hw <- datagrp$zip

data <- data %>%
  dplyr::filter(!zip %in% no_hw)

data <- data %>%
  mutate(sept_heatwave=ifelse(heatwave==1 & date %in% c("2016-09-23","2016-09-24","2016-09-25", "2016-09-26", "2016-09-27", "2016-09-28", "2016-09-29", "2016-09-30", 
                                                        "2017-09-23","2017-09-24","2017-09-25", "2017-09-26", "2017-09-27", "2017-09-28", "2017-09-29", "2017-09-30", 
                                                        "2018-09-23","2018-09-24","2018-09-25", "2018-09-26", "2018-09-27", "2018-09-28", "2018-09-29", "2018-09-30", 
                                                        "2019-09-23","2019-09-24","2019-09-25", "2019-09-26", "2019-09-27", "2019-09-28", "2019-09-29", "2019-09-30"), 1, 0))

data <- data %>%
  mutate(old_heatwave=heatwave)

data <- data %>%
  mutate(heatwave=ifelse(old_heatwave==1 & sept_heatwave==0, 1, 0))%>%
  dplyr::select(-old_heatwave, -sept_heatwave)

# Match days

dat<-data
head(dat)

zip_list <- unique(dat$zip) 
setorder(dat, zip, date)

gc()

#Match days

for (i in 1:length(zip_list)) {
  df <- subset(dat, zip == zip_list[i])
  
  # exclude the 3 days within any other heatwave
  df$time <- 1:nrow(df)
  cand_control <- unique(c(which(df$heatwave == 1), which(df$heatwave == 1) + 1, which(df$heatwave == 1) - 1))
  df$cand_control <- TRUE
  df$cand_control[cand_control[cand_control <= nrow(df)]] <- FALSE # Ensure cand_control indices are within df's range
  
  case_dates <- subset(df, heatwave == 1)
  control_dates <- subset(df, heatwave == 0)
  
  for (j in 1:nrow(case_dates)) {
    # choose lags (lagged 0 to 3)
    lag_dates <- case_dates[j, ]$date + 0:3
    lag_case <- subset(df, date %in% lag_dates)
    
    # choose 10 comparable unexposed days for each heatwave-exposed day
    control_range <- case_dates[j, ]$doy + -3:3
    control_subset <- subset(control_dates,
                             control_dates$year != case_dates[j, ]$year &
                               doy %in% control_range &
                               cand_control)
    controls <- dplyr::sample_n(control_subset, 3) #Changed from 10 to 3
    
    # choose lagged days for selected unexposed days
    la_con <- c(1:7)
    for (p in 1:length(la_con)) {
      lag_control_dates <- controls$date + la_con[p]  # Fixed index from j to p
      lag_control_each <- subset(df, date %in% lag_control_dates)
      
      if (p == 1) {
        lag_control <- lag_control_each
      } else {
        lag_control <- rbind(lag_control, lag_control_each)
      }
    }
    j_stratum <- rbind(lag_case, controls, lag_control)
    stratum <- paste("stratum", j, sep = ".")
    j_stratum$stratum <- stratum
    status <- c(rep("case", nrow(lag_case)), rep("control", nrow(controls)), rep("control", nrow(lag_control)))
    j_stratum$status <- status
    # Adjusted the length calculation for the lag column
    lag <- c(rep(0:3, length.out = nrow(lag_case)), rep(0, length.out = nrow(controls)), rep(c(1:3), length.out = nrow(lag_control)))
    j_stratum$lag <- lag
    
    if (j == 1) {
      new_df <- j_stratum
    } else {
      new_df <- rbind(new_df, j_stratum)
    }
  }
  if (i == 1) {
    matched_df <- new_df
  } else {
    matched_df <- rbind(matched_df, new_df)
  }
}

gc()

# We use "dlnm" package to generate the distributed lag function for "heatwave" 
library(lme4); library(dlnm); library(dplyr); library(stats)

for(i in 1:length(zip_list)){
  orig_dat <- subset(dat, zip == zip_list[i])
  match_dat <- subset(matched_df, zip == zip_list[i])
  
  orig_cb <- dlnm::crossbasis(orig_dat$heatwave, lag = c(0,3),
                              argvar = list(fun = "lin"),
                              arglag = list(fun = "integer"))
  obs_n <- nrow(orig_dat)
  orig_cb_matr <- as.data.frame(subset(orig_cb, nrow = obs_n))
  orig_cb_matr$date <- orig_dat$date
  matched_date <- match_dat %>% dplyr::select(date)
  matched_cb_matr <- orig_cb_matr %>%
    dplyr::right_join(matched_date, by = "date") %>%
    dplyr::select(-date) %>% as.matrix()
  
  if(i == 1){
    matched_cb_matrix <- matched_cb_matr
  }else{
    matched_cb_matrix <- rbind(matched_cb_matrix, matched_cb_matr)
  }
  # add attributes
  matched_dim <- dim(matched_cb_matrix)
  attr <- attributes(orig_cb)
  attr$dim <- matched_dim
  matched_cb <- matched_cb_matrix
  attributes(matched_cb) <- attr
}

gc()

# Check for incomplete strata

table(matched_df$lag)

restricted_df <- matched_df %>%
  group_by(zip, stratum)%>%
  summarize(n=n())
setorder(restricted_df, zip, stratum)

### GET CUMULATIVE LAG ###

outcomes <- c("gdm", "ptb", "hdp", "smi", "pmad", "mdp", "sub", "smm21")
outcome_data <- list()

for (outcome in outcomes) {
  formula <- as.formula(paste0(outcome, " ~ matched_cb + factor(dow) + year * zip"))
  
  fit <- gnm::gnm(formula,
                  eliminate = factor(zip), family = quasipoisson(link = "log"),
                  data = matched_df, offset = log(denom))
  
  pred <- dlnm::crosspred(matched_cb, fit, at = 1)
  
  over_rr <- sum(pred$matRRfit) / 4
  
  library(msm)
  estvar <- pred$vcov
  estmean <- c(pred$coefficients)
  
  over_rr_se <- msm::deltamethod(~ (exp(x1) + exp(x2) + exp(x3) + exp(x4)) / 4, 
                                 estmean, estvar)
  over_rr_low <- over_rr / exp(1.96 * over_rr_se)
  over_rr_high <- over_rr * exp(1.96 * over_rr_se)
  
  # Create a data frame for the current outcome
  outcome_df <- data.frame(outcome = outcome,
                           over_rr = over_rr,
                           over_rr_se = over_rr_se,
                           over_rr_low = over_rr_low,
                           over_rr_high = over_rr_high,
                           estvar = estvar,
                           estmean = estmean)
  
  # Store the data frame in the list
  outcome_data[[paste0(outcome, "_cumlag")]] <- outcome_df
}

data_frames <- unname(outcome_data)
merged_df <- do.call(rbind, data_frames)
print(merged_df)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Age > 35")
write.csv(merged_df, "Results/cumulative_lag_Age.csv")

#### MEDICAID ####

data <- read_parquet("Data/Heatwave_Metrics.parquet")%>%
  setDT()%>%
  filter(Zip >= 27006 & Zip <= 28909)%>%
  rename(date=Date,
         heatwave=Heatwave,
         zip=Zip)%>%
  mutate(X=row_number(),
         year=year(date),
         doy=yday(date),
         dow=wday(date))%>%
  dplyr::select(date, heatwave, year, doy, dow, zip)%>%
  mutate(heatwave=as.integer(heatwave)) %>%
  filter(date >= "2016-01-01" & date <= "2019-12-31")%>%
  filter(month(date) >= 5 & month(date) <= 9)

datagrp <- data %>%
  group_by(zip)%>%
  summarize(hw=sum(heatwave))

#DELIVERIES
deliv <- read_parquet("Data/Deliv_Heatwaves.parquet")%>%
  setDT()%>%
  mutate(date=as.Date(admitdt),
         zip=as.numeric(ptzip))%>%
  filter(zip >= 27006 & zip <= 28909)%>%
  filter(Insurance=="Medicaid") %>%
  dplyr::select(zip, date, SMM21, GDM, HDP, SMI, PTB, MDP, PMAD, SUB)%>%
  filter(month(date) >= 5 & month(date) <= 9)%>%
  filter(year(date) >= 2016 & year(date) <= 2019)

#POP
pop <- read.csv("Data/TotalPop_F_18_44_2016_2021.csv")%>%
  mutate(denom=as.numeric(denom))

#CROSSWALK
cw <- read.csv("Data/Crosswalk_NC.csv")%>%
  setDT()%>%
  rename(zip=ZIP)

# Merge data

deliv <- deliv %>%
  left_join(cw, by=c('zip'))%>%
  dplyr::select(-zip, ZIP_TYPE)%>%
  rename(zip=ZCTA)

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

data <- data %>%
  left_join(delivzip, by=c('zip', 'date'))

data[is.na(data)] <- 0

data <- data %>%
  left_join(pop, by=c('zip', 'year'))

table(is.na(data))

data <- data %>% # Remove zip codes where denom = NA
  dplyr::filter(zip != 28583 & zip != 28647 & zip != 28018) 

datagrp <- data %>%
  group_by(zip)%>%
  summarize(hw=sum(heatwave))

datagrp <- datagrp %>%
  filter(hw==0)
no_hw <- datagrp$zip

data <- data %>%
  dplyr::filter(!zip %in% no_hw)

data <- data %>%
  mutate(sept_heatwave=ifelse(heatwave==1 & date %in% c("2016-09-23","2016-09-24","2016-09-25", "2016-09-26", "2016-09-27", "2016-09-28", "2016-09-29", "2016-09-30", 
                                                        "2017-09-23","2017-09-24","2017-09-25", "2017-09-26", "2017-09-27", "2017-09-28", "2017-09-29", "2017-09-30", 
                                                        "2018-09-23","2018-09-24","2018-09-25", "2018-09-26", "2018-09-27", "2018-09-28", "2018-09-29", "2018-09-30", 
                                                        "2019-09-23","2019-09-24","2019-09-25", "2019-09-26", "2019-09-27", "2019-09-28", "2019-09-29", "2019-09-30"), 1, 0))

data <- data %>%
  mutate(old_heatwave=heatwave)

data <- data %>%
  mutate(heatwave=ifelse(old_heatwave==1 & sept_heatwave==0, 1, 0))%>%
  dplyr::select(-old_heatwave, -sept_heatwave)

# Match days

dat<-data
head(dat)

zip_list <- unique(dat$zip) 
setorder(dat, zip, date)

gc()

#Match days

for (i in 1:length(zip_list)) {
  df <- subset(dat, zip == zip_list[i])
  
  # exclude the 3 days within any other heatwave
  df$time <- 1:nrow(df)
  cand_control <- unique(c(which(df$heatwave == 1), which(df$heatwave == 1) + 1, which(df$heatwave == 1) - 1))
  df$cand_control <- TRUE
  df$cand_control[cand_control[cand_control <= nrow(df)]] <- FALSE # Ensure cand_control indices are within df's range
  
  case_dates <- subset(df, heatwave == 1)
  control_dates <- subset(df, heatwave == 0)
  
  for (j in 1:nrow(case_dates)) {
    # choose lags (lagged 0 to 3)
    lag_dates <- case_dates[j, ]$date + 0:3
    lag_case <- subset(df, date %in% lag_dates)
    
    # choose 10 comparable unexposed days for each heatwave-exposed day
    control_range <- case_dates[j, ]$doy + -3:3
    control_subset <- subset(control_dates,
                             control_dates$year != case_dates[j, ]$year &
                               doy %in% control_range &
                               cand_control)
    controls <- dplyr::sample_n(control_subset, 3) #Changed from 10 to 3
    
    # choose lagged days for selected unexposed days
    la_con <- c(1:7)
    for (p in 1:length(la_con)) {
      lag_control_dates <- controls$date + la_con[p]  # Fixed index from j to p
      lag_control_each <- subset(df, date %in% lag_control_dates)
      
      if (p == 1) {
        lag_control <- lag_control_each
      } else {
        lag_control <- rbind(lag_control, lag_control_each)
      }
    }
    j_stratum <- rbind(lag_case, controls, lag_control)
    stratum <- paste("stratum", j, sep = ".")
    j_stratum$stratum <- stratum
    status <- c(rep("case", nrow(lag_case)), rep("control", nrow(controls)), rep("control", nrow(lag_control)))
    j_stratum$status <- status
    # Adjusted the length calculation for the lag column
    lag <- c(rep(0:3, length.out = nrow(lag_case)), rep(0, length.out = nrow(controls)), rep(c(1:3), length.out = nrow(lag_control)))
    j_stratum$lag <- lag
    
    if (j == 1) {
      new_df <- j_stratum
    } else {
      new_df <- rbind(new_df, j_stratum)
    }
  }
  if (i == 1) {
    matched_df <- new_df
  } else {
    matched_df <- rbind(matched_df, new_df)
  }
}

gc()

# We use "dlnm" package to generate the distributed lag function for "heatwave" 
library(lme4); library(dlnm); library(dplyr); library(stats)

for(i in 1:length(zip_list)){
  orig_dat <- subset(dat, zip == zip_list[i])
  match_dat <- subset(matched_df, zip == zip_list[i])
  
  orig_cb <- dlnm::crossbasis(orig_dat$heatwave, lag = c(0,3),
                              argvar = list(fun = "lin"),
                              arglag = list(fun = "integer"))
  obs_n <- nrow(orig_dat)
  orig_cb_matr <- as.data.frame(subset(orig_cb, nrow = obs_n))
  orig_cb_matr$date <- orig_dat$date
  matched_date <- match_dat %>% dplyr::select(date)
  matched_cb_matr <- orig_cb_matr %>%
    dplyr::right_join(matched_date, by = "date") %>%
    dplyr::select(-date) %>% as.matrix()
  
  if(i == 1){
    matched_cb_matrix <- matched_cb_matr
  }else{
    matched_cb_matrix <- rbind(matched_cb_matrix, matched_cb_matr)
  }
  # add attributes
  matched_dim <- dim(matched_cb_matrix)
  attr <- attributes(orig_cb)
  attr$dim <- matched_dim
  matched_cb <- matched_cb_matrix
  attributes(matched_cb) <- attr
}

gc()

# Check for incomplete strata

table(matched_df$lag)

restricted_df <- matched_df %>%
  group_by(zip, stratum)%>%
  summarize(n=n())
setorder(restricted_df, zip, stratum)

### GET CUMULATIVE LAG ###

outcomes <- c("gdm", "ptb", "hdp", "smi", "pmad", "mdp", "sub", "smm21")
outcome_data <- list()

for (outcome in outcomes) {
  formula <- as.formula(paste0(outcome, " ~ matched_cb + factor(dow) + year * zip"))
  
  fit <- gnm::gnm(formula,
                  eliminate = factor(zip), family = quasipoisson(link = "log"),
                  data = matched_df, offset = log(denom))
  
  pred <- dlnm::crosspred(matched_cb, fit, at = 1)
  
  over_rr <- sum(pred$matRRfit) / 4
  
  library(msm)
  estvar <- pred$vcov
  estmean <- c(pred$coefficients)
  
  over_rr_se <- msm::deltamethod(~ (exp(x1) + exp(x2) + exp(x3) + exp(x4)) / 4, 
                                 estmean, estvar)
  over_rr_low <- over_rr / exp(1.96 * over_rr_se)
  over_rr_high <- over_rr * exp(1.96 * over_rr_se)
  
  # Create a data frame for the current outcome
  outcome_df <- data.frame(outcome = outcome,
                           over_rr = over_rr,
                           over_rr_se = over_rr_se,
                           over_rr_low = over_rr_low,
                           over_rr_high = over_rr_high,
                           estvar = estvar,
                           estmean = estmean)
  
  # Store the data frame in the list
  outcome_data[[paste0(outcome, "_cumlag")]] <- outcome_df
}

data_frames <- unname(outcome_data)
merged_df <- do.call(rbind, data_frames)
print(merged_df)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Medicaid")
write.csv(merged_df, "Results/cumulative_lag_Medicaid.csv")

#### COMMERCIAL ####

data <- read_parquet("Data/Heatwave_Metrics.parquet")%>%
  setDT()%>%
  filter(Zip >= 27006 & Zip <= 28909)%>%
  rename(date=Date,
         heatwave=Heatwave,
         zip=Zip)%>%
  mutate(X=row_number(),
         year=year(date),
         doy=yday(date),
         dow=wday(date))%>%
  dplyr::select(date, heatwave, year, doy, dow, zip)%>%
  mutate(heatwave=as.integer(heatwave)) %>%
  filter(date >= "2016-01-01" & date <= "2019-12-31")%>%
  filter(month(date) >= 5 & month(date) <= 9)

datagrp <- data %>%
  group_by(zip)%>%
  summarize(hw=sum(heatwave))

#DELIVERIES
deliv <- read_parquet("Data/Deliv_Heatwaves.parquet")%>%
  setDT()%>%
  mutate(date=as.Date(admitdt),
         zip=as.numeric(ptzip))%>%
  filter(zip >= 27006 & zip <= 28909)%>%
  filter(Insurance=="Commercial") %>%
  dplyr::select(zip, date, SMM21, GDM, HDP, SMI, PTB, MDP, PMAD, SUB)%>%
  filter(month(date) >= 5 & month(date) <= 9)%>%
  filter(year(date) >= 2016 & year(date) <= 2019)

#POP
pop <- read.csv("Data/TotalPop_F_18_44_2016_2021.csv")%>%
  mutate(denom=as.numeric(denom))

#CROSSWALK
cw <- read.csv("Data/Crosswalk_NC.csv")%>%
  setDT()%>%
  rename(zip=ZIP)

# Merge data

deliv <- deliv %>%
  left_join(cw, by=c('zip'))%>%
  dplyr::select(-zip, ZIP_TYPE)%>%
  rename(zip=ZCTA)

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

data <- data %>%
  left_join(delivzip, by=c('zip', 'date'))

data[is.na(data)] <- 0

data <- data %>%
  left_join(pop, by=c('zip', 'year'))

table(is.na(data))

data <- data %>% # Remove zip codes where denom = NA
  dplyr::filter(zip != 28583 & zip != 28647 & zip != 28018) 

datagrp <- data %>%
  group_by(zip)%>%
  summarize(hw=sum(heatwave))

datagrp <- datagrp %>%
  filter(hw==0)
no_hw <- datagrp$zip

data <- data %>%
  dplyr::filter(!zip %in% no_hw)

data <- data %>%
  mutate(sept_heatwave=ifelse(heatwave==1 & date %in% c("2016-09-23","2016-09-24","2016-09-25", "2016-09-26", "2016-09-27", "2016-09-28", "2016-09-29", "2016-09-30", 
                                                        "2017-09-23","2017-09-24","2017-09-25", "2017-09-26", "2017-09-27", "2017-09-28", "2017-09-29", "2017-09-30", 
                                                        "2018-09-23","2018-09-24","2018-09-25", "2018-09-26", "2018-09-27", "2018-09-28", "2018-09-29", "2018-09-30", 
                                                        "2019-09-23","2019-09-24","2019-09-25", "2019-09-26", "2019-09-27", "2019-09-28", "2019-09-29", "2019-09-30"), 1, 0))

data <- data %>%
  mutate(old_heatwave=heatwave)

data <- data %>%
  mutate(heatwave=ifelse(old_heatwave==1 & sept_heatwave==0, 1, 0))%>%
  dplyr::select(-old_heatwave, -sept_heatwave)

# Match days

dat<-data
head(dat)

zip_list <- unique(dat$zip) 
setorder(dat, zip, date)

gc()

#Match days

for (i in 1:length(zip_list)) {
  df <- subset(dat, zip == zip_list[i])
  
  # exclude the 3 days within any other heatwave
  df$time <- 1:nrow(df)
  cand_control <- unique(c(which(df$heatwave == 1), which(df$heatwave == 1) + 1, which(df$heatwave == 1) - 1))
  df$cand_control <- TRUE
  df$cand_control[cand_control[cand_control <= nrow(df)]] <- FALSE # Ensure cand_control indices are within df's range
  
  case_dates <- subset(df, heatwave == 1)
  control_dates <- subset(df, heatwave == 0)
  
  for (j in 1:nrow(case_dates)) {
    # choose lags (lagged 0 to 3)
    lag_dates <- case_dates[j, ]$date + 0:3
    lag_case <- subset(df, date %in% lag_dates)
    
    # choose 10 comparable unexposed days for each heatwave-exposed day
    control_range <- case_dates[j, ]$doy + -3:3
    control_subset <- subset(control_dates,
                             control_dates$year != case_dates[j, ]$year &
                               doy %in% control_range &
                               cand_control)
    controls <- dplyr::sample_n(control_subset, 3) #Changed from 10 to 3
    
    # choose lagged days for selected unexposed days
    la_con <- c(1:7)
    for (p in 1:length(la_con)) {
      lag_control_dates <- controls$date + la_con[p]  # Fixed index from j to p
      lag_control_each <- subset(df, date %in% lag_control_dates)
      
      if (p == 1) {
        lag_control <- lag_control_each
      } else {
        lag_control <- rbind(lag_control, lag_control_each)
      }
    }
    j_stratum <- rbind(lag_case, controls, lag_control)
    stratum <- paste("stratum", j, sep = ".")
    j_stratum$stratum <- stratum
    status <- c(rep("case", nrow(lag_case)), rep("control", nrow(controls)), rep("control", nrow(lag_control)))
    j_stratum$status <- status
    # Adjusted the length calculation for the lag column
    lag <- c(rep(0:3, length.out = nrow(lag_case)), rep(0, length.out = nrow(controls)), rep(c(1:3), length.out = nrow(lag_control)))
    j_stratum$lag <- lag
    
    if (j == 1) {
      new_df <- j_stratum
    } else {
      new_df <- rbind(new_df, j_stratum)
    }
  }
  if (i == 1) {
    matched_df <- new_df
  } else {
    matched_df <- rbind(matched_df, new_df)
  }
}

gc()

# We use "dlnm" package to generate the distributed lag function for "heatwave" 
library(lme4); library(dlnm); library(dplyr); library(stats)

for(i in 1:length(zip_list)){
  orig_dat <- subset(dat, zip == zip_list[i])
  match_dat <- subset(matched_df, zip == zip_list[i])
  
  orig_cb <- dlnm::crossbasis(orig_dat$heatwave, lag = c(0,3),
                              argvar = list(fun = "lin"),
                              arglag = list(fun = "integer"))
  obs_n <- nrow(orig_dat)
  orig_cb_matr <- as.data.frame(subset(orig_cb, nrow = obs_n))
  orig_cb_matr$date <- orig_dat$date
  matched_date <- match_dat %>% dplyr::select(date)
  matched_cb_matr <- orig_cb_matr %>%
    dplyr::right_join(matched_date, by = "date") %>%
    dplyr::select(-date) %>% as.matrix()
  
  if(i == 1){
    matched_cb_matrix <- matched_cb_matr
  }else{
    matched_cb_matrix <- rbind(matched_cb_matrix, matched_cb_matr)
  }
  # add attributes
  matched_dim <- dim(matched_cb_matrix)
  attr <- attributes(orig_cb)
  attr$dim <- matched_dim
  matched_cb <- matched_cb_matrix
  attributes(matched_cb) <- attr
}

gc()

# Check for incomplete strata

table(matched_df$lag)

restricted_df <- matched_df %>%
  group_by(zip, stratum)%>%
  summarize(n=n())
setorder(restricted_df, zip, stratum)

### GET CUMULATIVE LAG ###

outcomes <- c("gdm", "ptb", "hdp", "smi", "pmad", "mdp", "sub", "smm21")
outcome_data <- list()

for (outcome in outcomes) {
  formula <- as.formula(paste0(outcome, " ~ matched_cb + factor(dow) + year * zip"))
  
  fit <- gnm::gnm(formula,
                  eliminate = factor(zip), family = quasipoisson(link = "log"),
                  data = matched_df, offset = log(denom))
  
  pred <- dlnm::crosspred(matched_cb, fit, at = 1)
  
  over_rr <- sum(pred$matRRfit) / 4
  
  library(msm)
  estvar <- pred$vcov
  estmean <- c(pred$coefficients)
  
  over_rr_se <- msm::deltamethod(~ (exp(x1) + exp(x2) + exp(x3) + exp(x4)) / 4, 
                                 estmean, estvar)
  over_rr_low <- over_rr / exp(1.96 * over_rr_se)
  over_rr_high <- over_rr * exp(1.96 * over_rr_se)
  
  # Create a data frame for the current outcome
  outcome_df <- data.frame(outcome = outcome,
                           over_rr = over_rr,
                           over_rr_se = over_rr_se,
                           over_rr_low = over_rr_low,
                           over_rr_high = over_rr_high,
                           estvar = estvar,
                           estmean = estmean)
  
  # Store the data frame in the list
  outcome_data[[paste0(outcome, "_cumlag")]] <- outcome_df
}

data_frames <- unname(outcome_data)
merged_df <- do.call(rbind, data_frames)
print(merged_df)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Commercial")
write.csv(merged_df, "Results/cumulative_lag_Commercial.csv")

#### OTHER RACE ####

data <- read_parquet("Data/Heatwave_Metrics.parquet")%>%
  setDT()%>%
  filter(Zip >= 27006 & Zip <= 28909)%>%
  rename(date=Date,
         heatwave=Heatwave,
         zip=Zip)%>%
  mutate(X=row_number(),
         year=year(date),
         doy=yday(date),
         dow=wday(date))%>%
  dplyr::select(date, heatwave, year, doy, dow, zip)%>%
  mutate(heatwave=as.integer(heatwave)) %>%
  filter(date >= "2016-01-01" & date <= "2019-12-31")%>%
  filter(month(date) >= 5 & month(date) <= 9)

datagrp <- data %>%
  group_by(zip)%>%
  summarize(hw=sum(heatwave))

#DELIVERIES
deliv <- read_parquet("Data/Deliv_Heatwaves.parquet")%>%
  setDT()%>%
  mutate(date=as.Date(admitdt),
         zip=as.numeric(ptzip))%>%
  filter(zip >= 27006 & zip <= 28909)%>%
  filter(Race != "Black" & Race != "White") %>%
  dplyr::select(zip, date, SMM21, GDM, HDP, SMI, PTB, MDP, PMAD, SUB)%>%
  filter(month(date) >= 5 & month(date) <= 9)%>%
  filter(year(date) >= 2016 & year(date) <= 2019)

#POP
pop <- read.csv("Data/TotalPop_F_18_44_2016_2021.csv")%>%
  mutate(denom=as.numeric(denom))

#CROSSWALK
cw <- read.csv("Data/Crosswalk_NC.csv")%>%
  setDT()%>%
  rename(zip=ZIP)

# Merge data

deliv <- deliv %>%
  left_join(cw, by=c('zip'))%>%
  dplyr::select(-zip, ZIP_TYPE)%>%
  rename(zip=ZCTA)

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

data <- data %>%
  left_join(delivzip, by=c('zip', 'date'))

data[is.na(data)] <- 0

data <- data %>%
  left_join(pop, by=c('zip', 'year'))

table(is.na(data))

data <- data %>% # Remove zip codes where denom = NA
  dplyr::filter(zip != 28583 & zip != 28647 & zip != 28018) 

datagrp <- data %>%
  group_by(zip)%>%
  summarize(hw=sum(heatwave))

datagrp <- datagrp %>%
  filter(hw==0)
no_hw <- datagrp$zip

data <- data %>%
  dplyr::filter(!zip %in% no_hw)

data <- data %>%
  mutate(sept_heatwave=ifelse(heatwave==1 & date %in% c("2016-09-23","2016-09-24","2016-09-25", "2016-09-26", "2016-09-27", "2016-09-28", "2016-09-29", "2016-09-30", 
                                                        "2017-09-23","2017-09-24","2017-09-25", "2017-09-26", "2017-09-27", "2017-09-28", "2017-09-29", "2017-09-30", 
                                                        "2018-09-23","2018-09-24","2018-09-25", "2018-09-26", "2018-09-27", "2018-09-28", "2018-09-29", "2018-09-30", 
                                                        "2019-09-23","2019-09-24","2019-09-25", "2019-09-26", "2019-09-27", "2019-09-28", "2019-09-29", "2019-09-30"), 1, 0))

data <- data %>%
  mutate(old_heatwave=heatwave)

data <- data %>%
  mutate(heatwave=ifelse(old_heatwave==1 & sept_heatwave==0, 1, 0))%>%
  dplyr::select(-old_heatwave, -sept_heatwave)

# Match days

dat<-data
head(dat)

zip_list <- unique(dat$zip) 
setorder(dat, zip, date)

gc()

#Match days

for (i in 1:length(zip_list)) {
  df <- subset(dat, zip == zip_list[i])
  
  # exclude the 3 days within any other heatwave
  df$time <- 1:nrow(df)
  cand_control <- unique(c(which(df$heatwave == 1), which(df$heatwave == 1) + 1, which(df$heatwave == 1) - 1))
  df$cand_control <- TRUE
  df$cand_control[cand_control[cand_control <= nrow(df)]] <- FALSE # Ensure cand_control indices are within df's range
  
  case_dates <- subset(df, heatwave == 1)
  control_dates <- subset(df, heatwave == 0)
  
  for (j in 1:nrow(case_dates)) {
    # choose lags (lagged 0 to 3)
    lag_dates <- case_dates[j, ]$date + 0:3
    lag_case <- subset(df, date %in% lag_dates)
    
    # choose 10 comparable unexposed days for each heatwave-exposed day
    control_range <- case_dates[j, ]$doy + -3:3
    control_subset <- subset(control_dates,
                             control_dates$year != case_dates[j, ]$year &
                               doy %in% control_range &
                               cand_control)
    controls <- dplyr::sample_n(control_subset, 3) #Changed from 10 to 3
    
    # choose lagged days for selected unexposed days
    la_con <- c(1:7)
    for (p in 1:length(la_con)) {
      lag_control_dates <- controls$date + la_con[p]  # Fixed index from j to p
      lag_control_each <- subset(df, date %in% lag_control_dates)
      
      if (p == 1) {
        lag_control <- lag_control_each
      } else {
        lag_control <- rbind(lag_control, lag_control_each)
      }
    }
    j_stratum <- rbind(lag_case, controls, lag_control)
    stratum <- paste("stratum", j, sep = ".")
    j_stratum$stratum <- stratum
    status <- c(rep("case", nrow(lag_case)), rep("control", nrow(controls)), rep("control", nrow(lag_control)))
    j_stratum$status <- status
    # Adjusted the length calculation for the lag column
    lag <- c(rep(0:3, length.out = nrow(lag_case)), rep(0, length.out = nrow(controls)), rep(c(1:3), length.out = nrow(lag_control)))
    j_stratum$lag <- lag
    
    if (j == 1) {
      new_df <- j_stratum
    } else {
      new_df <- rbind(new_df, j_stratum)
    }
  }
  if (i == 1) {
    matched_df <- new_df
  } else {
    matched_df <- rbind(matched_df, new_df)
  }
}

gc()

# We use "dlnm" package to generate the distributed lag function for "heatwave" 
library(lme4); library(dlnm); library(dplyr); library(stats)

for(i in 1:length(zip_list)){
  orig_dat <- subset(dat, zip == zip_list[i])
  match_dat <- subset(matched_df, zip == zip_list[i])
  
  orig_cb <- dlnm::crossbasis(orig_dat$heatwave, lag = c(0,3),
                              argvar = list(fun = "lin"),
                              arglag = list(fun = "integer"))
  obs_n <- nrow(orig_dat)
  orig_cb_matr <- as.data.frame(subset(orig_cb, nrow = obs_n))
  orig_cb_matr$date <- orig_dat$date
  matched_date <- match_dat %>% dplyr::select(date)
  matched_cb_matr <- orig_cb_matr %>%
    dplyr::right_join(matched_date, by = "date") %>%
    dplyr::select(-date) %>% as.matrix()
  
  if(i == 1){
    matched_cb_matrix <- matched_cb_matr
  }else{
    matched_cb_matrix <- rbind(matched_cb_matrix, matched_cb_matr)
  }
  # add attributes
  matched_dim <- dim(matched_cb_matrix)
  attr <- attributes(orig_cb)
  attr$dim <- matched_dim
  matched_cb <- matched_cb_matrix
  attributes(matched_cb) <- attr
}

gc()

# Check for incomplete strata

table(matched_df$lag)

restricted_df <- matched_df %>%
  group_by(zip, stratum)%>%
  summarize(n=n())
setorder(restricted_df, zip, stratum)

### GET CUMULATIVE LAG ###

outcomes <- c("gdm", "ptb", "hdp", "smi", "pmad", "mdp", "sub", "smm21")
outcome_data <- list()

for (outcome in outcomes) {
  formula <- as.formula(paste0(outcome, " ~ matched_cb + factor(dow) + year * zip"))
  
  fit <- gnm::gnm(formula,
                  eliminate = factor(zip), family = quasipoisson(link = "log"),
                  data = matched_df, offset = log(denom))
  
  pred <- dlnm::crosspred(matched_cb, fit, at = 1)
  
  over_rr <- sum(pred$matRRfit) / 4
  
  library(msm)
  estvar <- pred$vcov
  estmean <- c(pred$coefficients)
  
  over_rr_se <- msm::deltamethod(~ (exp(x1) + exp(x2) + exp(x3) + exp(x4)) / 4, 
                                 estmean, estvar)
  over_rr_low <- over_rr / exp(1.96 * over_rr_se)
  over_rr_high <- over_rr * exp(1.96 * over_rr_se)
  
  # Create a data frame for the current outcome
  outcome_df <- data.frame(outcome = outcome,
                           over_rr = over_rr,
                           over_rr_se = over_rr_se,
                           over_rr_low = over_rr_low,
                           over_rr_high = over_rr_high,
                           estvar = estvar,
                           estmean = estmean)
  
  # Store the data frame in the list
  outcome_data[[paste0(outcome, "_cumlag")]] <- outcome_df
}

data_frames <- unname(outcome_data)
merged_df <- do.call(rbind, data_frames)
print(merged_df)

merged_df <- merged_df %>%
  dplyr::select(outcome, over_rr, over_rr_se, over_rr_low, over_rr_high, estmean) 
merged_df <- merged_df %>%
  dplyr::mutate(Subgroup="Other Race")
write.csv(merged_df, "Results/cumulative_lag_Other.csv")

other_cumlag <- read.csv("Results/cumulative_lag_Other.csv")
medicaid_cumlag <- read.csv("Results/cumulative_lag_Medicaid.csv")
commercial_cumlag <- read.csv("Results/cumulative_lag_Commercial.csv")
full_cumlag <- read.csv("Results/cumulative_lag_full.csv")
age_cumlag <- read.csv("Results/cumulative_lag_Age.csv")
black_cumlag <- read.csv("Results/cumulative_lag_Black.csv")
white_cumlag <- read.csv("Results/cumulative_lag_White.csv")
hispanic_cumlag <- read.csv("Results/cumulative_lag_Hispanic.csv")

forest_data <- rbind(full_cumlag, age_cumlag, black_cumlag, hispanic_cumlag, white_cumlag,  other_cumlag,  medicaid_cumlag, commercial_cumlag)

forest_data <- forest_data %>%
  dplyr::select(-X, -estmean)%>%
  distinct(over_rr, over_rr_se, over_rr_high, over_rr_low, outcome, Subgroup)%>%
  rename(Outcome=outcome)

outcomes <- c("gdm", "ptb", "hdp", "pmad", "mdp", "sub", "smm21")

for (outcome in outcomes) {
  dt <- forest_data %>%
    filter(Outcome == outcome)%>%
    dplyr::select(-Outcome)
  
  dt$`RR (95% CI)` <- ifelse(is.na(dt$over_rr_se), "",
                             sprintf("%.2f (%.2f to %.2f)",
                                     dt$over_rr, dt$over_rr_low, dt$over_rr_high))
  
  dt <- dt %>%
    dplyr::select('Subgroup', 'RR (95% CI)', 'over_rr', 'over_rr_low', 'over_rr_high', 'over_rr_se')
  
  dt$` ` <- paste(rep(" ", 40), collapse = " ")
  est = dt$over_rr
  lower = dt$over_rr_low
  upper = dt$over_rr_high
  sizes = dt$over_rr_se
  
  dt <- dt %>%
    dplyr::select(-over_rr, -over_rr_high, -over_rr_low, -over_rr_se)
  
  p <- forest(dt,
              est = est,
              lower = lower, 
              upper = upper,
              sizes = sizes,
              ci_column = 3,
              ref_line = 1,
              xlim = c(0.5, 1.5),
              ticks_at = c(0.5, 1, 1.5),
              title = paste0(outcome))
  
  png(paste0(outcome, "_plot_indiv.png"))
  plot(p)  # Replace with your plot code
  dev.off()
  
}

#plot smi separately because different xlims

dt <- forest_data %>%
  filter(Outcome=="smi")%>%
  dplyr::select(-Outcome)

dt$`RR (95% CI)` <- ifelse(is.na(dt$over_rr_se), "",
                           sprintf("%.2f (%.2f to %.2f)",
                                   dt$over_rr, dt$over_rr_low, dt$over_rr_high))

dt <- dt %>%
  dplyr::select('Subgroup', 'RR (95% CI)', 'over_rr', 'over_rr_low', 'over_rr_high', 'over_rr_se')

dt$` ` <- paste(rep(" ", 40), collapse = " ")
est = dt$over_rr
lower = dt$over_rr_low
upper = dt$over_rr_high
sizes = dt$over_rr_se

dt <- dt %>%
  dplyr::select(-over_rr, -over_rr_high, -over_rr_low, -over_rr_se)

p <- forest(dt,
            est = est,
            lower = lower, 
            upper = upper,
            sizes = sizes,
            ci_column = 3,
            ref_line = 1,
            xlim = c(0.5, 4),
            ticks_at = c(0.75, 1, 2, 4),
            title = "smi")
png("smi_plot_indiv.png")
plot(p)  # Replace with your plot code
dev.off()
