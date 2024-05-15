
## The following code is used to match heatwave days and lag periods to non-heatwave days and lag periods. 
## The matching procedure is based on Yan et al. (2021) https://pubmed.ncbi.nlm.nih.gov/33591048/#:~:text=Results%3A%20For%201999%2D2010%2C,elevated%20throughout%20the%20storm%20period.

## The data file used in this code should have the following columns or similar: 

## "zip": zip code or other geographic identifier
## "date": date of each day
## "year": year for each day
## "doy": day of year, used to find similar non-heatwave days in other years
## "dow": day of week
## "heatwave": an indicator variable for heatwave exposure, "heatwave" is 1 
##          if a day is exposed to storm and 0 for non-heatwave exposed days (or
##          other binary exposure of interest)
## "denom": relevant denominator for the health outcome of interest; this example uses the total female population
##          between age 18 and 44 for each ZCTA and year
## "births": total number of deliveries for each day for each ZCTA; unique to this analysis 
## "SMM21, GDM, HDP, SMI, PTB, MDP, PMAD, SUB": outcomes included in delivery data 


# Read in packages and set working directory
library(dplyr)
library(arrow)
library(lubridate)
library(data.table)
library(haven)
library(sjPlot)
library(MASS)

setwd("~/RStudio/EHF Matched Analysis/Example")
getwd()

#### MATCH HEATWAVE DAYS TO NON-EXPOSED DAYS ####

## Read in the parquet file created in DataPrep_MA.R

dat<- read_parquet("Data/Deliv_Heatwaves.parquet")
head(dat)

zip_list <- unique(dat$zip) 
setorder(dat, zip, date)

gc()

#### MATCHING PART ####
## This section will match every heatwave and the 7 days following each heatwave day with three non-heatwave, unexposed periods. 
## The matched periods will occur in the same ZCTA and during the same doy (+/- 3 days) but primarily during a different year

## This section creates our matched dataframe (matched_df)

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
    # choose lags (lagged 0 to 7)
    lag_dates <- case_dates[j, ]$date + 0:7
    lag_case <- subset(df, date %in% lag_dates)
    
    # choose 10 comparable unexposed days for each heatwave-exposed day
    control_range <- case_dates[j, ]$doy + -3:3
    control_subset <- subset(control_dates,
                             control_dates$year != case_dates[j, ]$year &
                               doy %in% control_range &
                               cand_control)
    controls <- dplyr::sample_n(control_subset, 10) # I used 3 in my thesis instead of 10
    
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
    lag <- c(rep(0:7, length.out = nrow(lag_case)), rep(0, length.out = nrow(controls)), rep(c(1:7), length.out = nrow(lag_control)))
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

## We use "dlnm" package to generate the distributed lag function for "heatwave" 
## This section creates the crossbasis (matched_cb) that will be used in our models

library(lme4); library(dlnm); library(dplyr); library(stats)

for(i in 1:length(zip_list)){
  orig_dat <- subset(dat, zip == zip_list[i])
  match_dat <- subset(matched_df, zip == zip_list[i])
  
  orig_cb <- dlnm::crossbasis(orig_dat$heatwave, lag = c(0,7),
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

## Check for incomplete strata to ensure the matching procedure worked correctly 
## Each stratum for each ZCTA consists of 8 case days (lag0 heatwave day plus lag1 to lag7) 
## and 80 control days (each case day matched to 10 control days) to create 88 observations for each stratum.
## These numbers will change if the lag periods or control parameters are adjusted

# Each lag should have an equal number of observations
table(matched_df$lag) 

# Each stratum should have an equal number of observations (n) 
restricted_df <- matched_df %>%
  group_by(zip, stratum)%>%
  summarize(n=n())
setorder(restricted_df, zip, stratum)
View(restricted_df)

#### MODELS ####

## Models use the matched dataframe (matched_df) and the crossbasis (matched_cb) 
## The original example in Yan et al. uses a mixed-effect model and includes additional alternative models as a sensitivity analysis
## In this section, you will have to manually replace for each outcome of interest (smm21, gdm, hdp, etc.)

# mixed-effect model 
fit <- lme4::glmer(smm21 ~ matched_cb + factor(year) + factor(dow) + (1|zip),
                   data = matched_df, offset=log(denom),
                   family = poisson(link = "log"),
                   control = glmerControl(optimizer = "bobyqa",
                                          optCtrl = list(maxfun = 2e5)))

pred <- dlnm::crosspred(matched_cb, fit, at = 1)

# calculate overdispersion parameter
blmeco::dispersion_glmer(fit)

# estimates of distributed relative risks（RR）for cardiovascular hospitalizations 
# on heatwave-exposed days compared with matched unexposed days 
# for all heatwaves and across the seven counties
pred$matRRfit

# lower and upper confidence intervals for distriuted RRs
pred$matRRlow; pred$matRRhigh

# estimate of overall RR for the entire heatwave-exposure period (i.e., ten days)
over_rr <- sum(pred$matRRfit)/10

# we used "delta method" to calculate confidence interval for the overall 
# RR for the entire heatwave period
library(msm)
estvar <- pred$vcov
estmean <- c(pred$coefficients)

over_rr_se <- msm::deltamethod(~ (exp(x1) + exp(x2) + exp(x3) + exp(x4) + 
                                    exp(x5) + exp(x6) + exp(x7) + 
                                    exp(x8) + exp(x9) + exp(x10))/10, 
                               estmean, estvar)
over_rr_low <- over_rr / exp(1.96*over_rr_se)
over_rr_high <- over_rr * exp(1.96*over_rr_se)

################# Sensitivity analyses ###################
# alternative model 1
fit_am1 <- glm(smm21 ~ matched_cb + factor(dow) + factor(year) + factor(zip),
               data = matched_df, offset = log(denom), 
               family = poisson(link = "log"))
pred_am1 <- dlnm::crosspred(matched_cb, fit_am1, at = 1)

# alternative model 2
matched_df$year <- as.numeric(matched_df$year)
matched_df$zip <- factor(matched_df$zip)
fit_am2 <- glm(smm21 ~ matched_cb + factor(dow) + year * zip,
               data = matched_df, offset = log(denom), 
               family = poisson(link = "log"))
pred_am2 <- dlnm::crosspred(matched_cb, fit_am2, at = 1)

# alternative model 3
matched_df$year <- as.numeric(matched_df$year)
fit_am3 <- lme4::glmer(smm21 ~ matched_cb + factor(dow) + (1 + year|zip) + year,
                       data = matched_df, offset = log(denom), family = poisson(link = "log"), 
                       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
pred_am3 <- dlnm::crosspred(matched_cb, fit_am3, at = 1)

# alternative model 4: conditional Poisson model
fit_am4 <- gnm::gnm(smm21 ~ matched_cb + factor(dow) + factor(year),
                    eliminate = factor(zip), family = poisson(link = "log"),
                    data = matched_df, offset=log(denom))
pred_am4 <- dlnm::crosspred(matched_cb, fit_am4, at = 1)


# alternative model 5: conditional quasiPoisson model
fit_am5 <- gnm::gnm(smm21 ~ matched_cb + factor(dow) + factor(year),
                    eliminate = factor(zip), family = quasipoisson(link = "log"),
                    data = matched_df, offset=log(denom))
pred_am5 <- dlnm::crosspred(matched_cb, fit_am5, at = 1)
summary(fit_am5)

# alternative 6: 
matched_df$obs <- factor(1:nrow(matched_df))
fit_am6 <- glm(smm21 ~ matched_cb + factor(year) + factor(dow) + (1|zip) + (1|obs),
               data = matched_df, offset = log(denom), 
               family = poisson(link = "log"),
               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
pred_am6 <- dlnm::crosspred(matched_cb, fit_am6, at = 1)

# alternative 7
fit_am7 <- MASS::glmmPQL(smm21 ~ matched_cb + factor(year) + factor(dow) + offset(log(denom)),
                         random = ~ 1|zip, data = matched_df, 
                         family = quasipoisson(link = "log"),
                         control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
pred_am7 <- dlnm::crosspred(matched_cb, fit_am7, at = 1)

pred_am7$matRRfit; pred_am7$matRRlow; pred_am7$matRRhigh

# Negative binomial model
nb_model <- glm.nb(smm21 ~ matched_cb
                   + factor(dow) + year * zip, data = matched_df)
summary(nb_model)

#### MODELS AS LOOP ####

## This section allows you to run all the selected maternal outcomes through different models as a for loop
## This saves a lot of time!
## This example includes a negative binomial, generalized nonlinear, and case-crossover model
## Outputs matRRfit, matRRlow, and matRRhigh values for each lag and each model as a table for each outcome
## We will use the values in these tables in our forest plots later

# Create list of outcomes to run through our loop
outcomes <- c("gdm", "ptb", "hdp", "smi", "pmad", "mdp", "sub", "smm21")

# Additional prep before the model loop 
datasets <- list()
combined_data_list <- list()

# Lines 276-285 are specific to the case crossover model
cb_crossover_heat <- crossbasis(matched_df$heatwave, lag = c(0, 7), argvar = list("lin"),
                                arglag = list("integer"))
matched_df$month <- month(matched_df$date)
matched_df$year <- year(matched_df$date)
matched_df$year <- factor(matched_df$year)
matched_df$month <- factor(matched_df$month)
matched_df$dow <- factor(matched_df$dow)
matched_df$stratum <- matched_df$year:matched_df$month:matched_df$dow
matched_df$year <- as.numeric(matched_df$year)

# Loop outcomes through the models
for (outcome in outcomes) {
  formula <- as.formula(paste0(outcome, " ~ matched_cb + factor(dow) + year * zip"))
  
  # Negative binomial
  nb_model <- glm.nb(formula, data = matched_df)
  pred_nb <- dlnm::crosspred(matched_cb, nb_model, at = 1)
 
  # GNM
  fit_am5 <- gnm::gnm(formula,
                      eliminate = factor(zip), family = quasipoisson(link = "log"),
                      data = matched_df, offset=log(denom))
  pred_am5 <- dlnm::crosspred(matched_cb, fit_am5, at = 1)
  
  # Case crossover
  ccformula <- as.formula(paste0(outcome, "~ cb_crossover_heat + stratum + year")) 
  fit<-glm(ccformula, 
           data = matched_df, offset = log(denom+1),
           family = quasipoisson(link = "log"), 
           control = glm.control(epsilon = 10E-8, maxit = 5000))
  
  pred_cc <- dlnm::crosspred(cb_crossover_heat , fit, at = 1)
 
  datasets[[paste0(outcome, "_nb")]] <- pred_nb
  datasets[[paste0(outcome, "_gnm")]] <- pred_am5
  datasets[[paste0(outcome, "_cc")]] <- pred_cc
  
  combined_data_name <- paste0(outcome, "_combined_data")
  
  # Create an empty data frame with dynamic column names
  max_columns <- max(sapply(datasets, function(x) ncol(x$matRRfit), simplify = TRUE),
                     sapply(datasets, function(x) ncol(x$matRRlow), simplify = TRUE),
                     sapply(datasets, function(x) ncol(x$matRRhigh), simplify = TRUE))
  
  combined_data <- data.frame(matrix(NA, nrow = 0, ncol = max_columns + 2))
  colnames(combined_data) <- c(names(datasets[[1]]$matRRfit), "Type", "Pred_Model")
  
  # Populate the combined data frame
  for (name in names(datasets)) {
    matRRfit <- datasets[[name]]$matRRfit
    matRRlow <- datasets[[name]]$matRRlow
    matRRhigh <- datasets[[name]]$matRRhigh
    
    combined_data <- rbind(combined_data, cbind(matRRfit, Type = "matRRfit", Pred_Model = name))
    combined_data <- rbind(combined_data, cbind(matRRlow, Type = "matRRlow", Pred_Model = name))
    combined_data <- rbind(combined_data, cbind(matRRhigh, Type = "matRRhigh", Pred_Model = name))
  }
  
  # Assign the combined data frame to the list with dynamic name
  combined_data_list[[combined_data_name]] <- combined_data
}

write.csv(combined_data, file = "Results/Combined_Matrices_Example_MA", row.names = FALSE)






















