
## This code pulls selected maternal health outcomes from ED delivery data ("deliveries.sas7bdat") using ICD-10 codes. 

## Call packages and set working directory 

library(tableone)
library(gmodels)
library(dplyr)
library(lubridate)
library(arrow)
library(data.table)
library(dplyr)
library(geepack)
library(haven)
library(lubridate)
library(lme4)
setwd("~/RStudio/EHF Matched Analysis/")

## Read in data files

# Hospital delivery data for North Carolina
# Already includes SMM20 and SMM21, so we do not need to code for those outcomes

df_wrong <- read_sas("Data/deliveries.sas7bdat")%>%
  mutate(Date=as.Date(admitdt))%>%
  filter(year(Date) >= 2016)
df_wrong[is.na(df_wrong)] <- 0 # Replace NAs with 0

# File with the ICD-10 codes that will be used to define each outcome of interest
codes <- read.csv("Data/ICD10_Codes_Maternal_Outcomes.csv")

## The following section creates a list of ICD-10 codes for each outcome of interest, and then
## checks the diagnosis codes for each delivery to see if they match. 
## Creates a column that indicates (0/1) whether the condtion is present for the delivery

#### Gestational Diabetes Mellitus (GDM) ####

# Define the codes for Gestational Diabetes Mellitus without decimal points
gdm_codes <- codes %>%
  filter(condition=="GDM")%>%
  dplyr::select(ICD10)
gdm_codes <- gdm_codes$ICD10

# Function to check if a diagnosis code is related to GDM
GDM_code <- function(code) {
  any(code %in% gdm_codes)
}

# Loop through diagnostic code columns and check for GDM codes
df_wrong$GDM <- apply(df_wrong[, c("diag1", "diag2", "diag3", "diag4", "diag5", "diag6", "diag7", "diag8", "diag9", "diag10",
                                   "diag11", "diag12", "diag13", "diag14", "diag15", "diag16", "diag17", "diag18", "diag19",
                                   "diag20", "diag21", "diag22", "diag23", "diag24", "diag25")], 1, function(x) as.numeric(any(GDM_code(gsub("\\.", "", x)))))

# Check the result
head(df_wrong$GDM)

# Subset cases where GDM is equal to 1
cases_with_gdm <- df_wrong[df_wrong$GDM == 1, ]

# Take a look at the sample
head(cases_with_gdm)
table(cases_with_gdm$GDM)

#### Hypertensive Disorders of Pregnancy (HDP) ####

# Define the codes for Hypertensive Disorders of Pregnancy without decimal points
hdp_codes <- codes %>%
  filter(condition=="HDD")%>%
  dplyr::select(ICD10)
hdp_codes <- hdp_codes$ICD10

# Function to check if a diagnosis code is related to hdd
HDP_code <- function(code) {
  any(code %in% hdp_codes)
}

# Loop through diagnostic code columns and check for hdd codes
df_wrong$HDP <- apply(df_wrong[, c("diag1", "diag2", "diag3", "diag4", "diag5", "diag6", "diag7", "diag8", "diag9", "diag10",
                                   "diag11", "diag12", "diag13", "diag14", "diag15", "diag16", "diag17", "diag18", "diag19",
                                   "diag20", "diag21", "diag22", "diag23", "diag24", "diag25")], 1, function(x) as.numeric(any(HDP_code(gsub("\\.", "", x)))))

# Check the result
head(df_wrong$HDP)

# Subset cases where HDP is equal to 1
cases_with_hdp <- df_wrong[df_wrong$HDP == 1, ]

# Take a look at the sample
head(cases_with_hdp)
table(cases_with_hdp$HDP)

#### Maternal Mental Disorders of Pregnancy (MDP) ####

# Define the codes for maternal mental health without decimal points
mdp_codes <- codes %>%
  filter(condition=="MMDP")%>%
  dplyr::select(ICD10)
mdp_codes <- mdp_codes$ICD10

# Function to check if a diagnosis code is related to mmdp
MDP_code <- function(code) {
  any(code %in% mdp_codes)
}

# Loop through diagnostic code columns and check for mmdp codes
df_wrong$MDP <- apply(df_wrong[, c("diag1", "diag2", "diag3", "diag4", "diag5", "diag6", "diag7", "diag8", "diag9", "diag10",
                                   "diag11", "diag12", "diag13", "diag14", "diag15", "diag16", "diag17", "diag18", "diag19",
                                   "diag20", "diag21", "diag22", "diag23", "diag24", "diag25")], 1, function(x) as.numeric(any(MDP_code(gsub("\\.", "", x)))))

# Check the result
head(df_wrong$MDP)

# Subset cases where Is_mmdp is equal to 1
cases_with_mdp <- df_wrong[df_wrong$MDP == 1, ]

# Take a look at the sample
head(cases_with_mdp)
table(cases_with_mdp$MDP)

#### Perinatal Mood and Anxiety Disorder (PMAD) ####

# Define the codes for perinatal mood and anxiety disorder without decimal points
pmad_codes <- codes %>%
  filter(condition=="PMAD")%>%
  dplyr::select(ICD10)
pmad_codes <- pmad_codes$ICD10

# Function to check if a diagnosis code is related to sga
PMAD_code <- function(code) {
  any(code %in% pmad_codes)
}

# Loop through diagnostic code columns and check for sga codes
df_wrong$PMAD <- apply(df_wrong[, c("diag1", "diag2", "diag3", "diag4", "diag5", "diag6", "diag7", "diag8", "diag9", "diag10",
                                    "diag11", "diag12", "diag13", "diag14", "diag15", "diag16", "diag17", "diag18", "diag19",
                                    "diag20", "diag21", "diag22", "diag23", "diag24", "diag25")], 1, function(x) as.numeric(any(PMAD_code(gsub("\\.", "", x)))))

# Check the result
head(df_wrong$PMAD)

# Subset cases where PMAD is equal to 1
cases_with_pmad <- df_wrong[df_wrong$PMAD == 1, ]

# Take a look at the sample
head(cases_with_pmad)
table(cases_with_pmad$PMAD)

#### Severe or Serious Mental Illness (SMI) ####

# Define the codes for severe mental illness without decimal points
smi_codes <- codes %>%
  filter(condition=="SMI")%>%
  dplyr::select(ICD10)
smi_codes <- smi_codes$ICD10

# Function to check if a diagnosis code is related to sga
SMI_code <- function(code) {
  any(code %in% smi_codes)
}

# Loop through diagnostic code columns and check for sga codes
df_wrong$SMI <- apply(df_wrong[, c("diag1", "diag2", "diag3", "diag4", "diag5", "diag6", "diag7", "diag8", "diag9", "diag10",
                                   "diag11", "diag12", "diag13", "diag14", "diag15", "diag16", "diag17", "diag18", "diag19",
                                   "diag20", "diag21", "diag22", "diag23", "diag24", "diag25")], 1, function(x) as.numeric(any(SMI_code(gsub("\\.", "", x)))))

# Check the result
head(df_wrong$SMI)

# Subset cases where SMI is equal to 1
cases_with_smi <- df_wrong[df_wrong$SMI == 1, ]

# Take a look at the sample
head(cases_with_smi)
table(cases_with_smi$SMI)

#### Pre-Term Birth (PTB) ####

# Define the codes for preterm birth without decimal points
ptb_codes <- codes %>%
  filter(condition=="PTB")%>%
  dplyr::select(ICD10)
ptb_codes <- ptb_codes$ICD10

# Function to check if a diagnosis code is related to ptb
PTB_code <- function(code) {
  any(code %in% ptb_codes)
}

# Loop through diagnostic code columns and check for ptb codes
df_wrong$PTB <- apply(df_wrong[, c("diag1", "diag2", "diag3", "diag4", "diag5", "diag6", "diag7", "diag8", "diag9", "diag10",
                                   "diag11", "diag12", "diag13", "diag14", "diag15", "diag16", "diag17", "diag18", "diag19",
                                   "diag20", "diag21", "diag22", "diag23", "diag24", "diag25")], 1, function(x) as.numeric(any(PTB_code(gsub("\\.", "", x)))))

# Check the result
head(df_wrong$PTB)

# Subset cases where PTB is equal to 1
cases_with_ptb <- df_wrong[df_wrong$PTB == 1, ]

# Take a look at the sample
head(cases_with_ptb)
table(cases_with_ptb$PTB)

#### Substance Abuse (SUB) ####

# Define the codes for preterm birth without decimal points
sub_codes <- codes %>%
  filter(condition=="SUB")%>%
  dplyr::select(ICD10)
sub_codes <- sub_codes$ICD10

# Function to check if a diagnosis code is related to sub
SUB_code <- function(code) {
  any(code %in% sub_codes)
}

# Loop through diagnostic code columns and check for sub codes
df_wrong$SUB <- apply(df_wrong[, c("diag1", "diag2", "diag3", "diag4", "diag5", "diag6", "diag7", "diag8", "diag9", "diag10",
                                   "diag11", "diag12", "diag13", "diag14", "diag15", "diag16", "diag17", "diag18", "diag19",
                                   "diag20", "diag21", "diag22", "diag23", "diag24", "diag25")], 1, function(x) as.numeric(any(SUB_code(gsub("\\.", "", x)))))

# Check the result
head(df_wrong$SUB)

# Subset cases where SUB is equal to 1
cases_with_sub <- df_wrong[df_wrong$SUB == 1, ]

# Take a look at the sample
head(cases_with_sub)
table(cases_with_sub$SUB)

#### Create Race, Ethnicity, Age, and Insurance columns ####

tabdat <- df_wrong

tabdat$Race <- fifelse(tabdat$Race == "5", "White",
                       fifelse(tabdat$Race == "3", "Black",
                               fifelse(tabdat$Race == "2", "Asian",
                                       fifelse(tabdat$Race=="1","Indigenous American",
                                               fifelse(tabdat$Race=="4","Native Hawaiian/PI",
                                                       fifelse(tabdat$Race=="6", "Other",
                                                               "Unknown"))))))
tabdat$Race <- as.factor(tabdat$Race)

tabdat$Ethnicity <- fifelse(tabdat$Ethnicity == "1", "Not Hispanic",
                            fifelse(tabdat$Ethnicity == "2", "Hispanic",
                                    "Unknown"))
tabdat$Ethnicity <- factor(tabdat$Ethnicity)

tabdat$Age <- as.factor(ifelse(tabdat$agey < 20, '18-19', 
                               ifelse(tabdat$agey<25, '20-24',
                                      ifelse(tabdat$agey <30, '25-29',
                                             ifelse(tabdat$agey<35,'30-34',
                                                    ifelse(tabdat$agey<40, '35-39',
                                                           ifelse(tabdat$agey <= 44, '40+', 0)))))))
tabdat$Age <- factor(tabdat$Age)


tabdat$Insurance <- as.factor(ifelse(tabdat$payer1 == "09", "Self-pay",
                                     ifelse(tabdat$payer1 == "MC", "Medicaid",
                                            ifelse(tabdat$payer1 %in% c("MA", "MB", "OF", "VA", "TV", "11", "16", "CH", "MM"), "Other gov't", 
                                                   ifelse(tabdat$payer1 %in% c("LM", "LI", "HM", "DS", "CI", "BL", "AM", "12", "13", "14", "15"),
                                                          "Commercial", "Other/Unknown")))))
tabdat$Insurance <- factor(tabdat$Insurance)

## Write off dataframe with outcomes and new demographic columns as a parquet file for use in DataPrep_MA.R

write_parquet(tabdat, "Data/Deliveries_Outcomes.parquet")
