
## This code reads in the file created in OutcomePrep_MA.R and creates tables for each outcome with demographics

## Call packages and set working directory 

library(tableone)
library(dplyr)
library(lubridate)
library(arrow)
library(data.table)

setwd("~/RStudio/EHF Matched Analysis/")

## Read in Deliveries_Outcomes.parquet file created in OutcomePrep_MA.R
## Should have select maternal health outcomes and demographic variables (race, ethnicity, age, insurance)

tabdat <- read_parquet("Data/Deliveries_Outcomes.parquet")

tabdat <- tabdat %>%
  dplyr::select(Shepsid, Date, Age, Race, Ethnicity, Insurance, SMM21, GDM, PTB, MDP, PMAD, SMI, SUB, HDP)

catVars <- c("Age", "Race", "Ethnicity", "Insurance")

full_tab <- CreateTableOne(vars = catVars, data = tabdat, factorVars = catVars)
full_tabmat <- print(full_tab, showAllLevels = TRUE, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, smd=FALSE)
write.csv(full_tabmat, file = "Results/tab1_full.csv")

GDM <- subset(tabdat, GDM=='1') #We are only interested in women who have PMAD/SMI/MDP
PTB <- subset(tabdat, PTB=='1')
MDP <- subset(tabdat, MDP=='1')
HDP <- subset(tabdat, HDP=='1')
SMM <- subset(tabdat, SMM21=='1')
PMAD <- subset(tabdat, PMAD=='1')
SMI <- subset(tabdat, SMI=='1')
SUB <- subset(tabdat, SUB=='1')

gdm_tab <- CreateTableOne(vars = catVars, data = GDM, factorVars = catVars)
ptb_tab <- CreateTableOne(vars = catVars, data = PTB, factorVars = catVars)
mdp_tab <- CreateTableOne(vars = catVars, data = MDP, factorVars = catVars)
hdp_tab <- CreateTableOne(vars = catVars, data = HDP, factorVars = catVars)
smm_tab<- CreateTableOne(vars = catVars, data = SMM, factorVars = catVars)
smi_tab<- CreateTableOne(vars = catVars, data = SMI, factorVars = catVars)
sub_tab<- CreateTableOne(vars = catVars, data = SUB, factorVars = catVars)
pmad_tab<- CreateTableOne(vars = catVars, data = PMAD, factorVars = catVars)

gdm_tabmat <- print(gdm_tab, showAllLevels = TRUE, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, smd=FALSE)
write.csv(gdm_tabmat, file = "Results/gdm_cases_tab1.csv")

ptb_tabmat <- print(ptb_tab, showAllLevels = TRUE, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, smd=FALSE)
write.csv(ptb_tabmat, file = "Results/ptb_cases_tab1.csv")

mdp_tabmat <- print(mdp_tab, showAllLevels = TRUE, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, smd=FALSE)
write.csv(mdp_tabmat, file = "Results/mdp_cases_tab1.csv")

hdp_tabmat <- print(hdp_tab, showAllLevels = TRUE, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, smd=FALSE)
write.csv(hdp_tabmat, file = "Results/hdp_cases_tab1.csv")

smm_tabmat <- print(smm_tab, showAllLevels = TRUE, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, smd=FALSE)
write.csv(smm_tabmat, file = "Results/smm_cases_tab1.csv")

smi_tabmat <- print(smi_tab, showAllLevels = TRUE, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, smd=FALSE)
write.csv(smi_tabmat, file = "Results/smi_cases_tab1.csv")

pmad_tabmat <- print(pmad_tab, showAllLevels = TRUE, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, smd=FALSE)
write.csv(pmad_tabmat, file = "Results/pmad_cases_tab1.csv")

sub_tabmat <- print(sub_tab, showAllLevels = TRUE, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, smd=FALSE)
write.csv(sub_tabmat, file = "Results/sub_cases_tab1.csv")

