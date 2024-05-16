# Matched-Analysis

- Data files *Deliveries_Outcomes.parquet*, *deliveries.sas7bdat* are not included in this repository but may be available upon request. 

- *Heatwave_Metrics.parquet* and heatwave documentation is available from *[Luke Wertis's Github](https://github.com/wertisml/Heatwave/blob/main/Data/Zip/Heatwave_Metrics.parquet)*

### Order as follows: 

## 1. OutcomePrep_MA.R
- Adds maternal outcome indicators to individual ED data using ICD-10 codes.
- **Input:** deliveries.sas7bdat, ICD10_Codes_Maternal_Outcomes.csv
- **Output:** Deliveries_Outcomes.parquet
## 2. (Optional) TableOne_MA.R
- Creates tables of outcomes with case demographics.
- **Input:** Deliveries_Outcomes.parquet
## 3. DataPrep_MA.R
- Aggregates delivery data to the daily ZCTA level and merges with heatwave data.
- **Input:** Heatwave_Metrics.parquet, Deliveries_Outcomes.parquet, Crosswalk_NC.csv, TotalPop_F_18_44_2016_2021.csv
- **Output:** Deliv_Heatwaves.parquet
## 4. Matched_Example.R
- Matches each heatwave day and lag period to unexposed non-heatwave control period in the same ZCTA.
- Creates matched dataframe and matched crossbasis for modeling.
- Outputs relative risk for selected models for each lag day.
- **Input:** Deliv_Heatwaves.parquet
- **Output:** Combined_Matrices_Example_MA.csv
## 5. ForestPlot3.R
- Plot cumulative lags with different subgroups (Black, Hispanic, Medicaid, etc.)
- **Input:** Forest_Plot_Results.csv (created externally)
## 6. (Optional) Create_Shapefile_MA.R
- Aggregates maternal outcomes and heatwave events at the ZCTA level and merges with a shapefile for mapping.
- **Input:** Deliv_Heatwaves.parquet, tl_2020_us_zcta510.shp
- **Output:** Outcomes_HW_ZCTA.shp
