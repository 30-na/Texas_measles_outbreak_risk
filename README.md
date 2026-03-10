# Measles Outbreak Risk Modeling in Texas Counties

This repository contains a framework for modeling the risk of measles outbreaks across counties in Texas. The approach uses data on vaccination rates, population sizes, and human mobility, and evaluates the impact of several vaccination improvement strategies.

## Objectives

- Estimate the probability that an outbreak in one county will trigger secondary outbreaks in others.
- Evaluate how vaccination improvements reduce outbreak risk.
- Visualize risk using geographic and statistical plots.

## Data Sources

#### Vaccination Coverage
- **Source**: Texas Department of State Health Services  
- **Link**: [MMR Vaccination Coverage (DSHS)](https://www.dshs.texas.gov/immunizations/data/school/coverage)  
- File used: `2023–2024_School_Vaccination_Coverage_Levels_Kindergarten.xlsx`

#### County Population
- **Probable Source**: Texas Legislative Council (2020 Census Redistricting Data)  
- **Link**: [Texas Capitol Data Portal](https://data.capitol.texas.gov/dataset/vtds)  
- File used: `Counties_Pop.txt`  
- Contains population data per county.

#### Mobility Flow Data
- **Source**: GeoDS COVID-19 US Flows Repository  
- **Link**: [https://github.com/GeoDS/COVID19USFlows](https://github.com/GeoDS/COVID19USFlows)  
- Weekly county-to-county flows were downloaded using:
 
```bash 

python download_weekly_data.py --start_year 2019 --start_month 1 --start_day 7 
--end_year 2019 \--end_month 12 --end_day 30 
--output_folder RawData/weekly_flows --county
```

#### Geometries Data
- **County Geometries**: TIGER/Line shapefiles via the `tigris` R package



## Vaccination Strategies

Four vaccination scenarios are analyzed:

- **Strategy 0 (Baseline)**: Uses the actual vaccination rates.
- **Strategy 1**: Sets all counties with MMR < 90% to exactly 90%.
- **Strategy 2**: Sets all counties with MMR < 92% to exactly 92%.
- **Strategy 3**: Increases all MMR rates by 5%, capped at 100%.

These strategies are implemented in `calculate_outbreak_probability.R`.

## Modeling Process

- **Susceptible Population**: Calculated from MMR coverage and vaccine efficacy (97%).
- **Outbreak Probability**: Probability of large outbreaks in each county.
- **Internal Infection Probability**: Solves an implicit equation for the attack rate.
- **Contact Matrices**:
  -  Gravity model based on population and distance
  -  Empirical 2019 mobility flow data

- **Transmission Probability**:
  \( p_{ij} = 1 - (1 - q \cdot P_j^I \cdot p_i^s \cdot p_i^{mo})^{C_{ij}} \)  
  where \( q = 0.9 \)

## Main Scripts

- `read_clean_merge_data.R`:  
  Prepares the cleaned datasets, merges vaccination and population info, computes distance matrices, and saves them.

- `calculate_outbreak_probability.R`:  
  Implements the epidemiological model and calculates transmission probabilities under each strategy and method.

- `plots.R`:  
  Generates maps and histograms for transmission probability and population at risk, under each strategy and method.

## Outputs

- **ProcessedData/**:
  - `map_county.rds`: County-level merged data
  - `distance_matrix_haversine_county.rds`: Distance matrix
  - `pij_M*.rds`: Transmission matrices for all methods and strategies

- **Figures/county_transmission/**:
  - Transmission risk maps per county
  - Histograms of risk values and population at risk


