# 01861979-math70076-assessment-2

## Project Status

Completed.

## Project Description

This project presents an [in-depth report](reports/landfall-frequency-report.pdf) analysing tropical cyclone landfall patterns across the
American continent. Our analysis uses Bayesian Poisson regression models and cyclone track
data from 2000 to 2024. We explore and forecast landfall frequency per month and year, and per country/territory from 2000-2029. 
Intended users include humanitarian aid organisations and government policy makers.

## Project Structure

This Project contains 6 main folders, each containing sub directories and files of their own.
To ensure code scripts run smoothly, it is crucial to run the code files in the *src/* folder first, followed
by the files in the *analyses/* folder (in increasing numerical order, 00 to 07).

```
├── 01861979-math70076-assessment-2.Rproj 
├── README.md
├── LICENSE.md
├── data
│   ├── derived
│   │   ├── hurricane-data.csv
│   │   ├── landfall-data.csv
│   │   ├── stan-landfall-monthly-freq.csv
│   │   └── stan-landfall-per-country-freq.csv
│   └── raw
│       ├── atlantic-data.txt
│       └── pacific-data.txt
├── src
│   ├── check-model-diagnostics.R
│   └── check-posterior-predictions.R
└── tests
    ├── testthat
    │   ├── test-check-model-diagnostics.R
    │   └── test-check-posterior-predictions.R
    └── testthat.R
├── analyses
│   ├── 00_wrangle-hurricane-data.R
│   ├── 01_create-landfalls-dataset.R
│   ├── 02_eda-hurricane-data.R
│   ├── 03_model-landfall-freq-hsgp.R
│   ├── 04_model-landfall-freq-no-year.R
│   ├── 05_find-best-landfall-freq-model.R
│   ├── 06_analyse-landfalls-per-country.R
│   └── 07_forecast-landfalls.R
├── outputs
│   ├── bayesian-analysis-country-freq
│   │   ├── country-model-diagnostics.html
│   │   ├── country-model-diagnostics.png
│   │   ├── country-post-pred-checks.pdf
│   │   ├── landfall-per-country-density-plots.pdf
│   │   ├── landfalls-per-country-forecasts.html
│   │   ├── landfalls-per-country-forecasts.png
│   │   ├── simple-landfalls-per-country-forecasts.html
│   │   └── simple-landfalls-per-country-forecasts.png
│   ├── bayesian-analysis-monthly-freq
│   │   ├── best-model-comp.html
│   │   ├── best-model-comp.png
│   │   ├── landfall-monthly-density-plots.pdf
│   │   ├── landfalls-monthly-forecasts.html
│   │   ├── landfalls-monthly-forecasts.png
│   │   ├── model-HSGP
│   │   │   ├── HSGP-model-diagnostics.html
│   │   │   ├── HSGP-model-diagnostics.png
│   │   │   └── HSGP-post-pred-checks.pdf
│   │   ├── model-no-year-effect
│   │   │   ├── no-year-model-diagnostics.html
│   │   │   ├── no-year-model-diagnostics.png
│   │   │   └── no-year-post-pred-checks.pdf
│   │   ├── simple-landfalls-monthly-forecasts.html
│   │   └── simple-landfalls-monthly-forecasts.png
│   ├── eda-hurricane-data
│   │   ├── avg-landfalls-per-month.pdf
│   │   ├── landfall-countries.pdf
│   │   ├── landfall-locations.pdf
│   │   ├── landfall-map.pdf
│   │   ├── landfalls-per-year.pdf
│   │   └── nbr-cyclone-landfalls.pdf
│   └── stan-models
│       ├── country-model-cmdstanr.rds
│       ├── hsgp-model-cmdstanr.rds
│       ├── model_19923362a52f46d4c2c225083ab891bf
│       ├── model_19923362a52f46d4c2c225083ab891bf.stan
│       ├── model_59891ce15e0ed7c11692029e3b681abd.stan
│       ├── model_6652b09b8928d3e5214f509e7780b627
│       ├── model_6652b09b8928d3e5214f509e7780b627.stan
│       ├── model_8594d4def8fb7faea2aac6bfce509198
│       ├── model_8594d4def8fb7faea2aac6bfce509198.stan
│       ├── model_cb326800f9c69b887ad71221288b89a0
│       ├── model_cb326800f9c69b887ad71221288b89a0.stan
│       ├── model_d9a2196c112e560002299851d393c465
│       ├── model_d9a2196c112e560002299851d393c465.stan
│       ├── model_dfbc62b3c4e43050cf00e271d4496af2
│       ├── model_dfbc62b3c4e43050cf00e271d4496af2.stan
│       └── no-year-effect-model-cmdstanr.rds
└── reports
    ├── harvard-imperial-college-london.csl
    ├── landfall-frequency-report.aux
    ├── landfall-frequency-report.pdf
    ├── landfall-frequency-report.Rmd
    ├── landfall-frequency-report.tex
    └── references.bib

```

## Packages and Dependencies

This project is developed in R and requires the installation of the packages listed below:

- here
- readr
- stringr
- rnaturalearth
- sf
- dplyr
- data.table
- magrittr
- tidyr
- ggplot2
- scales
- cowplot
- ggsci 
- hexbin 
- bayesplot 
- kableExtra
- cmdstanr 
- webshot2
- loo

For Stan and cmdstanr installation, please refer to the following link:

- https://mc-stan.org/cmdstanr/

## Methodology

We modeled the number of cyclone landfalls per month and year using Bayesian
Poisson regression models and the **cmdstanr package**, in order to explore any seasonality effects and 
yearly trends. Our first model, Model A, included both a monthly effect and a non-linear yearly effect,
due to the high fluctuations in cyclone landfall numbers per year. Our second model, Model B, 
only included the monthly effect, ignoring any long-term trends. 
We then used Leave-One-Out Cross Validation (LOO) to
select the model with the best performance. 

Our final model, Model C, then investigated where cyclones were most likely to land, amongst the top 
10 countries and territories which have experienced the most cyclone landfalls since
2000. These include Antigua and Barbuda, the Bahamas, Belize, Canada, Cuba, 
the Dominican Republic, Mexico, Nicaragua, Puerto Rico and the United States of 
America.

Using our best performing time-based and geographical models, 
we then simulated landfall counts for each month as well as for each of the ten most 
affected countries and territories over the next five years (2025 to 2029). 


## Data

### Data Source

We based our analysis on the Atlantic hurricane database (HURDAT2), 1851-2024, 
more recently updated on April 4th, 2025 to include the 2024 hurricane season, 
as well as the Northeast and North Central Pacific hurricane database (HURDAT2), 
1949-2024, most recently updated on March 17, 2025 to include the 2024 hurricane season. 
These data sets are provided by the [National Hurricane Center (NHC)](https://www.nhc.noaa.gov/data/), 
as part of the National Oceanic and Atmospheric Administration (NOAA). 

### Data Description

Both databases contain records of cyclone coordinates, maximum winds, central pressure, 
system status at 6 hour intervals, as well as additional records outside of these standard intervals 
if a cyclone makes landfall, unexpectedly changes status, intensity or reaches a peak in
terms of wind speed or pressure. The NHC provides detailed documentation of the 
[structure of the HURDAT2 databases](https://www.nhc.noaa.gov/data/hurdat/hurdat2-format-atl-1851-2021.pdf), 
including examples and terminology.

For the purpose of this report, the two databases were combined and the following 
features were extracted from cyclones in 2000 to 2024:

- storm ID;
- name, if available;
- day;
- month;
- year;
- time;
- record identifier;
- cyclone position (measured in latitude and longitude, converted to positive and negative decimals);
- system status;
- maximum sustained wind (in knots, converted to integers);
- minimum pressure (in millibars, converted to integers).

This data set was then used to identify all landfall events in the studied period, 
using the record identifier "L". We then combined the landfall data with geographical data, 
containing country outlines and coordinates (acquired using the **sf** and **rnaturalearth** packages).
In cases where landfall coordinates did not fall within a country or territory boundary, we mapped them to 
the nearest country. 



## Areas for Project Extension and Development

There are many potential directions for further development of this project, such as the inclusion of 
storm intensity, wind speed and pressure into the models as well as exploring differences between cyclones originating in the Atlantic and Pacific basins. 

## License

This project is licensed under the MIT License. \
Please refer to the [LICENSE](LICENSE.md) file for further details and copyright policies.

Copyright (c) 2025 Leyli Cowell

