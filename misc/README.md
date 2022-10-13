## Age-specific weekly deaths data reported by the CDC 
The age-specific deaths data across the US states in ```data/``` have been reported by the CDC [here](https://data.cdc.gov/NCHS/Provisional-COVID-19-Death-Counts-by-Sex-Age-and-S/9bhg-hcku), and continuously extracted and made publicly available [here](https://github.com/rearc-data/covid-19-death-counts-sex-age-state)

The data sets are gathered, cleaned and saved in ```process-CDC-data.R```.

## All-ages daily deaths reported by the JHU
The data set is downloaded, cleaned and saved in ```prepare_JHU_data.R```.

## Age-specific vaccination rates reported by the CDC
The age-specific vaccination rate data in ```data-vaccination/COVID-19_Vaccinations_in_the_United_States_Jurisdiction_220623.csv``` have been reported by the CDC [here](https://data.cdc.gov/Vaccinations/COVID-19-Vaccinations-in-the-United-States-Jurisdi/unsk-b7fc).

The data set is cleaned and saved in ```prepare_vaccination_data_byage.R```.

## All-ages vaccination rates reported by the CDC
The all-ages vaccination rate data in ```data-vaccination/COVID-19_Vaccination_Trends_in_the_United_States_National_and_Jurisdictional.csv``` have been reported by the CDC [here](https://data.cdc.gov/Vaccinations/COVID-19-Vaccination-Trends-in-the-United-States-N/rh2h-3yt2).

The data set is cleaned and saved in ```prepare_vaccination_data.R```.