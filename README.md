# International-COVID-IFR

Code and data to reproduce the analysis in "Age-specific mortality and immunity patterns of SARS-CoV-2 infection in 45 countries".

### **_Code_**
- ```RunModel.R``` is the main RScript to call all data, functions & model fitting.
- ```IFRInternational.stan``` contains the core model code.
- ```IFRInternational_singlesero.stan``` contains core model code for use when only a single seroprevalence datapoint is being included in the model likelihood.
- ```adjust-deaths.R``` estimates age-specific non-nursing home COVID-19 deaths aged 65+ for a subset of countries.



### **_Data_**
- ```deaths_age.csv``` contains reported age-specific COVID-19 death data from 45 countries.
- ```deaths_location.csv``` contains the number of COVID-19 deaths associated with nursing home/long-term care facilities for a subset of countries.
- ```deaths_age_adjusted.csv``` is the estimated age-specific number of COVID-19 deaths for a subset of countries where nursing home/long-term care deaths
could be excluded.
- ```time_series_covid19_deaths_global.csv``` is the daily time-series of COVID-19 deaths from the COVID-19 Data Repository 
by the Center for Systems Science and Engineering (CSSE) at Johns Hopkins University (https://github.com/CSSEGISandData/COVID-19).
- ```deathsT_region.csv``` contains daily time-series of COVID-19 deaths for a subset of locations, where unavailable from the Johns Hopkins repository.
- ```serostudies.csv``` are the results from each of the seroprevalence surveys included in the analysis.
- ```data_sources.csv``` contains URL links to sources of age-specific death data for each country.
- ```France_NursingHomePop.csv``` is the nursing home population demographics for France.
