# International-COVID-IFR

Code and data to reproduce the analysis in "Age-specific mortality and immunity patterns of SARS-CoV-2 infection in 45 countries".

### **_Code_**
- 'RunModel.R' is the main RScript to call all data, functions & model fitting.
- 'IFRInternational.stan' contains the core model code.
- 'IFRInternational_singlesero.stan' contains core model code for use when only a single seroprevalence datapoint is being included in the model likelihood.
- 'adjust-deaths.R' estimates age-specific non-nursing home COVID-19 deaths aged 65+ for a subset of countries.



### **_Data_**
- 'deaths_age.csv' contains reported age-specific COVID-19 death data from 45 countries.
- 'deaths_location.csv' contains the number of COVID-19 deaths associated with nursing home/long-term care facilities for a subset of countries.
- 'deaths_age_adjusted.csv' is the estimated age-specific number of COVID-19 deaths for a subset of countries where nursing home/long-term care deaths
could be excluded.
- 'deathsT_region.csv' contains daily time-series of COVID-19 deaths for a subset of locations, where unavailable from the COVID-19 Data Repository 
by the Center for Systems Science and Engineering (CSSE) at Johns Hopkins University.
- 'serostudies.csv' are the results from each seroprevalence surveys included in the analysis.
- 'CharlesDeGaulle.csv' and 'DiamondPrincess.csv' inlcude the reported number of age-specific infections detected from active surveillance campaigns on the Charles de Gaulle aircraft 
carrier and Diamond Princess cruise ship.
- 'France_NursingHomePop.csv' is the nursing home population demographics of France.
