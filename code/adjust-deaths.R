#===== Script to adjust for nursing home deaths =====#


# source functions
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR')
source('./code/FunctionsForIFR.R')

# death data
df <- read.csv('./data/deaths_age.csv')

# death location data
loc <- read.csv('./data/deaths_location.csv')
loc$pLTC <- loc$LTC_deaths/loc$total_deaths


#--- Adjustment to remove nursing home deaths
# assumes proportion of COVID-19 LTC deaths is constant in time
# assumes age-distribution of LTC deaths to be the same as all COVID-19 deaths
new <- list()
for(c in unique(loc$country)){
  
  # country-specific death data
  dfc <- df[df$country==c, ]
  
  # number of deaths associated with LTC settings
  LTC <- sum(dfc$deaths)*loc$pLTC[loc$country==c]
  
  if(dfc$sex[1]=='B'){
    
    W <- dfc$deaths[dfc$age_max>65]/sum(dfc$deaths[dfc$age_max>65]) # proportions
    dfc$deaths[dfc$age_max>65] <- round(W*(sum(dfc$deaths[dfc$age_max>65]) - LTC))
    
  }else{
    
    Wf <- dfc$deaths[dfc$age_max>65 & dfc$sex=='F']/sum(dfc$deaths[dfc$age_max>65]) # proportions F
    Wm <- dfc$deaths[dfc$age_max>65 & dfc$sex=='M']/sum(dfc$deaths[dfc$age_max>65]) # proportions M
    dfc$deaths[dfc$age_max>65 & dfc$sex=='F'] <- round(Wf*(sum(dfc$deaths[dfc$age_max>65 & dfc$sex=='F']) - LTC))
    dfc$deaths[dfc$age_max>65 & dfc$sex=='M'] <- round(Wm*(sum(dfc$deaths[dfc$age_max>65 & dfc$sex=='M']) - LTC))
    
  }
  
  new[[c]] <- dfc
}

# output adjusted death data csv
adjusted_deaths <- do.call('rbind', new)
adjusted_deaths <- rbind(adjusted_deaths, df[df$country %in% c('France'), ]) # add France hospital only data
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR/data')
write.csv(adjusted_deaths, 'deaths_age_adjusted.csv', row.names=F)
