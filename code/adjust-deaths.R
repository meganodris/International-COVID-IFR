#===== Script to adjust for nursing home deaths =====#
library(readxl)


# source functions
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR')
source('./code/FunctionsForIFR.R')

# death data
df <- read.csv('./data/deaths_age.csv')
countries <- sort(unique(df$country))
df <- agg_deaths80p(df, countries)

# death location data
loc <- read_excel('./data/deaths_location.xlsx', sheet='Sheet1')
loc$pLTC <- loc$LTC_deaths/loc$total_deaths
loc <- loc[!(loc$country %in% c('England','Wales','Geneva')),]

# nursing home population
NHpop <- list()
NH <- list.files('./data/population/nursing-home-population/')
for(k in unique(NH)) NHpop[[k]] <- read.csv(paste('./data/population/nursing-home-population/', k, sep=''))
for(k in 1:length(NHpop)) names(NHpop)[k] <- substr(names(NHpop)[k],1,nchar(names(NHpop)[k])-4)

# remove nursing home deaths
new <- list()
for(c in unique(loc$country)){
  
  # extract death & nursing home data
  dfc <- df[df$country==c, ]
  NHc <- NHpop[[c]]
  
  # number of deaths associated with LTC settings
  LTC <- sum(dfc$deaths)*loc$pLTC[loc$country==c]
  
  if(dfc$sex[1]=='B'){
    
    NHc$B <- NHc$M + NHc$`F` # sum M + F
    W <- NHc$B/sum(NHc$B) # weights
    
    for(k in which(dfc$age_max>65)){ # align age groupings
      minA <- which(dfc$age_min[k]==NHc$age_min)
      maxA <- which(dfc$age_max[k]==NHc$age_max)
      if(length(minA)==0) minA <- 1
      dfc$deaths[k] <- dfc$deaths[k] - round(sum(W[minA:maxA])*LTC)
    }
  }else{
    
    Wf <- NHc$`F`/sum(NHc$`F` + NHc$M) # weights F
    Wm <- NHc$M/sum(NHc$`F` + NHc$M) # weights M
    
    for(k in which(dfc$age_max>65)){ # align age groupings
      minA <- which(dfc$age_min[k]==NHc$age_min)
      maxA <- which(dfc$age_max[k]==NHc$age_max)
      if(length(minA)==0) minA <- 1
      if(dfc$sex[k]=='M') dfc$deaths[k] <- dfc$deaths[k] - round(sum(Wm[minA:maxA])*LTC)
      if(dfc$sex[k]=='F') dfc$deaths[k] <- dfc$deaths[k] - round(sum(Wf[minA:maxA])*LTC)
    }
  }
  
  new[[c]] <- dfc
}

# output adjusted death data csv
adjusted_deaths <- do.call('rbind', new)
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR/data')
write.csv(adjusted_deaths, 'deaths_age_adjusted.csv', row.names=F)
