#===== Script to remove nursing home population from general population counts =====#
rm(list=ls())

# population
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR/data/population')
poplist <- list()
countries <- list.files()
for(k in unique(countries)) poplist[[k]] <- read.csv(k)

# nursing home population
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR/data/NHpopulation')
NHpop <- list()
NH <- list.files()
for(k in unique(NH)) NHpop[[k]] <- read.csv(k)
countriesNH <- paste(countries, '')

c <- vector()
for(x in 1:length(countries)) c[x] <- substr(countries[x],1,nchar(countries[x])-4) 

cnh <- paste(c, '-nursing-homes.csv')

setdiff(cnh, NH)
