#===== Script to remove nursing home population from general population counts =====#

setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR/data/population')


# general population
poplist <- list()
countries <- list.files('./general-population/')
for(k in unique(countries)) poplist[[k]] <- read.csv(paste('./general-population/', k, sep=''))
for(k in 1:length(poplist)) names(poplist)[k] <- substr(names(poplist)[k],1,nchar(names(poplist)[k])-9)


# nursing home population
NHpop <- list()
NH <- list.files('./nursing-home-population/')
for(k in unique(NH)) NHpop[[k]] <- read.csv(paste('./nursing-home-population/', k, sep=''))
for(k in 1:length(NHpop)) names(NHpop)[k] <- substr(names(NHpop)[k],1,nchar(names(NHpop)[k])-23)


# check which countries we have data for
missing <- setdiff(names(poplist), names(NHpop))


# remove nursing home populations from general population data & output
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR/data/population/non-nursing-home-population')
poplist <- poplist[!(names(poplist) %in% missing)]
nonNH <- poplist
for(ci in names(poplist)){
  
  # nursing home populations
  nhp <- NHpop[[ci]] 
  nhM <- round(c(rep(0,13), nhp$M[1:6], c(0.95*nhp$M[7], 0.05*nhp$M[7])))
  nhF <- round(c(rep(0,13), nhp$`F`[1:6], c(0.95*nhp$`F`[7], 0.05*nhp$`F`[7])))
  
  # remove
  nonNH[[ci]]$M <- poplist[[ci]]$M - nhM
  nonNH[[ci]]$`F` <- poplist[[ci]]$`F` - nhF
  
  # check for negative values
  nonNH[[ci]]$M[nonNH[[ci]]$M<0] <- 0
  nonNH[[ci]]$`F`[nonNH[[ci]]$`F`<0] <- 0
  
  # output
  write.csv(nonNH[[ci]], paste(ci,'.csv',sep=''))
}
