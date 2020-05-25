rm(list=ls())

# source functions
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR')
source('./code/FunctionsForIFR.R')

# death data
df <- read.csv('./data/deaths_age_adjusted.csv')
countries <- sort(unique(df$country))

# population data
poplist <- list()
for(k in unique(countries)){
  poplist[[k]] <- read.csv(paste('./data/population/non-nursing-home-population/', k,'.csv', sep=''))
} 

# Diamond Princess data
dpraw <- read.csv('./data/DiamondPrincess.csv')
dpd <- adjust_DPdata(dpraw)


# align death and pop data
df$age_max[df$age_max==17] <- 19 
df$age_min[df$age_min==18] <- 20
df <- agg_deaths80p(df, countries)


# data inputs for model
Inputs <- get_agesex(df, countries) # age/sex by country
Inputs$NArea <- length(countries) # N countries
Cpop <- compile_pop(poplist, countries)
Inputs <- c(Inputs, Cpop) # population data
ddata <- compile_deathsA(df, countries)
Inputs <- c(Inputs, ddata) # death data
ageInd <- index_ages(df, countries)
Inputs <- c(Inputs, ageInd) # age indices
Inputs$relProbInfection <- rep(1,17) # rel prob infection
Inputs$DP_pos_m <- dpd$pos_m # Diamond Princess data
Inputs$DP_pos_f <- dpd$pos_f
Inputs$DP_deathsTot <- 18


# run model
library(rstan)
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR/code')
fit <- runStan(model='IFRinternational.stan', Inputs, N_iter=8000, N_chains=3, max_td=15, ad=0.9)

# check chains
traceplot(fit)
chains <- rstan::extract(fit)

# plot outputs
plot_IFR(fit)
plot_Pinfection(chains, countries)
plotFit <- fit_deaths(chains, Inputs, df, countries)
plotFit[8]

fit_DP(fit, Inputs)

