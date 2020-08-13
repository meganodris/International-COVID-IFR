library(rstan)

# Source functions
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR')
source('./code/FunctionsForIFR.R')


#----- Data Inputs -----#

# Age-specific death data
df <- read.csv('./data/deaths_age.csv')
df <- df[df$nat.reg=='nat', ]
countries <- sort(as.character(unique(df$country)))
df <- agg_deathsAp(df, countries, 80) # aggregate 80+ age groups
df <- agg_deathsu5(df, countries) # aggregate ages <4
df$age_max[df$age_max==17] <- 19 # where 17/18 is used as an age-band, switch to 19/20 (Greece & NYC only)
df$age_min[df$age_min==18] <- 20
dfU65 <- df[df$age_max<65, ]


# Adjusted death data 65+
df65p <- read.csv('./data/deaths_age_adjusted.csv')
df65p <- df65p[df65p$country=='England',]
df65p <- agg_deathsAp(df65p, 'England', 80)
df65p <- df65p[df65p$age_min>63, ]

# Demographic data
poplist <- list()
for(k in unique(countries)){
  poplist[[k]] <- read.csv(paste('./data/population/', k,'-2019.csv', sep=''))
} 

# Active Surveillance data
dpraw <- read.csv('./data/DiamondPrincess.csv')
dpd <- adjust_DPdata(dpraw)
cdg <- read.csv('./data/CharlesDeGaulle.csv')

# Serology data
sero <- read.csv('./data/serostudies.csv')
sero$region <- as.character(sero$region)
sero <- sero[sero$region %in% countries, ]

# Death time series from Johns Hopkins repo
setwd('C:/Users/Megan/Documents/GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series') 
deathsT <- read.csv('time_series_covid19_deaths_global.csv')
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR/data/')
deathsTR <- read.csv('deathsT_region.csv')
deathsT <- tidy_deathsT(deathsT, deathsTR, countries)

# Distribute deaths backwards for timing of seroconversion & infection
g <- delay_deaths(deathsT, countries)


# List of inputs for model
Inputs <- get_inputs(countries, poplist, poplist, dfU65, df65p, deathsT, sero, cdg, dpd, 13)
Inputs$deathsTinfec <- g$deathsTinfec
Inputs$deathsTsero <- g$deathsTsero
Inputs$relProbInfection <- rep(1,17) 



#----- Run model -----#
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR/code')
fit <- runStan(model='IFRinternational.stan', Inputs, N_iter=10000, N_chains=3, max_td=15, ad=0.95)


# Check convergence
chains <- rstan::extract(fit)
traceplot(fit, pars=names(fit)[1:10])
stan_dens(fit, pars=names(fit)[1:10])



#----- Model outputs -----#

# IFR estimates
ifrAge <- plot_IFR_age(chains, Inputs)
ifrPop <- plot_IFR_area(chains, Inputs, countries)
ifrAge$ifrAp
ifrPop$ifrCp

# Lambda estimates
lambdaFit <- serofit(chains, sero)
ggplot(lambdaFit, aes(seroprev,fit))+ geom_point(aes(col=region))+
  geom_line(aes(seroprev, seroprev))

# Other paramater estimates
pars <- extract_pars(chains)

# Fit to active surveillance data
dp <- fit_active(chains, Inputs)
dp

# Seroprevalence time series
seroT <- sero_time(chains, Inputs, countries, continent)
imm <- plot_immunity(chains, Inputs, countries, plotfit=T)
imm[[16]]

# Fit to age data
df <- agg_deathsAp(df, countries, 80)
plotFit <- fit_deaths_direct(chains, Inputs, df, countries)
plotFit$plots[[14]]



# output pars
write.csv(ifrAge$ifrA, 'IFRage.csv', row.names=F)
write.csv(ifrPop$ifrc, 'IFRPop.csv', row.names=F)
write.csv(pars, 'Params.csv', row.names=F)
write.csv(lambdaFit, 'Fitserology.csv', row.names=F)
write.csv(seroT$seroT, 'seroT.csv', row.names=F)
write.csv(seroT$infecT, 'infecT.csv', row.names=F)
write.csv(plotFit$ests, 'FitDeaths.csv', row.names=F)




