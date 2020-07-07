library(rstan)

# Source functions
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR')
source('./code/FunctionsForIFR.R')


#----- Data Inputs -----#

# Age-specific death data
df <- read.csv('./data/deaths_age.csv')
df <- df[df$nat.reg=='nat', ]
countries <- sort(unique(df$country))
df <- agg_deathsAp(df, countries, 80) # aggregate 80+ ages
df <- agg_deathsu5(df, countries) # aggregate ages <4
df$age_max[df$age_max==17] <- 19 
df$age_min[df$age_min==18] <- 20

# Demographic data
poplist <- list()
for(k in unique(countries)){
  poplist[[k]] <- read.csv(paste('./data/population/general-population/', k,'-2019.csv', sep=''))
} 

# Active Surveillance data
dpraw <- read.csv('./data/DiamondPrincess.csv')
dpd <- adjust_DPdata(dpraw)
cdg <- read.csv('./data/CharlesDeGaulle.csv')

# Serology data
sero <- read.csv('./data/seroprev_estimates.csv')
sero <- sero[sero$region %in% countries, ]

# Death time series
setwd('C:/Users/Megan/Documents/GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series')
deathsT <- read.csv('time_series_covid19_deaths_global.csv')
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR/data/')
deathsTR <- read.csv('deathsT_region.csv')
deathsT <- tidy_deathsT(deathsT, deathsTR, countries)

# Distribute deaths backwards
g <- delay_deaths(deathsT, countries)

# Hack the timings for now
df$asof <- as.Date(df$asof, format('%d/%m/%Y'))
colnames(deathsT)[ncol(deathsT)]
unique(df$country[df$asof>'2020-06-19'])
df$asof[df$asof>'2020-06-19'] <- '2020-06-19'
colnames(deathsT)[ncol(deathsT)-20]
sero$tmax <- as.Date(sero$tmax, format('%d/%m/%Y'))
unique(sero$region[sero$tmax>'2020-05-30'])
sero$tmax[sero$tmax>'2020-05-30'] <- '2020-05-30'
sero$tmin <- as.Date(sero$tmin, format('%d/%m/%Y'))
sero <- sero[!sero$tmin>max(dates), ]


# List of inputs for model
dfU65 <- df[df$age_max<65, ]
countries <- sort(as.character(countries))
continent <- vector()
for(i in 1:length(countries)) continent[i] <- paste(df$continent[df$country==countries[i]][1])
Inputs <- get_inputs(countries, poplist, poplist, dfU65, cdg, dpd)
Inputs <- c(Inputs, get_deathsT(deathsT, countries, dfU65))
sero <- sero[!is.na(sero$tmax), ]
Inputs <- c(Inputs, get_sero(sero, countries, Inputs$Ndays))
Inputs$deathsTinfec <- g$deathsTinfec
Inputs$deathsTsero <- g$deathsTsero
Inputs$relProbInfection <- rep(1,17) 


#----- Run model -----#
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR/code')
fit <- runStan(model='IFRinternational_main.stan', Inputs, N_iter=6000, N_chains=3, max_td=15, ad=0.7)


# Check convergence
chains <- rstan::extract(fit)
traceplot(fit, pars=names(fit)[1:10])
stan_dens(fit, pars=names(fit)[1:10])



#----- Model outputs -----#

# IFR estimates
ifrAge <- plot_IFR_age(chains, Inputs)
ifrPop <- plot_IFR_area(chains, Inputs, countries)

# Lambda estimates
sero <- sero[!is.na(sero$tmax), ]
lambdaFit <- serofit(chains, sero)

# Other paramater estimates
pars <- extract_pars(chains, linmod=T)

# Fit to active surveillance data
dp <- fit_active(chains, Inputs)

# Seroprevalence time series
seroT <- sero_time(chains, Inputs, countries, continent)
imm <- plot_immunity(chains, Inputs, countries, plotfit=T)

# Fit to age data
plotFit <- fit_deaths(chains, Inputs, df, countries, 65)


# Output path
opath <- 'C:/Users/Megan/OneDrive - University of Cambridge/COVID/IFRModel/international/NationalOnly'
setwd(opath)
folder <- 'Lin-u65-sero'
dir.create(folder, showWarnings=F)
setwd(paste(opath, folder, sep='/'))

# output pars
write.csv(ifrAge$ifrA, 'IFRage.csv', row.names=F)
write.csv(ifrPop$ifrc, 'IFRPop.csv', row.names=F)
write.csv(pars, 'Params.csv', row.names=F)
write.csv(lambdaFit, 'Fitserology.csv', row.names=F)
write.csv(seroT$seroT, 'seroT.csv', row.names=F)
write.csv(seroT$infecT, 'infecT.csv', row.names=F)
write.csv(plotFit$ests, 'FitDeaths.csv', row.names=F)

# plots
library(gridExtra)
png(filename='ImmunityT1.png', width=20, height=15, res=400, units='cm')
grid.arrange(grobs=imm[1:6], ncol=3, left='Proportion infected', bottom='date')
dev.off()
png(filename='ImmunityT2.png', width=20, height=15, res=400, units='cm')
grid.arrange(grobs=imm[7:12], ncol=3, left='Proportion infected', bottom='date')
dev.off()
png(filename='ImmunityT3.png', width=20, height=15, res=400, units='cm')
grid.arrange(grobs=imm[13:18], ncol=3, left='Proportion infected', bottom='date')
dev.off()
png(filename='ImmunityT4.png', width=20, height=15, res=400, units='cm')
grid.arrange(grobs=imm[19:24], ncol=3, left='Proportion infected', bottom='date')
dev.off()

png(filename='IFRlog.png', width=10, height=10, res=400, units='cm')
plot(xx)
dev.off()
png(filename='IFR.png', width=10, height=10, res=400, units='cm')
plot(ifrA$ifrAp)
dev.off()
png(filename='Active_fit.png', width=10, height=10, res=400, units='cm')
plot(dp)
dev.off()
png(filename='ProbInfec.png', width=18, height=10, res=400, units='cm')
plot(pinf$plotPI)
dev.off()
png(filename='ifrW.png', width=18, height=10, res=400, units='cm')
plot(ifrC$ifrCp)
dev.off()

library(gridExtra)
png(filename='FitDeaths1.png', width=20, height=15, res=400, units='cm')
grid.arrange(grobs=plotFit[1:6], ncol=3, left='deaths', bottom='age')
dev.off()
png(filename='FitDeaths2.png', width=20, height=15, res=400, units='cm')
grid.arrange(grobs=plotFit[7:12], ncol=3, left='deaths', bottom='age')
dev.off()
png(filename='FitDeaths3.png', width=20, height=15, res=400, units='cm')
grid.arrange(grobs=plotFit[13:18], ncol=3, left='deaths', bottom='age')
dev.off()
png(filename='FitDeaths4.png', width=20, height=15, res=400, units='cm')
grid.arrange(grobs=plotFit[19:24], ncol=3, left='deaths', bottom='age')
dev.off()
png(filename='FitDeaths5.png', width=20, height=15, res=400, units='cm')
grid.arrange(grobs=plotFit[25:30], ncol=3, left='deaths', bottom='age')
dev.off()
png(filename='FitDeaths6.png', width=20, height=15, res=400, units='cm')
grid.arrange(grobs=plotFit[30:36], ncol=3, left='deaths', bottom='age')
dev.off()
png(filename='FitDeaths7.png', width=20, height=15, res=400, units='cm')
grid.arrange(grobs=plotFit[37:42], ncol=3, left='deaths', bottom='age')
dev.off()
png(filename='FitDeaths8.png', width=20, height=15, res=400, units='cm')
grid.arrange(grobs=plotFit[43:44], ncol=3, left='deaths', bottom='age')
dev.off()
saveRDS(fit, 'Fit.RDS')
saveRDS(Inputs, 'Inputs.RDS')

