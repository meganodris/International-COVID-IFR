rm(list=ls())
library(rstan)

# Source functions
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR')
source('./code/FunctionsForIFR.R')


#### Read in data ####

# Age-specific death data 
df <- read.csv('./data/deaths_age.csv')
countries <- sort(unique(df$country))
countries <- countries[!countries %in% c('Catalonia','Rio de Janeiro','Canada')]
df <- df[df$country %in% countries,]
df <- agg_deathsAp(df, countries, 80) # aggregate ages 80+


# Population data
poplist <- list()
for(k in unique(countries)){
  poplist[[k]] <- read.csv(paste('./data/population/general-population/', k,'-2019.csv', sep=''))
} 


# Diamond Princess & Charles de Gaulle data
dpraw <- read.csv('./data/DiamondPrincess.csv')
dpd <- adjust_DPdata(dpraw)
cdg <- read.csv('./data/CharlesDeGaulle.csv')


# Serology data
sero <- read.csv('./data/seroprev_estimates.csv')
sero2 <- sero[sero$region %in% countries,]


# align death and pop data
dff <- df[df$age_max<65,]
dff$age_max[dff$age_max==17] <- 19 
dff$age_min[dff$age_min==18] <- 20




# death time series
setwd('C:/Users/Megan/Documents/GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series')
deathsT <- read.csv('time_series_covid19_deaths_global.csv')
deathsT$Country.Region <- as.character(deathsT$Country.Region)
deathsT$Country.Region[deathsT$Country.Region=='Korea, South'] <- 'Republic of Korea'
deathsT$Country.Region[deathsT$Country.Region=='Czechia'] <- 'Czech Republic'
deathsT$Country.Region[deathsT$Country.Region=='US'] <- 'United States of America'
deathsT <- deathsT[deathsT$Country.Region %in% countries, ]
deathsT <- deathsT[deathsT$Province.State=='', ]
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR/data/')
deathsTR <- read.csv('deathsT_region.csv')
deathsTR <- deathsTR[deathsTR[,1] %in% countries, ]
deathsTR$ï..region <- as.character(deathsTR$ï..region)
deathsTR <- deathsTR[,1:132]
for(c in 1:nrow(deathsT)){
  tt <- deathsT[c,5:ncol(deathsT)]
  deathsTR[nrow(deathsTR)+1, ] <- c(paste(deathsT$Country.Region[c]), tt[1:ncol(deathsTR)-1])
}

for(j in 1:nrow(deathsTR)){
  for(k in 7:(ncol(deathsTR)-5)) deathsTR[j,k] <- round(mean(as.numeric(deathsTR[j,(k-5):(k+5)]), na.rm=T))
}   

# data inputs for model
countries <- sort(as.character(countries))
Inputs <- get_inputs(countries, poplist, poplist, dff, cdg, dpd)
Inputs <- c(Inputs, get_deathsT(deathsTR, countries, dff, infectodeath=25, infectosero=15))
Inputs <- c(Inputs, get_sero(sero, countries, Inputs$Ndays, Inputs$serotodeath))
Inputs$relProbInfection <- rep(1,17) # rel prob infection


# Run model
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR/code')
fit <- runStan(model='IFRinternationalLinear.stan', Inputs, N_iter=8000, N_chains=3, max_td=15, ad=0.7)


# Check convergence
chains <- rstan::extract(fit)
traceplot(fit, pars=names(fit)[1:10])
stan_dens(fit, pars=names(fit)[1:10])


#=== plot outputs
imm <- plot_immunity(chains, Inputs, countries, plotfit=T)
imm[[24]]

# IFR estimates by age
ifrA <- plot_IFR_age(chains, Inputs)
ifrA$ifrAp
x <- ifrA$ifrA
x$age_mid <- c(2,7,12,17,22,27,32,37,42,47,52,57,62,67,72,77,85)
xx <- ggplot(x, aes(reorder(ageG,age), log(mean), col=sex))+ geom_point()+
  theme_minimal()+ theme(axis.text.x=element_text(angle=60, hjust=1))+
  xlab('')+ ylab('log IFR')
xx
# IFR estimates by country
ifrC <- plot_IFR_area(chains, Inputs, countries)
ifrC$ifrCp
ifrC$ifrc

# fit to active surveillance campaigns
dp <- fit_active(chains, Inputs)
dp

# Prob infection by country
pinf <- plot_Pinfection(chains, countries)
pinf$plotPI
pp <- pinf$estsPI

# fits
df <- agg_deathsAp(df, countries, 80)
plotFit <- fit_deaths(chains, Inputs, df, countries, 65)
plotFit[27]
library(cowplot)
plot_grid(plotlist=plotFit[1:4])

# outputs
opath <- 'C:/Users/Megan/OneDrive - University of Cambridge/COVID/IFRModel/international/test_fits'
setwd(opath)
folder <- 'age-constant-lin-u65'
dir.create(folder, showWarnings=F)
setwd(paste(opath, folder, sep='/'))

# ifrs
write.csv(ifrA$ifrA, 'IFRage.csv')
write.csv(ifrC$ifrc, 'IFRW.csv')
write.csv(pinf$estsPI, 'SeroEsts.csv')
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

