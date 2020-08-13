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


# Adjusted data 65+
df2 <- read.csv('./data/deaths_age_adjusted2.csv')
df2 <- df2[df2$country=='England',]
df2 <- agg_deathsAp(df2, 'England', 80)
df2 <- df2[df2$age_min>63, ]

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
sero$region <- as.character(sero$region)
sero <- sero[sero$region %in% countries, ]
sero <- sero[!sero$study %in% c('France-2'), ]

# Death time series
setwd('C:/Users/Megan/Documents/GitHub/COVID-19/csse_covid_19_data/csse_covid_19_time_series')
deathsT <- read.csv('time_series_covid19_deaths_global.csv')
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR/data/')
deathsTR <- read.csv('deathsT_region.csv')
deathsT <- tidy_deathsT(deathsT, deathsTR, countries)

# Distribute deaths backwards
g <- delay_deaths(deathsT, countries)

# Hack the timings for now
#df$asof <- as.Date(df$asof, format('%d/%m/%Y'))
#colnames(deathsT)[ncol(deathsT)]
#unique(df$country[df$asof>'2020-07-10'])
#df$asof[df$asof>'2020-06-19'] <- '2020-06-19'
#colnames(deathsT)[ncol(deathsT)-20]
#sero$tmax <- as.Date(sero$tmax, format('%d/%m/%Y'))
#unique(sero$region[sero$tmax>'2020-05-30'])
#sero$tmax[sero$tmax>'2020-05-30'] <- '2020-05-30'
#sero$tmin <- as.Date(sero$tmin, format('%d/%m/%Y'))
#sero <- sero[!sero$tmin>'2020-05-30', ]
#sero <- sero[!is.na(sero$study),]
#for(i in unique(sero$study)){
 # nr <- nrow(sero[sero$study==i, ])
  #for(n in 1:nr){
   # sero$n[sero$study==i][n] <- floor(5000/nr)
  #}
  #sero$n_pos[sero$study==i] <- round(sero$seroprev[sero$study==i]*sero$n[sero$study==i])
#}



# List of inputs for model
dfU65 <- df[df$age_max<65, ]
countries <- sort(as.character(countries))
continent <- vector()
for(i in 1:length(countries)) continent[i] <- paste(df$continent[df$country==countries[i]][1])
Inputs <- get_inputs(countries, poplist, poplist, dfU65, cdg, dpd, 13)
Inputs <- c(Inputs, get_deathsT(deathsT, countries, dfU65))
sero <- sero[!is.na(sero$tmax), ]
Inputs <- c(Inputs, get_sero(sero, countries, Inputs$Ndays))
Inputs$deathsTinfec <- g$deathsTinfec
Inputs$deathsTsero <- g$deathsTsero
Inputs$relProbInfection <- rep(1,17) 
deaths65p_m <- df2$deaths[df2$sex=='M']
deaths65p_f <- df2$deaths[df2$sex=='F']
Inputs$indArea65p <- which(countries=='England')
Inputs$deaths65p_m <- deaths65p_m
Inputs$deaths65p_f <- deaths65p_f
Inputs$relProbInfection[14:17] <- 0.7

#----- Run model -----#
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR/code')
fit <- runStan(model='IFRinternational_mainnew.stan', Inputs, N_iter=10000, N_chains=3, max_td=15, ad=0.95)


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
ggplot(ifrAge$ifrA[ifrAge$ifrA$sex %in% c('M','F'), ], aes(age, log(mean)))+
  geom_point()

# Lambda estimates
sero <- sero[!sero$study=='Denmark-3', ]
lambdaFit <- serofit(chains, sero)
ggplot(lambdaFit, aes(seroprev,fit,col=region))+ geom_point()+
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


ifro <- data.frame(country=co2, sex=c(rep('M',length(co2)),rep('F',length(co2))),
                   mean=NA, ciL=NA, ciU=NA)
for(i in 1:length(co2)){
  ifro[i,3:5] <- quantile(chains$ifrm65p[,i], c(0.5,0.025,0.975))
  ifro[i+length(co2),3:5] <- quantile(chains$ifrf65p[,i], c(0.5,0.025,0.975))
}         
enm <- quantile(chains$ifr_m[,14], c(0.5,0.025,0.975))
enf <- quantile(chains$ifr_f[,14], c(0.5,0.025,0.975))

ifro$country <- as.character(ifro$country)
ifro$country[ifro$country=='United States of America'] <- 'USA'
ifrOld <- ggplot(ifro, aes(reorder(country,mean),mean))+ 
  geom_hline(aes(yintercept=enm[1]), col='blue')+
  geom_hline(aes(yintercept=enf[1]), col='red')+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=enm[2],ymax=enm[3]), fill='lightblue', alpha=0.1)+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=enf[2],ymax=enf[3]), fill='lightcoral', alpha=0.1)+
  geom_point(position=position_dodge(width=0.3), aes(col=sex))+
  geom_linerange(aes(ymin=ciL, ymax=ciU,col=sex),position=position_dodge(width=0.3))+
  theme_minimal()+ theme(axis.title.x=element_blank(), 
                         axis.text.x=element_text(angle=60,hjust=1))+
  ylab('IFR 65+')
ifrOld
# Output path
opath <- 'C:/Users/Megan/OneDrive - University of Cambridge/COVID/IFRModel/international/NationalOnly'
setwd(opath)
folder <- 'Main-sero-new-07'
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
saveRDS(Inputs, 'Inputs.RDS')
saveRDS(chains, 'Chains.RDS')
png(filename='IFR65p.png', width=20, height=15, res=400, units='cm')
plot(ifrOld)
dev.off()

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

