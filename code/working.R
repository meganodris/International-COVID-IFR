rm(list=ls())

# source functions
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR')
source('./code/FunctionsForIFR.R')

# death data
df <- read.csv('./data/deaths_age.csv')
countries <- sort(unique(df$country))

# population data
poplist <- list()
for(k in unique(countries)){
  poplist[[k]] <- read.csv(paste('./data/population/', k,'-2019.csv', sep=''))
} 

# Diamond Princess data
dpraw <- read.csv('./data/DiamondPrincess.csv')
dpd <- adjust_DPdata(dpraw)

# select countries
countries <- c('Argentina','Colombia','Ecuador', 'New York City',
               'Mexico','Panama','Peru','Philippines','South Africa')
countries <- countries[!countries %in% c('England','Wales')]
df <- df[df$country %in% countries, ]
poplist <- poplist[names(poplist) %in% countries]

# align death and pop data
df$age_max[df$age_max==17] <- 19
df$age_min[df$age_min==18] <- 20
df <- agg_deaths80p(df, countries)
new <- list()
for(c in unique(df$country)){
  
  dc <- df[df$country==c, ]
  new[[c]] <- align_ages(data=dc, popdata=poplist[[c]])
}
countries

# data inputs for model
NAges <- vector()
for(c in 1:length(countries)) NAges[c] <- length(unique(df$age_min[df$country==countries[c]]))

gender <- vector()
for(c in 1:length(countries)){
  cc <- df[df$country==countries[c], ]
  if(cc$sex[1] %in% c('B')) gender[c] <- 1
  if(cc$sex[1] %in% c('M','F')) gender[c] <- 2
}

Inputs <- list()
Inputs$NArea <- length(countries)
Cpop <- compile_pop(poplist, countries)
Inputs <- c(Inputs, Cpop)
Inputs$NAges <- NAges
Inputs$gender <- gender

ddata <- compile_deathsA(df, countries)
Inputs <- c(Inputs, ddata)


ageInd <- index_ages(df, countries)
Inputs <- c(Inputs, ageInd)
Inputs$relProbInfection <- rep(1,17)
Inputs$DP_pos_m <- dpd$pos_m
Inputs$DP_pos_f <- dpd$pos_f
Inputs$DP_deathsTot <- 18

library(rstan)
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR/code')
fit <- runStan(model='IFRinternational.stan', Inputs, N_iter=6000, N_chains=3, max_td=15, ad=0.9)

traceplot(fit)
chains <- rstan::extract(fit)
mean(chains$relifrsex)
infec <- data.frame(co=countries, mean=NA, ciL=NA, ciU=NA)
for(c in 1:length(countries)){
  infec[c,2:4] <- quantile(chains$probInfec[,c], c(0.5,0.025,0.975))
}

pinf <- ggplot(infec, aes(reorder(co, mean),mean))+ geom_point(col='purple')+
  geom_linerange(aes(ymin=ciL, ymax=ciU),col='purple')+ theme_minimal()+
  theme(axis.text.x=element_text(angle=60, hjust=1))+ ylab('Prob Infection')+
  xlab('')#+ ylim(0,0.25)
getwd()
png(filename='EG_ProbInf.png', width=15, height=15, res=400, units='cm')
plot(pinf)
dev.off()


plot_IFR(fit)


getwd()
png(filename='IFRests.png', width=15, height=15, res=400, units='cm')
plot(ifrp)
dev.off()

library(gridExtra)
fp <- fit_deaths(chains, Inputs, df, countries)
fp[9]
countries
library(gridExtra)
getwd()
png(filename='EG_fits.png', width=18, height=27, res=400, units='cm')
grid.arrange(grobs=fp[1:8], ncol=2, left='deaths', bottom='age')
dev.off()
fp[8]
fit_DP(fit, Inputs)

