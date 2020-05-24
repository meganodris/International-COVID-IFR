rm(list=ls())

setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR')
source('./code/FunctionsForIFR.R')

# death data
df <- read.csv('./data/deaths_age.csv')
cds <- unique(df$country)

# population data
files <- list.files('./data/population')
countries <- vector()
for(c in 1:length(files)) countries[c] <- substr(files[c],1,nchar(files[c])-9)
setdiff(cds, countries)
poplist <- list()
for(k in 1:length(files)) poplist[[countries[k]]] <- read.csv(paste('./data/population',files[k],sep='/'))

# Diamond Princess data
dpraw <- read.csv('./data/DiamondPrincess.csv',sep=';')
dpd <- tidy_DPdata(dpraw)


# align death and pop data
poplist <- poplist[!(names(poplist) %in% c('England & Wales', 'Geneva'))]
countries <- countries[!(countries %in% c('England & Wales', 'Geneva'))]
df$age_max[df$age_max==17] <- 19
df$age_min[df$age_min==18] <- 20
new <- list()
for(c in names(poplist)){
  
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
Cpop <- compile_pop(poplist)
Inputs <- c(Inputs, Cpop)
Inputs$NAges <- NAges
Inputs$gender <- gender
Inputs$NGroups <- 8
deaths_m <- matrix(-999, nrow=21, ncol=length(countries))
deaths_f <- matrix(-999, nrow=21, ncol=length(countries))
deaths_b <- matrix(-999, nrow=21, ncol=length(countries))

for(c in 1:length(countries)){
  
  dfc <- df[df$country==countries[c], ]
  if(dfc$sex[1]=='B'){
    deaths_b[1:nrow(dfc),c] <- dfc$deaths
  }else{
    deaths_m[1:(nrow(dfc)/2),c] <- dfc$deaths[dfc$sex=='M']
    deaths_f[1:(nrow(dfc)/2),c] <- dfc$deaths[dfc$sex=='F']
  }
}

Inputs$deaths_b <- deaths_b
Inputs$deaths_m <- deaths_m
Inputs$deaths_f <- deaths_f

ageInd <- index_ages(df, countries)
Inputs <- c(Inputs, ageInd)
Inputs$relProbInfection <- rep(1,21)
Inputs$DP_pos_m <- dpd$pos_m
Inputs$DP_pos_f <- dpd$pos_f
Inputs$DP_deathsTot <- 18

library(rstan)
setwd('C:/Users/Megan/Documents/GitHub/International-COVID-IFR/code')
fit <- runStan(model='IFRinternational.stan', Inputs, N_iter=6000, N_chains=3, max_td=15, ad=0.9)

traceplot(fit)
chains <- rstan::extract(fit)
mean(chains$relprobDeathSex)
infec <- data.frame(co=countries, mean=NA, ciL=NA, ciU=NA)
for(c in 1:length(countries)){
  infec[c,2:4] <- quantile(chains$probInfec[,c], c(0.5,0.025,0.975))
}

ggplot(infec, aes(reorder(co, mean),mean))+ geom_point()+
  geom_linerange(aes(ymin=ciL, ymax=ciU))+
  theme(axis.text.x=element_text(angle=60, hjust=1))
ifr_f <- data.frame(mean=NA, ciL=NA, ciU=NA)
ifr_m <- data.frame(mean=NA, ciL=NA, ciU=NA)
for(a in 1:21){
  ifr_f[a,] <- quantile(chains$ifr_f[,a], c(0.5,0.025,0.975))
  ifr_m[a,] <- quantile(chains$ifr_m[,a], c(0.5,0.025,0.975))
}
ifr <- rbind(ifr_m, ifr_f)
ifr$sex <- c(rep('M',21),rep('F',21))
ifr$age <- rep(seq(1:21),2)

ggplot(ifr, aes(age, mean, col=sex))+ geom_point(position=position_dodge(width=0.4))+
  geom_linerange(aes(ymin=ciL, ymax=ciU, col=sex),position=position_dodge(width=0.4))





# data plots
plots <- list()
df$age_mid <- df$age_min+((df$age_max-df$age_min)+1)/2
df$ageG <- paste(df$age_min, df$age_max, sep='-')
df$ageG[df$age_max==110] <- paste(df$age_min[df$age_max==110], '+', sep='')
for(c in 1:length(countries)){
  
  dfc <- df[df$country==countries[c], ]
  
  if(dfc$sex[1]=='B'){
    plots[[c]] <- ggplot(dfc, aes(reorder(ageG,age_min), deaths))+ geom_col(fill='seagreen2')+
      theme_minimal()+ ylab('')+ labs(subtitle=paste(dfc$country[1], ', N=', sum(dfc$deaths), sep=''))+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_text(angle=60, hjust=1))
  }else{
    plots[[c]] <- ggplot(dfc, aes(reorder(ageG,age_min), deaths, fill=sex))+ 
      geom_col(position=position_dodge())+ theme_minimal()+ ylab('')+ 
      labs(subtitle=paste(dfc$country[1], ', N=', sum(dfc$deaths), sep=''))+
      theme(axis.title.x=element_blank(), legend.position='none',
            axis.text.x=element_text(angle=60, hjust=1))+
      scale_fill_manual(values=c('indianred1','royalblue1'), labels=c('female','male'))
  
    }
}
countries
plots[[16]]
library(gridExtra)
setwd('C:/Users/Megan/OneDrive - University of Cambridge/COVID/IFRModel/international/DataPlots')
png(filename='AgeDeaths_1.png', width=18, height=27, res=400, units='cm')
grid.arrange(grobs=plots[1:10], ncol=2, left='deaths', bottom='age')
dev.off()
png(filename='AgeDeaths_2.png', width=18, height=27, res=400, units='cm')
grid.arrange(grobs=plots[11:20], ncol=2, left='deaths', bottom='age')
dev.off()
png(filename='AgeDeaths_3.png', width=18, height=27, res=400, units='cm')
grid.arrange(grobs=plots[21:30], ncol=2, left='deaths', bottom='age')
dev.off()
png(filename='AgeDeaths_4.png', width=18, height=13, res=400, units='cm')
grid.arrange(grobs=plots[31:35], ncol=2, left='deaths', bottom='age')
dev.off()


# age-specific RR
yay <- do.call('rbind', new)
yay$agemid <- yay$age_min+((yay$age_max-yay$age_min)+1)/2
yay$ageG <- paste(yay$age_min, yay$age_max, sep='-')
yay$ageG[yay$age_max==110] <- paste(yay$age_min[yay$age_max==110], '+', sep='')
library(ggplot2)
RRplots <- list()
for(c in 1:length(countries)){
  
  dfc <- yay[yay$country==countries[c], ]
  dfc$probD <- dfc$deaths/dfc$pop
  dfc$ref <- 0
  dfc$RR <- 0
  if(dfc$sex[1]=='B'){
    if(any(dfc$age_min==60)){
      dfc$ref[dfc$age_min==60] <- 1
    }else{
      dfc$ref[dfc$age_min==65] <- 1
    } 
    dfc$RR <- dfc$probD/dfc$probD[dfc$ref==1]
    RRplots[[c]] <- ggplot(dfc, aes(reorder(ageG,age_min),RR))+ geom_point(col='seagreen2')+
      theme_minimal()+ geom_hline(aes(yintercept=1), linetype='dashed')+
      geom_line(aes(reorder(ageG,age_min),RR,group=1),col='seagreen2')+
      theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=60,hjust=1),
            legend.position='none')+ labs(subtitle=paste(dfc$country[1]))+ ylab('')
  }else{
    if(any(dfc$age_min==60)){
      dfc$ref[dfc$age_min==60] <- 1
    }else{
      dfc$ref[dfc$age_min==65] <- 1
    } 
    dfc$RR[dfc$sex=='M'] <- dfc$probD[dfc$sex=='M']/dfc$probD[which(dfc$sex=='M' & dfc$ref==1)]
    dfc$RR[dfc$sex=='F'] <- dfc$probD[dfc$sex=='F']/dfc$probD[which(dfc$sex=='F' & dfc$ref==1)]
    RRplots[[c]] <- ggplot(dfc, aes(reorder(ageG,age_min),RR,col=sex,group=sex))+ geom_point(position=position_dodge(width=0.3))+
      theme_minimal()+ geom_hline(aes(yintercept=1), linetype='dashed')+
      geom_line(aes(reorder(ageG,age_min),RR,col=sex),position=position_dodge(width=0.3))+
      theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=60,hjust=1),
            legend.position='none')+ labs(subtitle=paste(dfc$country[1]))+
      scale_colour_manual(values=c('indianred1','royalblue1'))+ ylab('')
  }
}

RRplots[[10]]
library(gridExtra)
setwd('C:/Users/Megan/OneDrive - University of Cambridge/COVID/IFRModel/international/DataPlots')
png(filename='RRDeaths_1.png', width=18, height=27, res=400, units='cm')
grid.arrange(grobs=RRplots[1:10], ncol=2, left='RR', bottom='age')
dev.off()
png(filename='RRDeaths_2.png', width=18, height=27, res=400, units='cm')
grid.arrange(grobs=RRplots[11:20], ncol=2, left='RR', bottom='age')
dev.off()
png(filename='RRDeaths_3.png', width=18, height=27, res=400, units='cm')
grid.arrange(grobs=RRplots[21:30], ncol=2, left='RR', bottom='age')
dev.off()
png(filename='RRDeaths_4.png', width=18, height=13, res=400, units='cm')
grid.arrange(grobs=RRplots[31:35], ncol=2, left='RR', bottom='age')
dev.off()