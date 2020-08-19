#===== Functions for International IFR model analysis =====#
library(ggplot2)
library(epitools)


#----- 1. Functions for compiling model inputs -----#


# Function to compile demographic data for model
compile_pop <- function(poplist, countries){
  
  # matrices for model
  pop_m <- matrix(NA, nrow=17, ncol=length(countries)) 
  pop_f <- matrix(NA, nrow=17, ncol=length(countries)) 
  pop_b <- matrix(NA, nrow=17, ncol=length(countries)) 
  
  # compile for model
  for(c in 1:length(countries)){
    ci <- countries[c]
    pop_m[ ,c] <- poplist[[ci]]$M[1:17]
    pop_f[ ,c] <- poplist[[ci]]$`F`[1:17]
    pop_b[ ,c] <- poplist[[ci]]$M[1:17] + poplist[[ci]]$`F`[1:17]
    
    # add ages 90+
    pop_m[17,c] <- sum(poplist[[ci]]$M[17:21])
    pop_f[17,c] <- sum(poplist[[ci]]$`F`[17:21])
    pop_b[17,c] <- sum(poplist[[ci]]$M[17:21]) + sum(poplist[[c]]$`F`[17:21])
  }
  
  
  # return list of matrices
  return(list(pop_m=pop_m, pop_f=pop_f, pop_b=pop_b))
}


# Function to compile death data for model 
compile_deathsA <- function(deathsA, countries, nAges){
  
  # matrices for model
  deaths_m <- matrix(-999, nrow=nAges, ncol=length(countries))
  deaths_f <- matrix(-999, nrow=nAges, ncol=length(countries))
  deaths_b <- matrix(-999, nrow=nAges, ncol=length(countries))
  
  # compile for model
  for(c in 1:length(countries)){
    
    dfc <- deathsA[deathsA$country==countries[c], ]
    dfc <- dfc[order(dfc$age_min), ]
    
    if(dfc$sex[1]=='B'){
      deaths_b[1:nrow(dfc),c] <- dfc$deaths
    }else{
      deaths_m[1:(nrow(dfc)/2),c] <- dfc$deaths[dfc$sex=='M']
      deaths_f[1:(nrow(dfc)/2),c] <- dfc$deaths[dfc$sex=='F']
    }
  }
  
  return(list(deaths_m=deaths_m, deaths_f=deaths_f, deaths_b=deaths_b))
}


# Function to aggregate deaths A+ (e.g. 80+)
agg_deathsAp <- function(deathsA, countries, A){
  
  # list for aggregated/tidied dataframes
  agg <- list()
  
  # for each country
  for(c in 1:length(countries)){
    
    dfc <- deathsA[deathsA$country==countries[c], ]
    
    # if multiple age groups 80+, aggregate
    if(any(dfc$age_min>A)){ 
      new <- dfc[dfc$age_max<A, ]
      for(s in unique(new$sex)){
        d80p <- sum(dfc$deaths[dfc$age_max>A & dfc$sex==s])
        minA <- min(dfc$age_min[dfc$age_max>A])
        n80p <- data.frame(country=countries[c], sex=s, age_min=minA, age_max=110,
                           deaths=as.numeric(d80p), asof=new$asof[1], nat.reg=new$nat.reg[1], continent=new$continent[1])
        new <- rbind(new,n80p)
      }
      agg[[c]] <- new[order(new$age_min), ]
    }else{
      agg[[c]] <- dfc[order(dfc$age_min), ]
    }
  }

  # return updated dataframe
  deathsAge <- do.call('rbind', agg)
  return(deathsAge)
}

# Function to aggregate deaths under 5 years
agg_deathsu5 <- function(deathsA, countries){
  
  # list for aggregated/tidied dataframes
  agg <- list()
  
  # for each country
  for(c in 1:length(countries)){
    
    dfc <- deathsA[deathsA$country==countries[c], ]
    
    # if multiple age groups 80+, aggregate
    if(any(dfc$age_max<4)){ 
      new <- dfc[dfc$age_min>4, ]
      for(s in unique(new$sex)){
        u5_s <- sum(dfc$deaths[dfc$sex==s & dfc$age_min<=4]) 
        maxA <- max(dfc$age_max[dfc$age_min<=4])
        nu5s <- data.frame(country=countries[c], sex=s, age_min=0, age_max=maxA,
                           deaths=as.numeric(u5_s), asof=new$asof[1], nat.reg=new$nat.reg[1], continent=new$continent[1])
        new <- rbind(new,nu5s)
      }
      agg[[c]] <- new[order(new$age_min), ]
    }else{
      agg[[c]] <- dfc[order(dfc$age_min), ]
    }
  }

  # return updated dataframe
  deathsAge <- do.call('rbind', agg)
  return(deathsAge)
}


# Align population and death data age groupings
align_ages <- function(data, popdata){
  
  # min & max bounds of population age groups
  popdata$age_min <- seq(0,100,5)
  popdata$age_max <- c(seq(4,99,5),110)
  
  # ensure age groups are alignable
  if(length(setdiff(c(data$age_min, data$age_max), c(popdata$age_min, popdata$age_max)))>0){
    stop('Non-alignable age groups in data')
  }
  
  # align
  data$pop <- NA
  for(a in 1:nrow(data)){
    minA <- which(popdata$age_min==data$age_min[a])
    maxA <- which(popdata$age_max==data$age_max[a])
    if(data$sex[a]=='M') data$pop[a] <- sum(popdata$M[minA:maxA])
    if(data$sex[a]=='F') data$pop[a] <- sum(popdata$`F`[minA:maxA])
    if(data$sex[a]=='B') data$pop[a] <- sum(popdata$M[minA:maxA]+popdata$`F`[minA:maxA])
  }
  
  # return 
  return(data)
}


# Index age groups for model
index_ages <- function(data, countries, NAges){
  
  # min & max bounds of age groupings
  Amin <- seq(0,80,5)
  Amax <- c(seq(4,79,5),110)
  
  # matrices for storage
  ageG_min <- matrix(-999, nrow=NAges, ncol=length(countries))
  ageG_max <- matrix(-999, nrow=NAges, ncol=length(countries))
  
  # get indices
  for(c in 1:length(countries)){
    
    dfc <- data[data$country==countries[c], ]
    dfc <- dfc[order(dfc$age_min), ]
    
    if(dfc$sex[1]=='B'){
      for(a in 1:nrow(dfc)){
        ageG_min[a,c] <- which(Amin==dfc$age_min[a])
        ageG_max[a,c] <- which(Amax==dfc$age_max[a])
      }
    }else{
      dfc <- dfc[dfc$sex=='M', ]
      for(a in 1:nrow(dfc)){
        ageG_min[a,c] <- which(Amin==dfc$age_min[a])
        ageG_max[a,c] <- which(Amax==dfc$age_max[a])
      }
    }
  }
  
  # return list of matrices 
  return(list(ageG_min=ageG_min, ageG_max=ageG_max))
}


# Tidy Diamond Princess data
adjust_DPdata <- function(DP_raw){
  
  # age-sex breakdown for 634/712 cases
  pos_m <- DP_raw$pos_m
  pos_f <- DP_raw$pos_f
  
  AdditionalPositives <- 712 - sum(pos_m+pos_f) # cases with no age/sex
  wts <- cbind(pos_m, pos_f)/(sum(pos_m+pos_f)) # assume same age/sex dist of new positives
  DP_raw$pos_m <- round(pos_m + (AdditionalPositives*wts)[,1])
  DP_raw$pos_f <- round(pos_f + (AdditionalPositives*wts)[,2])
  
  return(DP_raw) # return adjusted data
}

# Compile death time series data
tidy_deathsT <- function(deathsT, deathsTR, countries){
  
  # align names 
  deathsT$Country.Region <- as.character(deathsT$Country.Region)
  deathsT$Country.Region[deathsT$Country.Region=='Korea, South'] <- 'Republic of Korea'
  deathsT$Country.Region[deathsT$Country.Region=='Czechia'] <- 'Czech Republic'
  deathsT$Country.Region[deathsT$Country.Region=='US'] <- 'United States of America'
  
  # sum Canadian & Chinese regions
  ca <- colSums(deathsT[deathsT$Country.Region=='Canada', 5:ncol(deathsT)])
  ci <- colSums(deathsT[deathsT$Country.Region=='China', 5:ncol(deathsT)])
  deathsT[nrow(deathsT)+1, ] <- c('','Canada',NA,NA,ca)
  deathsT[nrow(deathsT)+1, ] <- c('','China',NA,NA,ci)
  
  # subset for inlcuded countries 
  deathsT <- deathsT[deathsT$Country.Region %in% countries, ]
  deathsT <- deathsT[deathsT$Province.State=='', ]
  
  # regional time series
  deathsTR <- deathsTR[deathsTR$region %in% countries, ]
  deathsTR$region <- as.character(deathsTR$region)
  
  # merge datasets
  latestdate <- max(which(!is.na(colSums(deathsTR[,3:ncol(deathsTR)]))))+2
  deathsTR <- deathsTR[,2:latestdate]
  for(c in 1:nrow(deathsT)){
    tt <- deathsT[c,5:ncol(deathsT)]
    deathsTR[nrow(deathsTR)+1, ] <- c(paste(deathsT$Country.Region[c]), tt[1:ncol(deathsTR)-1])
  }
  
  # return tidy dataset
  return(deathsT=deathsTR)
}


# Function to distribute deaths back to date of onset & seroconversion
delay_deaths <- function(Cdeaths, countries){
  
  # daily deaths
  Ndays <- ncol(Cdeaths)-1
  ddeathsT <- matrix(0, nrow=Ndays, ncol=length(countries))
  for(c in 1:length(countries)){
    dt <- as.numeric(as.vector(t(Cdeaths[Cdeaths$region==countries[c],2:(Ndays+1)])))
    for(i in 2:length(dt)){
      ddeathsT[i,c] <- dt[i] - dt[i-1]
    }
  }
  
  # matrices for storage
  deathsTinfec <- matrix(0, nrow=Ndays, ncol=length(countries))
  deathsTsero <- matrix(0, nrow=Ndays, ncol=length(countries))
  deathsTonset <- matrix(0, nrow=Ndays, ncol=length(countries))
  
  
  # infection to death
  id <- rgamma(1000,shape=(22.9/12.4)^2, rate=22.9/(12.4)^2)
  for(c in 1:length(countries)){
    for(i in 1:nrow(deathsTinfec)){
      if(ddeathsT[i,c]>0){
        for(k in 1:ddeathsT[i,c]){
          di <- floor(id[runif(1,1,length(id))])
          deathsTinfec[(i-di),c] <- deathsTinfec[(i-di),c]+1
        }
      }
    }
    deathsTinfec[,c] <- cumsum(deathsTinfec[,c])
  }
  
  # onset to seroconversion
  od <- rgamma(1000,shape=(17.8/8)^2, rate=17.8/(8)^2)
  os <- rgamma(1000,shape=(10/4)^2, rate=10/(4)^2)
  for(c in 1:length(countries)){
    for(i in 1:nrow(deathsTsero)){
      if(ddeathsT[i,c]>0){
        for(k in 1:ddeathsT[i,c]){
          di <- floor(od[runif(1,1,length(od))])
          di2 <- floor(os[runif(1,1,length(os))])
          deathsTonset[(i-di),c] <- deathsTonset[(i-di),c]+1
          if((i-di+di2)<nrow(deathsTsero)) deathsTsero[(i-di+di2),c] <- deathsTsero[(i-di+di2),c]+1
        }
      }
    }
    deathsTonset[,c] <- cumsum(deathsTonset[,c])
    deathsTsero[,c] <- cumsum(deathsTsero[,c])
  }
  
  # return matrices
  return(list(deathsTinfec=deathsTinfec[1:(Ndays-20), ], 
              deathsTonset=deathsTonset[1:(Ndays-20), ], 
              deathsTsero=deathsTsero[1:(Ndays-20), ]))
}


# Function to return N age groups & sex by country
get_agesex <- function(data, countries){
  
  # vectors to index Nages & sex by country
  NAges <- vector()
  gender <- vector()
  
  for(c in 1:length(countries)) NAges[c] <- length(unique(data$age_min[data$country==countries[c]]))
  for(c in 1:length(countries)){
    cc <- data[data$country==countries[c], ]
    if(cc$sex[1] %in% c('B')) gender[c] <- 1
    if(cc$sex[1] %in% c('M','F')) gender[c] <- 2
  }
  
  return(list(NAges=NAges, gender=gender))
}

# List of data inputs for model
get_inputs <- function(countries, poplist, dataA, df65p, deathsT, sero, cdg, dpd, NAges){
  
  Inputs <- list()
  
  # Death and population data
  Inputs <- get_agesex(dataA, countries) # age/sex by country
  Inputs$NArea <- length(countries) # N countries
  Inputs <- c(Inputs, compile_pop(poplist, countries)) # population data
  Inputs <- c(Inputs, compile_deathsA(dataA, countries, NAges)) # death data
  Inputs <- c(Inputs, index_ages(dataA, countries, NAges)) # age indices
  Inputs <- c(Inputs, get_deathsT(deathsT, countries, dataA)) # death time series
  Inputs <- c(Inputs, get_sero(sero, countries, Inputs$Ndays)) # serology data
  Inputs$indexArea65p <- which(countries=='England') # adjusted death data 65+
  Inputs$deaths65p_m <- df65p$deaths[df65p$sex=='M']
  Inputs$deaths65p_f <- df65p$deaths[df65p$sex=='F']
  
  # Active surveillance data
  Inputs$CDG_pos_m <- cdg$pos_m
  Inputs$CDG_pos_f <- cdg$pos_f
  Inputs$CDGamin <- c(5,6,8,10)
  Inputs$CDGamax <- c(5,7,9,12)
  Inputs$agemid <- c(2,7,12,17,22,27,32,37,42,47,52,57,62,67,72,77,85)
  Inputs$DP_pos_m <- dpd$pos_m 
  Inputs$DP_pos_f <- dpd$pos_f
  Inputs$DPamin <- c(1,5,7,9,11,13,15,17)
  Inputs$DPamax <- c(4,6,8,10,12,14,16,17)

  return(Inputs)
}


# Function to get model inputs for time series of death data
get_deathsT <- function(deathsT, countries, deathsA){
  
  # N days since 22/01/2020
  Ndays <- ncol(deathsT)-1
  
  # matrix of deaths over time
  dt <- matrix(NA, ncol=length(countries), nrow=Ndays)
  for(c in 1:length(countries)) dt[,c] <- as.numeric(deathsT[deathsT$region==countries[c],2:(Ndays+1)])
  
  # index for date of age-specific deaths
  TdeathsA <- vector()
  dates <- seq.Date(as.Date('2020-01-22'),by=1, length.out=Ndays)
  deathsA$asof <- as.Date(deathsA$asof, format('%d/%m/%Y'))
  for(c in 1:length(countries)) TdeathsA[c] <- which(dates==deathsA$asof[deathsA$country==countries[c]][1])
  
  # return list of inputs for model
  return(list(Ndays=Ndays-20, deathsT=dt, TdeathsA=TdeathsA))
}


# Function to get model inputs for serology data
get_sero <- function(sero, countries, Ndays){
  
  # countries for model
  sero$region <- as.character(sero$region)
  sero <- sero[sero$region %in% countries, ]
  
  # format dates
  sero$tmin <- as.Date(sero$tmin, format('%d/%m/%Y'))
  sero$tmax <- as.Date(sero$tmax, format('%d/%m/%Y'))
  sero <- sero[!(is.na(sero$tmin)), ]
  dates <- seq.Date(as.Date('2020-01-22'),by=1, length.out=Ndays)
  sero <- sero[!sero$tmin>max(dates), ]
  
  # inputs
  NSero <- nrow(sero)
  SeroAreaIndex <- vector()
  tmin <- vector()
  tmax <- vector()
  for(i in 1:NSero){
    SeroAreaIndex[i] <- which(countries==sero$region[i])
    tmin[i]  <- which(dates==sero$tmin[i])
    tmax[i]  <- which(dates==sero$tmax[i])
  }
  
  # return inputs for model
  return(list(NSero=NSero, NSamples=sero$n, NPos=sero$n_pos, SeroAreaIndex=SeroAreaIndex, tmin=tmin, tmax=tmax))
}


# Run stan model 
runStan <- function(model, inputs, N_iter, N_chains, max_td, ad){
  f <- stan(file=model,data=inputs,iter=N_iter,
            chains=N_chains,control=list(max_treedepth=max_td, adapt_delta=ad))
  beepr::beep('coin')
  return(f)
}



#----- 2. Functions for model outputs -----#

# Plot age/sex IFRs
plot_IFR_age <- function(chains, inputs){
  
  # min & max bounds of age groupings
  Amin <- seq(0,80,5)
  Amax <- c(seq(4,79,5),110)
  
  # dataframe for ests
  ifr <- data.frame(mean=NA, ciL=NA, ciU=NA, 
                    sex=c(rep('M',17),rep('F',17), rep('B',17), rep('RR',17)))
  for(a in 1:17){
    ifr[a,1:3] <- quantile(chains$ifr_m[,a], c(0.5,0.025,0.975))
    ifr[a+17,1:3] <- quantile(chains$ifr_f[,a], c(0.5,0.025,0.975))
    ifr[a+34,1:3] <- quantile(chains$ifr_b[,a], c(0.5,0.025,0.975))
    ifr[a+51,1:3] <- quantile(chains$ifr_RR[,a], c(0.5,0.025,0.975))
  }
  
  # plot
  ifr$age <- rep(seq(1:17),2)
  ifr$ageG <- paste(Amin,Amax, sep='-')
  ifr$ageG[ifr$ageG=='65-110'] <- '65+'
  ifrAp <- ggplot(ifr[ifr$sex %in% c('M','F'),], aes(reorder(ageG,age), mean*100, col=sex))+ 
    geom_point(position=position_dodge(width=0.4))+
    geom_linerange(aes(ymin=ciL*100, ymax=ciU*100, col=sex), position=position_dodge(width=0.4))+
    theme_minimal()+ ylab('IFR (%)')+ xlab('')+ labs(col='')+
    theme(axis.text.x=element_text(angle=60,hjust=1), legend.position=c(0.1,0.85))+
    scale_colour_manual(values=c('indianred1','royalblue1'), labels=c('female','male'))
  
  # return plot & data
  return(list(ifrA=ifr, ifrAp=ifrAp))
}


# Plot fit to Diamond Princess & Charles de Gaulle data
fit_active <- function(chains, inputs){
  
  # model estimates
  active_fit <- data.frame(outbreak=c('Diamond Princess','Charles de Gaulle'), 
                           N=c(15,0), mean=NA, ciL=NA, ciU=NA)
  active_fit[1,3:5] <- quantile(chains$estDPdeaths, c(0.5,0.025,0.975))
  active_fit[2,3:5] <- quantile(chains$estCDGdeaths, c(0.5,0.025,0.975))
  
  m_fit <- ggplot(active_fit, aes(outbreak, N))+ geom_col(fill='royalblue')+ 
    geom_point(aes(outbreak,mean))+ xlab('')+ ylab('N deaths')+
    geom_linerange(aes(outbreak,ymin=ciL,ymax=ciU))+ theme_minimal()+ 
    labs(subtitle='Predicted deaths from active surveillance campaigns')
  
  return(m_fit)
}


# Seroprevalence point fits
serofit <- function(chains, sero){
  
  sfit <- sero[ ,1:12]
  sfit[c('fit','fciL','fciU')] <- NA
  for(i in 1:nrow(sfit)){
    sfit[i,c('fit','fciL','fciU')] <- quantile(chains$serofit[,i], c(0.5,0.025,0.975))
  } 
  return(sfit)
}


# Seroprevance time series
sero_time <- function(chains, Inputs, countries, continent){
  
  # matrices for samples
  seroT <- matrix(ncol=5, nrow=Inputs$Ndays*length(countries))
  infecT <- matrix(ncol=5, nrow=Inputs$Ndays*length(countries))
  seroT[ ,1] <- sort(rep(countries,Inputs$Ndays))
  infecT[ ,1] <- sort(rep(countries,Inputs$Ndays))

  # extract estimates
  ind <- 0
  for(c in 1:length(countries)){
    for(i in 1:Inputs$Ndays){
      ind <- ind+1
      seroT[ind,3:5] <- quantile(chains$seroT[,i,c], c(0.5,0.025,0.975))
      infecT[ind,3:5] <- quantile(chains$infecT[,i,c], c(0.5,0.025,0.975))
      seroT[ind,2] <- continent[c]
      infecT[ind,2] <- continent[c]
    }
  }
  
  # to dataframe
  seroT <- as.data.frame(seroT)
  infecT <- as.data.frame(infecT)
  colnames(seroT) <- c('country','continent','mean','ciL','ciU')
  colnames(infecT) <- c('country','continent','mean','ciL','ciU')
  
  return(list(seroT=seroT, infecT=infecT))
}


# Plot model fits
fit_deaths <- function(chains, inputs, data, countries){
  
  # age groups for plots
  data$ageG <- paste(data$age_min, data$age_max, sep='-')
  data$ageG[data$age_max==110] <- paste(data$age_min[data$age_max==110], '+', sep='')
  indexA <- index_ages(data, countries, NAges=17)
  
  # list to store plots & estimates
  plotD <- list()
  fitD <- list()
  
  for(c in 1:length(countries)){
    
    dfc <- data[data$country==countries[c], ]
    nages <- length(unique(dfc$age_min))
    
    if(dfc$sex[1]=='B'){
      
      # extract samples
      fit_b <- t(apply(chains$natDeath_b[,,c], 2, quantile, probs=c(0.5,0.025,0.975)))
      
      # aggregate to data age groups
      fit_ba <- matrix(NA, nrow=nages, ncol=3)
      for(a in 1:nages){
        for(ci in 1:3){
          fit_ba[a,ci] <- sum(fit_b[indexA$ageG_min[a,c]:indexA$ageG_max[a,c],ci])
        }
      }
      
      fit <- data.frame(country=countries[c], continent=dfc$continent[1], 
                        n=dfc$deaths, sex=rep('B',nages), age=dfc$ageG, age_min=dfc$age_min,
                        age_max=dfc$age_max, fit=fit_ba[,1], ciL=fit_ba[,2], ciU=fit_ba[,3])
      
      
      fitD[[c]] <- fit
      plotD[[c]] <- ggplot(fit, aes(reorder(age,age_min), n))+ geom_col(fill='seagreen2')+
        xlab('')+ ylab('')+ theme_minimal()+ labs(subtitle=paste(inputs$county[c]))+
        theme(axis.text.x=element_text(angle=60, hjust=1),legend.position=c(0.15,0.9))+
        labs(fill='')+ geom_point(aes(age, fit))+
        geom_linerange(aes(age, ymin=ciL, ymax=ciU))+
        labs(subtitle=paste(dfc$country[1]))
      
    }else{
      
      # extract samples
      fit_m <- t(apply(chains$natDeath_m[,,c], 2, quantile, probs=c(0.5,0.025,0.975)))
      fit_f <- t(apply(chains$natDeath_f[,,c], 2, quantile, probs=c(0.5,0.025,0.975)))
      
      # aggregate to data age groups
      fit_ma <- matrix(NA, nrow=nages, ncol=3)
      fit_fa <- matrix(NA, nrow=nages, ncol=3)
      for(a in 1:nages){
        for(ci in 1:3){
          fit_ma[a,ci] <- sum(fit_m[indexA$ageG_min[a,c]:indexA$ageG_max[a,c], ci])
          fit_fa[a,ci] <- sum(fit_f[indexA$ageG_min[a,c]:indexA$ageG_max[a,c], ci])
        }
      }
      
      dfc$sex <- factor(dfc$sex,levels=c('M','F'))
      dfc <- with(dfc, dfc[order(sex),])
      fit <- data.frame(country=countries[c], continent=dfc$continent[1], 
                        n=dfc$deaths, sex=dfc$sex, age_min=dfc$age_min, 
                        age_max=dfc$age_max, age=dfc$ageG,fit=c(fit_ma[,1], fit_fa[,1]), 
                        ciL=c(fit_ma[,2], fit_fa[,2]),ciU=c(fit_ma[,3], fit_fa[,3]))
      
      fit$sex <- factor(fit$sex, levels=c('F','M'))
      fitD[[c]] <- fit
      plotD[[c]] <- ggplot(fit, aes(reorder(age,age_min), n))+ xlab('')+ ylab('')+ theme_minimal()+ 
        geom_col(aes(fill=sex), position=position_dodge())+
        labs(subtitle=paste(inputs$county[c]))+ labs(fill='')+ 
        theme(axis.text.x=element_text(angle=60, hjust=1),legend.position=c(0.15,0.9))+
        geom_point(aes(age, fit, group=sex), position=position_dodge(width=0.9))+
        geom_linerange(aes(age, ymin=ciL, ymax=ciU, group=sex), position=position_dodge(width=0.9))+
        scale_fill_manual(values=c('indianred1','royalblue1'), labels=c('female','male'))+
        labs(subtitle=paste(dfc$country[1]))
      
    }
  }
  fitD <- do.call('rbind',fitD)
  return(list(plots=plotD, ests=fitD))
}


# Plot country-specific IFRs
plot_IFR_area <- function(chains, inputs, countries){
  
  # country-specific IFRs
  ifrc <- data.frame(country=countries, ifr=NA, ciL=NA, ciU=NA)
  for(c in 1:length(countries)){
    ifrc[c,2:4] <- quantile(chains$ifr_C[,c], c(0.5,0.025,0.975))
  }
  
  # plot
  ifrCp <- ggplot(ifrc, aes(reorder(country,ifr), ifr))+ 
    geom_point(col='royalblue')+
    geom_linerange(aes(ymin=ciL, ymax=ciU), col='royalblue')+
    theme_minimal()+ ylab('IFR')+ xlab('')+
    theme(axis.text.x=element_text(angle=60,hjust=1))
  
  # return plot & data
  return(list(ifrc=ifrc, ifrCp=ifrCp))
}


# Function to plot immunity over time fits
plot_immunity <- function(chains, inputs, countries, plotfit=F){
  
  # dates for plots
  dates <- seq.Date(as.Date('2020-01-22'),by=1, length.out=inputs$Ndays)
  
  # plots
  plots <- list()
  for(c in 1:length(countries)){
    
    # extract values
    dts <- t(apply(chains$seroT[,,c], 2, quantile, probs=c(0.5,0.025,0.975, 0.25,0.75)))
    dts <- as.data.frame(dts)
    colnames(dts) <- c('mean','ciL','ciU','mciL','mciU')
    
    # plot
    dts$date <- dates[1:inputs$Ndays]
    plots[[c]] <- ggplot()+ geom_line(data=dts, aes(date, mean),col='green4')+
      theme_minimal()+ geom_ribbon(data=dts,aes(date, ymin=ciL, ymax=ciU), fill='lightgreen',alpha=0.4)+
      geom_ribbon(data=dts, aes(date, ymin=mciL, ymax=mciU), fill='lightgreen',alpha=0.7)+ xlab('')+
      ylab('')+ labs(subtitle=paste(countries[c]))
    
    # plot fit if serology data used
    if(plotfit==T & any(inputs$SeroAreaInd==c)){
      seroc <- data.frame(NPos=inputs$NPos[inputs$SeroAreaInd==c],NSamples=inputs$NSamples[inputs$SeroAreaInd==c],
                          tmin=inputs$tmin[inputs$SeroAreaInd==c],tmax=inputs$tmax[inputs$SeroAreaInd==c])
      seroc[,c('prop','sciL','sciU')] <- binom.exact(seroc$NPos, seroc$NSamples)[,3:5]
      seroc$dmin <- dates[seroc$tmin]
      seroc$dmax <- dates[seroc$tmax]
      seroc$dmean <- (seroc$dmax - seroc$dmin)/2 + seroc$dmin
      
      plots[[c]] <- plots[[c]] + geom_point(data=seroc, aes(dmean, prop))+ 
        geom_linerange(data=seroc, aes(x=dmean, ymin=sciL, ymax=sciU))
    }
  }
  return(plots)
}

