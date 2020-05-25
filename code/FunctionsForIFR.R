#===== Functions for International IFR model analysis =====#


# Function to compile demographic data for model
compile_pop <- function(poplist, countries){
  
  # matrices for model
  pop_m <- matrix(NA, nrow=17, ncol=length(poplist)) 
  pop_f <- matrix(NA, nrow=17, ncol=length(poplist)) 
  pop_b <- matrix(NA, nrow=17, ncol=length(poplist)) 
  
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
compile_deathsA <- function(deathsA, countries){
  
  # matrices for model
  deaths_m <- matrix(-999, nrow=17, ncol=length(countries))
  deaths_f <- matrix(-999, nrow=17, ncol=length(countries))
  deaths_b <- matrix(-999, nrow=17, ncol=length(countries))
  
  # compile for model
  for(c in 1:length(countries)){
    
    dfc <- deathsA[deathsA$country==countries[c], ]
    if(dfc$sex[1]=='B'){
      deaths_b[1:nrow(dfc),c] <- dfc$deaths
    }else{
      deaths_m[1:(nrow(dfc)/2),c] <- dfc$deaths[dfc$sex=='M']
      deaths_f[1:(nrow(dfc)/2),c] <- dfc$deaths[dfc$sex=='F']
    }
  }
  
  return(list(deaths_m=deaths_m, deaths_f=deaths_f, deaths_b=deaths_b))
}


# Function to aggregate deaths 80+ (if broken down)
agg_deaths80p <- function(deathsA, countries){
  
  # list for aggregated/tidied dataframes
  agg <- list()
  
  # for each country
  for(c in 1:length(countries)){
    
    dfc <- deathsA[deathsA$country==countries[c], ]
    
    # if multiple age groups 80+, aggregate
    if(any(dfc$age_min>80)){ 
      new <- dfc[dfc$age_max<80, ]
      for(s in unique(new$sex)){
        d80p <- sum(dfc$deaths[dfc$age_max>80 & dfc$sex==s])
        minA <- min(dfc$age_min[dfc$age_max>80])
        n80p <- data.frame(country=countries[c], sex=s, age_min=minA, age_max=110,
                           deaths=as.numeric(d80p), asof=new$asof[1], nat.reg=new$nat.reg[1])
        new <- rbind(new,n80p)
      }
      agg[[c]] <- new
    }else{
      agg[[c]] <- dfc
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
index_ages <- function(data, countries){
  
  # min & max bounds of age groupings
  Amin <- seq(0,80,5)
  Amax <- c(seq(4,79,5),110)
  
  # matrices for storage
  ageG_min <- matrix(-999, nrow=17, ncol=length(countries))
  ageG_max <- matrix(-999, nrow=17, ncol=length(countries))
  
  # get indices
  for(c in 1:length(countries)){
    
    dfc <- data[data$country==countries[c], ]
    
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


# Run stan model 
runStan <- function(model, inputs, N_iter, N_chains, max_td, ad){
  f <- stan(file=model,data=inputs,iter=N_iter,
            chains=N_chains,control=list(max_treedepth=max_td, adapt_delta=ad))
  beepr::beep('coin')
  return(f)
}


# Plot fit to Diamond Princess data
fit_DP <- function(fit, inputs){
  
  fit_m <- t(apply(chains$ifr_m, 2, quantile, probs=c(0.5,0.025,0.975)))
  fit_f <- t(apply(chains$ifr_f, 2, quantile, probs=c(0.5,0.025,0.975)))
  
  df <- data.frame(Obs=inputs$DP_deaths, expec=sum(fit_m[,1]*inputs$DP_pos_m + fit_f[,1]*inputs$DP_pos_f),
                   ciL=sum(fit_m[,2]*inputs$DP_pos_m + fit_f[,2]*inputs$DP_pos_f), 
                   ciU=sum(fit_m[,3]*inputs$DP_pos_m + fit_f[,3]*inputs$DP_pos_f))
  
  DP_fit <- ggplot(df, aes(x=1, Obs))+ geom_col(fill='royalblue')+ geom_point(aes(x=1,expec))+
    geom_linerange(aes(x=1,ymin=ciL,ymax=ciU))+ theme_minimal()+ 
    theme(axis.text.x=element_blank())+ xlab('')+ ylab('N')+
    labs(subtitle='Diamond Princess deaths')
  
  png(filename='DP_Fit.png', width=10, height=10, res=400, units='cm')
  plot(DP_fit)
  dev.off()
}


# Plot model fits
fit_deaths <- function(chains, inputs, data, countries){
  
  # age mid & age groups for plots
  data$age_mid <- data$age_min+((data$age_max-data$age_min)+1)/2
  data$ageG <- paste(data$age_min, data$age_max, sep='-')
  data$ageG[data$age_max==110] <- paste(data$age_min[data$age_max==110], '+', sep='')
  
  # list to store plots
  plotD <- list()
  
  for(c in 1:length(countries)){
    
    dfc <- data[data$country==countries[c], ]
    nages <- inputs$NAges[c]
    
    if(dfc$sex[1]=='B'){
      
      # extract samples
      fit_b <- t(apply(chains$natDeath_b[,,c], 2, quantile, probs=c(0.5,0.025,0.975)))
      
      # aggregate to data age groups
      fit_ba <- matrix(NA, nrow=nages, ncol=3)
      for(a in 1:nages){
        for(ci in 1:3){
          fit_ba[a,ci] <- sum(fit_b[inputs$ageG_min[a,c]:inputs$ageG_max[a,c],ci])
        }
      }
      
      fit <- data.frame(n=inputs$deaths_b[1:nages,c], sex=rep('B',nages), amin=dfc$age_min,
                        age=dfc$ageG, fit=fit_ba[,1], ciL=fit_ba[,2], ciU=fit_ba[,3])
      
      plotD[[c]] <- ggplot(fit, aes(reorder(age,amin), n))+ geom_col(fill='seagreen2')+
        xlab('')+ ylab('')+ theme_minimal()+ labs(subtitle=paste(inputs$county[c]))+
        theme(axis.text.x=element_text(angle=60, hjust=1),legend.position=c(0.15,0.9))+
        labs(fill='')+ geom_point(aes(age, fit))+geom_linerange(aes(age, ymin=ciL, ymax=ciU))+
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
          fit_ma[a,ci] <- sum(fit_m[inputs$ageG_min[a,c]:inputs$ageG_max[a,c], ci])
          fit_fa[a,ci] <- sum(fit_f[inputs$ageG_min[a,c]:inputs$ageG_max[a,c], ci])
        }
      }
      
      fit <- data.frame(n=c(inputs$deaths_m[1:nages,c],inputs$deaths_f[1:nages,c]),
                        sex=c(rep('M',nages),rep('F',nages)), amin=sort(unique(dfc$age_min)), age=NA,
                        fit=c(fit_ma[,1], fit_fa[,1]), ciL=c(fit_ma[,2], fit_fa[,2]),
                        ciU=c(fit_ma[,3], fit_fa[,3]))
      
      for(a in 1:nrow(fit)) fit$age[a] <- dfc$ageG[dfc$age_min==fit$amin[a]][1]
      fit$sex <- factor(fit$sex, levels=c('F','M'))
      plotD[[c]] <- ggplot(fit, aes(reorder(age,amin), n))+ geom_col(aes(fill=sex), position=position_dodge())+
        xlab('')+ ylab('')+ theme_minimal()+ labs(subtitle=paste(inputs$county[c]))+
        theme(axis.text.x=element_text(angle=60, hjust=1),legend.position=c(0.15,0.9))+
        labs(fill='')+ geom_point(aes(age, fit, group=sex), position=position_dodge(width=0.9))+
        geom_linerange(aes(age, ymin=ciL, ymax=ciU, group=sex), position=position_dodge(width=0.9))+
        scale_fill_manual(values=c('indianred1','royalblue1'), labels=c('female','male'))+
        labs(subtitle=paste(dfc$country[1]))
    }
  }
  N <- ceiling(length(countries)/10)
  for(p in 1:N){
    u <- p*10
    l <- u-9
    png(filename=paste('ModelFit',p,'.png',sep=''), width=18, height=27, res=400, units='cm')
    grid.arrange(grobs=plotD[l:u], ncol=2, left='deaths', bottom='age')
    dev.off()
  }
  return(plotD)
}


# Plot age/sex IFRs
plot_IFR <- function(fit){
  
  # extract samples
  chains <- rstan::extract(fit)
  
  # min & max bounds of age groupings
  Amin <- seq(0,80,5)
  Amax <- c(seq(4,79,5),110)
  
  # dataframe for ests
  ifr_f <- data.frame(mean=NA, ciL=NA, ciU=NA)
  ifr_m <- data.frame(mean=NA, ciL=NA, ciU=NA)
  for(a in 1:17){
    ifr_f[a,] <- quantile(chains$ifr_f[,a], c(0.5,0.025,0.975))
    ifr_m[a,] <- quantile(chains$ifr_m[,a], c(0.5,0.025,0.975))
  }
  
  ifr <- rbind(ifr_m, ifr_f)
  ifr$sex <- c(rep('M',17),rep('F',17))
  ifr$age <- rep(seq(1:17),2)
  ifr$ageG <- paste(Amin,Amax, sep='-')
  ifr$ageG[ifr$ageG=='80-110'] <- '80+'
  ifrp <- ggplot(ifr, aes(reorder(ageG,age), mean, col=sex))+ 
    geom_point(position=position_dodge(width=0.4))+
    geom_linerange(aes(ymin=ciL, ymax=ciU, col=sex), position=position_dodge(width=0.4))+
    theme_minimal()+ ylab('IFR')+ xlab('')+ labs(col='')+
    theme(axis.text.x=element_text(angle=60,hjust=1), legend.position=c(0.1,0.85))+
    scale_colour_manual(values=c('indianred1','royalblue1'), labels=c('female','male'))
  
  return(ifrp)
}



# Function to distribute cumulative immunity over time
immune_time <- function(deathsT, probInfec){
  
  
}

