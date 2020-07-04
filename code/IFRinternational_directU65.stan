functions {
  
  real[] alignSUM(real[] d, int Nage, int[] AMin, int[] AMax){
    
    real Da[Nage];
    for(a in 1:Nage) Da[a] = sum(d[AMin[a]:AMax[a]]);
    return(Da);
  }
  
  real[] alignMEAN(real[] d, int Nage, int[] AMin, int[] AMax){
    
    real Da[Nage];
    for(a in 1:Nage) Da[a] = mean(d[AMin[a]:AMax[a]]);
    return(Da);
  }
}



data {
   
  // Number of areas, age-groups & gender by area
  int <lower=1> NArea; 
  int <lower=1> NAges[NArea];
  int gender[NArea]; 
  
  // Population by 5 year age groups, excluding LTC
  matrix[13,NArea] pop_m;
  matrix[13,NArea] pop_f;
  int pop_b[13,NArea];
  
  // Population by 5 year age groups, total
  matrix[13,NArea] Tpop_m;
  matrix[13,NArea] Tpop_f;
  matrix[13,NArea] Tpop_b;
  
  // Age-specific death data
  int deaths_m[13,NArea];
  int deaths_f[13,NArea];
  int deaths_b[13,NArea];
  
  // Min & max of age bands for indexing
  int ageG_min[13,NArea];
  int ageG_max[13,NArea];
  real agemid[13];
  
  // Age-specific probabilities of infection
  real relProbInfection[13];
  
  // Charles de Gaulle data 
  int <lower=0> CDG_pos_m[4]; 
  int <lower=0> CDG_pos_f[4]; 
  int <lower=0> CDGamin[4];
  int <lower=0> CDGamax[4];
  int <lower=0> CDG_deathsTot;
  
  // Serology data
  int NSero;
  int SeroAreaInd[NSero];
  int NSamples[NSero];
  int NPos[NSero];
  int tmin[NSero];
  int tmax[NSero];
 
  // Time series of deaths
  int Ndays;
  int deathsTinfec[Ndays,NArea];
  int deathsTsero[Ndays,NArea];
  int deathsT[Ndays,NArea];
  int TdeathsA[NArea];
  
}

parameters {
  
  real <lower=-50, upper=-0.001> log_probInfec[NArea];  
  real <lower=-50, upper=-0.001> log_ifr_m[13];
  real <lower=-50, upper=-0.001> log_ifr_f[13];

}

transformed parameters {
  
  // cumulative prob infection
  real probInfec[NArea];

  // infection fatality ratios
  real ifr_f[13];
  real ifr_m[13];
  
  // estimated deaths by age sex & area
  real natDeath_m[13,NArea];
  real natDeath_f[13,NArea];
  real natDeath_b[13,NArea];
  
  // time series of infections & seroprevalence
  matrix[Ndays,NArea] seroT;
  matrix[Ndays,NArea] infecT;
  
  // total deaths <65 years
  real u65deaths[NArea];
  
  // expected seroprevalence at location & time of serosuvey
  real serofit[NSero];
  
  // transformed parameters
  for(c in 1:NArea) probInfec[c] = exp(log_probInfec[c]);
  for(a in 1:13) ifr_m[a] = exp(log_ifr_m[a]);
  for(a in 1:13) ifr_f[a] = exp(log_ifr_f[a]);
  
  // age-specific attack rates
  for(a in 1:13){
    ifr_m[a] = ifr_m[a]/relProbInfection[a];
    ifr_f[a] = ifr_f[a]/relProbInfection[a];
  }
  
  // probs by age, sex & region
  for (a in 1:13){
    for(c in 1:NArea){
      natDeath_m[a,c] = pop_m[a,c]*probInfec[c]*ifr_m[a];
      natDeath_f[a,c] = pop_f[a,c]*probInfec[c]*ifr_f[a];
      natDeath_b[a,c] = natDeath_f[a,c] + natDeath_m[a,c];
    }
  }
  
  // total expected deaths <65
  for(c in 1:NArea){
    u65deaths[c] = sum(alignSUM(natDeath_b[,c], NAges[c], ageG_min[,c], ageG_max[,c]));
  } 
  
  // distribute immunity over time
  for(c in 1:NArea){
    seroT[1:Ndays,c] = (probInfec[c]/deathsT[TdeathsA[c],c])*to_vector(deathsTsero[1:Ndays,c]);
    infecT[1:Ndays,c] = (probInfec[c]/deathsT[TdeathsA[c],c])*to_vector(deathsTinfec[1:Ndays,c]);
  }
  
  // expected seroprevalence at time of serosurveys
  for(i in 1:NSero) serofit[i] = mean(seroT[tmin[i]:tmax[i],SeroAreaInd[i]]);
}


model {

  real estDeaths_b[13,NArea];
  real estDeaths_m[13,NArea];
  real estDeaths_f[13,NArea];

  // Priors
  for(c in 1:NArea) log_probInfec[c] ~ uniform(-50,-0.001);
  for(a in 1:13) log_ifr_m[a] ~ uniform(-50,-0.001);
  for(a in 1:13) log_ifr_f[a] ~ uniform(-50,-0.001);
  
  // fit to age & sex-specific data
  for(c in 1:NArea){
    
    if(gender[c]==1){
      estDeaths_b[1:NAges[c],c] = alignSUM(natDeath_b[,c], NAges[c], ageG_min[,c], ageG_max[,c]);
      deaths_b[1:NAges[c],c] ~ poisson(estDeaths_b[1:NAges[c],c]);
    }
    if(gender[c]==2){
      estDeaths_m[1:NAges[c],c] = alignSUM(natDeath_m[,c], NAges[c], ageG_min[,c], ageG_max[,c]);
      estDeaths_f[1:NAges[c],c] = alignSUM(natDeath_f[,c], NAges[c], ageG_min[,c], ageG_max[,c]);
      deaths_m[1:NAges[c],c] ~ poisson(estDeaths_m[1:NAges[c],c]);
      deaths_f[1:NAges[c],c] ~ poisson(estDeaths_f[1:NAges[c],c]);
    }
  }
  
  // Likelihood
  for(i in 1:NSero){
    NPos[i] ~ binomial(NSamples[i], mean(seroT[tmin[i]:tmax[i],SeroAreaInd[i]]));
  }
}
 
generated quantities {
  
  real ifr_C[NArea]; 
  real ifr_b[13];
  real ifr_RR[13];
  real estCDGdeaths;
  real cdgifr_f[4];
  real cdgifr_m[4];
  
  
  // population-weighted IFRs
  for(c in 1:NArea){
    ifr_C[c] = sum(to_vector(ifr_m).*to_vector(Tpop_m[,c]) + to_vector(ifr_f).*to_vector(Tpop_f[,c]))/sum(Tpop_b[,c]);
  }
  
  // average IFR for males and females
  for(a in 1:13) ifr_b[a] = (ifr_m[a]+ifr_m[a])/2;
  
  // IFR relative to 55-59 group
  for(a in 1:13) ifr_RR[a] = ifr_b[a]/ifr_b[12];
  
  // align to CDG age groups
  cdgifr_m = alignMEAN(ifr_m, 4, CDGamin, CDGamax);
  cdgifr_f = alignMEAN(ifr_f, 4, CDGamin, CDGamax);
 
  // sum expected deaths across age groups for CDG
  estCDGdeaths=0;
  for(j in 1:4){ 
    estCDGdeaths=estCDGdeaths+(cdgifr_f[j]*CDG_pos_f[j]+cdgifr_m[j]*CDG_pos_m[j]);
  }

}