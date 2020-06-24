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
  
  // population by 5 year age groups, excluding LTC
  matrix[17,NArea] pop_m;
  matrix[17,NArea] pop_f;
  int pop_b[17,NArea];
  
  // population by 5 year age groups, total
  matrix[17,NArea] Tpop_m;
  matrix[17,NArea] Tpop_f;
  matrix[17,NArea] Tpop_b;
  
  // age-specific death data
  int deaths_m[17,NArea];
  int deaths_f[17,NArea];
  int deaths_b[17,NArea];
  
  // min & max of age bands for indexing
  int ageG_min[17,NArea];
  int ageG_max[17,NArea];
  
  // age-specific probabilities of infection
  real relProbInfection[17];
  real agemid[17];
  
  // Diamond Princess data for likelihood
  int <lower=0> DP_pos_m[8]; 
  int <lower=0> DP_pos_f[8]; 
  int <lower=0> DPamin[8];
  int <lower=0> DPamax[8];
  int <lower=0> DP_deathsTot;
  
  // CDG data for likelihood
  int <lower=0> CDG_pos_m[4]; 
  int <lower=0> CDG_pos_f[4]; 
  int <lower=0> CDGamin[4];
  int <lower=0> CDGamax[4];
  int <lower=0> CDG_deathsTot;
  
  // serology data
  int NSero;
  int SeroAreaInd[NSero];
  int NSamples[NSero];
  int NPos[NSero];
  int tmin[NSero];
  int tmax[NSero];
 
  // time series of deaths
  int Ndays;
  int deathsTinfec[Ndays,NArea];
  int deathsTsero[Ndays,NArea];
  int deathsT[Ndays,NArea];
  int TdeathsA[NArea];
  
}

parameters {
  
  real <lower=-50, upper=-0.001> log_probInfec[NArea];  
  real <lower=0> slope_m;
  real intercept; 
  real <lower=-50, upper=-0.001> log_ifr_mY[6];
  real <lower=-50, upper=-0.001> log_ifr_fY[10];
  real log_relifrsex; 
  //real phi; 
  
}

transformed parameters {
  
  // cumulative prob infection
  real probInfec[NArea];

  // infection fatality ratios
  real ifr_f[17];
  real ifr_m[17];
  real log_ifr_m[17];
  real relifrsex;
  
  // estimated deaths by age sex & area
  real natDeath_m[17,NArea];
  real natDeath_f[17,NArea];
  real natDeath_b[17,NArea];
  
  // time series of infections & seroprevalence
  matrix[Ndays,NArea] seroT;
  matrix[Ndays,NArea] infecT;
  
  // transformed parameters
  for(c in 1:NArea) probInfec[c] = exp(log_probInfec[c]);
  relifrsex = exp(log_relifrsex);
  for(a in 1:6) log_ifr_m[a] = log_ifr_mY[a];
  for(a in 1:10) ifr_f[a] = exp(log_ifr_fY[a]);
  for(a in 7:17) log_ifr_m[a] = slope_m*(agemid[a]) + intercept;
  for(a in 1:17) ifr_m[a] = exp(log_ifr_m[a]);
  for(a in 11:17) ifr_f[a] = exp(log_ifr_m[a])*relifrsex;

  
  // probs by age, sex & region
  for (a in 1:17){
    for(c in 1:NArea){
      natDeath_m[a,c] = pop_m[a,c]*probInfec[c]*relProbInfection[a]*ifr_m[a];
      natDeath_f[a,c] = pop_f[a,c]*probInfec[c]*relProbInfection[a]*ifr_f[a];
      natDeath_b[a,c] = natDeath_f[a,c] + natDeath_m[a,c];
    }
  }
  
  // distribute immunity over time
  for(c in 1:NArea){
    seroT[1:Ndays,c] = (probInfec[c]/deathsT[TdeathsA[c],c])*to_vector(deathsTsero[1:Ndays,c]);
    infecT[1:Ndays,c] = (probInfec[c]/deathsT[TdeathsA[c],c])*to_vector(deathsTinfec[1:Ndays,c]);
  }
}


model {

  real estDeaths_b[17,NArea];
  real estDeaths_m[17,NArea];
  real estDeaths_f[17,NArea];
  real estDPdeaths;
  real estCDGdeaths;
  real dpifr_f[8];
  real dpifr_m[8];
  real cdgifr_f[4];
  real cdgifr_m[4];

  // Priors
  for(c in 1:NArea) log_probInfec[c] ~ uniform(-50,-0.001);
  for(a in 1:6) log_ifr_mY[a] ~ uniform(-50,-0.001);
  for(a in 1:10) log_ifr_fY[a] ~ uniform(-50,-0.001);
  intercept ~ uniform(-50,-0.001);
  log_relifrsex ~ uniform(-5,5);
  slope_m ~ normal(0,1);
  
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

  // align IFRs to DP & CDG age groups
  dpifr_m = alignMEAN(ifr_m, 8, DPamin, DPamax);
  dpifr_f = alignMEAN(ifr_f, 8, DPamin, DPamax);
  cdgifr_m = alignMEAN(ifr_m, 4, CDGamin, CDGamax);
  cdgifr_f = alignMEAN(ifr_f, 4, CDGamin, CDGamax);
 
  // Sum expected deaths across age groups
  estDPdeaths=0;
  estCDGdeaths=0;
  for(j in 1:8){ 
    estDPdeaths=estDPdeaths+(dpifr_f[j]*DP_pos_f[j]+dpifr_m[j]*DP_pos_m[j]);
  }
  for(j in 1:4){ 
    estCDGdeaths=estCDGdeaths+(cdgifr_f[j]*CDG_pos_f[j]+cdgifr_m[j]*CDG_pos_m[j]);
  }
  
  // Likelihood
  DP_deathsTot ~ poisson(estDPdeaths);
  CDG_deathsTot ~ poisson(estCDGdeaths);
  for(c in 1:NSero){
    NPos[c] ~ binomial(NSamples[c], mean(seroT[tmin[c]:tmax[c],SeroAreaInd[c]]));
  }
}
 
generated quantities {
  
  real ifr_C[NArea]; 
  real estDPdeaths;
  real estCDGdeaths;
  real dpifr_f[8];
  real dpifr_m[8];
  real cdgifr_f[4];
  real cdgifr_m[4];
  
  
  // population-weighted IFRs
  for(c in 1:NArea){
    ifr_C[c] = sum(to_vector(ifr_m).*to_vector(Tpop_m[,c]) + to_vector(ifr_f).*to_vector(Tpop_f[,c]))/sum(Tpop_b[,c]);
  }
  
  // align to DP & CDG age groups
  dpifr_m = alignMEAN(ifr_m, 8, DPamin, DPamax);
  dpifr_f = alignMEAN(ifr_f, 8, DPamin, DPamax);
  cdgifr_m = alignMEAN(ifr_m, 4, CDGamin, CDGamax);
  cdgifr_f = alignMEAN(ifr_f, 4, CDGamin, CDGamax);
 
  // sum expected deaths across age groups
  estDPdeaths=0;
  estCDGdeaths=0;
  for(j in 1:8){ 
    estDPdeaths=estDPdeaths+(dpifr_f[j]*DP_pos_f[j]+dpifr_m[j]*DP_pos_m[j]);
  }
  for(j in 1:4){ 
    estCDGdeaths=estCDGdeaths+(cdgifr_f[j]*CDG_pos_f[j]+cdgifr_m[j]*CDG_pos_m[j]);
  }


}


