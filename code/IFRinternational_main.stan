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
  
  // Population by 5 year age groups
  matrix[17,NArea] pop_m;
  matrix[17,NArea] pop_f;
  int pop_b[17,NArea];
  
  // Age-specific death data <65 years
  int deaths_m[13,NArea];
  int deaths_f[13,NArea];
  int deaths_b[13,NArea];
  
  // Min & max of age bands for indexing
  int ageG_min[13,NArea];
  int ageG_max[13,NArea];

  // Age-specific probabilities of infection
  real relProbInfection[17];
  
  // Diamond Princess data 
  int <lower=0> DP_pos_m[8]; 
  int <lower=0> DP_pos_f[8]; 
  int <lower=0> DPamin[8];
  int <lower=0> DPamax[8];
  int <lower=0> DP_deathsTot;
  
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
  int deathsT[Ndays+20,NArea];
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

  // estimated deaths by age, sex & area
  real natDeath_m[17,NArea];
  real natDeath_f[17,NArea];
  real natDeath_b[17,NArea];
  
  // time series of infections & seroprevalence
  matrix[Ndays,NArea] seroT;
  matrix[Ndays,NArea] infecT;
  
  // total deaths <65 years
  real u65deaths[NArea];
  
  // expected seroprevalence at location & time of serosurvey
  real serofit[NSero];
  
  // mean increase in log IFR estimates by age (20-64)
  real diff_ifr_m[8]; 
  real diff_ifr_f[8];
  real mean_increase_ifr;
  
  // IFRs
  real log_ifr_fO[4]; 
  real log_ifr_mO[4];
  real ifr_m[17];
  real ifr_f[17];
  real ifr_b[17];
  
  // mean increase in log IFRs 20-64
  for(i in 1:8){
    diff_ifr_m[i] = log_ifr_m[i+5] - log_ifr_m[i+4];
    diff_ifr_f[i] = log_ifr_f[i+5] - log_ifr_f[i+4];
  }
  mean_increase_ifr = exp((mean(diff_ifr_m)+mean(diff_ifr_f))/2);
  
  // infer IFRs >64
  log_ifr_mO[1] = log_ifr_m[13]+ mean(diff_ifr_m);
  log_ifr_fO[1] = log_ifr_f[13]+ mean(diff_ifr_f);
  for(a in 2:4){
    if(a<4){
      log_ifr_mO[a] = log_ifr_mO[a-1] + mean(diff_ifr_m);
      log_ifr_fO[a] = log_ifr_fO[a-1] + mean(diff_ifr_f);
    }else{
      log_ifr_mO[a] = log_ifr_mO[a-1] + mean(diff_ifr_m)*1.5; // assume 85 mid-point for 80+
      log_ifr_fO[a] = log_ifr_fO[a-1] + mean(diff_ifr_f)*1.5;
    }
  }

  // IFRs all ages
  for(a in 1:13) ifr_m[a] = exp(log_ifr_m[a]);
  for(a in 1:13) ifr_f[a] = exp(log_ifr_f[a]);
  for(a in 14:17) ifr_m[a] = exp(log_ifr_mO[a-13])/relProbInfection[a];
  for(a in 14:17) ifr_f[a] = exp(log_ifr_fO[a-13])/relProbInfection[a];
  for(a in 1:17) ifr_b[a] = (ifr_m[a]+ifr_f[a])/2;
  
  // transformed parameters
  for(c in 1:NArea) probInfec[c] = exp(log_probInfec[c]);

  // probs by age, sex & region
  for (a in 1:17){
    for(c in 1:NArea){
      natDeath_m[a,c] = pop_m[a,c]*probInfec[c]*ifr_m[a]*relProbInfection[a];
      natDeath_f[a,c] = pop_f[a,c]*probInfec[c]*ifr_f[a]*relProbInfection[a];
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
  real ifr_RR[17];
  real estDPdeaths;
  real estCDGdeaths;
  real dpifr_f[8];
  real dpifr_m[8];
  real cdgifr_f[4];
  real cdgifr_m[4];
  
  // IFR relative to 55-59 group
  for(a in 1:17) ifr_RR[a] = ifr_b[a]/ifr_b[12];
  
  // population-weighted IFRs
  for(c in 1:NArea){
    ifr_C[c] = sum(to_vector(ifr_m).*to_vector(pop_m[,c]) + to_vector(ifr_f).*to_vector(pop_f[,c]))/sum(pop_b[,c]);
  }
  
  // align to DP & CDG age groups
  dpifr_m = alignMEAN(ifr_m, 8, DPamin, DPamax);
  dpifr_f = alignMEAN(ifr_f, 8, DPamin, DPamax);
  cdgifr_m = alignMEAN(ifr_m, 4, CDGamin, CDGamax);
  cdgifr_f = alignMEAN(ifr_f, 4, CDGamin, CDGamax);
 
  // sum expected deaths across age groups for DP and CDG
  estDPdeaths=0;
  estCDGdeaths=0;
  for(j in 1:8){ 
    estDPdeaths=estDPdeaths+(dpifr_f[j]*DP_pos_f[j]+dpifr_m[j]*DP_pos_m[j]);
  }
  for(j in 1:4){ 
    estCDGdeaths=estCDGdeaths+(cdgifr_f[j]*CDG_pos_f[j]+cdgifr_m[j]*CDG_pos_m[j]);
  }

}


