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
  matrix[14,NArea] pop_m;
  matrix[14,NArea] pop_f;
  int pop_b[14,NArea];
  
  // Age-specific death data <65 years
  int deaths_m[13,NArea];
  int deaths_f[13,NArea];
  int deaths_b[13,NArea];
  //int deathso_m[17,NArea];
  //int deathso_f[17,NArea];
  //int deathso_b[17,NArea];
  
  // Min & max of age bands for indexing
  int ageG_min[13,NArea];
  int ageG_max[13,NArea];
  
  //int NAreao65;
  //int deaths_mo65[4,NAreao65];
  //int deaths_fo65[4,NAreao65];
  //int indNAreao65[NAreao65];

  // Age-specific probabilities of infection
  real relProbInfection[14];
  
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
  real <lower=-50, upper=-0.001> log_ifr_m[14];
  real <lower=-50, upper=-0.001> log_ifr_f[14];
  //real <lower=0.001> relAR65p;
}

transformed parameters {
  
  // cumulative prob infection
  real probInfec[NArea];

  // estimated deaths by age, sex & area
  real natDeath_m[14,NArea];
  real natDeath_f[14,NArea];
  real natDeath_b[14,NArea];
  
  // time series of infections & seroprevalence
  matrix[Ndays,NArea] seroT;
  matrix[Ndays,NArea] infecT;
  
  // total deaths <65 years
  real u65deaths[NArea];
  real o65deaths[NArea];
  
  // expected seroprevalence at location & time of serosurvey
  real serofit[NSero];
  
  // ifr 65+
  //real ifro65f[4,NAreao65];
  //real ifro65m[4,NAreao65];
  
  // mean increase in IFR estimates
  //real diff_ifr_b[13]; // 20+
  //real mean_increase_ifr;
  
  // IFRs
  real ifr_m[14];
  real ifr_f[14];
  real ifr_b[14];
  
  // transformed parameters
  for(c in 1:NArea) probInfec[c] = exp(log_probInfec[c]);
  
  // IFRs all ages
  for(a in 1:14) ifr_m[a] = exp(log_ifr_m[a]);
  for(a in 1:14) ifr_f[a] = exp(log_ifr_f[a]);
  for(a in 1:14) ifr_b[a] = (ifr_m[a]+ifr_f[a])/2;

  // probs by age, sex & region
  for (a in 1:14){
    for(c in 1:NArea){
      natDeath_m[a,c] = pop_m[a,c]*probInfec[c]*ifr_m[a]*relProbInfection[a];
      natDeath_f[a,c] = pop_f[a,c]*probInfec[c]*ifr_f[a]*relProbInfection[a];
      natDeath_b[a,c] = natDeath_f[a,c] + natDeath_m[a,c];
    }
  }
  
  // mean increase in log IFRs 20-64
  //for(i in 1:12){
    //diff_ifr_b[i] = ifr_b[i+5] - ifr_b[i+4];
  //}
  //mean_increase_ifr = mean(diff_ifr_b);
  
  // total expected deaths <65
  for(c in 1:NArea){
    u65deaths[c] = sum(natDeath_b[1:13,c]);
    o65deaths[c] = natDeath_b[14,c];
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
  real estDPdeaths;
  real dpifr_f[8];
  real dpifr_m[8];

  // Priors
  for(c in 1:NArea) log_probInfec[c] ~ uniform(-50,-0.001);
  for(a in 1:14) log_ifr_m[a] ~ uniform(-50,-0.001);
  for(a in 1:14) log_ifr_f[a] ~ uniform(-50,-0.001);
  //relAR65p ~ uniform(0,5);
  
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
  
  // align to DP & CDG age groups
  dpifr_m = alignMEAN(ifr_m, 8, DPamin, DPamax);
  dpifr_f = alignMEAN(ifr_f, 8, DPamin, DPamax);
 
  // sum expected deaths across age groups for DP and CDG
  estDPdeaths=0;
  for(j in 1:8){ 
    estDPdeaths=estDPdeaths+(dpifr_f[j]*DP_pos_f[j]+dpifr_m[j]*DP_pos_m[j]);
  }
  DP_deathsTot ~ poisson(estDPdeaths);
  
}
 
generated quantities {
  
  real ifr_C[NArea]; 
  real ifr_RR[14];
  real estDPdeaths;
  real estCDGdeaths;
  real dpifr_f[8];
  real dpifr_m[8];
  real cdgifr_f[4];
  real cdgifr_m[4];
  
  // IFR relative to 55-59 group
  for(a in 1:14) ifr_RR[a] = ifr_b[a]/ifr_b[12];
  
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


