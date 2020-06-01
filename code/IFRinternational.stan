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
  matrix[17,NArea] pop_b;
  
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
  
}

parameters {
  
  real logit_probInfec[NArea]; // cumulative probability of infection  
  real logit_ifr_m[17]; // prob death given infection (males)
  real log_relifrsex[2]; // relative prob death (females)
  real phi; // over-dispersion parameter
  
}

transformed parameters {
  
  real probInfec[NArea];
  real relifrsex[2];

  // infection fatality ratios
  real ifr_f[17];
  real ifr_m[17];
  
  // estimated deaths by age sex & area
  real natDeath_m[17,NArea];
  real natDeath_f[17,NArea];
  real natDeath_b[17,NArea];
  
  // transformed parameters
  for(s in 1:2) relifrsex[s] = exp(log_relifrsex[s]);
  for(c in 1:NArea) probInfec[c] = inv_logit(logit_probInfec[c]);
  for(a in 1:17) ifr_m[a] = inv_logit(logit_ifr_m[a]);
  for(a in 1:12) ifr_f[a] = ifr_m[a]*relifrsex[1];
  for(a in 13:17) ifr_f[a] = ifr_m[a]*relifrsex[2];
  
  // Probs by age, sex & region
  for (a in 1:17){
    for(c in 1:NArea){
      natDeath_m[a,c] = pop_m[a,c]*probInfec[c]*relProbInfection[a]*ifr_m[a];
      natDeath_f[a,c] = pop_f[a,c]*probInfec[c]*relProbInfection[a]*ifr_f[a];
      natDeath_b[a,c] = natDeath_f[a,c] + natDeath_m[a,c];
    }
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
  for(s in 1:2) log_relifrsex[s] ~ normal(0,1);
  for(c in 1:NArea) logit_probInfec[c] ~ cauchy(0,1);
  for(a in 1:17) logit_ifr_m[a] ~ cauchy(0,1);
  
  // fit to age & sex-specific data
  for(c in 1:NArea){
    
    if(gender[c]==1){
      estDeaths_b[1:NAges[c],c] = alignSUM(natDeath_b[,c], NAges[c], ageG_min[,c], ageG_max[,c]);
      deaths_b[1:NAges[c],c] ~ neg_binomial_2(estDeaths_b[1:NAges[c],c], phi);
    }
    if(gender[c]==2){
      estDeaths_m[1:NAges[c],c] = alignSUM(natDeath_m[,c], NAges[c], ageG_min[,c], ageG_max[,c]);
      estDeaths_f[1:NAges[c],c] = alignSUM(natDeath_f[,c], NAges[c], ageG_min[,c], ageG_max[,c]);
      deaths_m[1:NAges[c],c] ~ neg_binomial_2(estDeaths_m[1:NAges[c],c], phi);
      deaths_f[1:NAges[c],c] ~ neg_binomial_2(estDeaths_f[1:NAges[c],c], phi);
    }
  }
  

  // align to DP & CDG age groups
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
  DP_deathsTot ~ neg_binomial_2(estDPdeaths, phi);
  CDG_deathsTot ~ neg_binomial_2(estCDGdeaths, phi);

}
 
generated quantities {
  
  real ifr_C[NArea]; // population weighted IFRs
  
  for(c in 1:NArea){
    ifr_C[c] = sum(to_vector(ifr_m).*to_vector(Tpop_m[,c]) + to_vector(ifr_f).*to_vector(Tpop_f[,c]))/sum(Tpop_b[,c]);
  }


}

