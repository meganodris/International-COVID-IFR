functions {
  
  real[] align(real[] d, int Nage, int[] AMin, int[] AMax){
    
    real Da[Nage];
    for(a in 1:Nage) Da[a] = sum(d[AMin[a]:AMax[a]]);
    return(Da);
  }
}



data {
   
  // Number of areas, age-groups & gender by area
  int <lower=1> NArea; 
  int <lower=1> NAges[NArea];
  int gender[NArea]; 
  
  // population by 5 year age groups
  matrix[17,NArea] pop_m;
  matrix[17,NArea] pop_f;
  matrix[17,NArea] pop_b;
  
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
  int <lower=0> DP_deathsTot;
  
}

parameters {
  
  real logit_probInfec[NArea]; // cumulative probability of infection  
  real logit_ifr_m[17]; // prob death given infection (males)
  real log_relifrsex; // relative prob death (females)
  real phi; // over-dispersion parameter
  
}

transformed parameters {
  
  real probInfec[NArea];
  real relifrsex;

  // infection fatality ratios
  real ifr_f[17];
  real ifr_m[17];
  
  // estimated deaths by age sex & area
  real natDeath_m[17,NArea];
  real natDeath_f[17,NArea];
  real natDeath_b[17,NArea];
  
  // transformed parameters
  relifrsex = exp(log_relifrsex);
  for(c in 1:NArea) probInfec[c] = inv_logit(logit_probInfec[c]);
  for(a in 1:17) ifr_m[a] = inv_logit(logit_ifr_m[a]);
  for(a in 1:17) ifr_f[a] = ifr_m[a]*relifrsex;
  
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

  int ObservedDeaths;
  real estDPdeaths;
  real estDeaths_b[17,NArea];
  real estDeaths_m[17,NArea];
  real estDeaths_f[17,NArea];
  real difr_f[8];
  real difr_m[8];

  // Priors
  log_relifrsex ~ normal(0.,0.5);
  for(c in 1:NArea) logit_probInfec[c] ~ normal(0.,2);
  for(a in 1:17) logit_ifr_m[a] ~ normal(0.,2);
  
  // fit to age & sex-specific data
  for(c in 1:NArea){
    
    if(gender[c]==1){
      estDeaths_b[1:NAges[c],c] = align(natDeath_b[,c], NAges[c], ageG_min[,c], ageG_max[,c]);
      deaths_b[1:NAges[c],c] ~ neg_binomial_2(estDeaths_b[1:NAges[c],c], phi);
    }
    if(gender[c]==2){
      estDeaths_m[1:NAges[c],c] = align(natDeath_m[,c], NAges[c], ageG_min[,c], ageG_max[,c]);
      estDeaths_f[1:NAges[c],c] = align(natDeath_f[,c], NAges[c], ageG_min[,c], ageG_max[,c]);
      deaths_m[1:NAges[c],c] ~ neg_binomial_2(estDeaths_m[1:NAges[c],c], phi);
      deaths_f[1:NAges[c],c] ~ neg_binomial_2(estDeaths_f[1:NAges[c],c], phi);
    }
  }
  

  estDPdeaths=0;
  difr_m[1] = mean(ifr_m[1:4]);
  difr_m[2] = mean(ifr_m[5:6]);
  difr_m[3] = mean(ifr_m[7:8]);
  difr_m[4] = mean(ifr_m[9:10]);
  difr_m[5] = mean(ifr_m[11:12]);
  difr_m[6] = mean(ifr_m[13:14]);
  difr_m[7] = mean(ifr_m[15:16]);
  difr_m[8] = ifr_m[17];
  difr_f[1] = mean(ifr_f[1:4]);
  difr_f[2] = mean(ifr_f[5:6]);
  difr_f[3] = mean(ifr_f[7:8]);
  difr_f[4] = mean(ifr_f[9:10]);
  difr_f[5] = mean(ifr_f[11:12]);
  difr_f[6] = mean(ifr_f[13:14]);
  difr_f[7] = mean(ifr_f[15:16]);
  difr_f[8] = ifr_f[17];
  // Sum expected deaths across age groups
  for (j in 1:8){ 
    estDPdeaths=estDPdeaths+(difr_f[j]*DP_pos_f[j]+difr_m[j]*DP_pos_m[j]);
  }
  
  // Likelihood
  DP_deathsTot ~ poisson(estDPdeaths);

}
 
generated quantities {


}

