functions {
  
  real[] align(real[] d, int Nage, int[] AMin, int[] AMax){
    
    real Da[Nage];
    for(a in 1:Nage) Da[a] = sum(d[AMin[a]:AMax[a]]);
    return(Da);
  }
}



data {
   
  int <lower=1> NArea; 
  int <lower=1> NAges[NArea];
  int gender[NArea]; // 1=both, 2=male & female
  
  // population by 5 year age groups
  matrix[19,NArea] pop_m;
  matrix[19,NArea] pop_f;
  matrix[19,NArea] pop_b;
  
  // age-specific death data
  int deaths_m[19,NArea];
  int deaths_f[19,NArea];
  int deaths_b[19,NArea];
  
  // min & max of age bands for indexing
  int ageG_min[19,NArea];
  int ageG_max[19,NArea];
  int NGroups;
 
  int <lower=0> DP_pos_m[NGroups]; 
  int <lower=0> DP_pos_f[NGroups]; 
  int <lower=0> DP_deathsTot;
  real relProbInfection[19];
  
}

parameters {
  
  real logit_probInfec[NArea];  
  real logit_probDeathMi[19];
  real log_relprobDeathSex;
  real phi;
  
}

transformed parameters {
  
  real probInfec[NArea]; 
  real probDeath_m[19,NArea]; 
  real probDeath_f[19,NArea];

  // probs given infection
  real probDeathMi[19];
  real relprobDeathSex;
  real ifr_f[19];
  real ifr_m[19];
  
  // national age sex estimates
  real natDeath_m[19,NArea];
  real natDeath_f[19,NArea];
  real natDeath_b[19,NArea];
  
  // transformed parameters
  relprobDeathSex = exp(log_relprobDeathSex);
  for(c in 1:NArea) probInfec[c] = inv_logit(logit_probInfec[c]);
  for(a in 1:19) probDeathMi[a] = inv_logit(logit_probDeathMi[a]);
  
  // Probs by age, sex & region
  for (a in 1:19){
    for(c in 1:NArea){
      probDeath_m[a,c] = probInfec[c]*relProbInfection[a]*probDeathMi[a];
      probDeath_f[a,c] = probInfec[c]*relProbInfection[a]*probDeathMi[a]*relprobDeathSex;
      natDeath_m[a,c] = probDeath_m[a,c]*pop_m[a,c];
      natDeath_f[a,c] = probDeath_f[a,c]*pop_f[a,c];
      natDeath_b[a,c] = probDeath_m[a,c]*pop_m[a,c] + probDeath_f[a,c]*pop_f[a,c];
    }
    ifr_m[a] = probDeathMi[a]; 
    ifr_f[a] = probDeathMi[a]*relprobDeathSex;
  }
  
}


model {

  int ObservedDeaths;
  real estDPdeaths;
  real estDeaths_b[19,NArea];
  real estDeaths_m[19,NArea];
  real estDeaths_f[19,NArea];
  real difr_f[NGroups];
  real difr_m[NGroups];

  // Priors
  log_relprobDeathSex ~ normal(0.,0.5);
  for(c in 1:NArea) logit_probInfec[c] ~ normal(0.,2);
  for(a in 1:19) logit_probDeathMi[a] ~ normal(0.,2);
  
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
  difr_m[8] = mean(ifr_m[17:19]);
  difr_f[1] = mean(ifr_f[1:4]);
  difr_f[2] = mean(ifr_f[5:6]);
  difr_f[3] = mean(ifr_f[7:8]);
  difr_f[4] = mean(ifr_f[9:10]);
  difr_f[5] = mean(ifr_f[11:12]);
  difr_f[6] = mean(ifr_f[13:14]);
  difr_f[7] = mean(ifr_f[15:16]);
  difr_f[8] = mean(ifr_f[17:19]);
  // Sum expected deaths across age groups
  for (j in 1:NGroups){ 
    estDPdeaths=estDPdeaths+(difr_f[j]*DP_pos_f[j]+difr_m[j]*DP_pos_m[j]);
  }
  
  // Likelihood
  DP_deathsTot ~ poisson(estDPdeaths);

}
 
generated quantities {

  real probDeathFi[19];
  
  // probabilties of death and ICU for females 
  for(a in 1:19) probDeathFi[a] = probDeathMi[a]*relprobDeathSex; 

}

