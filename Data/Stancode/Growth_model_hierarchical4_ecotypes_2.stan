data {
  int<lower=0> N;//Number of samples
  
  //data vectors
  vector<lower=0, upper=13> [N]Agei;
  vector<lower=0> [N]Radmm;
  
  //groups
  int<lower=1> n_fish;
  int<lower=1> n_eco;
  int<lower=1, upper=n_fish> fish_id [N];
  int<lower=1, upper=n_eco> eco_id[N];
}

parameters {
  vector<lower=0> [n_fish] LInf_i;
  vector<lower=0> [n_fish] k_i;
  vector<upper=0> [n_fish] tZero_i;
  vector<lower=0> [n_eco] LInf_eco; 
  vector<lower=0> [n_eco] k_eco;
  vector<upper=0> [n_eco] tZero_eco;
  real<lower=0> sigma;
  //hyperprameters
  real<lower=0> LInf;
  real<lower=0> k;
  real<upper=0> tZero;
  vector<lower=0> [n_fish] rate_LInf_fish;
  vector<lower=0> [n_eco] rate_LInf_eco;
  vector<lower=0> [n_fish] rate_k_fish;
  vector<lower=0> [n_eco] rate_k_eco;
  }
  
model {
  //storage
  vector[N] mu;
  //priors
  LInf ~ gamma(0.001,0.001);
  rate_LInf_fish ~ gamma(0.001,0.001);
  rate_LInf_eco ~ gamma(0.001,0.001);
  k ~ gamma(0.001,0.001);
  rate_k_fish ~ gamma(0.001,0.001);
  rate_k_eco ~ gamma(0.001,0.001);
  tZero ~ cauchy(0,5);
  sigma ~ uniform(0, 100);
  
  //VBGM likelyhood
  for(i in 1:N){
	    mu[i] = LInf_i[fish_id[i]]*(1-exp(-k_i[fish_id[i]]*(Agei[i]-tZero_i[fish_id[i]])));
	    target += gamma_lpdf(Radmm[i]|mu[i], sigma); //likelihood
	  }
	for(i in 1:n_fish){
	  LInf_i[i] ~ gamma(LInf_eco[eco_id[i]]*rate_LInf_fish[eco_id[i]],rate_LInf_fish[eco_id[i]]);
	  k_i[i] ~ gamma(k_eco[eco_id[i]]*rate_k_fish[eco_id[i]],rate_k_fish[eco_id[i]]);
	  tZero_i[i] ~ cauchy(0,5);
	}
	for(i in 1:n_eco){
	  LInf_eco[i] ~ gamma(LInf*rate_LInf_eco,rate_LInf_eco);
	  k_eco[i] ~ gamma(k*rate_k_eco, rate_k_eco);
	  tZero_eco[i] ~ cauchy(0,5);
	}
  }
