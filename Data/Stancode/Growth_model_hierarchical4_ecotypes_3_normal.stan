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
  real<lower=0> phi;
  //hyperprameters
  real<lower=0> LInf;
  real<lower=0> LInf_sig;
  real<lower=0> k;
  real<lower=0> k_sig;
  real<upper=0> tZero;
  real<lower=0> tZero_sig;
  vector<lower=0> [n_eco] sig_LInf_eco;
  vector<lower=0> [n_eco] sig_k_eco;
  vector<lower=0> [n_eco] tZero_sig_eco;
  }
  
transformed parameters {
    //storage
  vector[N] mu;
  vector <lower = 0> [N] alpha;
	vector <lower = 0> [N] beta;
   //VBGM likelyhood
  for(i in 1:N){
	    mu[i] = LInf_i[fish_id[i]]*(1-exp(-k_i[fish_id[i]]*(Agei[i]-tZero_i[fish_id[i]])));
	  }
	  
	alpha = square(mu)/phi;
	beta = mu/phi;
}
  
model {
  Radmm ~ gamma(alpha,beta); //likelihood
  
  //hierarchical priors
  for(i in 1:n_fish){
	  LInf_i[i] ~ normal(LInf_eco[eco_id[i]],sig_LInf_eco[eco_id[i]]);
	  k_i[i] ~ normal(k_eco[eco_id[i]],sig_k_eco[eco_id[i]]);
	  tZero_i[i] ~ normal(tZero_eco[eco_id[i]],tZero_sig_eco[eco_id[i]]);
	}
	for(i in 1:n_eco){
	  LInf_eco[i] ~ normal(LInf,LInf_sig);
	  k_eco[i] ~ normal(k,k_sig);
	  tZero_eco[i] ~ normal(tZero,tZero_sig);
	}
  
  //hyperpriors
  LInf ~ normal(0,2.5);
  LInf_sig ~ gamma(0.001, 0.001);
  k ~ uniform(0,1);
  k_sig ~ gamma(0.001,0.001);
  tZero ~ cauchy(0,5);
  tZero_sig ~ gamma(0.001,0.001);
  sig_LInf_eco ~ gamma(0.001,0.001);
  sig_k_eco ~ gamma(0.001,0.001);
  tZero_sig_eco ~ gamma(0.001,0.001);
  phi ~ gamma(0.001, 0.001);
  }
