
data {
  int<lower=0> N;//Number of samples
  
  //data vectors
  vector<lower=0, upper=13> [N]Agei;
  vector<lower=0> [N]Radmm;
  real<upper=0> tZero_mu;
  real<lower=0> tZero_sig;
  
  //groups
  int<lower=1> n_groupmax;
  int<lower=1> n_group1;
  int<lower=1> n_group2;
  int<lower=1> n_group3;
  int<lower=1> n_group4;
  int<lower=1> n_group5;
  int<lower=1> n_group6;
  int<lower=1, upper=n_groupmax> group_id [N];
  int<lower=0, upper=6> ecotype[N];
}

parameters {
  vector<lower=0> [n_group1] LInf1;
  vector<lower=0> [n_group2] LInf2;
  vector<lower=0> [n_group3] LInf3;
  vector<lower=0> [n_group4] LInf4;
  vector<lower=0> [n_group5] LInf5;
  vector<lower=0> [n_group6] LInf6;
  vector<lower=0> [n_group1] k1;
  vector<lower=0> [n_group2] k2;
  vector<lower=0> [n_group3] k3;
  vector<lower=0> [n_group4] k4;
  vector<lower=0> [n_group5] k5;
  vector<lower=0> [n_group6] k6;
  vector <upper=0> [n_group1] tZero1; 
  vector <upper=0> [n_group2] tZero2;
  vector <upper=0> [n_group3] tZero3;
  vector <upper=0> [n_group4] tZero4;
  vector <upper=0> [n_group5] tZero5;
  vector <upper=0> [n_group6] tZero6;
  real<lower=0> phi;
  //hyperprameters
  real<lower=0> k1_mu;
  real<lower=0> k1_sig;
  real<lower=0> k2_mu;
  real<lower=0> k2_sig;
  real<lower=0> k3_mu;
  real<lower=0> k3_sig;
  real<lower=0> k4_mu;
  real<lower=0> k4_sig;
  real<lower=0> k5_mu;
  real<lower=0> k5_sig;
  real<lower=0> k6_mu;
  real<lower=0> k6_sig;
  real<lower=0> LInf1_mu;
  real<lower=0> LInf1_sig;
  real<lower=0> LInf2_mu;
  real<lower=0> LInf2_sig;
  real<lower=0> LInf3_mu;
  real<lower=0> LInf3_sig;
  real<lower=0> LInf4_mu;
  real<lower=0> LInf4_sig;
  real<lower=0> LInf5_mu;
  real<lower=0> LInf5_sig;
  real<lower=0> LInf6_mu;
  real<lower=0> LInf6_sig;
  }
transformed parameters{
  vector <lower = 0> [N] mu;
	vector <lower = 0> [N] alpha;
	vector <lower = 0> [N] beta;
	
	for(i in 1:N){
	  if(ecotype[i]==1){
	    mu[i] = LInf1[group_id[i]]*(1-exp(-k1[group_id[i]]*(Agei[i]-tZero1[group_id[i]])));
	  }
	  else{
	    if(ecotype[i]==2){
	      mu[i] = LInf2[group_id[i]]*(1-exp(-k2[group_id[i]]*(Agei[i]-tZero2[group_id[i]])));
	    }
	    else{
	    if(ecotype[i]==3){
	      mu[i] = LInf3[group_id[i]]*(1-exp(-k3[group_id[i]]*(Agei[i]-tZero3[group_id[i]])));
	    }
	    else{
	    if(ecotype[i]==4){
	      mu[i] = LInf4[group_id[i]]*(1-exp(-k4[group_id[i]]*(Agei[i]-tZero4[group_id[i]])));
	    }
	    else{
	    if(ecotype[i]==5){
	      mu[i] = LInf5[group_id[i]]*(1-exp(-k5[group_id[i]]*(Agei[i]-tZero5[group_id[i]])));
	    }
	    else{
	    if(ecotype[i]==6){
	      mu[i] = LInf6[group_id[i]]*(1-exp(-k6[group_id[i]]*(Agei[i]-tZero6[group_id[i]])));
	    }
	    }
	    }
	    }
	    }
	  }
	}

	alpha = square(mu)/phi;
	beta = mu/phi;
}

model {
  Radmm ~ gamma(alpha, beta);
  LInf1 ~ normal(LInf1_mu, LInf1_sig); 
  LInf2 ~ normal(LInf2_mu, LInf2_sig);
  LInf3 ~ normal(LInf3_mu, LInf3_sig);
  LInf4 ~ normal(LInf4_mu, LInf4_sig);
  LInf5 ~ normal(LInf5_mu, LInf5_sig);
  LInf6 ~ normal(LInf6_mu, LInf6_sig);
  k1 ~ normal(k1_mu, k1_sig); 
  k2 ~ normal(k2_mu, k2_sig);
  k3 ~ normal(k3_mu, k3_sig);
  k4 ~ normal(k4_mu, k4_sig);
  k5 ~ normal(k5_mu, k5_sig);
  k6 ~ normal(k6_mu, k6_sig);
  tZero1 ~ normal(tZero_mu, tZero_sig);
  tZero2 ~ normal(tZero_mu, tZero_sig);
  tZero3 ~ normal(tZero_mu, tZero_sig);
  tZero4 ~ normal(tZero_mu, tZero_sig);
  tZero5 ~ normal(tZero_mu, tZero_sig);
  tZero6 ~ normal(tZero_mu, tZero_sig);
  phi~gamma(0.001, 0.001);
  //hyperpriors
  k1_mu ~ gamma(0.001, 0.001);
  k1_sig ~ gamma(0.001, 0.001);
  k2_mu ~ gamma(0.001, 0.001);
  k2_sig ~ gamma(0.001, 0.001);
  k3_mu ~ gamma(0.001, 0.001);
  k3_sig ~ gamma(0.001, 0.001);
  k4_mu ~ gamma(0.001, 0.001);
  k4_sig ~ gamma(0.001, 0.001);
  k5_mu ~ gamma(0.001, 0.001);
  k5_sig ~ gamma(0.001, 0.001);
  k6_mu ~ gamma(0.001, 0.001);
  k6_sig ~ gamma(0.001, 0.001);
  LInf1_mu ~ gamma(0.001, 0.001);
  LInf1_sig ~ gamma(0.001, 0.001);
  LInf2_mu ~ gamma(0.001, 0.001);
  LInf2_sig ~ gamma(0.001, 0.001);
  LInf3_mu ~ gamma(0.001, 0.001);
  LInf3_sig ~ gamma(0.001, 0.001);
  LInf4_mu ~ gamma(0.001, 0.001);
  LInf4_sig ~ gamma(0.001, 0.001);
  LInf5_mu ~ gamma(0.001, 0.001);
  LInf5_sig ~ gamma(0.001, 0.001);
  LInf6_mu ~ gamma(0.001, 0.001);
  LInf6_sig ~ gamma(0.001, 0.001);
  }

