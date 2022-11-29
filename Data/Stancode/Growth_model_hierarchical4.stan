
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
  int<lower=1, upper=n_groupmax> group_id [N];
  int<lower=0, upper=1> Sexbin[N];
}

parameters {
  vector<lower=0> [n_group1] LInf1;
  vector<lower=0> [n_group1] LInf2;
  vector<lower=0> [n_group1] k1;
  vector<lower=0> [n_group2] k2;
  vector <upper=0> [n_group1] tZero1; 
  vector <upper=0> [n_group2] tZero2;
  real<lower=0> phi;
  //hyperprameters
  real<lower=0> k1_mu;
  real<lower=0> k1_sig;
  real<lower=0> k2_mu;
  real<lower=0> k2_sig;
  real<lower=0> LInf1_mu;
  real<lower=0> LInf1_sig;
  real<lower=0> LInf2_mu;
  real<lower=0> LInf2_sig;
  }
transformed parameters{
  vector <lower = 0> [N] mu;
	vector <lower = 0> [N] alpha;
	vector <lower = 0> [N] beta;
	
	for(i in 1:N){
	  if(Sexbin[i]==0){
	    mu[i] = LInf1[group_id[i]]*(1-exp(-k1[group_id[i]]*(Agei[i]-tZero1[group_id[i]])));
	  }
	  else{mu[i] = LInf2[group_id[i]]*(1-exp(-k2[group_id[i]]*(Agei[i]-tZero2[group_id[i]])));
	  }
	}

	alpha = square(mu)/phi;
	beta = mu/phi;
}

model {
  Radmm ~ gamma(alpha, beta);
  LInf1 ~ normal(LInf1_mu, LInf1_sig); 
  LInf2 ~ normal(LInf2_mu, LInf2_sig); 
  k1 ~ normal(k1_mu, k1_sig); 
  k2 ~ normal(k2_mu, k2_sig);
  tZero1 ~ normal(tZero_mu, tZero_sig);
  tZero2 ~ normal(tZero_mu, tZero_sig);
  phi~gamma(0.001, 0.001);
  //hyperpriors
  k1_mu ~ gamma(0.001, 0.001);
  k1_sig ~ gamma(0.001, 0.001);
  k2_mu ~ gamma(0.001, 0.001);
  k2_sig ~ gamma(0.001, 0.001);
  LInf1_mu ~ gamma(0.001, 0.001);
  LInf1_sig ~ gamma(0.001, 0.001);
  LInf2_mu ~ gamma(0.001, 0.001);
  LInf2_sig ~ gamma(0.001, 0.001);
  }

