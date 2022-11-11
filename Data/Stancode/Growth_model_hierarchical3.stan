
data {
  int<lower=0> N;//Number of samples
  
  //data vectors
  vector<lower=0, upper=13> [N]Agei;
  vector<lower=0> [N]Radmm;
  real<lower=0> alinf;
  real<lower=0> blinf;
  
  //groups
  int<lower=1> n_groupmax;
  int<lower=1> n_group1;
  int<lower=1> n_group2;
  int<lower=1, upper=n_groupmax> group_id [N];
  int<lower=0, upper=1> Sexbin[N];
}

parameters {
  real<lower=0> LInf1;
  real<lower=0> LInf2;
  real<lower=0> phi;
  vector<lower=0> [n_group1] k1;
  vector<lower=0> [n_group2] k2;
  real <upper=0>tZero1; 
  real <upper=0>tZero2;
  //hyperprameters
  real<lower=0> k1_mu;
  real<lower=0> k1_sig;
  real<lower=0> k2_mu;
  real<lower=0> k2_sig;
  }
transformed parameters{
  vector <lower = 0> [N] mu;
	vector <lower = 0> [N] alpha;
	vector <lower = 0> [N] beta;
	
	for(i in 1:N){
	  if(Sexbin[i]==0){
	    mu[i] = LInf1*(1-exp(-k1[group_id[i]]*(Agei[i]-tZero1)));
	  }
	  else{mu[i] = LInf2*(1-exp(-k2[group_id[i]]*(Agei[i]-tZero2)));
	  }
	}

	alpha = square(mu)/phi;
	beta = mu/phi;
}

model {
  Radmm ~ gamma(alpha, beta);
  LInf1 ~ gamma(alinf, blinf); 
  LInf2 ~ gamma(alinf, blinf); 
  k1 ~ normal(k1_mu, k1_sig); 
  k2 ~ normal(k2_mu, k2_sig);
  tZero1 ~ cauchy(0,5);
  tZero2 ~ cauchy(0,5);
  phi~gamma(0.01, 0.01);
  //hyperpriors
  k1_mu ~ gamma(0.01, 0.01);
  k1_sig ~ gamma(0.01, 0.01);
  k2_mu ~ gamma(0.01, 0.01);
  k2_sig ~ gamma(0.01, 0.01);
  }

