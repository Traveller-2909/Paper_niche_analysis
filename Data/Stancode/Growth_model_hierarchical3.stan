
data {
  int<lower=0> N;//Number of samples
  
  //data vectors
  vector<lower=0, upper=13> [N]Agei;
  vector<lower=0> [N]Lengthi;
  real<lower=0> alinf;
  real<lower=0> blinf;
  
  //groups
  int<lower=1> n_groupmax;
  int<lower=1> n_group;
  int<lower=1> n_group1;
  int<lower=1, upper=n_groupmax> group_id [N];
  int<lower=0, upper=1> Sexbin[N];
}

parameters {
  real<lower=0> LInf;
  real<lower=0> LInf1;
  real<lower=0> phi;
  vector<lower=0> [n_group] k;
  vector<lower=0> [n_group1] k1;
  real <upper=0>tZero; 
  real <upper=0>tZero1;
  //hyperprameters
  real<lower=0> k_mu;
  real<lower=0> k_sig;
  real<lower=0> k1_mu;
  real<lower=0> k1_sig;
  }
transformed parameters{
  vector <lower = 0> [N] mu;
	vector <lower = 0> [N] alpha;
	vector <lower = 0> [N] beta;
	
	for(i in 1:N){
	  if(Sexbin[i]==0){
	    mu[i] = LInf*(1-exp(-k[group_id[i]]*(Agei[i]-tZero)));
	  }
	  else{mu[i] = LInf1*(1-exp(-k1[group_id[i]]*(Agei[i]-tZero1)));
	  }
	}

	alpha = square(mu)/phi;
	beta = mu/phi;
}

model {
  Lengthi ~ gamma(alpha, beta);
  LInf ~ gamma(alinf, blinf); 
  LInf1 ~ gamma(alinf, blinf); 
  k ~ normal(k_mu, k_sig); 
  k1 ~ normal(k1_mu, k1_sig);
  tZero ~ cauchy(0,5);
  tZero1 ~ cauchy(0,5);
  phi~gamma(0.01, 0.01);
  //hyperpriors
  k_mu ~ gamma(0.01, 0.01);
  k_sig ~ gamma(0.01, 0.01);
  k1_mu ~ gamma(0.01, 0.01);
  k1_sig ~ gamma(0.01, 0.01);
  }

