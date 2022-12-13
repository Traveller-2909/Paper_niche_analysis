//Lifelong growth dynamics of pike ecotypes
//Hierarchical design with three levels:
//   (1) multiple data (otolith annual growth marks) per fish 
//   (2) several fish per ecotype
//   (3) several ecotypes
// Main objective: Looking for differences in Linf and k between the ecotypes

data {
  int<lower=0> N;//Number of observations
  
  //data vectors
  vector<lower=0, upper=13> [N]Agei; //Age in years
  vector<lower=0> [N]Radmm; //otolith growth increments (in mm)
  
  //groups
  int<lower=1> n_fish; //number of individual fish
  int<lower=1> n_eco; //number of ecotypes
  int<lower=1, upper=n_fish> fish_id [N];
  int<lower=1, upper=n_eco> eco_id[N];
}

parameters {
  vector<lower=0> [n_fish] LInf_i;  //Mean maximum adult size at fish level
  vector<lower=0> [n_fish] k_i;     //Growth completion parameter at fish level
  vector<upper=0> [n_fish] tZero_i; //Age at wich size = 0 at fish level (size can´t be zero, so has to be negative)
  vector<lower=0> [n_eco] LInf_eco; //Mean maximum adult size at ecotype level
  vector<lower=0> [n_eco] k_eco;    //Growth completion parameter at ecotype level
  vector<upper=0> [n_eco] tZero_eco;//Age at wich size = 0 at ecotype level
  real<lower=0> phi;
  //hyperprameters
  real<lower=0> LInf;     //Baseline expectation for mean maximum adult size
  real<lower=0> LInf_sig; //Baseline SD for mean maximum adult size
  real<lower=0> k;        //Baseline expectation for growth completion parameter
  real<lower=0> k_sig;    //Baseline SD for growth completion parameter
  real<upper=0> tZero;    //Baseline expectation for age at wich size = 0 
  real<lower=0> tZero_sig;//Baseline SD for age at wich size = 0 
  vector<lower=0> [n_eco] LInf_eco_sig; //SD for LInf at ecotype level
  vector<lower=0> [n_eco] k_eco_sig;    //SD for k at ecotype level
  vector<lower=0> [n_eco] tZero_eco_sig;//SD for tZero at ecotype level
  }
  
transformed parameters {
    //storage vectors
  vector[N] mu;
  vector <lower = 0> [N] alpha;
	vector <lower = 0> [N] beta;
	
   //VBGM Expectation for Length at age (fish level)
  for(i in 1:N){
	    mu[i] = LInf_i[fish_id[i]]*(1-exp(-k_i[fish_id[i]]*(Agei[i]-tZero_i[fish_id[i]])));
	  }
	  
	alpha = square(mu)/phi;
	beta = mu/phi;
}
  
model {
  Radmm ~ gamma(alpha,beta);
  
  //priors for fish level
  for(i in 1:n_fish){
	  LInf_i[i] ~ normal(LInf_eco[eco_id[i]],LInf_eco_sig[eco_id[i]]);
	  k_i[i] ~ normal(k_eco[eco_id[i]],k_eco_sig[eco_id[i]]);
	  tZero_i[i] ~ normal(tZero_eco[eco_id[i]],tZero_eco_sig[eco_id[i]]);
	}
	//priors for ecotype level
	for(i in 1:n_eco){
	  LInf_eco[i] ~ normal(LInf,LInf_sig);
	  k_eco[i] ~ normal(k,k_sig);
	  tZero_eco[i] ~ normal(tZero,tZero_sig);
	}
  
  //hyperpriors
  LInf ~ normal(0,2.5); //2.3 mm was the largest otolith radius in sample
  LInf_sig ~ gamma(0.001, 0.001);
  k ~ uniform(0,1); //let the model choose k (literature suggestion)
  k_sig ~ gamma(0.001,0.001);
  tZero ~ cauchy(0,5); //commonly used for this param
  tZero_sig ~ gamma(0.001,0.001);
  LInf_eco_sig ~ gamma(0.001,0.001);
  k_eco_sig ~ gamma(0.001,0.001);
  tZero_eco_sig ~ gamma(0.001,0.001);
  phi ~ gamma(0.001, 0.001);
  }
