// Lactin-2 model

// set up a function to calculate the growth rate from the parameters

functions{

  real lactin2(real T, real rho, real Tmax, real delta, real lambda){

    real r;
    r = exp(rho*T) - exp(rho*Tmax - (Tmax - T)/delta) + lambda;
    return r;

  }

}

// data block

data{
  int N;
  array[N] real growth;
  array[N] real temperature;
}

// parameters block

parameters{
  real<lower = 0> rho;
  real Tmax;
  real<lower = 0> delta_t;
  real lambda;
  real<lower = 0> sigma;
}

// model block

model{
  
  // priors
  rho ~ normal(0,0.1);
  Tmax ~ normal(40,5);
  delta_t ~ normal(0,5);
  lambda ~ normal(0,2.5);
  sigma ~ normal(0,2.5);
  
  // likelihood
  
  for(i in 1:N){
  growth[i] ~ normal(lactin2(temperature[i], rho, Tmax, delta_t, lambda), sigma);
  }
  
  
}
  












